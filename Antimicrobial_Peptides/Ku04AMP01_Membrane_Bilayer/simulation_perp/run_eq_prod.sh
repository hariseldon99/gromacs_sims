#!/bin/bash

set -euo pipefail

# Return true if a file exists at the exact path OR as a GROMACS -noappend .partXXXX variant
# Usage: part_exists "step7_2.gro"  ->  also matches step7_2.part*.gro
part_exists() { local b="${1%.*}" e="${1##*.}"; test -f "$1" || ls "${b}".part*."${e}" > /dev/null 2>&1; }

#--------------------------------------------------
# USER SETTINGS
#--------------------------------------------------

# Toggle equilibration stages (6.1 -> 6.6)
RUN_EQUIL=false
RUN_PROD=false

# Production chunk range: 10 ns each chunk
# For resume after completed step7_6, set to START_CHUNK=7, else 1
START_CHUNK=1
MAX_CHUNK=10
CHUNK_TIME=10000 # In picoseconds

#--------------------------------------------------
# INITIALIZATION
#--------------------------------------------------


LOGFILE="wrapper_$(date +%Y%m%d_%H%M%S).log"

# Persistent unbuffered logging
exec > >(stdbuf -oL tee -a "$LOGFILE")
exec 2>&1

echo "================================================="
echo "JOB START"
echo "================================================="
echo "Date        : $(date)"
echo "Hostname    : $(hostname)"
echo "Working dir : $(pwd)"
echo "================================================="

#--------------------------------------------------
# MODULES
#--------------------------------------------------


export NCPUS="47"

echo "[INFO] Detected NCPUS=${NCPUS}"

#--------------------------------------------------
# START TIMER
#--------------------------------------------------

start=$(date +%s.%N)

#--------------------------------------------------
# EQUILIBRATION
#--------------------------------------------------

if [ "$RUN_EQUIL" = true ]; then

    echo
    echo "================================================="
    echo "RUNNING EQUILIBRATION STAGES"
    echo "================================================="

    for step in 6.1 6.2 6.3 6.4 6.5 6.6; do

        prev=$(echo $step | awk -F. '{
            if ($2==1)
                print "em";
            else
                printf "step%s.%s", $1, $2-1
        }')

        echo
        echo "-------------------------------------------------"
        echo "[$(date)] Starting equilibration step${step}"
        echo "Previous structure: ${prev}"
        echo "-------------------------------------------------"

        test -f "${prev}.gro"

        echo "[INFO] Running grompp for step${step}"

        bash -c "
            gmx grompp \
                -f step${step}_equilibration.mdp \
                -c ${prev}.gro \
                -r system.gro \
                -p topol.top \
                -n index.ndx \
                -o step${step}.tpr \
                -maxwarn 5
        "

        test -f "step${step}.tpr"

        echo "[INFO] Starting mdrun for step${step}"

        bash -c "
            gmx mdrun \
                -nobackup \
                -nocopyright \
                -v \
                -deffnm step${step} \
                -pin on \
                -nb gpu \
                -pme gpu \
                -bonded gpu \
                -nt ${NCPUS} \
                -gpu_id 0
        "

        echo "[INFO] mdrun completed for step${step}"

        test -f "step${step}.gro"
        test -f "step${step}.cpt"
        test -f "step${step}.log"

        echo ">>> Completed equilibration step${step}"
        echo ">>> Checkpoint: step${step}.cpt"

    done

    echo
    echo "================================================="
    echo "ALL EQUILIBRATION STAGES COMPLETE"
    echo "================================================="

else

    echo
    echo "================================================="
    echo "SKIPPING EQUILIBRATION STAGES"
    echo "RUN_EQUIL=false"
    echo "================================================="

fi

if [ "$RUN_PROD" = true ]; then

    #--------------------------------------------------
    # PRODUCTION SAFETY CHECK
    #--------------------------------------------------

    if [ "$START_CHUNK" -gt 1 ]; then

        prev_chunk=$((START_CHUNK - 1))

        echo
        echo "[INFO] Resume mode detected"
        echo "[INFO] Expecting checkpoint from step7_${prev_chunk}"

        part_exists "step7_${prev_chunk}.gro" || { echo "[ERROR] No .gro found for step7_${prev_chunk} (checked plain and .partXXXX)"; exit 1; }
        test -f "step7_${prev_chunk}.cpt"
        test -f "step7_${prev_chunk}.tpr"

    fi

    #--------------------------------------------------
    # PRODUCTION MD
    #--------------------------------------------------

    cnt_start=$START_CHUNK
    set_cntmax=$MAX_CHUNK
    chunk_time=$CHUNK_TIME

    echo
    echo "================================================="
    echo "STARTING PRODUCTION RUNS"
    echo "Chunks ${cnt_start} -> ${set_cntmax}"
    echo "================================================="

    for cnt in $(seq ${cnt_start} ${set_cntmax}); do

        pcnt=$((cnt - 1))
        istep="step7_${cnt}"

        echo
        echo "-------------------------------------------------"
        echo "[$(date)] Starting production chunk ${cnt}/${set_cntmax}"
        echo "Current chunk : ${istep}"
        echo "-------------------------------------------------"

        if [ ${cnt} -eq 1 ]; then

            pstep="step6.6"

            test -f "${pstep}.gro"

            echo "[INFO] Running grompp for ${istep}"

            bash -c "
                gmx grompp \
                    -f step7_production.mdp \
                    -c ${pstep}.gro \
                    -p topol.top \
                    -n index.ndx \
                    -o ${istep}.tpr \
                    -maxwarn 2
            "

        else

            pstep="step7_${pcnt}"

            part_exists "${pstep}.gro" || { echo "[ERROR] No .gro found for ${pstep} (checked plain and .partXXXX)"; exit 1; }
            test -f "${pstep}.cpt"
            test -f "${pstep}.tpr"

            # Explicitly setting target time stops convert-tpr compounding steps incorrectly
            target_ps=$(( cnt * chunk_time ))

            echo "[INFO] Extending TPR for ${istep} until ${target_ps} ps"

            bash -c "
                gmx convert-tpr \
                    -s ${pstep}.tpr \
                    -until ${target_ps} \
                    -o ${istep}.tpr
            "
        fi

        test -f "${istep}.tpr"
        
        # Determine correct -cpi flag and whether to append
        if [ -f "${istep}.cpt" ]; then
            # Resume within this chunk (job was killed mid-chunk) — safe to append
            cpi_flag="-cpi ${istep}.cpt"
            noappend_flag=""
        elif [ ${cnt} -gt 1 ]; then
            # Continue from end of previous chunk — output files differ, must use -noappend
            cpi_flag="-cpi ${pstep}.cpt"
            noappend_flag="-noappend"
        else
            # Fresh start for chunk 1
            cpi_flag=""
            noappend_flag=""
        fi

        echo "[INFO] Starting mdrun for ${istep} (cpi: ${cpi_flag:-none}, noappend: ${noappend_flag:-no})"

        bash -c "
            gmx mdrun \
                -nobackup \
                -nocopyright \
                -v \
                -s ${istep}.tpr \
                -deffnm ${istep} \
                ${cpi_flag} \
                ${noappend_flag} \
                -pin on \
                -nb gpu \
                -pme gpu \
                -bonded gpu \
                -nt ${NCPUS} \
                -gpu_id 0
        "

        echo "[INFO] mdrun completed for ${istep}"

        part_exists "${istep}.gro" || { echo "[ERROR] No .gro found for ${istep} (checked plain and .partXXXX)"; exit 1; }
        test -f "${istep}.cpt"
        part_exists "${istep}.log" || { echo "[ERROR] No .log found for ${istep} (checked plain and .partXXXX)"; exit 1; }

        echo ">>> Completed production chunk ${cnt}/${set_cntmax}"
        echo ">>> Checkpoint: ${istep}.cpt"

    done

    echo
    echo "================================================="
    echo "ALL PRODUCTION CHUNKS COMPLETE"
    echo "================================================="

else

    echo
    echo "================================================="
    echo "SKIPPING PRODUCTION MD"
    echo "RUN_PROD=false"
    echo "================================================="

fi

#--------------------------------------------------
# POST-PROCESSING
#--------------------------------------------------

echo
echo "================================================="
echo "POST-PROCESSING"
echo "================================================="

rm -f traj_prod.xtc
rm -f traj_whole.xtc
rm -f traj_nojump.xtc
rm -f traj_dry.xtc

xtc_list=$(ls -v step7_*.xtc 2>/dev/null | tr '\n' ' ')

if [ -z "${xtc_list}" ]; then

    echo "[WARNING] No production XTC files found."
    echo "[WARNING] Skipping post-processing."

else

    echo "[INFO] Concatenating trajectories sequentially"

    # Removed -cat because the chunks now possess native, sequential time markers
    bash -c "
        gmx trjcat \
            -f ${xtc_list} \
            -o traj_prod.xtc
    "

    test -f traj_prod.xtc

    echo ">>> Concatenated trajectory created: traj_prod.xtc"

    echo "System" > trjconv_whole.in

    echo "[INFO] Running whole-PBC repair"

    bash -c "
        gmx trjconv \
            -s step7_1.tpr \
            -f traj_prod.xtc \
            -pbc whole \
            -o traj_whole.xtc \
            < trjconv_whole.in
    "

    test -f traj_whole.xtc

    echo ">>> Whole trajectory created: traj_whole.xtc"

    echo "System" > trjconv_nojump.in

    echo "[INFO] Running nojump to remove PBC-crossing artifacts"

    bash -c "
        gmx trjconv \
            -s step7_1.tpr \
            -f traj_whole.xtc \
            -pbc nojump \
            -o traj_nojump.xtc \
            < trjconv_nojump.in
    "

    test -f traj_nojump.xtc

    echo ">>> Nojump trajectory created: traj_nojump.xtc"

    printf "MEMB\nMEMB\nMEMB\n" > trjconv_dry.in

    echo "[INFO] Creating dry clustered trajectory"

    bash -c "
        gmx trjconv \
            -s step7_1.tpr \
            -f traj_nojump.xtc \
            -n index.ndx \
            -center \
            -pbc cluster \
            -ur compact \
            -o traj_dry.xtc \
            < trjconv_dry.in
    "

    test -f traj_dry.xtc

    echo ">>> Dry trajectory created: traj_dry.xtc"

    echo "[INFO] Creating dry TPR file matching dry clustered trajectory selections"

    printf "MEMB\nMEMB\n" > convert_tpr_dry.in

    bash -c "
        gmx convert-tpr \
            -s step7_1.tpr \
            -n index.ndx \
            -o traj_dry.tpr \
            < convert_tpr_dry.in
    "

    test -f traj_dry.tpr

    echo ">>> Dry TPR file created: traj_dry.tpr"

    rm -f traj_whole.xtc
    rm -f trjconv_whole.in
    rm -f trjconv_nojump.in
    rm -f trjconv_dry.in
    rm -f convert_tpr_dry.in
fi
