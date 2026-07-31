#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=ku04-gim
#SBATCH --partition=GPU
#SBATCH --gres=gpu
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
# This allocates the walltime to 10  days. The program will not run for longer.
#SBATCH --time=240:00:00
# This sets the quality of service to 'elevated'
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703
#----------------------------------------------------------

#--------------------------------------------------
# USER SETTINGS
#--------------------------------------------------

# Toggle equilibration 
RUN_EQUIL=true

# Toggle Production run
RUN_PROD=true #Set to false if only PP is desired

RESUME_MD=false #Set to true if resuming from chkpt

EM_OUT="em" # Output prefix from energy minimization step
TOPOL="ku04_gim_multi" # Topology file prefix
EQUI_MDP="npt" # MDP file prefix for equilibration
MD_MDP="md" # MDP file prefix for production run


# Singularity container
SIF_PATH="${HOME}/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif"

#--------------------------------------------------
# INITIALIZATION
#--------------------------------------------------


LOGFILE="wrapper_$(date +%Y%m%d_%H%M%S).log"

# Persistent unbuffered logging
exec > >(stdbuf -oL tee -a "$LOGFILE")
exec 2>&1

echo "================================================="
echo "SLURM JOB START"
echo "================================================="
echo "Date        : $(date)"
echo "Hostname    : $(hostname)"
echo "Working dir : $(pwd)"
echo "SLURM_JOBID   : ${SLURM_JOB_ID}"
echo "SLURM_JOBNAME : ${SLURM_JOB_NAME}"
echo "================================================="

export SLURM_NCPUS=$SLURM_CPUS_PER_TASK

echo "[INFO] Detected SLURM_NCPUS=${SLURM_NCPUS}"

#--------------------------------------------------
# SINGULARITY WRAPPER
#--------------------------------------------------

SINGULARITY() {
    singularity exec \
        --nv \
        -B "${PWD}:/host_pwd" \
        --pwd /host_pwd \
        "$SIF_PATH" "$@"
}

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
    echo "RUNNING EQUILIBRATION STAGE"
    echo "================================================="


    echo "[INFO] Running grompp for equilibration"

    SINGULARITY bash -c "
        gmx grompp \
            -f ${EQUI_MDP}.mdp \
            -c ${EM_OUT}.gro \
            -r ${EM_OUT}.gro \
            -p ${TOPOL}.top \
            -n index.ndx \
            -o npt.tpr \
            -maxwarn 5
    "

    test -f "npt.tpr"

    echo "[INFO] Starting mdrun for equilibration"

    SINGULARITY bash -c "
        gmx mdrun \
            -nobackup \
            -nocopyright \
            -v \
            -deffnm npt \
            -pin on \
            -nb gpu \
            -ntmpi 1 \
            -ntomp ${SLURM_NCPUS} \
            -gpu_id 0
    "

    echo "[INFO] mdrun completed for equilibration"

    test -f "npt.gro"
    test -f "npt.cpt"
    test -f "npt.log"

    echo ">>> Completed equilibration step"
    echo ">>> Checkpoint: npt.cpt"

    echo
    echo "================================================="
    echo "EQUILIBRATION  COMPLETE"
    echo "================================================="

else

    echo
    echo "================================================="
    echo "SKIPPING EQUILIBRATION"
    echo "RUN_EQUIL=false"
    echo "================================================="

fi
#--------------------------------------------------
# RESUME CHECKPOINT OR START NEW PRODUCTION
#--------------------------------------------------
if [ "$RESUME_MD" = true ]; then

    echo
    echo "================================================="
    echo "RESUMING PRODUCTION FROM CHECKPOINT"
    echo "RESUME_MD=true"
    echo "================================================="

    if [ ! -f "md.cpt" ]; then
        echo "[ERROR] Checkpoint file md.cpt not found. Cannot resume."
        exit 1
    fi
    echo "[INFO] resuming from checkpoint md.cpt"

    SINGULARITY bash -c "
        gmx mdrun \
            -nobackup \
            -nocopyright \
            -v \
            -deffnm md \
            -pin on \
            -cpi md.cpt \
            -ntmpi 1 \
            -ntomp ${SLURM_NCPUS} \
            -nb gpu \
            -gpu_id 0
    "
else
    if [ "$RUN_PROD" = true ]; then
        echo
        echo "================================================="
        echo "STARTING NEW PRODUCTION RUN"
        echo "RESUME_MD=false"
        echo "================================================="
        
        echo "[INFO] Running grompp for production"

        SINGULARITY bash -c "
            gmx grompp \
                -f ${MD_MDP}.mdp \
                -c npt.gro \
                -p ${TOPOL}.top \
                -n index.ndx \
                -o md.tpr \
                -t npt.cpt \
                -maxwarn 2
        "

        test -f "md.tpr"

        echo "[INFO] Starting mdrun for production"

        SINGULARITY bash -c "
            gmx mdrun \
                -nobackup \
                -nocopyright \
                -v \
                -deffnm md \
                -pin on \
                -nb gpu \
                -ntmpi 1 \
                -ntomp ${SLURM_NCPUS} \
                -gpu_id 0
        "
        echo "[INFO] mdrun completed for production"

        test -f "md.gro"
        test -f "md.cpt"
        test -f "md.log"

        echo ">>> Checkpoint: md.cpt"

        echo
        echo "================================================="
        echo "PRODUCTION COMPLETE"
        echo "================================================="
    else
        echo
        echo "================================================="
        echo "SKIPPING PRODUCTION RUN"
        echo "RUN_PROD=false"
        echo "================================================="
    fi        
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

xtc_list=$(ls -v md*.xtc 2>/dev/null | tr '\n' ' ')

if [ -z "${xtc_list}" ]; then

    echo "[WARNING] No production XTC files found."
    echo "[WARNING] Skipping post-processing."

else

    echo "[INFO] Concatenating trajectories sequentially"

    # Removed -cat because the chunks should possess native, sequential time markers
    SINGULARITY bash -c "
        gmx trjcat \
            -f ${xtc_list} \
            -o traj_prod.xtc
    "

    test -f traj_prod.xtc

    echo ">>> Concatenated trajectory created: traj_prod.xtc"

    #--------------------------------------------------
    # STEP 1: Make molecules whole
    #--------------------------------------------------
    echo "System" > trjconv_whole.in
    echo "[INFO] Running whole-PBC repair"
    SINGULARITY bash -c "
        gmx trjconv \
            -s md.tpr \
            -f traj_prod.xtc \
            -pbc whole \
            -o traj_whole.xtc \
            < trjconv_whole.in
    "
    test -f traj_whole.xtc

    #--------------------------------------------------
    # STEP 2: Center the membrane in the box (NEW STEP)
    #--------------------------------------------------
    # We need to provide TWO inputs for centering:
    # 1. The group to center on (e.g., Membrane or POPE_POPG...)
    # 2. The group to output (System)
    # Change "Membrane" below to match your actual index group name if needed
    echo "Membrane\nSystem" > trjconv_center.in
    echo "[INFO] Centering the bilayer sheet"
    SINGULARITY bash -c "
        gmx trjconv \
            -s md.tpr \
            -f traj_whole.xtc \
            -n index.ndx \
            -center \
            -pbc mol \
            -ur compact \
            -boxcenter rect \
            -o traj_centered.xtc \
            < trjconv_center.in
    "
    test -f traj_centered.xtc

    #--------------------------------------------------
    # STEP 3: Remove frame-skipping jumps
    #--------------------------------------------------
    echo "System" > trjconv_nojump.in
    echo "[INFO] Running nojump to remove PBC-crossing artifacts"
    SINGULARITY bash -c "
        gmx trjconv \
            -s md.tpr \
            -f traj_centered.xtc \
            -pbc nojump \
            -o traj_nojump.xtc \
            < trjconv_nojump.in
    "
    test -f traj_nojump.xtc

    # Clean up intermediate files
    rm -f traj_whole.xtc traj_centered.xtc
    rm -f trjconv_whole.in trjconv_center.in trjconv_nojump.in
fi

#--------------------------------------------------
# END TIMER
#--------------------------------------------------

end=$(date +%s.%N)

RUNTIME=$(echo "$end - $start" | bc -l)

echo
echo "================================================="
echo "JOB COMPLETE"
echo "================================================="
echo "Runtime: ${RUNTIME} sec"
echo "================================================="
