#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=ku04-gim-ext
#SBATCH --partition=GPU
#SBATCH --gres=gpu
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# This allocates the walltime to 10  days. The program will not run for longer.
#SBATCH --time=240:00:00
# This sets the quality of service to 'elevated'
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703
#----------------------------------------------------------

set -euo pipefail

#--------------------------------------------------
# USER SETTINGS
#--------------------------------------------------

# Toggle Simulation Extension
EXTEND_MD=false

# Checkpoint and Topology settings
OLD_TPR="md.tpr"
EXT_TPR="md_extend.tpr"
CPT_FILE="md.cpt"
INDEX_FILE="index.ndx"
DRY_INDEX_FILE="index_dry.ndx"
EXTEND_TIME=4000000 # extend by time in ps (4 microseconds)

SIF_PATH="${HOME}/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif"

#--------------------------------------------------
# INITIALIZATION
#--------------------------------------------------

LOGFILE="wrapper_extend_$(date +%Y%m%d_%H%M%S).log"

# Persistent unbuffered logging
exec > >(stdbuf -oL tee -a "$LOGFILE")
exec 2>&1

echo "================================================="
echo "SLURM EXTENSION JOB START"
echo "================================================="
echo "Date        : $(date)"
echo "Hostname    : $(hostname)"
echo "Working dir : $(pwd)"
echo "SLURM_JOBID : ${SLURM_JOB_ID}"
echo "SLURM_JOBNAME : ${SLURM_JOB_NAME}"
echo "================================================="


echo "[INFO] Detected SLURM_NCPUS=${SLURM_CPUS_PER_TASK}"

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
# EXTEND PRODUCTION RUN
#--------------------------------------------------
if [ "$EXTEND_MD" = true ]; then
    echo
    echo "================================================="
    echo "EXTENDING MD SIMULATION BY ${EXTEND_TIME} ps (${EXTEND_TIME} microseconds)"
    echo "================================================="

    # Validation check for required base files
    test -f "$OLD_TPR"
    test -f "$CPT_FILE"

    echo "[INFO] Modifying the run input file to extend time..."
    SINGULARITY bash -c "
        gmx convert-tpr \
            -s ${OLD_TPR} \
            -extend ${EXTEND_TIME} \
            -o ${EXT_TPR}
    "

    test -f "$EXT_TPR"
    echo "[INFO] Successfully generated extended topology: ${EXT_TPR}"

    echo "[INFO] Launching mdrun with GPU offloading configuration..."
    # Configured with explicit GPU offloading parameters required by this cluster partition
    SINGULARITY bash -c "
        gmx mdrun \
            -nobackup \
            -nocopyright \
            -v \
            -s ${EXT_TPR} \
            -cpi ${CPT_FILE} \
            -deffnm md \
            -pin on \
            -nb gpu \
            -ntmpi 1 \
            -ntomp ${SLURM_CPUS_PER_TASK} \
            -gpu_id 0
    "
    echo "[INFO] mdrun extension completed successfully"

    test -f "md.gro"
    test -f "md.cpt"
    test -f "md.log"

    echo "================================================="
    echo "EXTENSION RUN COMPLETE"
    echo "================================================="
fi

#--------------------------------------------------
# POST-PROCESSING
#--------------------------------------------------
echo
echo "================================================="
echo "POST-PROCESSING"
echo "================================================="

# Clean past post-processed outputs to ensure data consistency
rm -f traj_prod.xtc traj_whole.xtc traj_centered.xtc traj_nojump.xtc traj_dry.xtc

xtc_list=$(ls -v md*.xtc 2>/dev/null | tr '\n' ' ')

if [ -z "${xtc_list}" ]; then
    echo "[WARNING] No production XTC files found. Skipping post-processing."
else
    echo "[INFO] Concatenating trajectories sequentially"
    SINGULARITY bash -c "
        gmx trjcat \
            -f ${xtc_list} \
            -o traj_prod.xtc
    "
    test -f traj_prod.xtc
    echo ">>> Concatenated trajectory created: traj_prod.xtc"

    # STEP 1: Make molecules whole
    echo "System" > trjconv_whole.in
    echo "[INFO] Running whole-PBC repair"
    SINGULARITY bash -c "
        gmx trjconv \
            -s ${EXT_TPR} \
            -f traj_prod.xtc \
            -pbc whole \
            -o traj_whole.xtc \
            < trjconv_whole.in
    "
    test -f traj_whole.xtc

    # STEP 2: Center the membrane in the box
    # Note: Swap "Membrane" to "MEMB" if your index group mirrors file 2
    echo "Membrane\nSystem" > trjconv_center.in
    echo "[INFO] Centering the bilayer sheet"
    SINGULARITY bash -c "
        gmx trjconv \
            -s ${EXT_TPR} \
            -f traj_whole.xtc \
            -n ${INDEX_FILE} \
            -center \
            -pbc mol \
            -ur compact \
            -boxcenter rect \
            -o traj_centered.xtc \
            < trjconv_center.in
    "
    test -f traj_centered.xtc

    # STEP 3: Remove frame-skipping jumps
    echo "System" > trjconv_nojump.in
    echo "[INFO] Running nojump to remove PBC-crossing artifacts"
    SINGULARITY bash -c "
        gmx trjconv \
            -s ${EXT_TPR} \
            -f traj_centered.xtc \
            -pbc nojump \
            -o traj_nojump.xtc \
            < trjconv_nojump.in
    "
    test -f traj_nojump.xtc
    # STEP 4: Remove water and ions for dry trajectory
    echo "DRY" > trjconv_dry.in
    echo "[INFO] Running trjconv to remove water and ions"
    SINGULARITY bash -c "
        gmx trjconv \
            -s ${EXT_TPR} \
            -f traj_nojump.xtc \
            -n ${DRY_INDEX_FILE} \
            -o traj_dry.xtc \
            < trjconv_dry.in
    "
    # Clean up intermediate input files
    rm -f trjconv_whole.in trjconv_center.in trjconv_nojump.in trjconv_dry.in
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
