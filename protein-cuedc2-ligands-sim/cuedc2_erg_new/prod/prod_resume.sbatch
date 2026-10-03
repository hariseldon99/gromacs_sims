#!/bin/bash
#SBATCH --job-name=cuedc2-erg
#SBATCH --partition=GPU
#SBATCH --gres=gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --sockets-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --qos=elevated
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703
# =============================================================================
# SLURM job script — CUEDC2 + ERGOTAMINE protein-ligand MD
# Runs equilibration (5.5 ns total) + production (100 ns) inside Singularity.
#
# USAGE: Submit from the simulations/ directory on HPC:
#   cd /path/to/simulations
#   sbatch 09_hpc_submit.sbatch
#
# PREREQUISITES:
#   - checkpoint.em_YYYYMMDD.pycpt must exist in simulations/
#   - complex_solv/index.ndx must exist
#   - All input files from Steps 1-7 must be present
#
# To adapt for a different ligand, change LIGAND_NAME and ensure:
#   - The corresponding checkpoint exists
#   - The complex topology and coordinates reference that ligand
# =============================================================================

# ----------------------------------------------------------
# USER SETTINGS — adjust before submitting
# ----------------------------------------------------------
# Singularity image on HPC
SIF_HPC="${HOME}/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif"

# Python script to run

# ----------------------------------------------------------
# LOGGING
# ----------------------------------------------------------
LOGFILE="hpc_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1

echo "================================================="
echo "  SLURM JOB START — CUEDC2-ERG Simulation"
echo "================================================="
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  WorkDir   : $(pwd)"
echo "  SLURM_JOB : ${SLURM_JOB_ID}"
echo "  CPUs      : ${SLURM_CPUS_PER_TASK}"
echo "  SIF       : ${SIF_HPC}"
echo "================================================="

# ----------------------------------------------------------
# TIMER START
# ----------------------------------------------------------
start=$(date +%s.%N)

# ----------------------------------------------------------
# RUN PYTHON SCRIPT INSIDE SINGULARITY
# ----------------------------------------------------------
echo ""
echo "Starting Singularity container …"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
srun --cpu-bind=cores singularity exec \
    --nv \
    -B "${PWD}:/host_pwd" \
    --pwd /host_pwd \
    "$SIF_HPC" \
    gmx mdrun -pin off -ntomp ${SLURM_CPUS_PER_TASK} -ntmpi 1 -gpu_id 0 -cpi prod_cuedc2_erg.cpt -deffnm prod_cuedc2_erg -nb gpu -pme gpu -bonded gpu

EXIT_CODE=$?

# ----------------------------------------------------------
# TIMER END
# ----------------------------------------------------------
end=$(date +%s.%N)
RUNTIME=$(echo "$end - $start" | bc -l)

echo ""
echo "================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "  JOB COMPLETE (exit 0)"
else
    echo "  JOB FAILED (exit ${EXIT_CODE})"
fi
echo "  Runtime: ${RUNTIME} s"
echo "================================================="

exit $EXIT_CODE
