#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --partition=GPU
#SBATCH --gres=mps:50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
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
export TC_GRPS_CMPLX="Protein_RIF"  # Change this to match the ligand name in the complex

# Python script to run
HPC_SCRIPT="08_hpc_equi_prod.py"

# ----------------------------------------------------------
# LOGGING
# ----------------------------------------------------------
LOGFILE="hpc_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1

echo "================================================="
echo "  SLURM JOB START — ${SLURM_JOB_NAME} Simulation"
echo "================================================="
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  WorkDir   : $(pwd)"
echo "  SLURM_JOB : ${SLURM_JOB_ID}"
echo "  CPUs      : ${SLURM_CPUS_PER_TASK}"
echo "  SIF       : ${SIF_HPC}"
echo "================================================="

# ----------------------------------------------------------
# SANITY CHECKS
# ----------------------------------------------------------
if [ ! -f "$SIF_HPC" ]; then
    echo "[ERROR] Singularity image not found: $SIF_HPC"
    exit 1
fi

if ! ls checkpoint.em_*.pycpt 1>/dev/null 2>&1; then
    echo "[ERROR] No EM checkpoint found. Run Steps 1-7 locally first."
    exit 1
fi

if [ ! -f "complex_solv/index.ndx" ]; then
    echo "[ERROR] Index file not found: complex_solv/index.ndx"
    echo "        Run Step 6 to create the index file."
    exit 1
fi

if [ ! -f "$HPC_SCRIPT" ]; then
    echo "[ERROR] HPC script not found: $HPC_SCRIPT"
    exit 1
fi

# ----------------------------------------------------------
# TIMER START
# ----------------------------------------------------------
start=$(date +%s.%N)

# ----------------------------------------------------------
# RUN PYTHON SCRIPT INSIDE SINGULARITY
# ----------------------------------------------------------
echo ""
echo "Starting Singularity container …"

singularity exec \
    --nv \
    -B "${PWD}:/host_pwd" \
    --pwd /host_pwd \
    "$SIF_HPC" \
    python3 "$HPC_SCRIPT"

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
