#!/bin/bash

# =============================================================================
# job script — CUEDC2 + ERGOTAMINE protein-ligand MD
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

export TC_GRPS_CMPLX="Protein_MID"  # Change this to match the ligand name in the complex

# Python script to run
HPC_SCRIPT="prod_resume.py"

# Set manually
export PBS_NCPUS=48


# ----------------------------------------------------------
# LOGGING
# ----------------------------------------------------------
LOGFILE="hpc_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(stdbuf -oL tee -a "$LOGFILE") 2>&1

echo "================================================="
echo " JOB START — Simulation"
echo "================================================="
echo "  Date      : $(date)"
echo "  Host      : $(hostname)"
echo "  WorkDir   : $(pwd)"
echo "================================================="

# ----------------------------------------------------------
# SANITY CHECKS
# ----------------------------------------------------------

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
