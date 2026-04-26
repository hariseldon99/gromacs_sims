#!/bin/bash
#SBATCH --job-name=g16_smopt
#SBATCH --output=logs/g16_sm_%j.out
#SBATCH --error=logs/g16_sm_%j.err
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --partition=CPU
#SBATCH --qos=elevated
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=telegram:-1001215703472

mkdir -p logs
echo "Starting Gaussian jobs"
# 2. Gaussian Environment (Source without strict mode)
source /etc/g16setup

# 3. Scratch Directory Setup
export GAUSS_SCRDIR=$(pwd)/scratch_g16
mkdir -p "$GAUSS_SCRDIR"

# 4. Correct Filenames (Ensuring they match your 'ls')
SM_OPT="PDL_small_opt.gjf"
CPUS=$(taskset -cp $$ | awk -F':' '{print $2}')

echo "[$(date)] Starting Gaussian Job"
echo "Using File: $SM_OPT"

# 5. Patch Files in-place
sed -i -e "/^.*%cpu/I s/=.*$/=$CPUS/" "$SM_OPT"
sed -i "s|%Mem=.*|%Mem=32GB|" "$SM_OPT"
echo "Patched $SM_OPT"


echo "[$(date)] Small model opt restart — MaxStep=5, MaxCycles=500"
echo "  Checkpoint: PDL_small_opt.chk"
echo "  GAUSS_SCRDIR: $GAUSS_SCRDIR"

# ── Sanity checks ────────────────────────────────────────────────────────────
if [[ ! -f PDL_small_opt.chk ]]; then
    echo "[WARNING] PDL_small_opt.chk not found."
    echo "        This job must continue from the existing checkpoint."
fi
if [[ ! -f PDL_small_opt.gjf ]]; then
    echo "[WARNING] PDL_small_opt.gjf not found."
    exit 1
fi

# ── Run optimisation ─────────────────────────────────────────────────────────
echo "[$(date)] Running g16 < PDL_small_opt.gjf ..."
g16 < PDL_small_opt.gjf > logs/PDL_small_opt2.log 2>&1

# ── Check result ─────────────────────────────────────────────────────────────
echo ""
if grep -q "Normal termination" logs/PDL_small_opt2.log; then
    echo "[$(date)] [OK] Optimisation converged — Normal termination."
    FINAL_STEP=$(grep "Step number" logs/PDL_small_opt2.log | tail -1)
    FINAL_E=$(grep "SCF Done" logs/PDL_small_opt2.log | tail -1 | awk '{print $5}')
    echo "  $FINAL_STEP"
    echo "  Final SCF energy: $FINAL_E Ha"
else
    LAST_STEP=$(grep "Step number" logs/PDL_small_opt2.log | tail -1)
    LAST_CONV=$(grep -A5 "Item.*Value.*Threshold" logs/PDL_small_opt2.log | tail -8)
    echo "[$(date)] [WARN] Did not reach Normal termination."
    echo "  Last: $LAST_STEP"
    echo "  Last convergence table:"
    echo "$LAST_CONV"
    echo ""
    echo "  If oscillating: reduce MaxStep further (try MaxStep=3 in .gjf)"
    echo "  If still descending: increase walltime and resubmit from checkpoint"
fi

rm -rf "$GAUSS_SCRDIR"
echo "[$(date)] Done."
