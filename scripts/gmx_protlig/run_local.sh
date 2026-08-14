#!/bin/bash
# =============================================================================
# run_local.sh — wrapper to execute any Python (or shell) script inside the
#                Singularity container on the local workstation.
#
# Usage (from the simulations/ directory):
#   bash run_local.sh <script.py> [args...]
#   bash run_local.sh <script.sh> [args...]
#
# The script mounts the current directory as /host_pwd inside the container
# so all relative paths used by gromacs_py remain consistent.
# =============================================================================

set -euo pipefail

SIF_LOCAL="${HOME}/images/gromacs_2024.5.1-GPU_20260718.sif"

if [ ! -f "$SIF_LOCAL" ]; then
    echo "[ERROR] Singularity image not found at: $SIF_LOCAL"
    echo "        Please update SIF_LOCAL in this script."
    exit 1
fi

if [ $# -lt 1 ]; then
    echo "Usage: bash run_local.sh <script.py|script.sh> [args...]"
    exit 1
fi

SCRIPT="$1"
shift

# Determine how to run the script
case "$SCRIPT" in
    *.py)
        RUNNER="python3"
        ;;
    *.sh)
        RUNNER="bash"
        ;;
    *)
        RUNNER=""
        ;;
esac

echo "================================================="
echo "  Singularity local runner"
echo "  Image : $SIF_LOCAL"
echo "  Script: $SCRIPT"
echo "  CWD   : $(pwd)"
echo "================================================="

singularity exec \
    --nv \
    -B "${PWD}:/host_pwd" \
    --pwd /host_pwd \
    "$SIF_LOCAL" \
    $RUNNER "$SCRIPT" "$@"

echo "================================================="
echo "  Done: $SCRIPT"
echo "================================================="
