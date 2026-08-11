#!/bin/bash
# =============================================================================
# Step 3 — In-vacuo ligand simulation (optional but recommended).
#
# Purpose: Verify the acpype AMBER topology by running a short vacuum MD
# on the ligand alone and visualising the trajectory in VMD.
# A failed or unphysical vacuum MD indicates a problem with the topology
# (e.g. wrong charges, missing parameters) that must be fixed BEFORE
# running the full protein-ligand simulation.
#
# This script:
#   1. Runs a quick energy minimisation of the ligand in vacuo.
#   2. Runs a short (100 ps) vacuum MD using GPU (non-bonded) locally.
#
# Run from the simulations/ directory (inside Singularity):
#   bash run_local.sh 03_ligand_invacuo.sh --residue-name ERG
#   bash run_local.sh 03_ligand_invacuo.sh --residue-name DIG
#
# Arguments:
#   --residue-name RESNAME   3-letter residue name (default: ERG)
#   --ligand-dir DIR         Ligand directory (default: ligand)
#   --invacuo-dir DIR        Output directory (default: invacuo)
#
# After completion, open invacuo/md_vac.trr in VMD to inspect.
# =============================================================================

set -euo pipefail

# Parse arguments
RESIDUE_NAME="ERG"
LIGAND_DIR="ligand"
INVACUO_DIR="invacuo"

while [[ $# -gt 0 ]]; do
    case $1 in
        --residue-name)
            RESIDUE_NAME="$2"
            shift 2
            ;;
        --ligand-dir)
            LIGAND_DIR="$2"
            shift 2
            ;;
        --invacuo-dir)
            INVACUO_DIR="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

RESIDUE_NAME=$(echo "$RESIDUE_NAME" | tr '[:lower:]' '[:upper:]')
NTOMP=$(nproc)

ACPYPE_DIR="${LIGAND_DIR}/${RESIDUE_NAME}.acpype"

# ─── Sanity checks ──────────────────────────────────────────────────────────
for f in "${RESIDUE_NAME}_GMX.gro" "${RESIDUE_NAME}_GMX.top"; do
    if [ ! -f "${ACPYPE_DIR}/${f}" ]; then
        echo "[ERROR] ${ACPYPE_DIR}/${f} not found. Run Step 2 first."
        exit 1
    fi
done

mkdir -p "$INVACUO_DIR"
cp "${ACPYPE_DIR}/${RESIDUE_NAME}_GMX.gro" "${INVACUO_DIR}/${RESIDUE_NAME}.gro"
cp "${ACPYPE_DIR}/${RESIDUE_NAME}_GMX.top" "${INVACUO_DIR}/${RESIDUE_NAME}.top"
# Copy any .itp files referenced by the topology
cp "${ACPYPE_DIR}"/*.itp "${INVACUO_DIR}/" 2>/dev/null || true

# ─── Write EM MDP ────────────────────────────────────────────────────────────
cat > "${INVACUO_DIR}/em_vac.mdp" << 'EOF'
; In-vacuo energy minimisation
integrator              = steep
nsteps                  = 10000
emtol                   = 100.0
emstep                  = 0.01
nstlist                 = 1
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = cutoff
rcoulomb                = 1.2
rvdw                    = 1.2
pbc                     = no
EOF

# ─── Write vacuum MD MDP ─────────────────────────────────────────────────────
cat > "${INVACUO_DIR}/md_vac.mdp" << 'EOF'
; In-vacuo MD — 100 ps to check topology
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 500
nstvout                 = 500
nstfout                 = 500
nstlog                  = 500
nstenergy               = 500
nstlist                 = 10
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = cutoff
rcoulomb                = 1.2
rvdw                    = 1.2
pbc                     = no
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
constraints             = h-bonds
constraint_algorithm    = lincs
EOF

cd "${INVACUO_DIR}"

echo "================================================="
echo "  In-vacuo EM"
echo "================================================="

gmx grompp \
    -f em_vac.mdp \
    -c "${RESIDUE_NAME}.gro" \
    -p "${RESIDUE_NAME}.top" \
    -o em_vac.tpr \
    -maxwarn 5

gmx mdrun \
    -nobackup -nocopyright -v \
    -deffnm em_vac \
    -ntmpi 1 \
    -ntomp "${NTOMP}" \
    -nb gpu \
    -gpu_id 0 \
    -pme cpu

echo "================================================="
echo "  In-vacuo MD (100 ps)"
echo "================================================="

gmx grompp \
    -f md_vac.mdp \
    -c em_vac.gro \
    -p "${RESIDUE_NAME}.top" \
    -o md_vac.tpr \
    -maxwarn 5

gmx mdrun \
    -nobackup -nocopyright -v \
    -deffnm md_vac \
    -ntmpi 1 \
    -ntomp "${NTOMP}" \
    -nb gpu \
    -gpu_id 0 \
    -pme cpu

cd ..

echo ""
echo "================================================="
echo "  In-vacuo simulation complete."
echo "  Open invacuo/md_vac.trr in VMD to inspect."
echo "================================================="
echo ""
echo "Next step: bash run_local.sh 04_prepare_protein.py"
