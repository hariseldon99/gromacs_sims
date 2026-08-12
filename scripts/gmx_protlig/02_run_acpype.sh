#!/bin/bash
# =============================================================================
# Step 2 — Run acpype on the protonated ligand PDB to generate AMBER
#           GROMACS topology (ITP + GRO + position restraints).
#
# Run from the simulations/ directory:
#   bash run_local.sh 02_run_acpype.sh --residue-name ERG
#   bash run_local.sh 02_run_acpype.sh --residue-name DIG
#
# Arguments:
#   --residue-name RESNAME   3-letter residue name (default: ERG)
#   --ligand-dir DIR         Ligand directory (default: ligand)
#
# Outputs (inside ligand/<RESNAME>.acpype/):
#   <RESNAME>_GMX.gro       — ligand coordinates
#   <RESNAME>_GMX.top       — complete AMBER topology (includes atomtypes)
#   <RESNAME>_GMX.itp       — molecule-type section only (for inclusion)
#   posre_<RESNAME>.itp     — position restraint file
#   em.mdp, nvt.mdp         — template MDP files (used for in-vacuo test)
# =============================================================================

set -euo pipefail

# Parse arguments
RESIDUE_NAME="ERG"
LIGAND_DIR="ligand"

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
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

RESIDUE_NAME=$(echo "$RESIDUE_NAME" | tr '[:lower:]' '[:upper:]')

LIGAND_PDB="${LIGAND_DIR}/${RESIDUE_NAME}.pdb"

if [ ! -f "$LIGAND_PDB" ]; then
    echo "[ERROR] $LIGAND_PDB not found. Please run Step 1 first."
    exit 1
fi

echo "================================================="
echo "  Running acpype on $LIGAND_PDB"
echo "  Residue: $RESIDUE_NAME"
echo "================================================="

# Change to the ligand directory so acpype output stays there.
# We are already inside the Singularity container when run via run_local.sh,
# so 'acpype' is available directly.
cd "$LIGAND_DIR"

acpype \
    -i "${RESIDUE_NAME}.pdb" \
    -n 0 \
    -d 
    #-c bcc \
    #-a amber \
    #-o gmx \

cd ..

# Verify key output files exist
ACPYPE_DIR="${LIGAND_DIR}/${RESIDUE_NAME}.acpype"
for f in "${RESIDUE_NAME}_GMX.gro" "${RESIDUE_NAME}_GMX.top" "${RESIDUE_NAME}_GMX.itp"; do
    if [ ! -f "${ACPYPE_DIR}/${f}" ]; then
        echo "[ERROR] Expected output not found: ${ACPYPE_DIR}/${f}"
        echo "  Check acpype output above for errors."
        exit 1
    fi
done

echo ""
echo "================================================="
echo "  acpype completed successfully."
echo "  Output files in: $ACPYPE_DIR"
echo "================================================="
ls -lh "$ACPYPE_DIR"
echo ""
echo "Next: run Step 3 (in-vacuo ligand test):"
echo "  bash run_local.sh 03_ligand_invacuo.sh --residue-name $RESIDUE_NAME"
echo ""
echo "Then run Step 4 (prepare protein topology):"
echo "  bash run_local.sh 04_prepare_protein.py --protein-pdb <name-of-protein-pdb> --residue-name $RESIDUE_NAME"
