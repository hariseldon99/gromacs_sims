#!/usr/bin/env python3
"""
Step 4 — Prepare the protein topology with pdb2gmx (via gromacs_py),
then combine protein and ligand coordinates into a single GRO file.

Run from the simulations/ directory inside Singularity:

    bash run_local.sh 04_prepare_protein.py --residue-name ERG
    bash run_local.sh 04_prepare_protein.py --residue-name DIG

Inputs
------
../CUEDC2_corrected.pdb                           — protein structure
ligand/<RESIDUE_NAME>.pdb                         — protonated ligand (from Step 1)
ligand/<RESIDUE_NAME>.acpype/<RESIDUE_NAME>_GMX.gro — ligand coordinates from acpype

Outputs
-------
protein/                                          — gromacs_py protein topology folder
complex_raw/
  complex_raw.gro                                 — protein + ligand combined GRO
checkpoint.prepare_protein_YYYYMMDD.pycpt         — pickle checkpoint
"""

import os, subprocess
import sys
import argparse
import pickle
import datetime

from gromacs_py import gmx

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Prepare protein topology and combine with ligand coordinates."
)
parser.add_argument(
    "--residue-name",
    type=str,
    default="ERG",
    help="3-letter GROMACS residue name for the ligand (default: ERG)",
)
parser.add_argument(
    "--ligand-dir",
    type=str,
    default="ligand",
    help="Ligand directory (default: ligand)",
)
parser.add_argument(
    "--protein-pdb",
    type=str,
    default="../CUEDC2_corrected.pdb",
    help="Path to protein PDB file (default: ../CUEDC2_corrected.pdb)",
)

args = parser.parse_args()

RESIDUE_NAME = args.residue_name.upper()
LIGAND_DIR = args.ligand_dir
PROTEIN_PDB = args.protein_pdb
LIGAND_GRO = os.path.join(LIGAND_DIR, f"{RESIDUE_NAME}.acpype", f"{RESIDUE_NAME}_GMX.gro")
PROTEIN_OUTDIR = "protein"
COMPLEX_OUTDIR = "complex_raw"

for p in (PROTEIN_PDB, LIGAND_GRO):
    if not os.path.exists(p):
        print(f"[ERROR] Required file not found: {p}")
        sys.exit(1)

os.makedirs(COMPLEX_OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Create protein topology with pdb2gmx
# ---------------------------------------------------------------------------
print("=" * 60)
print("  [Step 4a] Running pdb2gmx on CUEDC2 protein …")
print("=" * 60)

protein_sys = gmx.GmxSys(name="CUEDC2", coor_file=PROTEIN_PDB)

# amber99sb-ildn force field, TIP3P water
# ignh=True: re-add hydrogen atoms via pdb2gmx (recommended when the PDB
#             already has H from pdbfixer, to avoid clashes with FF H-names)
# ter selection: AMBER default = charged termini (NH3+ N-term, COO- C-term)
#protein_sys.prepare_top(
#    out_folder=PROTEIN_OUTDIR,
#    ff="amber99sb-ildn",
#    #ignh=True,
#    #vsite='hydrogens',
#    #maxwarn=5,
#)

protein_sys.add_top(
    out_folder=PROTEIN_OUTDIR,
    ff="amber99sb-ildn",
    pdb2gmx_option_dict={'ignh':None}
        )

# --- ADD THIS CONVERSION STEP ---
prot_gro_path = os.path.join(PROTEIN_OUTDIR, "CUEDC2_pdb2gmx.gro")
print(f"\n[INFO] Converting {protein_sys.coor_file} to {prot_gro_path}...")
subprocess.run(["gmx", "editconf", "-f", protein_sys.coor_file, "-o", prot_gro_path], check=True)
# --------------------------------

print(f"\n[OK] Protein topology created in: {PROTEIN_OUTDIR}")
print(f"     coor_file : {prot_gro_path}")
print(f"     top_file  : {protein_sys.top_file}")

# ---------------------------------------------------------------------------
# 2. Combine protein GRO + ligand GRO into a single GRO file
#    GRO format:
#      line 1:  title
#      line 2:  total number of atoms
#      lines 3…N+2: residue_nr residue_name atom_name  atom_nr  x  y  z
#      last line:  box vectors
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("  [Step 4b] Combining protein + ligand GRO files …")
print("=" * 60)


def read_gro(path):
    """Parse a GRO file into (title, atom_lines, box_line)."""
    with open(path) as fh:
        lines = fh.readlines()
    title    = lines[0]
    n_atoms  = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms]
    return title, atom_lines, box_line


def write_combined_gro(prot_gro, lig_gro, out_gro):
    """Write a combined GRO of protein + ligand."""
    _, prot_atoms, prot_box = read_gro(prot_gro)
    _, lig_atoms,  _        = read_gro(lig_gro)

    total = len(prot_atoms) + len(lig_atoms)

    # GRO atom lines: columns 0-4 = resnum (5 chars), 5-9 = resname (5),
    # 10-14 = atomname (5), 15-19 = atomnum (5), then x y z
    # We just renumber atom indices sequentially.
    combined_atoms = []
    atom_idx = 1
    for line in prot_atoms + lig_atoms:
        # Keep everything except atom index (cols 15-19)
        new_line = line[:15] + f"{atom_idx % 100000:5d}" + line[20:]
        combined_atoms.append(new_line)
        atom_idx += 1

    with open(out_gro, "w") as fh:
        fh.write("CUEDC2-ERG complex\n")
        fh.write(f"{total:5d}\n")
        fh.writelines(combined_atoms)
        fh.write(prot_box)

    return total


combined_gro = os.path.join(COMPLEX_OUTDIR, "complex_raw.gro")
n_total = write_combined_gro("protein/CUEDC2_pdb2gmx.gro", LIGAND_GRO, combined_gro)

print(f"[OK] Combined GRO written: {combined_gro}  ({n_total} atoms)")

# ---------------------------------------------------------------------------
# 3. Checkpoint
# ---------------------------------------------------------------------------
chkpt = f"checkpoint.prepare_protein_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
with open(chkpt, "wb") as fh:
    pickle.dump(protein_sys, fh)
print(f"\n[Checkpoint] Saved: {chkpt}")

print("""
╔══════════════════════════════════════════════════════════╗
║  Step 4 complete.                                        ║
║                                                          ║
║  MANUAL STEP REQUIRED before continuing:                 ║
║  Run Step 5 (topology merge):                            ║
║    bash run_local.sh 05_merge_topologies.py              ║
╚══════════════════════════════════════════════════════════╝
""")
