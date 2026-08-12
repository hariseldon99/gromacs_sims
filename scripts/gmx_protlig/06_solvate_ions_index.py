#!/usr/bin/env python3
"""
Step 6 — Box creation, solvation, ion addition, and index group creation.

Usage:
    python 06_solvate_ions_index.py --residue-name ERG
    python 06_solvate_ions_index.py --residue-name DIG --ion-conc 0.15
"""

import os
import sys
import glob
import pickle
import argparse
import datetime
import subprocess

from gromacs_py import gmx

# ---------------------------------------------------------------------------
# Command Line Argument Parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Step 6: Solvate, add ions, and create robust GROMACS index groups."
)

parser.add_argument(
    "-r", "--residue-name",
    type=str,
    default="UNL",
    help="3-letter residue name of the ligand (default: UNL). Needed for indexx group creation."
)
parser.add_argument(
    "-c", "--ion-conc",
    type=float,
    default=0.15,
    help="Ion concentration in mol/L (default: 0.15)"
)
parser.add_argument(
    "-b", "--box-dist",
    type=float,
    default=1.5,
    help="Minimum distance from protein to box edge in nm (default: 1.5)"
)
args = parser.parse_args()

RESIDUE_NAME      = args.residue_name.upper()
ION_CONCENTRATION = args.ion_conc
BOX_DIST_NM       = args.box_dist
BOX_TYPE          = "dodecahedron"
SOLV_OUTDIR       = "complex_solv"
MAXWARN           = 5

# ---------------------------------------------------------------------------
# Helper: Pure Python Index Group Generator
# ---------------------------------------------------------------------------
def generate_custom_index_groups(gro_file, ndx_file, res_name):
    """
    Parses .gro file and appends required custom groups to index.ndx:
      - Protein_<RESNAME>
      - Dry
      - Water_and_ions
    """
    # 1. Generate default GROMACS index file first
    subprocess.run(
        ["gmx", "make_ndx", "-f", gro_file, "-o", ndx_file],
        input="q\n",
        text=True,
        capture_output=True,
    )

    # 2. Parse GRO file structure line-by-line
    water_ion_resnames = {"SOL", "NA", "CL", "NA+", "CL-", "WAT", "HOH", "ION"}
    
    protein_atoms = []
    ligand_atoms = []
    water_ion_atoms = []

    with open(gro_file, "r") as fh:
        lines = fh.readlines()

    # Atom lines start at line 3 and end before the box vector line (last line)
    atom_lines = lines[2:-1]

    for idx, line in enumerate(atom_lines, start=1):
        if len(line) < 10:
            continue
        res_name_in_gro = line[5:10].strip().upper()

        if res_name_in_gro == res_name:
            ligand_atoms.append(idx)
        elif res_name_in_gro in water_ion_resnames:
            water_ion_atoms.append(idx)
        else:
            protein_atoms.append(idx)

    protein_ligand_atoms = protein_atoms + ligand_atoms

    # Helper to format atom lists into standard GROMACS 15-column index blocks
    def format_group(group_name, atom_indices):
        out = [f"\n[ {group_name} ]\n"]
        for i in range(0, len(atom_indices), 15):
            chunk = atom_indices[i : i + 15]
            out.append(" ".join(f"{atom:6d}" for atom in chunk) + "\n")
        return "".join(out)

    # 3. Append custom groups
    with open(ndx_file, "a") as fh:
        fh.write(format_group(f"Protein_{res_name}", protein_ligand_atoms))
        fh.write(format_group("Dry", protein_ligand_atoms))
        fh.write(format_group("Water_and_ions", water_ion_atoms))

    print(f"[OK] Successfully appended groups to {ndx_file}:")
    print(f"     - Protein_{res_name} ({len(protein_ligand_atoms)} atoms)")
    print(f"     - Dry ({len(protein_ligand_atoms)} atoms)")
    print(f"     - Water_and_ions ({len(water_ion_atoms)} atoms)")

# ---------------------------------------------------------------------------
# Load latest topology-merge checkpoint
# ---------------------------------------------------------------------------
chkpts = sorted(glob.glob("checkpoint.merge_topologies_*.pycpt"))
if not chkpts:
    print("[ERROR] No merge-topology checkpoint found. Run Step 5 first.")
    sys.exit(1)
latest = chkpts[-1]
print(f"[INFO] Loading checkpoint: {latest}")
with open(latest, "rb") as fh:
    complex_sys = pickle.load(fh)

if complex_sys.name is not None:
    print(f"[INFO] Loaded system name: {complex_sys.name}")
    SYS_NAME = complex_sys.name
else:
    SYS_NAME = f"complex_{RESIDUE_NAME.lower()}"
    print(f"[INFO] No system name found in checkpoint. Using default: {SYS_NAME}")


print(f"[INFO] Target Residue Name : {RESIDUE_NAME}")
print(f"[INFO] coor_file            : {complex_sys.coor_file}")
print(f"[INFO] top_file             : {complex_sys.top_file}")

# Sanity check
for p in (complex_sys.coor_file, complex_sys.top_file):
    if not os.path.exists(p):
        print(f"[ERROR] File missing: {p}")
        sys.exit(1)

# ---------------------------------------------------------------------------
# 1. Create dodecahedral simulation box
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("  [Step 6a] Creating simulation box …")
print("=" * 60)
complex_sys.create_box(
    dist=BOX_DIST_NM,
    box_type=BOX_TYPE,
    check_file_out=True,
)
print(f"[OK] Box created. coor_file: {complex_sys.coor_file}")

# ---------------------------------------------------------------------------
# 2. Solvate and add ions
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("  [Step 6b] Solvating and adding ions …")
print("=" * 60)
complex_sys.solvate_add_ions(
    out_folder=SOLV_OUTDIR,
    name=SYS_NAME,
    ion_C=ION_CONCENTRATION,
    create_box_flag=False,
    maxwarn=MAXWARN,
)
print(f"[OK] Solvated system: {complex_sys.coor_file}")
print(f"     top_file       : {complex_sys.top_file}")

# ---------------------------------------------------------------------------
# 3. Create GROMACS index file using Python parser
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("  [Step 6c] Creating custom index groups …")
print("=" * 60)

solv_gro = complex_sys.coor_file
ndx_file = os.path.join(SOLV_OUTDIR, "index.ndx")

generate_custom_index_groups(solv_gro, ndx_file, RESIDUE_NAME)
complex_sys.ndx = ndx_file

# ---------------------------------------------------------------------------
# 4. Checkpoint
# ---------------------------------------------------------------------------
chkpt = f"checkpoint.solvate_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
with open(chkpt, "wb") as fh:
    pickle.dump(complex_sys, fh)
print(f"\n[Checkpoint] Saved: {chkpt}")

print("""
╔══════════════════════════════════════════════════════════╗
║  Step 6 complete.                                        ║
║  Next: python 07_energy_minimization.py                  ║
╚══════════════════════════════════════════════════════════╝
""")
