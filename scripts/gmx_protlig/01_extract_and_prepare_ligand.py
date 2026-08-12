#!/usr/bin/env python3
"""
Step 1 — Extract the best-docked pose from a PDBQT file, convert it to PDB
format, and protonate it at a specified pH with OpenBabel.

Run from the simulations/ directory inside Singularity:

    bash run_local.sh 01_extract_and_prepare_ligand.py \\
        --ligand-name ERGOTAMINE \\
        --residue-name ERG \\
        --input-pdbqt ../CHEMBL442_ERGOTAMINE_out.pdbqt \\
        [--ph 7.4] [--ligand-dir ligand]

Or with defaults (ERGOTAMINE):

    bash run_local.sh 01_extract_and_prepare_ligand.py

Outputs (in ligand_dir/)
-------
<ligand_name>_pose1.pdbqt   — best-docked pose (MODEL 1)
<ligand_name>_raw.pdb       — converted PDB without hydrogens
<RESNAME>.pdb               — protonated PDB at specified pH
"""

import os
import argparse
import subprocess
import sys

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Extract and prepare a ligand from docked PDBQT for gromacs_py."
)
parser.add_argument(
    "--ligand-name",
    type=str,
    default="UnknownLigand",
    help="Name of the ligand (default: UnknownLigand)",
)
parser.add_argument(
    "--residue-name",
    type=str,
    default="UNL",
    help="3-letter GROMACS residue name for the ligand (default: UNL)",
)
parser.add_argument(
    "--input-pdbqt",
    type=str,
    help="Path to the input PDBQT file containing docked poses",
)

parser.add_argument(
    "--ph",
    type=float,
    default=None,
    help="pH for protonation, if None, no pH adjustment is made",
)

parser.add_argument(
    "--full_protonate",
    action="store_false", 
    help="Fully protonate the ligand (default: False - protonates only polar hydrogens)",
)

args = parser.parse_args()

# Validate inputs
if not os.path.exists(args.input_pdbqt):
    print(f"[ERROR] Input PDBQT not found: {args.input_pdbqt}")
    sys.exit(1)

if len(args.residue_name) != 3:
    print(f"[ERROR] Residue name must be exactly 3 letters: {args.residue_name}")
    sys.exit(1)

LIGAND_NAME = args.ligand_name
RESIDUE_NAME = args.residue_name.upper()
PDBQT_FILE = args.input_pdbqt
PH = args.ph
LIGAND_DIR = args.ligand_dir

os.makedirs(LIGAND_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Read all models and confirm MODEL 1 is the best (most-negative score)
# ---------------------------------------------------------------------------
print(f"[Step 1] Reading docked poses from {LIGAND_NAME} PDBQT …")
print(f"         File: {PDBQT_FILE}")

with open(PDBQT_FILE) as fh:
    raw = fh.read()

blocks = {}
current_model = None
current_lines = []

for line in raw.splitlines(keepends=True):
    if line.startswith("MODEL "):
        current_model = int(line.split()[1])
        current_lines = [line]
    elif line.startswith("ENDMDL") and current_model is not None:
        current_lines.append(line)
        blocks[current_model] = current_lines
        current_model = None
        current_lines = []
    elif current_model is not None:
        current_lines.append(line)

scores = {}
for mid, lines in blocks.items():
    for l in lines:
        if "VINA RESULT:" in l:
            scores[mid] = float(l.split()[3])
            break

print("  Docked poses found:")
for mid in sorted(scores):
    flag = " ← BEST" if mid == min(scores, key=scores.get) else ""
    print(f"    MODEL {mid}: {scores[mid]:.3f} kcal/mol{flag}")

best_model = min(scores, key=scores.get)
if best_model != 1:
    print(f"\n[WARNING] Best pose is MODEL {best_model}, not MODEL 1!")
    print(f"          Proceeding with MODEL {best_model}.")
else:
    print(f"\n[OK] Confirmed: MODEL 1 is the best pose ({scores[1]:.3f} kcal/mol).")

# ---------------------------------------------------------------------------
# 2. Write the best pose to a new PDBQT
# ---------------------------------------------------------------------------
best_pdbqt = os.path.join(LIGAND_DIR, f"{LIGAND_NAME.lower()}_pose1.pdbqt")
with open(best_pdbqt, "w") as fh:
    fh.writelines(blocks[best_model])
print(f"\n[Step 2] Best-pose PDBQT written to: {best_pdbqt}")

# ---------------------------------------------------------------------------
# 3. Convert PDBQT → PDB with OpenBabel
# ---------------------------------------------------------------------------
raw_pdb = os.path.join(LIGAND_DIR, f"{LIGAND_NAME.lower()}_raw.pdb")
print(f"\n[Step 3] Converting PDBQT → PDB …")
ret = subprocess.run(
    [
        "obabel",
        "-ipdbqt", best_pdbqt,
        "-opdb", "-O", raw_pdb,
        "--partialcharge", "none",
    ],
    capture_output=True, text=True
)
if ret.returncode != 0:
    print(f"[ERROR] obabel failed:\n{ret.stderr}")
    sys.exit(1)
print(f"  Written: {raw_pdb}")

# ---------------------------------------------------------------------------
# 4. Protonate at specified pH and save with residue name
#    Comment: Disabled pH setting due to wrong charge assignments
#    The filename uses RESIDUE_NAME so acpype picks it up correctly.
# ---------------------------------------------------------------------------
prot_pdb = os.path.join(LIGAND_DIR, f"{RESIDUE_NAME}.pdb")
cmd =  [
        "obabel",
        "-ipdb", raw_pdb,
        "-opdb", "-O", prot_pdb,
        ]

if args.full_protonate:
    cmd.append("-h")  # Add hydrogens to all atoms
else:
    cmd.append("--AddPolarH")  # Protonate only polar hydrogens    

if PH is not None:
    print(f"\n[Step 4] Protonating at pH {PH} with OpenBabel …")
    cmd.append("-p")
    cmd.append(str(PH))

ret = subprocess.run(cmd, capture_output=True, text=True)
if ret.returncode != 0:
    print(f"[ERROR] obabel protonation failed:\n{ret.stderr}")
    sys.exit(1)

# Rename residue to RESIDUE_NAME in the PDB file (obabel may leave it as UNL)
with open(prot_pdb) as fh:
    content = fh.read()
content = content.replace(" UNL ", f" {RESIDUE_NAME} ")
with open(prot_pdb, "w") as fh:
    fh.write(content)

print(f"  Written: {prot_pdb}")

print(f"""
╔══════════════════════════════════════════════════════════╗
║  Step 1 complete — {LIGAND_NAME} prepared.               ║
║                                                          ║
║  Output: {LIGAND_DIR}/{RESIDUE_NAME}.pdb                ║
║                                                          ║
║  Next: run Step 2 to generate AMBER topology with acpype ║
║    bash run_local.sh 02_run_acpype.sh --residue-name {RESIDUE_NAME} ║
╚══════════════════════════════════════════════════════════╝
""")
