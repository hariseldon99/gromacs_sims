#!/usr/bin/env python3
"""
Example script demonstrating Pythonic batch processing of protein-ligand complexes
using the gmx_protlig API.

Run inside Singularity:
    bash run_local.sh example_batch_run.py
"""
import os
from pathlib import Path
PROD_OUTDIR = "prod"
LIG = "UNL"
PROTLIG_SELECT = f"Protein_{LIG}"

HOME = os.getenv("HOME", ".")
prod_checkpoint_path = Path("./")
chkpt_file = Path(prod_checkpoint_path) / "checkpoint.prod_20260821.pycpt"

import sys, subprocess
import logging
from protlig_api import ProteinLigandPipeline, BatchPipelineRunner, step1_extract_and_prepare_ligand
import pickle, re

# Configure logging to show progress in terminal
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
    datefmt="%H:%M:%S"
)


os.chdir(prod_checkpoint_path)
print(f"Using production checkpoint path: {chkpt_file}")

with open(chkpt_file, 'rb') as py_cpt:
    complex_sys = pickle.load(py_cpt)


# Match any leading '../' sequence followed by 'root/<anything>/'
ROOT_PATTERN = re.compile(r"^(?:\.\./)+root/[^/]+/")

base_dir = Path(prod_checkpoint_path)

# Use getattr / setattr to support properties, slots, or normal instance dicts
attrs_to_check = getattr(complex_sys, "__slots__", None) or vars(complex_sys).keys()

for attr in list(attrs_to_check):
    val = getattr(complex_sys, attr, None)
    if val is None:
        continue
    
    val_str = str(val)
    if "/root/" in val_str:
        # Splits on '/root/<subfolder>/' and keeps everything after
        rel_path = val_str.split("/root/", 1)[1].split("/", 1)[1]
        new_path = base_dir / rel_path
        
        setattr(complex_sys, attr, new_path if isinstance(val, Path) else str(new_path))

# Verify updated attributes
print("\nUpdated state:")
complex_sys.display()
# Center trajectory
complex_sys.center_mol_box(traj=True)
complex_sys.convert_trj(select=f'Protein\nProtein | {LIG}', fit='rot+trans', pbc='none')

# Extract dry trajectory (protein + ligand only, no water/ions)
# Find the production TPR
prod_tpr = complex_sys.tpr
fitted_xtc = complex_sys.xtc
ndx_file = complex_sys.ndx

dry_xtc = os.path.join(PROD_OUTDIR, "traj_dry.xtc")
dry_tpr = os.path.join(PROD_OUTDIR, "traj_dry.tpr")

if prod_tpr and fitted_xtc and os.path.exists(ndx_file):
    # trjconv to extract dry
    ret = subprocess.run(
        [
            "gmx", "trjconv",
            "-s", prod_tpr,
            "-f", fitted_xtc,
            "-n", ndx_file,
            "-o", dry_xtc,
        ],
        input=f"Protein | {LIG}\n",
        text=True,
        capture_output=True,
    )
    if ret.returncode != 0:
        print(f"[WARNING] trjconv for dry trajectory failed:\n{ret.stderr}")
    else:
        print(f"[OK] Dry trajectory: {dry_xtc}")

    # Also extract a dry TPR (for VMD loading)
    ret2 = subprocess.run(
        [
            "gmx", "convert-tpr",
            "-s", prod_tpr,
            "-n", ndx_file,
            "-o", dry_tpr,
        ],
        input=f"Protein | {LIG}\n",
        text=True,
        capture_output=True,
    )
    if ret2.returncode == 0:
        print(f"[OK] Dry TPR: {dry_tpr}")
else:
    print("[WARNING] Could not auto-extract dry trajectory.")
    print(f"  prod_tpr   : {prod_tpr}")
    print(f"  fitted_xtc : {fitted_xtc}")
    print(f"  index_file : {ndx_file}")
    print("  Please run trjconv manually to extract the Dry group.")

print(f"""
╔═══════════════════════════════════════════════════════════════╗
║  All simulation stages complete.                              ║
║                                                               ║
║  Key outputs:                                                 ║
║    Production    : {PROD_OUTDIR}/                             ║
║    Dry trajectory: {dry_xtc}                 ║
║    Dry TPR       : {dry_tpr}                 ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
""")