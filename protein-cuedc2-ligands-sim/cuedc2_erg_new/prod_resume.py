#!/usr/bin/env python3
"""
Step 8 — Equilibration and 100 ns production MD (HPC).

This script is run INSIDE Singularity on the HPC cluster, launched by
the SLURM job script 09_hpc_submit.sbatch.

Do NOT run this script directly — use sbatch 09_hpc_submit.sbatch instead.

Workflow
--------
1. Load the equi-minimisation checkpoint (from Step 8).
2. Run 100 ns production MD.
3. Post-process: make whole → center → nojump → extract dry trajectory.
4. Save checkpoints after equilibration and production.

Temperature coupling groups  : Protein_ERG  Water_and_ions
Reference temperature        : 310 K (physiological)
Pressure coupling            : Parrinello-Rahman (NPT)
"""

import os
import sys
import glob
import pickle
import datetime
import subprocess
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gromacs_py import gmx

# ---------------------------------------------------------------------------
# Parallelisation — on HPC, SLURM sets SLURM_CPUS_PER_TASK
# ---------------------------------------------------------------------------
N_CORES = int(os.environ.get("SLURM_CPUS_PER_TASK") or 
            os.environ.get("PBS_NCPUS") or 
            os.environ.get("OMP_NUM_THREADS") or
            multiprocessing.cpu_count())
print(f"[INFO] Using {N_CORES} OpenMP threads.")

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------
VSITE       = None
DT          = 0.002       # ps
DT_HA       = 0.001       # ps for heavy-atom constrained equi

# Equilibration times (ns)
HA_TIME     = 0.5         # ns — heavy atom restraints
CA_TIME     = 1.0         # ns — Cα restraints
CA_LOW_TIME = 4.0         # ns — low Cα restraints

HA_STEP     = int(1000 * HA_TIME     / DT_HA)
CA_STEP     = int(1000 * CA_TIME     / DT)
CA_LOW_STEP = int(1000 * CA_LOW_TIME / DT)

# Production (100 ns)
PROD_TIME   = 100.0       # ns
PROD_STEPS  = int(1000 * PROD_TIME / DT)

TC_GRPS_CMPLX = os.environ.get("TC_GRPS_CMPLX") or "Protein"

TC_GRPS     = f"{TC_GRPS_CMPLX} Water_and_ions"
TAU_T       = "0.1 0.1"
REF_T       = "310 310"

EQUI_OUTDIR = "equi"
PROD_OUTDIR = "prod"
MAXWARN     = 10
INDEX_FILE  = "complex_solv/index.ndx"
DRY_INDEX   = "complex_solv/index.ndx"   # same file, group name = Dry

# ---------------------------------------------------------------------------
# Load equilibration checkpoint
# ---------------------------------------------------------------------------
chkpts = sorted(glob.glob("checkpoint.equi_*.pycpt"))
if not chkpts:
    print("[ERROR] No equilibration checkpoint found.")
    print("        Make sure checkpoint.equi_YYYYMMDD.pycpt is in simulations/")
    sys.exit(1)
latest = chkpts[-1]
print(f"[INFO] Loading equilibration checkpoint: {latest}")
with open(latest, "rb") as fh:
    complex_sys = pickle.load(fh)

print(f"[INFO] coor_file : {complex_sys.coor_file}")
print(f"[INFO] top_file  : {complex_sys.top_file}")

# Attach index file
if os.path.exists(INDEX_FILE):
    complex_sys.ndx = INDEX_FILE
    print(f"[INFO] ndx file  : {INDEX_FILE}")

# Configure parallelisation
complex_sys.nt     = N_CORES
complex_sys.ntmpi  = 1
complex_sys.gpu_id = "0"

# ===========================================================================
# STAGE B — Production MD (100 ns)
# ===========================================================================
print("\n" + "=" * 60)
print(f"  [Stage B] Production MD — {PROD_TIME} ns ({PROD_STEPS} steps) …")
print("=" * 60)

complex_sys.production(
    out_folder=PROD_OUTDIR,
    nsteps=PROD_STEPS,
    tc_grps=TC_GRPS,
    tau_t=TAU_T,
    ref_t=REF_T,
    dt=DT,
    vsite=VSITE,
    maxwarn=1,
    nstlist=200,
)

# Checkpoint after production
prod_chkpt = f"checkpoint.prod_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
with open(prod_chkpt, "wb") as fh:
    pickle.dump(complex_sys, fh)
print(f"[Checkpoint] Saved: {prod_chkpt}")

# ===========================================================================
# STAGE C — Post-processing
# ===========================================================================
print("\n" + "=" * 60)
print("  [Stage C] Post-processing trajectory …")
print("=" * 60)

# 1. Center on Protein_ERG (wrapping molecules)
print("[PP1] Centering trajectory (protein + ligand) …")
complex_sys.center_mol_box(traj=True)

# 2. rot+trans fit, remove PBC jumps
print("[PP2] Fitting and removing PBC jumps …")
complex_sys.convert_trj(
    select="Protein_ERG\nSystem",
    fit="rot+trans",
    pbc="none",
)

# 3. Extract dry trajectory (protein + ligand only, no water/ions)
print("[PP3] Extracting dry trajectory …")

# Find the production TPR
prod_tpr = complex_sys.tpr if hasattr(complex_sys, "tpr") else None
if prod_tpr is None:
    # Try to find it
    tpr_candidates = sorted(glob.glob(os.path.join(PROD_OUTDIR, "*.tpr")))
    if tpr_candidates:
        prod_tpr = tpr_candidates[-1]

# The fitted XTC is the current coor_file's xtc counterpart
fitted_xtc = complex_sys.xtc if hasattr(complex_sys, "xtc") else None
if fitted_xtc is None:
    xtc_candidates = sorted(glob.glob(os.path.join(PROD_OUTDIR, "*_fit.xtc")))
    if not xtc_candidates:
        xtc_candidates = sorted(glob.glob(os.path.join(PROD_OUTDIR, "*.xtc")))
    if xtc_candidates:
        fitted_xtc = xtc_candidates[-1]

dry_xtc = os.path.join(PROD_OUTDIR, "traj_dry.xtc")
dry_tpr = os.path.join(PROD_OUTDIR, "traj_dry.tpr")

if prod_tpr and fitted_xtc and os.path.exists(INDEX_FILE):
    # trjconv to extract dry
    ret = subprocess.run(
        [
            "gmx", "trjconv",
            "-s", prod_tpr,
            "-f", fitted_xtc,
            "-n", INDEX_FILE,
            "-o", dry_xtc,
        ],
        input="Dry\n",
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
            "-n", INDEX_FILE,
            "-o", dry_tpr,
        ],
        input="Dry\n",
        text=True,
        capture_output=True,
    )
    if ret2.returncode == 0:
        print(f"[OK] Dry TPR: {dry_tpr}")
else:
    print("[WARNING] Could not auto-extract dry trajectory.")
    print(f"  prod_tpr   : {prod_tpr}")
    print(f"  fitted_xtc : {fitted_xtc}")
    print(f"  index_file : {INDEX_FILE}")
    print("  Please run trjconv manually to extract the Dry group.")

print(f"""
╔═══════════════════════════════════════════════════════════════╗
║  All simulation stages complete.                              ║
║                                                               ║
║  Key outputs:                                                 ║
║    Equilibration : {EQUI_OUTDIR}/                             ║
║    Production    : {PROD_OUTDIR}/                             ║
║    Dry trajectory: {dry_xtc}                 ║
║    Dry TPR       : {dry_tpr}                 ║
║                                                               ║
║  Checkpoints:                                                 ║
║    {equi_chkpt}                                               ║
║    {prod_chkpt}                                               ║
╚═══════════════════════════════════════════════════════════════╝
""")
