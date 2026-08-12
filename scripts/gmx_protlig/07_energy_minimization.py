#!/usr/bin/env python3
"""
Step 7 — Energy minimisation (two-step: without bond constraints, then with).

This is the last step run locally before transferring to HPC.

Run from the simulations/ directory inside Singularity:

    bash run_local.sh 07_energy_minimization.py

Inputs
------
checkpoint.solvate_YYYYMMDD.pycpt   — from Step 6

Outputs
-------
em/                                  — energy minimisation output
checkpoint.em_YYYYMMDD.pycpt         — pickle for HPC resumption
"""

import os
import sys
import glob
import pickle
import datetime
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # headless
import matplotlib.pyplot as plt
import argparse
from gromacs_py import gmx

# ---------------------------------------------------------------------------
# Parallelization (use all available local cores, no GPU needed for EM)
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Step 7: Energy minimisation (two-step: without bond constraints, then with)."
)
parser.add_argument(
    "--nt",
    type=int,
    default=1,
    help="Number of threads to use for energy minimisation."
)
parser.add_argument(
    "--gpu",
    action="store_true",
    help="Use GPU for energy minimisation (default: False)."
)

args = parser.parse_args()

N_CORES = args.nt 
GPU_FLAG = args.gpu
# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
EM_OUTDIR      = "em"
EM_STEPS       = 10000   # steps per EM stage
EMTOL          = 10.0    # kJ mol⁻¹ nm⁻¹ — convergence criterion
EMSTEP         = 0.01    # nm — initial step size
MAXWARN        = 5

# ---------------------------------------------------------------------------
# Load latest solvation checkpoint
# ---------------------------------------------------------------------------
chkpts = sorted(glob.glob("checkpoint.solvate_*.pycpt"))
if not chkpts:
    print("[ERROR] No solvation checkpoint found. Run Step 6 first.")
    sys.exit(1)
latest = chkpts[-1]
print(f"[INFO] Loading checkpoint: {latest}")
with open(latest, "rb") as fh:
    complex_sys = pickle.load(fh)

print(f"[INFO] coor_file : {complex_sys.coor_file}")
print(f"[INFO] top_file  : {complex_sys.top_file}")

# ---------------------------------------------------------------------------
# Configure gromacs_py parallelization
# ---------------------------------------------------------------------------
complex_sys.nt     = N_CORES
complex_sys.ntmpi  = 1
complex_sys.gpu_id = None if not GPU_FLAG else '0'

print(f"[INFO] Using {N_CORES} CPU threads for energy minimisation.")

# ---------------------------------------------------------------------------
# Run two-step energy minimisation
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("  [Step 7] Running two-step energy minimisation …")
print("=" * 60)

complex_sys.em_2_steps(
    out_folder=EM_OUTDIR,
    no_constr_nsteps=EM_STEPS,
    constr_nsteps=EM_STEPS,
    posres=None,               # no position restraints during EM
    create_box_flag=False,
    emtol=EMTOL,
    emstep=EMSTEP,
    maxwarn=MAXWARN,
)

print(f"\n[OK] Energy minimisation complete.")
print(f"     coor_file: {complex_sys.coor_file}")

# ---------------------------------------------------------------------------
# Plot potential energy profile and save to PNG
# ---------------------------------------------------------------------------
try:
    ener_1 = complex_sys.sys_history[-1].get_ener(selection_list=["Potential"])
    ener_2 = complex_sys.get_ener(selection_list=["Potential"])

    ener_1["label"] = "no bond constraints"
    ener_2["label"] = "bond constraints"
    ener_2["Time (ps)"] = (
        ener_2["Time (ps)"].values
        + ener_1["Time (ps)"].max()
    )

    ener_all = pd.concat([ener_1, ener_2])

    fig, ax = plt.subplots(figsize=(10, 4))
    for label, grp in ener_all.groupby("label"):
        ax.plot(grp["Time (ps)"], grp["Potential"], label=label)
    ax.set_xlabel("Step")
    ax.set_ylabel("Potential energy (kJ mol⁻¹)")
    ax.set_title("Energy Minimisation")
    ax.legend()
    ax.grid(alpha=0.4)
    plt.tight_layout()
    fig.savefig(os.path.join(EM_OUTDIR, "em_energy.png"), dpi=150)
    plt.close(fig)
    print(f"[OK] Energy plot saved: {EM_OUTDIR}/em_energy.png")
except Exception as exc:
    print(f"[WARNING] Could not plot energy: {exc}")

# ---------------------------------------------------------------------------
# Checkpoint — this is the file transferred to HPC
# ---------------------------------------------------------------------------
chkpt = f"checkpoint.em_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
with open(chkpt, "wb") as fh:
    pickle.dump(complex_sys, fh)
print(f"\n[Checkpoint] Saved: {chkpt}")

print(f"""
╔═══════════════════════════════════════════════════════════════╗
║  Step 7 complete.                                             ║
║                                                               ║
║  Energy-minimised structure: {EM_OUTDIR}/                     ║
║  Checkpoint for HPC         : {chkpt}                         ║
║                                                               ║
║  Transfer the entire simulations/ directory to the HPC        ║
║  cluster, then submit the SLURM job:                          ║
║    sbatch 09_hpc_submit.sbatch                                ║
╚═══════════════════════════════════════════════════════════════╝
""")
