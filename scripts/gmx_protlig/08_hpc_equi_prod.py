#!/usr/bin/env python3
"""
Step 8 — Equilibration and 100 ns production MD (HPC).
CLI wrapper around gmx_protlig.step8_hpc_equi_prod.
"""

import os
import sys
import argparse

try:
    from gmx_protlig import step8_hpc_equi_prod
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gmx_protlig import step8_hpc_equi_prod

def main():
    parser = argparse.ArgumentParser(description="Step 8: Equilibration and Production MD.")
    parser.add_argument("--n-cores", "--n_cores", type=int, default=None, help="Number of CPU cores")
    parser.add_argument("--gpu-id", "--gpu_id", type=str, default="0", help="GPU ID (default: 0)")
    parser.add_argument("--prod-time", "--prod_time", type=float, default=100.0, help="Production time in ns (default: 100.0)")
    parser.add_argument("--tc_grps", "--tc-grps", type=str, default="Protein", help="Temperature coupling groups (default: Protein)")
    parser.add_argument("--ligand-resname", "--ligand_resname", type=str, default="UNL", help="Ligand residue name (default: UNL)")
    parser.add_argument("--pin", type=str, default="on", help="Thread pinning option for mdrun (default: on)")
    parser.add_argument("--gpu-pme", "--gpu_pme", type=str, default="gpu", help="PME calculation offload (default: gpu)")
    parser.add_argument("--gpu-bonded", "--gpu_bonded", type=str, default="gpu", help="Bonded interactions offload (default: gpu)")
    parser.add_argument("--gpu-update", "--gpu_update", type=str, default="gpu", help="Coordinate update offload (default: gpu)")
    parser.add_argument("--no-gpu-hook", action="store_true", help="Disable automatic GPU flags injection hook")

    args = parser.parse_args()

    step8_hpc_equi_prod(
        n_cores=args.n_cores,
        gpu_id=args.gpu_id,
        prod_time_ns=args.prod_time,
        tc_grps_cmplx=args.tc_grps,
        ligand_resname=args.ligand_resname,
        pin=args.pin,
        gpu_pme=args.gpu_pme,
        gpu_bonded=args.gpu_bonded,
        gpu_update=args.gpu_update,
        enable_gpu_hook=not args.no_gpu_hook,
    )

if __name__ == "__main__":
    main()

