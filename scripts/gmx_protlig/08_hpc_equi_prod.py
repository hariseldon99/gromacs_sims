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

    args = parser.parse_args()

    step8_hpc_equi_prod(
        n_cores=args.n_cores,
        gpu_id=args.gpu_id,
        prod_time_ns=args.prod_time,
        tc_grps_cmplx=args.tc_grps,
        ligand_resname=args.ligand_resname
    )

if __name__ == "__main__":
    main()
