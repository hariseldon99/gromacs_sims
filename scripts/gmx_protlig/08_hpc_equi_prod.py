#!/usr/bin/env python3
"""
Step 8 — Equilibration and 100 ns production MD (HPC).
CLI wrapper around protlig_api.step8_hpc_equi_prod.
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protlig_api import step8_hpc_equi_prod

def main():
    parser = argparse.ArgumentParser(description="Step 8: Equilibration and Production MD.")
    parser.add_argument("--n-cores", type=int, default=None, help="Number of CPU cores")
    parser.add_argument("--gpu-id", type=str, default="0", help="GPU ID (default: 0)")
    parser.add_argument("--prod-time", type=float, default=100.0, help="Production time in ns (default: 100.0)")
    parser.add_argument("--tc_grps", type=str, default="Protein", help="Temperature coupling groups (default: Protein)")

    args = parser.parse_args()

    step8_hpc_equi_prod(
        n_cores=args.n_cores,
        gpu_id=args.gpu_id,
        prod_time_ns=args.prod_time,
        tc_grps_cmplx=args.tc_grps
    )

if __name__ == "__main__":
    main()
