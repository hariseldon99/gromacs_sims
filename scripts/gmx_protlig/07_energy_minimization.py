#!/usr/bin/env python3
"""
Step 7 — Energy minimisation (two-step).
CLI wrapper around gmx_protlig.step7_energy_minimization.
"""

import os
import sys
import argparse

try:
    from gmx_protlig import step7_energy_minimization
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gmx_protlig import step7_energy_minimization

def main():
    parser = argparse.ArgumentParser(description="Step 7: Energy minimisation (two-step).")
    parser.add_argument("--nt", type=int, default=1, help="Number of CPU threads to use.")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for energy minimisation.")

    args = parser.parse_args()

    step7_energy_minimization(
        nt=args.nt,
        gpu=args.gpu,
    )

if __name__ == "__main__":
    main()
