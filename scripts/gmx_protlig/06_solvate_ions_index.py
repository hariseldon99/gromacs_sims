#!/usr/bin/env python3
"""
Step 6 — Solvation, ion addition, and index group creation.
CLI wrapper around protlig_api.step6_solvate_ions_index.
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protlig_api import step6_solvate_ions_index

def main():
    parser = argparse.ArgumentParser(description="Step 6: Solvate, add ions, and create custom index groups.")
    parser.add_argument("-r", "--residue-name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("-c", "--ion-conc", type=float, default=0.15, help="Ion concentration in mol/L (default: 0.15)")
    parser.add_argument("-b", "--box-dist", type=float, default=1.5, help="Distance from protein to box edge in nm (default: 1.5)")

    args = parser.parse_args()

    step6_solvate_ions_index(
        residue_name=args.residue_name,
        ion_conc=args.ion_conc,
        box_dist=args.box_dist,
    )

if __name__ == "__main__":
    main()
