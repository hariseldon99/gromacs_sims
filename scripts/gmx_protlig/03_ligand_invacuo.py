#!/usr/bin/env python3
"""
Step 3 — In-vacuo ligand simulation (energy minimisation + short MD).
CLI wrapper around gmx_protlig.step3_ligand_invacuo.
"""

import os
import sys
import argparse

try:
    from gmx_protlig import step3_ligand_invacuo
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gmx_protlig import step3_ligand_invacuo

def main():
    parser = argparse.ArgumentParser(description="In-vacuo ligand simulation test.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory (default: ligand)")
    parser.add_argument("--invacuo-dir", "--invacuo_dir", type=str, default="invacuo", help="Output directory (default: invacuo)")
    parser.add_argument("--ntomp", type=int, default=None, help="Number of OpenMP threads")

    args = parser.parse_args()

    step3_ligand_invacuo(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        invacuo_dir=args.invacuo_dir,
        ntomp=args.ntomp,
    )

if __name__ == "__main__":
    main()
