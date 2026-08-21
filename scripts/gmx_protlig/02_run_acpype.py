#!/usr/bin/env python3
"""
Step 2 — Run acpype on protonated ligand PDB to generate AMBER GROMACS topology.
CLI wrapper around protlig_api.step2_run_acpype.
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protlig_api import step2_run_acpype

def main():
    parser = argparse.ArgumentParser(description="Run ACPYPE for ligand topology generation.")
    parser.add_argument("--residue-name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("--ligand-dir", type=str, default="ligand", help="Ligand directory (default: ligand)")
    parser.add_argument("--net-charge", type=int, default=0, help="Net charge of ligand (default: 0)")
    parser.add_argument("--charge-method", type=str, default="bcc", help="Charge method (default: bcc)")

    args = parser.parse_args()

    step2_run_acpype(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        net_charge=args.net_charge,
        charge_method=args.charge_method,
    )

if __name__ == "__main__":
    main()
