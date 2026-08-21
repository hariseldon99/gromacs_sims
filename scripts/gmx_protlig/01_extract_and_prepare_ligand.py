#!/usr/bin/env python3
"""
Step 1 — Extract the best-docked pose from a PDBQT file, convert it to PDB
format, and protonate it at a specified pH with OpenBabel.
CLI wrapper around protlig_api.step1_extract_and_prepare_ligand.
"""

import os
import sys
import argparse

# Allow importing protlig_api when run from script directory or simulations/
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protlig_api import step1_extract_and_prepare_ligand

def main():
    parser = argparse.ArgumentParser(
        description="Extract and prepare a ligand from docked PDBQT for GROMACS MD."
    )
    parser.add_argument(
        "--ligand-name",
        type=str,
        default="UNKNOWN",
        help="Name of the ligand (default: UNKNOWN)",
    )
    parser.add_argument(
        "--ligand-dir",
        type=str,
        default="ligand",
        help="Directory to store ligand files (default: ligand)",
    )
    parser.add_argument(
        "--residue-name",
        type=str,
        default="UNL",
        help="3-letter GROMACS residue name for the ligand (default: UNL)",
    )
    parser.add_argument(
        "--input-pdbqt",
        type=str,
        help="Path to input PDBQT file",
    )

    parser.add_argument(
        "--ph",
        type=float,
        default=None,
        help="pH for protonation",
    )
    parser.add_argument(
        "--full_protonate",
        action="store_true", 
        help="Fully protonate ligand (default: polar H only)",
    )

    args = parser.parse_args()

    step1_extract_and_prepare_ligand(
        input_pdbqt=args.input_pdbqt,
        ligand_name=args.ligand_name,
        residue_name=args.residue_name,
        ph=args.ph,
        full_protonate=args.full_protonate,
        ligand_dir=args.ligand_dir,
    )

if __name__ == "__main__":
    main()
