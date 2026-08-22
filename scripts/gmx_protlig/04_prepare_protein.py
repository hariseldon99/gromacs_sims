#!/usr/bin/env python3
"""
Step 4 — Prepare the protein topology with pdb2gmx,
then combine protein and ligand coordinates into a single GRO file.
CLI wrapper around gmx_protlig.step4_prepare_protein.
"""

import os
import sys
import argparse

try:
    from gmx_protlig import step4_prepare_protein
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gmx_protlig import step4_prepare_protein

def main():
    parser = argparse.ArgumentParser(description="Prepare protein topology and combine with ligand coordinates.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory (default: ligand)")
    parser.add_argument("--protein-pdb", "--protein_pdb", type=str, required=True, help="Path to protein PDB file")
    parser.add_argument("--prepare-top", "--prepare_top", action="store_true", help="Run pdb2gmx to prepare protein topology")
    parser.add_argument("--ff", type=str, default="amber99sb-ildn", help="Force field (default: amber99sb-ildn)")

    args = parser.parse_args()

    step4_prepare_protein(
        protein_pdb=args.protein_pdb,
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        prepare_top=args.prepare_top,
        ff=args.ff,
    )

if __name__ == "__main__":
    main()
