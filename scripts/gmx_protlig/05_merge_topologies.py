#!/usr/bin/env python3
"""
Step 5 — Merge the protein topology with the ligand topology.
CLI wrapper around gmx_protlig.step5_merge_topologies.
"""

import os
import sys
import argparse

try:
    from gmx_protlig import step5_merge_topologies
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gmx_protlig import step5_merge_topologies

def main():
    parser = argparse.ArgumentParser(description="Merge protein and ligand topologies.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory (default: ligand)")
    parser.add_argument("--protein-top", "--protein_top", type=str, default=None, help="Protein topology path")
    parser.add_argument("--protein-pdb", "--protein_pdb", type=str, default=None, help="Path to protein PDB file")

    args = parser.parse_args()

    step5_merge_topologies(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        protein_top=args.protein_top,
        protein_pdb=args.protein_pdb,
    )

if __name__ == "__main__":
    main()
