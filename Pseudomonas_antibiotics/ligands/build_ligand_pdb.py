#!/usr/bin/env python
import os
from gromacs_py import gmx
import argparse
lig_name = 'phz'
lig_smiles = 'C1=CC=C(C=C1)NN'
parser = argparse.ArgumentParser(description="Build ligand PDB file from SMILES string.")
parser.add_argument("--lig_name", required=True, help="Ligand name")
parser.add_argument("--lig_smiles", required=True, help="Ligand SMILES string")
    
args = parser.parse_args()

lig_name = args.lig_name
lig_smiles = args.lig_smiles
current_directory = os.getcwd()
if __name__ == "__main__":
    lig_resname = lig_name.upper()
    gmx.gmxsys.ambertools.smile_to_pdb(lig_smiles,os.path.join(current_directory,f"{lig_name}.pdb"),lig_resname)
