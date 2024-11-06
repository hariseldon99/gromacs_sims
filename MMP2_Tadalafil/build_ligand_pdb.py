#!/usr/bin/env python

import sys
import os
import numpy as np
from gromacs_py import gmx
lig_name = 'tdl'
lig_smiles = 'CN1CC(=O)N2[C@@H](C1=O)CC3=C([C@H]2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36'
DATA_OUT="./input"

if __name__ == "__main__":
    os.makedirs(DATA_OUT, exist_ok = True)
    lig_resname = lig_name.upper()
    gmx.gmxsys.ambertools.smile_to_pdb(lig_smiles,os.path.join(DATA_OUT,f"{lig_name}.pdb"),lig_resname)
