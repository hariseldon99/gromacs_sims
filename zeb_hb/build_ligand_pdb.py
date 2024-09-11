#!/usr/bin/env python

import sys
import os
import numpy as np
from gromacs_py import gmx
lig_name = 'phz'
lig_smiles = 'C1=CC=C(C=C1)NN'
DATA_OUT="./zeb_hb_phz_docking/input"

if __name__ == "__main__":
    os.makedirs(DATA_OUT, exist_ok = True)
    lig_resname = lig_name.upper()
    gmx.gmxsys.ambertools.smile_to_pdb(lig_smiles,os.path.join(DATA_OUT,f"{lig_name}.pdb"),lig_resname)
