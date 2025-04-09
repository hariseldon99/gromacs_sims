#!/usr/bin/env python

import sys
import os
import numpy as np
from gromacs_py import gmx
lig_name = 'lfx'
lig_smiles = 'C[C@H]1COC2=C3N1C=C(C(=O)C3=CC(=C2N4CCN(CC4)C)F)C(=O)O'
DATA_OUT="./D0E7M2_PSEAI_Levofloxacin_docking/input"

if __name__ == "__main__":
    os.makedirs(DATA_OUT, exist_ok = True)
    lig_resname = lig_name.upper()
    gmx.gmxsys.ambertools.smile_to_pdb(lig_smiles,os.path.join(DATA_OUT,f"{lig_name}.pdb"),lig_resname)
