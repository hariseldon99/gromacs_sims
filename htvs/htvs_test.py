#/usr/bin/env python
import importmonkey
import os
home = os.environ["HOME"]

importmonkey.add_path(os.path.join(home, "gitrepos/gromacs_sims/scripts"))
from htvs_amdock import run_htvs

PROTEIN_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1hew_nag/experiment/1HEW_expt.pdb")
LIGAND_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1hew_nag/experiment/nag.pdb")

testjob = [{
            "protein_name": "1hew_protein",
            "ligand_name": "nag",
            "protein_pdb": PROTEIN_PDB,
            "ligand_pdb": LIGAND_PDB,
            "working_dir": os.path.join(home, "gitrepos/gromacs_sims/htvs/1hew_nag/temp_test"),
        }]

# Uncomment this to rerun the docking
run_htvs(testjob, verbose=True)
