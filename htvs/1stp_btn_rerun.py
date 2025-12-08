#/usr/bin/env python
import importmonkey
import os
home = os.environ["HOME"]

importmonkey.add_path(os.path.join(home, "gitrepos/gromacs_sims/scripts"))
from htvs_amdock import run_htvs


PROTEIN_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn/1STP_expt.pdb")
LIGAND_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn/btn.pdb")

complex_1stp_btn_rerun = [{
            "protein_name": "1stp_protein",
            "ligand_name": "btn",
            "protein_pdb": PROTEIN_PDB,
            "ligand_pdb": LIGAND_PDB,
            "working_dir": os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn_rerun/job_1stp_btn_script_20251208"),
        }]


# Uncomment this to rerun the docking
run_htvs(complex_1stp_btn_rerun)