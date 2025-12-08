#/usr/bin/env python
import importmonkey
import os
home = os.environ["HOME"]

importmonkey.add_path(os.path.join(home, "gitrepos/gromacs_sims/scripts"))
from htvs_amdock import run_htvs

PROTEIN_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1fkb_rap/1FKB_expt.pdb")
LIGAND_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1fkb_rap/rap.pdb")

docking_jobs = [{
            "protein_name": "1fkb_protein",
            "ligand_name": "rap",
            "protein_pdb": PROTEIN_PDB,
            "ligand_pdb": LIGAND_PDB,
            "working_dir": os.path.join(home, "gitrepos/gromacs_sims/htvs/1fkb_rap/job_1fkb_rap_script_20251208"),
        }]

PROTEIN_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn/1STP_expt.pdb")
LIGAND_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn/btn.pdb")

docking_jobs.append({
            "protein_name": "1stp_protein",
            "ligand_name": "btn",
            "protein_pdb": PROTEIN_PDB,
            "ligand_pdb": LIGAND_PDB,
            "working_dir": os.path.join(home, "gitrepos/gromacs_sims/htvs/1stp_btn/job_1stp_btn_script_20251208"),
        })


PROTEIN_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1afk_pap/1afk_expt.pdb")
LIGAND_PDB = os.path.join(home, "gitrepos/gromacs_sims/htvs/1afk_pap/pap.pdb")

docking_jobs.append({
            "protein_name": "1afk_protein",
            "ligand_name": "pap",
            "protein_pdb": PROTEIN_PDB,
            "ligand_pdb": LIGAND_PDB,
            "working_dir": os.path.join(home, "gitrepos/gromacs_sims/htvs/1afk_pap/job_1afk_pap_script_20251208"),
        })

# Uncomment this to rerun the docking
run_htvs(docking_jobs)
