"""
ProteinLigandPipeline orchestrator for single-complex workflow management.
"""

import os
import logging
from typing import Optional, Dict, Any

from .steps import (
    step1_extract_and_prepare_ligand,
    step2_run_acpype,
    step3_ligand_invacuo,
    step4_prepare_protein,
    step5_merge_topologies,
    step6_solvate_ions_index,
    step7_energy_minimization,
    step8_hpc_equi_prod,
)

logger = logging.getLogger("gmx_protlig.pipeline")


class ProteinLigandPipeline:
    """
    Object-oriented Python interface to run protein-ligand MD simulation steps
    for a single complex.
    """

    def __init__(
        self,
        protein_pdb: str,
        ligand_pdbqt: str,
        ligand_name: str = "Ligand",
        residue_name: str = "LIG",
        work_dir: str = ".",
        ph: Optional[float] = None,
        full_protonate: bool = False,
        prepare_top: bool = False,
        ion_conc: float = 0.15,
        box_dist: float = 1.5,
        box_type: str = "dodecahedron",
        prod_time_ns: float = 100.0,
    ):
        self.protein_pdb = os.path.abspath(protein_pdb)
        self.ligand_pdbqt = os.path.abspath(ligand_pdbqt)
        self.ligand_name = ligand_name
        self.residue_name = residue_name.upper()
        self.work_dir = os.path.abspath(work_dir)
        self.ph = ph
        self.full_protonate = full_protonate
        self.prepare_top = prepare_top
        self.ion_conc = ion_conc
        self.box_dist = box_dist
        self.box_type = box_type
        self.prod_time_ns = prod_time_ns

        self.results: Dict[str, Any] = {}

    def run_step1(self) -> Dict[str, Any]:
        """Step 1: Extract & prepare ligand."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 1...")
        res = step1_extract_and_prepare_ligand(
            input_pdbqt=self.ligand_pdbqt,
            ligand_name=self.ligand_name,
            residue_name=self.residue_name,
            ph=self.ph,
            full_protonate=self.full_protonate,
            work_dir=self.work_dir,
        )
        self.results["step1"] = res
        return res

    def run_step2(self, net_charge: int = 0, charge_method: str = "bcc") -> Dict[str, Any]:
        """Step 2: Run ACPYPE topology generation."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 2...")
        res = step2_run_acpype(
            residue_name=self.residue_name,
            net_charge=net_charge,
            charge_method=charge_method,
            work_dir=self.work_dir,
        )
        self.results["step2"] = res
        return res

    def run_step3(self, ntomp: Optional[int] = None) -> Dict[str, Any]:
        """Step 3: Optional in-vacuo ligand test."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 3...")
        res = step3_ligand_invacuo(
            residue_name=self.residue_name,
            ntomp=ntomp,
            work_dir=self.work_dir,
        )
        self.results["step3"] = res
        return res

    def run_step4(self) -> Dict[str, Any]:
        """Step 4: Prepare protein topology and combine coordinates."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 4...")
        res = step4_prepare_protein(
            protein_pdb=self.protein_pdb,
            residue_name=self.residue_name,
            prepare_top=self.prepare_top,
            work_dir=self.work_dir,
        )
        self.results["step4"] = res
        return res

    def run_step5(self) -> Dict[str, Any]:
        """Step 5: Merge topologies."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 5...")
        res = step5_merge_topologies(
            residue_name=self.residue_name,
            protein_pdb=self.protein_pdb,
            work_dir=self.work_dir,
        )
        self.results["step5"] = res
        return res

    def run_step6(self) -> Dict[str, Any]:
        """Step 6: Solvate, add ions, create custom index groups."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 6...")
        checkpoint_in = self.results.get("step5", {}).get("checkpoint_path")
        res = step6_solvate_ions_index(
            residue_name=self.residue_name,
            ion_conc=self.ion_conc,
            box_dist=self.box_dist,
            box_type=self.box_type,
            checkpoint_in=checkpoint_in,
            work_dir=self.work_dir,
        )
        self.results["step6"] = res
        return res

    def run_step7(self, nt: int = 1, gpu: bool = False) -> Dict[str, Any]:
        """Step 7: Energy minimisation."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 7...")
        checkpoint_in = self.results.get("step6", {}).get("checkpoint_path")
        res = step7_energy_minimization(
            nt=nt,
            gpu=gpu,
            checkpoint_in=checkpoint_in,
            work_dir=self.work_dir,
        )
        self.results["step7"] = res
        return res

    def run_step8(self, n_cores: Optional[int] = None, gpu_id: str = "0") -> Dict[str, Any]:
        """Step 8: Equilibration and production MD."""
        logger.info(f"Pipeline [{self.ligand_name}]: Starting Step 8...")
        checkpoint_in = self.results.get("step7", {}).get("checkpoint_path")
        res = step8_hpc_equi_prod(
            n_cores=n_cores,
            gpu_id=gpu_id,
            prod_time_ns=self.prod_time_ns,
            checkpoint_in=checkpoint_in,
            work_dir=self.work_dir,
        )
        self.results["step8"] = res
        return res

    def run_all_local_prep(self, run_vacuo_test: bool = False, nt_em: int = 1, gpu_em: bool = False) -> Dict[str, Any]:
        """
        Run local preparation steps (1 through 7).
        """
        logger.info(f"Starting local prep for {self.protein_pdb} + {self.ligand_name} in {self.work_dir}")
        self.run_step1()
        self.run_step2()
        if run_vacuo_test:
            self.run_step3()
        self.run_step4()
        self.run_step5()
        self.run_step6()
        self.run_step7(nt=nt_em, gpu=gpu_em)
        logger.info(f"Local preparation complete for {self.ligand_name}!")
        return self.results

    def run_full_pipeline(self, run_vacuo_test: bool = False, nt_em: int = 1, gpu_em: bool = False, n_cores_md: Optional[int] = None, gpu_id_md: str = "0") -> Dict[str, Any]:
        """
        Run entire pipeline (steps 1 through 8).
        """
        self.run_all_local_prep(run_vacuo_test=run_vacuo_test, nt_em=nt_em, gpu_em=gpu_em)
        self.run_step8(n_cores=n_cores_md, gpu_id=gpu_id_md)
        return self.results
