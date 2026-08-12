"""
BatchPipelineRunner for iterative processing of multiple protein-ligand complexes.
"""

import os
import csv
import json
import time
import logging
import traceback
from typing import List, Dict, Any, Union, Optional

from .pipeline import ProteinLigandPipeline

logger = logging.getLogger("gmx_protlig.batch")


class BatchPipelineRunner:
    """
    Batch runner to execute MD pipeline across large numbers of protein-ligand pairs.
    Supports CSV/JSON manifest files and provides summary reports.
    """

    def __init__(self, manifest: Optional[Union[str, List[Dict[str, Any]]]] = None):
        self.complexes: List[Dict[str, Any]] = []
        if manifest is not None:
            if isinstance(manifest, str):
                self.load_manifest(manifest)
            elif isinstance(manifest, list):
                self.complexes = manifest

    def load_manifest(self, filepath: str):
        """Load manifest from CSV or JSON file."""
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Manifest file not found: {filepath}")

        ext = os.path.splitext(filepath)[1].lower()
        if ext == ".csv":
            with open(filepath, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                self.complexes = list(reader)
        elif ext == ".json":
            with open(filepath, "r", encoding="utf-8") as fh:
                data = json.load(fh)
                if isinstance(data, list):
                    self.complexes = data
                elif isinstance(data, dict) and "complexes" in data:
                    self.complexes = data["complexes"]
                else:
                    raise ValueError("JSON manifest must be a list or object with 'complexes' key.")
        else:
            raise ValueError(f"Unsupported manifest extension '{ext}'. Use .csv or .json")

        logger.info(f"Loaded {len(self.complexes)} complex definitions from {filepath}")

    def run_batch(
        self,
        up_to_step: int = 7,
        continue_on_error: bool = True,
        output_summary_dir: str = ".",
    ) -> List[Dict[str, Any]]:
        """
        Run the pipeline iteratively over all complexes.

        :param up_to_step: Step number to run up to (e.g. 7 for local prep, 8 for full pipeline)
        :param continue_on_error: If True, log error and continue next complex on failure
        :param output_summary_dir: Directory to save batch summary report files
        :return: List of summary result dicts
        """
        results_summary = []
        total_start = time.time()

        logger.info(f"Starting batch execution for {len(self.complexes)} complexes (Up to step {up_to_step})...")

        for idx, comp in enumerate(self.complexes, start=1):
            comp_name = comp.get("ligand_name", f"complex_{idx}")
            work_dir = comp.get("work_dir", f"simulations/{comp_name}")
            start_time = time.time()

            summary_item = {
                "index": idx,
                "protein_pdb": comp.get("protein_pdb"),
                "ligand_pdbqt": comp.get("ligand_pdbqt"),
                "ligand_name": comp_name,
                "residue_name": comp.get("residue_name", "LIG"),
                "work_dir": work_dir,
                "status": "PENDING",
                "step_reached": 0,
                "elapsed_sec": 0.0,
                "error": None,
            }

            logger.info(f"\n==================================================")
            logger.info(f"  Batch [{idx}/{len(self.complexes)}]: Processing {comp_name}")
            logger.info(f"  Workdir: {work_dir}")
            logger.info(f"==================================================")

            try:
                pipeline = ProteinLigandPipeline(
                    protein_pdb=comp["protein_pdb"],
                    ligand_pdbqt=comp["ligand_pdbqt"],
                    ligand_name=comp_name,
                    residue_name=comp.get("residue_name", "LIG"),
                    work_dir=work_dir,
                    ph=float(comp["ph"]) if comp.get("ph") else None,
                    full_protonate=str(comp.get("full_protonate", "")).lower() in ("true", "1", "yes"),
                    prepare_top=str(comp.get("prepare_top", "")).lower() in ("true", "1", "yes"),
                    ion_conc=float(comp.get("ion_conc", 0.15)),
                    box_dist=float(comp.get("box_dist", 1.5)),
                    box_type=comp.get("box_type", "dodecahedron"),
                    prod_time_ns=float(comp.get("prod_time_ns", 100.0)),
                )

                # Run steps pythonically
                pipeline.run_step1()
                summary_item["step_reached"] = 1

                pipeline.run_step2(
                    net_charge=int(comp.get("net_charge", 0)),
                    charge_method=comp.get("charge_method", "bcc"),
                )
                summary_item["step_reached"] = 2

                if str(comp.get("run_invacuo", "")).lower() in ("true", "1", "yes"):
                    pipeline.run_step3()
                    summary_item["step_reached"] = 3

                if up_to_step >= 4:
                    pipeline.run_step4()
                    summary_item["step_reached"] = 4

                if up_to_step >= 5:
                    pipeline.run_step5()
                    summary_item["step_reached"] = 5

                if up_to_step >= 6:
                    pipeline.run_step6()
                    summary_item["step_reached"] = 6

                if up_to_step >= 7:
                    pipeline.run_step7(
                        nt=int(comp.get("nt_em", 1)),
                        gpu=str(comp.get("gpu_em", "")).lower() in ("true", "1", "yes"),
                    )
                    summary_item["step_reached"] = 7

                if up_to_step >= 8:
                    pipeline.run_step8(
                        n_cores=int(comp["n_cores_md"]) if comp.get("n_cores_md") else None,
                        gpu_id=comp.get("gpu_id_md", "0"),
                    )
                    summary_item["step_reached"] = 8

                summary_item["status"] = "SUCCESS"
                logger.info(f"[SUCCESS] {comp_name} completed up to step {summary_item['step_reached']}.")

            except Exception as exc:
                err_msg = f"{type(exc).__name__}: {str(exc)}"
                logger.error(f"[FAILURE] Complex {comp_name} failed at step {summary_item['step_reached'] + 1}: {err_msg}")
                logger.debug(traceback.format_exc())
                summary_item["status"] = "FAILED"
                summary_item["error"] = err_msg

                if not continue_on_error:
                    summary_item["elapsed_sec"] = round(time.time() - start_time, 2)
                    results_summary.append(summary_item)
                    raise exc
            finally:
                summary_item["elapsed_sec"] = round(time.time() - start_time, 2)
                results_summary.append(summary_item)

        total_elapsed = round(time.time() - total_start, 2)
        logger.info(f"\n==================================================")
        logger.info(f"  Batch Run Finished in {total_elapsed} seconds.")
        logger.info(f"==================================================")

        # Write batch report files
        os.makedirs(output_summary_dir, exist_ok=True)
        summary_json = os.path.join(output_summary_dir, "batch_summary.json")
        summary_csv = os.path.join(output_summary_dir, "batch_summary.csv")

        with open(summary_json, "w", encoding="utf-8") as fh:
            json.dump({
                "total_complexes": len(self.complexes),
                "total_elapsed_sec": total_elapsed,
                "results": results_summary,
            }, fh, indent=2)

        fieldnames = ["index", "ligand_name", "residue_name", "protein_pdb", "ligand_pdbqt", "work_dir", "status", "step_reached", "elapsed_sec", "error"]
        with open(summary_csv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for r in results_summary:
                writer.writerow({k: r.get(k) for k in fieldnames})

        logger.info(f"Saved batch summaries to {summary_json} and {summary_csv}")
        return results_summary
