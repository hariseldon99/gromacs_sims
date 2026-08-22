#!/usr/bin/env python3
"""
Example script demonstrating Pythonic batch processing of protein-ligand complexes
using the gmx_protlig API.

Run inside Singularity:
    bash run_local.sh example_batch_run.py
"""

import logging
try:
    from gmx_protlig import ProteinLigandPipeline, BatchPipelineRunner, step1_extract_and_prepare_ligand
except ImportError:
    from protlig_api import ProteinLigandPipeline, BatchPipelineRunner, step1_extract_and_prepare_ligand

# Configure logging to show progress in terminal
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
    datefmt="%H:%M:%S"
)

def example_1_single_complex_pipeline():
    """Example 1: Run local preparation for a single complex using Python object API."""
    print("\n--- Example 1: Single Complex Pipeline ---")
    pipeline = ProteinLigandPipeline(
        protein_pdb="protein.pdb",
        ligand_pdbqt="ligand_docked.pdbqt",
        ligand_name="Ligand1",
        residue_name="UNL",
        work_dir="simulations/complex_1",
    )
    
    # Run steps 1 to 7 (local preparation & EM)
    # results = pipeline.run_all_local_prep(nt_em=4)
    # print("Step 7 completed! EM checkpoint:", results["step7"]["checkpoint_path"])


def example_2_batch_from_manifest():
    """Example 2: Run batch pipeline from a CSV or JSON manifest."""
    print("\n--- Example 2: Batch Processing from Manifest ---")
    runner = BatchPipelineRunner(manifest="manifest_example.csv")
    
    # Run batch up to step 7 (local prep before HPC submission)
    # results = runner.run_batch(up_to_step=7, continue_on_error=True)
    # print("Batch run summary:", results)


def example_3_direct_step_functions():
    """Example 3: Call step functions directly for fine-grained script control."""
    print("\n--- Example 3: Direct Step Functions ---")
    # lig_info = step1_extract_and_prepare_ligand(
    #     input_pdbqt="ligand_docked.pdbqt",
    #     ligand_name="Ligand1",
    #     residue_name="UNL",
    #     work_dir="simulations/custom_dir",
    # )
    # print("Protonated ligand PDB created at:", lig_info["protonated_pdb"])


if __name__ == "__main__":
    print("=================================================")
    print("  gmx_protlig Modular Python API Examples")
    print("=================================================")
    example_1_single_complex_pipeline()
    example_2_batch_from_manifest()
    example_3_direct_step_functions()
