"""
gmx_protlig - Modular Python Package for GROMACS Protein-Ligand MD Simulations
"""

__version__ = "1.0.0"
__author__ = "Antigravity / GROMACS MD Team"

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

from .pipeline import ProteinLigandPipeline
from .batch import BatchPipelineRunner
from .utils import (
    working_directory,
    run_command,
    save_checkpoint,
    load_latest_checkpoint,
    relocate_checkpoint_paths,
    read_gro,
    write_combined_gro,
    generate_custom_index_groups,
    rewrite_include_paths,
    extract_atomtypes_block,
    extract_nonbonded_params_block,
    strip_global_directives,
)

__all__ = [
    "__version__",
    "step1_extract_and_prepare_ligand",
    "step2_run_acpype",
    "step3_ligand_invacuo",
    "step4_prepare_protein",
    "step5_merge_topologies",
    "step6_solvate_ions_index",
    "step7_energy_minimization",
    "step8_hpc_equi_prod",
    "ProteinLigandPipeline",
    "BatchPipelineRunner",
    "working_directory",
    "run_command",
    "save_checkpoint",
    "load_latest_checkpoint",
    "relocate_checkpoint_paths",
    "read_gro",
    "write_combined_gro",
    "generate_custom_index_groups",
    "rewrite_include_paths",
    "extract_atomtypes_block",
    "extract_nonbonded_params_block",
    "strip_global_directives",
]
