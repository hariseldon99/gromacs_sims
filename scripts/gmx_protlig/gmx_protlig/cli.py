"""
Unified CLI command dispatcher for gmx_protlig.
"""

import sys
import argparse
import logging
from typing import List, Optional

from .steps import (
    step1_cli,
    step2_cli,
    step3_cli,
    step4_cli,
    step5_cli,
    step6_cli,
    step7_cli,
    step8_cli,
)
from .batch import batch_cli
from .pipeline import ProteinLigandPipeline


def run_pipeline_cli():
    parser = argparse.ArgumentParser(description="Run complete protein-ligand simulation pipeline for a single complex.")
    parser.add_argument("--protein-pdb", "--protein_pdb", type=str, required=True, help="Path to protein PDB")
    parser.add_argument("--ligand-pdbqt", "--ligand_pdbqt", type=str, required=True, help="Path to docked ligand PDBQT")
    parser.add_argument("--ligand-name", "--ligand_name", type=str, default="UnknownLigand", help="Ligand name")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name (default: UNL)")
    parser.add_argument("--work-dir", "--work_dir", type=str, default=".", help="Working directory")
    parser.add_argument("--ph", type=float, default=None, help="pH for protonation")
    parser.add_argument("--full-protonate", action="store_true", help="Fully protonate ligand")
    parser.add_argument("--prepare-top", action="store_true", help="Run pdb2gmx to prepare topology")
    parser.add_argument("--ff", type=str, default="amber99sb-ildn", help="Force field")
    parser.add_argument("--net-charge", type=int, default=0, help="Net charge")
    parser.add_argument("--charge-method", type=str, default="bcc", help="Charge method")
    parser.add_argument("--ion-conc", type=float, default=0.15, help="Ion concentration (mol/L)")
    parser.add_argument("--box-dist", type=float, default=1.5, help="Box distance (nm)")
    parser.add_argument("--nt-em", type=int, default=1, help="CPU threads for EM")
    parser.add_argument("--gpu-em", action="store_true", help="Use GPU for EM")
    parser.add_argument("--up-to-step", type=int, default=7, help="Step to run up to (default: 7)")
    parser.add_argument("--prod-time", type=float, default=100.0, help="Production time in ns (if running step 8)")
    parser.add_argument("--n-cores-md", type=int, default=None, help="Cores for MD (step 8)")
    parser.add_argument("--gpu-id-md", type=str, default="0", help="GPU ID for MD (step 8)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%H:%M:%S"
    )

    pipeline = ProteinLigandPipeline(
        protein_pdb=args.protein_pdb,
        ligand_pdbqt=args.ligand_pdbqt,
        ligand_name=args.ligand_name,
        residue_name=args.residue_name,
        work_dir=args.work_dir,
        ph=args.ph,
        full_protonate=args.full_protonate,
        prepare_top=args.prepare_top,
        ff=args.ff,
        ion_conc=args.ion_conc,
        box_dist=args.box_dist,
        prod_time_ns=args.prod_time,
    )

    if args.up_to_step >= 7 and args.up_to_step < 8:
        pipeline.run_all_local_prep(
            nt_em=args.nt_em,
            gpu_em=args.gpu_em,
            net_charge=args.net_charge,
            charge_method=args.charge_method,
        )
    elif args.up_to_step >= 8:
        pipeline.run_full_pipeline(
            nt_em=args.nt_em,
            gpu_em=args.gpu_em,
            net_charge=args.net_charge,
            charge_method=args.charge_method,
            n_cores_md=args.n_cores_md,
            gpu_id_md=args.gpu_id_md,
        )
    else:
        pipeline.resume_from_step(start_step=1, up_to_step=args.up_to_step)


def main():
    """Main CLI dispatcher."""
    subcommands = {
        "step1": ("Extract and prepare ligand from PDBQT", step1_cli),
        "step2": ("Run ACPYPE topology generation", step2_cli),
        "step3": ("Run in-vacuo ligand test simulation", step3_cli),
        "step4": ("Prepare protein topology and combined coordinates", step4_cli),
        "step5": ("Merge protein and ligand topologies", step5_cli),
        "step6": ("Solvate, add ions, and create index groups", step6_cli),
        "step7": ("Run energy minimisation", step7_cli),
        "step8": ("Run equilibration and production MD (HPC)", step8_cli),
        "run": ("Run single complex pipeline", run_pipeline_cli),
        "batch": ("Run batch pipeline from CSV/JSON manifest", batch_cli),
    }

    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: gmx-protlig <subcommand> [options]\n")
        print("Available subcommands:")
        for name, (desc, _) in subcommands.items():
            print(f"  {name:10s}  {desc}")
        print("\nRun 'gmx-protlig <subcommand> --help' for options of a specific subcommand.")
        sys.exit(0)

    cmd = sys.argv[1].lower()
    if cmd not in subcommands:
        print(f"Error: Unknown subcommand '{cmd}'. Run 'gmx-protlig --help' for available commands.")
        sys.exit(1)

    # Shift sys.argv so subparser handles arguments correctly
    sys.argv.pop(1)
    _, func = subcommands[cmd]
    func()


if __name__ == "__main__":
    main()
