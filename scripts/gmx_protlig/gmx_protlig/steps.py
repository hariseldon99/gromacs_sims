"""
Modular step functions for GROMACS Protein-Ligand MD simulation pipeline.
Each step function can be called pythonically or via CLI.
"""

import os
import sys
import glob
import pickle
import datetime
import shutil
import logging
import subprocess
import multiprocessing
import re
import argparse
from typing import Dict, Any, Optional, Union

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    np = None
    pd = None
    matplotlib = None
    plt = None

try:
    from gromacs_py import gmx
except ImportError:
    gmx = None

from .utils import (
    working_directory,
    run_command,
    write_combined_gro,
    extract_atomtypes_block,
    extract_nonbonded_params_block,
    strip_global_directives,
    rewrite_include_paths,
    generate_custom_index_groups,
    save_checkpoint,
    load_latest_checkpoint,
    relocate_checkpoint_paths,
)

logger = logging.getLogger("gmx_protlig.steps")


def _ensure_gromacs_py():
    """Ensure gromacs_py is available when executing simulation steps."""
    if gmx is None:
        raise ImportError(
            "The 'gromacs_py' package is required for this step. "
            "Please run inside the Singularity container (e.g. via run_local.sh) "
            "or ensure gromacs_py is installed in your Python environment."
        )


def step1_extract_and_prepare_ligand(
    input_pdbqt: str,
    ligand_name: str = "UnknownLigand",
    residue_name: str = "UNL",
    ph: Optional[float] = None,
    full_protonate: bool = False,
    ligand_dir: str = "ligand",
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 1: Extract best pose from PDBQT, convert to PDB, and protonate.
    """
    with working_directory(work_dir):
        residue_name = residue_name.strip().upper()
        if len(residue_name) != 3:
            raise ValueError(f"Residue name must be exactly 3 letters: '{residue_name}'")

        abs_pdbqt = os.path.abspath(input_pdbqt)
        if not os.path.exists(abs_pdbqt):
            raise FileNotFoundError(f"Input PDBQT not found: {abs_pdbqt}")

        os.makedirs(ligand_dir, exist_ok=True)
        logger.info(f"[Step 1] Reading docked poses for {ligand_name} from {abs_pdbqt}...")

        with open(abs_pdbqt, "r", encoding="utf-8", errors="replace") as fh:
            raw = fh.read()

        blocks = {}
        current_model = None
        current_lines = []

        for line in raw.splitlines(keepends=True):
            if line.startswith("MODEL "):
                if current_model is not None and current_lines:
                    blocks[current_model] = current_lines
                try:
                    current_model = int(line.split()[1])
                except (IndexError, ValueError):
                    current_model = len(blocks) + 1
                current_lines = [line]
            elif line.startswith("ENDMDL") and current_model is not None:
                current_lines.append(line)
                blocks[current_model] = current_lines
                current_model = None
                current_lines = []
            elif current_model is not None:
                current_lines.append(line)

        if current_model is not None and current_lines:
            blocks[current_model] = current_lines

        if not blocks:
            # Single model PDBQT without MODEL tags
            blocks[1] = raw.splitlines(keepends=True)

        scores = {}
        for mid, lines in blocks.items():
            for l in lines:
                if "VINA RESULT:" in l or "REMARK VINA RESULT:" in l:
                    try:
                        parts = l.split()
                        # Find the first float after RESULT:
                        idx_res = next(i for i, p in enumerate(parts) if "RESULT" in p)
                        scores[mid] = float(parts[idx_res + 1])
                        break
                    except (StopIteration, IndexError, ValueError):
                        pass

        if scores:
            best_model = min(scores, key=scores.get)
            best_score = scores[best_model]
        else:
            best_model = min(blocks.keys())
            best_score = 0.0

        best_pdbqt = os.path.join(ligand_dir, f"{ligand_name.lower()}_pose1.pdbqt")
        with open(best_pdbqt, "w", encoding="utf-8") as fh:
            fh.writelines(blocks[best_model])

        raw_pdb = os.path.join(ligand_dir, f"{ligand_name.lower()}_raw.pdb")
        run_command([
            "obabel",
            "-ipdbqt", best_pdbqt,
            "-opdb", "-O", raw_pdb,
            "--partialcharge", "none",
        ])

        prot_pdb = os.path.join(ligand_dir, f"{residue_name}.pdb")
        cmd = ["obabel", "-ipdb", raw_pdb, "-opdb", "-O", prot_pdb]
        if full_protonate:
            cmd.append("-h")
        else:
            cmd.append("--AddPolarH")

        if ph is not None:
            cmd.extend(["-p", str(ph)])

        run_command(cmd)

        # Standardize residue name in output PDB
        with open(prot_pdb, "r", encoding="utf-8") as fh:
            pdb_lines = fh.readlines()

        new_pdb_lines = []
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                # Residue name in columns 18-20 (1-based: 17:20 in 0-based indexing)
                line = line[:17] + f"{residue_name:3s}" + line[20:]
            new_pdb_lines.append(line)

        with open(prot_pdb, "w", encoding="utf-8") as fh:
            fh.writelines(new_pdb_lines)

        logger.info(f"[Step 1 Complete] Prepared ligand PDB: {prot_pdb} (Score: {best_score})")
        return {
            "ligand_dir": ligand_dir,
            "best_pdbqt": best_pdbqt,
            "raw_pdb": raw_pdb,
            "protonated_pdb": prot_pdb,
            "best_model": best_model,
            "best_score": best_score,
        }


def step2_run_acpype(
    residue_name: str = "UNL",
    ligand_dir: str = "ligand",
    net_charge: int = 0,
    charge_method: str = "bcc",
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 2: Run acpype on protonated ligand PDB.
    """
    with working_directory(work_dir):
        resname = residue_name.strip().upper()
        ligand_pdb = os.path.join(ligand_dir, f"{resname}.pdb")
        if not os.path.exists(ligand_pdb):
            raise FileNotFoundError(f"Ligand PDB not found: {ligand_pdb}. Run Step 1 first.")

        logger.info(f"[Step 2] Running acpype on {ligand_pdb} (net charge: {net_charge})...")

        with working_directory(ligand_dir):
            cmd = ["acpype", "-i", f"{resname}.pdb", "-n", str(net_charge), "-d"]
            if charge_method:
                cmd.extend(["-c", charge_method])
            run_command(cmd)

        acpype_dir = os.path.join(ligand_dir, f"{resname}.acpype")
        gro_file = os.path.join(acpype_dir, f"{resname}_GMX.gro")
        top_file = os.path.join(acpype_dir, f"{resname}_GMX.top")
        itp_file = os.path.join(acpype_dir, f"{resname}_GMX.itp")
        posre_file = os.path.join(acpype_dir, f"posre_{resname}.itp")

        for f in (gro_file, top_file, itp_file):
            if not os.path.exists(f):
                raise FileNotFoundError(f"ACPYPE expected output missing: {f}")

        logger.info(f"[Step 2 Complete] ACPYPE output in: {acpype_dir}")
        return {
            "acpype_dir": acpype_dir,
            "gro": gro_file,
            "top": top_file,
            "itp": itp_file,
            "posre_itp": posre_file,
        }


def step3_ligand_invacuo(
    residue_name: str = "UNL",
    ligand_dir: str = "ligand",
    invacuo_dir: str = "invacuo",
    ntomp: Optional[int] = None,
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 3: In-vacuo ligand simulation (energy minimisation + short MD).
    """
    with working_directory(work_dir):
        resname = residue_name.strip().upper()
        if ntomp is None:
            ntomp = multiprocessing.cpu_count()

        acpype_dir = os.path.join(ligand_dir, f"{resname}.acpype")
        for f in (f"{resname}_GMX.gro", f"{resname}_GMX.top"):
            if not os.path.exists(os.path.join(acpype_dir, f)):
                raise FileNotFoundError(f"Missing acpype file in {acpype_dir}: {f}")

        os.makedirs(invacuo_dir, exist_ok=True)
        shutil.copy2(os.path.join(acpype_dir, f"{resname}_GMX.gro"), os.path.join(invacuo_dir, f"{resname}.gro"))
        shutil.copy2(os.path.join(acpype_dir, f"{resname}_GMX.top"), os.path.join(invacuo_dir, f"{resname}.top"))

        for itp in glob.glob(os.path.join(acpype_dir, "*.itp")):
            shutil.copy2(itp, invacuo_dir)

        em_mdp_content = """\
; In-vacuo energy minimisation
integrator              = steep
nsteps                  = 10000
emtol                   = 100.0
emstep                  = 0.01
nstlist                 = 1
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = cutoff
rcoulomb                = 1.2
rvdw                    = 1.2
pbc                     = no
"""
        md_mdp_content = """\
; In-vacuo MD — 100 ps
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 500
nstvout                 = 500
nstfout                 = 500
nstlog                  = 500
nstenergy               = 500
nstlist                 = 10
cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = cutoff
rcoulomb                = 1.2
rvdw                    = 1.2
pbc                     = no
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
constraints             = h-bonds
constraint_algorithm    = lincs
"""

        with open(os.path.join(invacuo_dir, "em_vac.mdp"), "w", encoding="utf-8") as fh:
            fh.write(em_mdp_content)
        with open(os.path.join(invacuo_dir, "md_vac.mdp"), "w", encoding="utf-8") as fh:
            fh.write(md_mdp_content)

        with working_directory(invacuo_dir):
            run_command([
                "gmx", "grompp",
                "-f", "em_vac.mdp",
                "-c", f"{resname}.gro",
                "-p", f"{resname}.top",
                "-o", "em_vac.tpr",
                "-maxwarn", "5",
            ])

            run_command([
                "gmx", "mdrun",
                "-nobackup", "-nocopyright", "-v",
                "-deffnm", "em_vac",
                "-ntmpi", "1",
                "-ntomp", str(ntomp),
                "-nb", "gpu",
                "-gpu_id", "0",
                "-pme", "cpu",
            ])

            run_command([
                "gmx", "grompp",
                "-f", "md_vac.mdp",
                "-c", "em_vac.gro",
                "-p", f"{resname}.top",
                "-o", "md_vac.tpr",
                "-maxwarn", "5",
            ])

            run_command([
                "gmx", "mdrun",
                "-nobackup", "-nocopyright", "-v",
                "-deffnm", "md_vac",
                "-ntmpi", "1",
                "-ntomp", str(ntomp),
                "-nb", "gpu",
                "-gpu_id", "0",
                "-pme", "cpu",
            ])

        logger.info(f"[Step 3 Complete] In-vacuo simulation finished in: {invacuo_dir}")
        return {
            "invacuo_dir": invacuo_dir,
            "em_tpr": os.path.join(invacuo_dir, "em_vac.tpr"),
            "md_tpr": os.path.join(invacuo_dir, "md_vac.tpr"),
            "md_trr": os.path.join(invacuo_dir, "md_vac.trr"),
            "md_gro": os.path.join(invacuo_dir, "md_vac.gro"),
        }


def step4_prepare_protein(
    protein_pdb: str,
    residue_name: str = "UNL",
    ligand_dir: str = "ligand",
    prepare_top: bool = False,
    ff: str = "amber99sb-ildn",
    protein_outdir: str = "protein",
    complex_outdir: str = "complex_raw",
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 4: Prepare protein topology and combine protein + ligand coordinates.
    """
    _ensure_gromacs_py()

    with working_directory(work_dir):
        resname = residue_name.strip().upper()
        abs_protein_pdb = os.path.abspath(protein_pdb)
        protein_name = os.path.splitext(os.path.basename(abs_protein_pdb))[0]

        ligand_gro = os.path.join(ligand_dir, f"{resname}.acpype", f"{resname}_GMX.gro")
        if not os.path.exists(abs_protein_pdb):
            raise FileNotFoundError(f"Protein PDB not found: {abs_protein_pdb}")
        if not os.path.exists(ligand_gro):
            raise FileNotFoundError(f"Ligand GRO missing: {ligand_gro}")

        os.makedirs(complex_outdir, exist_ok=True)
        os.makedirs(protein_outdir, exist_ok=True)
        logger.info(f"[Step 4] Preparing protein topology for {protein_name}...")

        protein_sys = gmx.GmxSys(name=protein_name, coor_file=abs_protein_pdb)

        if prepare_top:
            protein_sys.prepare_top(
                out_folder=protein_outdir,
                ff=ff,
                vsite="hydrogens",
                maxwarn=5,
            )
        else:
            protein_sys.add_top(
                out_folder=protein_outdir,
                ff=ff,
                pdb2gmx_option_dict={"ignh": None},
            )

        # Fix path bug: resolve coor_file and top_file relative to protein_outdir
        raw_coor = protein_sys.coor_file
        if not os.path.exists(raw_coor):
            raw_coor = os.path.join(protein_outdir, os.path.basename(raw_coor))

        raw_top = protein_sys.top_file
        if not os.path.exists(raw_top):
            raw_top = os.path.join(protein_outdir, os.path.basename(raw_top))

        prot_gro_path = os.path.join(protein_outdir, f"{protein_name}_pdb2gmx.gro")
        run_command(["gmx", "editconf", "-f", raw_coor, "-o", prot_gro_path])

        protein_sys.coor_file = prot_gro_path
        protein_sys.top_file = raw_top

        combined_gro = os.path.join(complex_outdir, "complex_raw.gro")
        n_atoms = write_combined_gro(
            prot_gro_path,
            ligand_gro,
            combined_gro,
            title=f"{protein_name}-{resname} complex",
        )

        chkpt_path = save_checkpoint(protein_sys, "prepare_protein")
        logger.info(f"[Step 4 Complete] Combined GRO: {combined_gro} ({n_atoms} atoms)")

        return {
            "protein_sys": protein_sys,
            "prot_gro": prot_gro_path,
            "prot_top": protein_sys.top_file,
            "complex_gro": combined_gro,
            "checkpoint_path": chkpt_path,
        }


def step5_merge_topologies(
    residue_name: str = "UNL",
    ligand_dir: str = "ligand",
    protein_top: Optional[str] = None,
    protein_pdb: Optional[str] = None,
    complex_dir: str = "complex",
    complex_gro_src: str = "complex_raw/complex_raw.gro",
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 5: Merge protein and ligand topologies.
    """
    _ensure_gromacs_py()

    with working_directory(work_dir):
        resname = residue_name.strip().upper()
        acpype_dir = os.path.join(ligand_dir, f"{resname}.acpype")

        if not os.path.exists(acpype_dir):
            raise FileNotFoundError(f"ACPYPE directory not found: {acpype_dir}")
        if not os.path.exists(complex_gro_src):
            raise FileNotFoundError(f"Combined GRO file not found: {complex_gro_src}")

        os.makedirs(complex_dir, exist_ok=True)

        prot_top_path = None
        if protein_top and os.path.exists(protein_top):
            prot_top_path = protein_top
        else:
            # Look for *.top in protein directory
            search_dir = os.path.dirname(protein_top) if protein_top else "protein"
            if not os.path.exists(search_dir):
                search_dir = "protein"
            if os.path.exists(search_dir):
                candidates = sorted(glob.glob(os.path.join(search_dir, "*.top")))
                # Prefer _pdb2gmx.top if present
                for cand in candidates:
                    if "_pdb2gmx" in cand:
                        prot_top_path = cand
                        break
                if not prot_top_path and candidates:
                    prot_top_path = candidates[0]

        if not prot_top_path or not os.path.exists(prot_top_path):
            raise FileNotFoundError(f"Cannot locate protein topology: {protein_top}")

        lig_top = os.path.join(acpype_dir, f"{resname}_GMX.top")
        lig_itp = os.path.join(acpype_dir, f"{resname}_GMX.itp")

        with open(prot_top_path, "r", encoding="utf-8") as fh:
            prot_text = fh.read()
        with open(lig_top, "r", encoding="utf-8") as fh:
            lig_top_text = fh.read()
        with open(lig_itp, "r", encoding="utf-8") as fh:
            lig_itp_text = fh.read()

        lig_atomtypes = extract_atomtypes_block(lig_top_text) or extract_atomtypes_block(lig_itp_text)
        lig_nonbonded = extract_nonbonded_params_block(lig_top_text) or extract_nonbonded_params_block(lig_itp_text)

        dst_itp = os.path.join(complex_dir, f"{resname}_GMX.itp")
        dst_posre = os.path.join(complex_dir, f"posre_{resname}.itp")

        clean_itp = strip_global_directives(lig_itp_text)
        with open(dst_itp, "w", encoding="utf-8") as fh:
            fh.write(clean_itp)

        src_posre = os.path.join(acpype_dir, f"posre_{resname}.itp")
        if os.path.exists(src_posre):
            shutil.copy2(src_posre, dst_posre)
        else:
            with open(dst_posre, "w", encoding="utf-8") as fh:
                fh.write(f"; Position restraints for {resname}\n[ position_restraints ]\n")

        for f in os.listdir(acpype_dir):
            if f.endswith(".itp") and f not in (f"{resname}_GMX.itp", f"posre_{resname}.itp"):
                shutil.copy2(os.path.join(acpype_dir, f), os.path.join(complex_dir, f))

        prot_text_rw = rewrite_include_paths(prot_text, prot_top_path, complex_dir)

        # 1. Insert atomtypes and nonbonded params after forcefield.itp
        ff_include_re = re.compile(r'(#include\s+"[^"]*(?:forcefield|ff)[^"]*\.itp"[^\n]*\n)', re.IGNORECASE)
        m = ff_include_re.search(prot_text_rw)
        if m:
            insert_after = m.end()
            block = f"\n; ── Ligand {resname} atomtypes ──────────────────────────\n" + lig_atomtypes
            if lig_nonbonded:
                block += "\n" + lig_nonbonded
            block += "; ─────────────────────────────────────────────────────────────\n\n"
            prot_text_rw = prot_text_rw[:insert_after] + block + prot_text_rw[insert_after:]
        elif lig_atomtypes:
            # Fallback: insert at beginning of topology
            prot_text_rw = f"\n; ── Ligand {resname} atomtypes ───\n{lig_atomtypes}\n" + prot_text_rw

        # 2. Insert ligand ITP include before water/ions include
        water_include_re = re.compile(r'(#include\s+"[^"]*(?:tip3p|tip4p|spc|spce|water|ions)[^"]*\.itp"[^\n]*\n)', re.IGNORECASE)
        lig_inc = f'#include "{resname}_GMX.itp"\n'
        m2 = water_include_re.search(prot_text_rw)
        if m2:
            prot_text_rw = prot_text_rw[:m2.start()] + f"; ── Ligand {resname} molecule type ───\n" + lig_inc + "\n" + prot_text_rw[m2.start():]
        else:
            # Fallback: insert before [ system ] or [ molecules ]
            sec_re = re.compile(r'(\[\s*(?:system|molecules)\s*\])', re.IGNORECASE)
            m_sec = sec_re.search(prot_text_rw)
            if m_sec:
                prot_text_rw = prot_text_rw[:m_sec.start()] + lig_inc + "\n" + prot_text_rw[m_sec.start():]
            else:
                prot_text_rw += "\n" + lig_inc

        # 3. Add ligand entry under [ molecules ]
        mol_sec_re = re.compile(r'\[\s*molecules\s*\][^\n]*\n', re.IGNORECASE)
        mol_m = mol_sec_re.search(prot_text_rw)
        if mol_m:
            mol_text = prot_text_rw[mol_m.end():]
            # Check if already present
            if not re.search(rf'^\s*{resname}\s+\d+', mol_text, re.MULTILINE):
                # Look for Protein line to append after, else append at end of [ molecules ]
                prot_line_re = re.compile(r'(Protein\S*\s+\d+\s*\n)', re.IGNORECASE)
                prot_matches = list(prot_line_re.finditer(prot_text_rw))
                if prot_matches:
                    last_match = prot_matches[-1]
                    prot_text_rw = prot_text_rw[:last_match.end()] + f"{resname:3s}            1\n" + prot_text_rw[last_match.end():]
                else:
                    prot_text_rw = prot_text_rw[:mol_m.end()] + f"{resname:3s}            1\n" + prot_text_rw[mol_m.end():]

        complex_top = os.path.join(complex_dir, "complex.top")
        with open(complex_top, "w", encoding="utf-8") as fh:
            fh.write(prot_text_rw)

        complex_gro = os.path.join(complex_dir, "complex.gro")
        shutil.copy2(complex_gro_src, complex_gro)

        protein_name = os.path.splitext(os.path.basename(protein_pdb))[0] if protein_pdb else "protein"
        complex_name = f"complex_{protein_name.lower()}_{resname.lower()}"
        complex_sys = gmx.GmxSys(name=complex_name, coor_file=complex_gro)
        complex_sys.top_file = complex_top

        chkpt_path = save_checkpoint(complex_sys, "merge_topologies")
        logger.info(f"[Step 5 Complete] Merged topology: {complex_top}")

        return {
            "complex_sys": complex_sys,
            "complex_top": complex_top,
            "complex_gro": complex_gro,
            "dst_itp": dst_itp,
            "dst_posre": dst_posre,
            "checkpoint_path": chkpt_path,
        }


def step6_solvate_ions_index(
    residue_name: str = "UNL",
    ion_conc: float = 0.15,
    box_dist: float = 1.5,
    box_type: str = "dodecahedron",
    solv_outdir: str = "complex_solv",
    maxwarn: int = 5,
    checkpoint_in: Optional[Any] = None,
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 6: Box creation, solvation, ion addition, and index group generation.
    """
    _ensure_gromacs_py()

    with working_directory(work_dir):
        resname = residue_name.strip().upper()

        if isinstance(checkpoint_in, gmx.GmxSys):
            complex_sys = checkpoint_in
        elif isinstance(checkpoint_in, str) and os.path.exists(checkpoint_in):
            with open(checkpoint_in, "rb") as fh:
                complex_sys = pickle.load(fh)
        else:
            complex_sys, _ = load_latest_checkpoint("merge_topologies")

        relocate_checkpoint_paths(complex_sys, ".")

        sys_name = complex_sys.name if complex_sys.name else f"complex_{resname.lower()}"

        logger.info(f"[Step 6] Creating box and solvating {sys_name}...")
        complex_sys.create_box(dist=box_dist, box_type=box_type, check_file_out=True)

        complex_sys.solvate_add_ions(
            out_folder=solv_outdir,
            name=sys_name,
            ion_C=ion_conc,
            create_box_flag=False,
            maxwarn=maxwarn,
        )

        # Fix relative path returned by solvate_add_ions
        solv_gro = complex_sys.coor_file
        if not os.path.exists(solv_gro):
            solv_gro = os.path.join(solv_outdir, os.path.basename(solv_gro))
        complex_sys.coor_file = solv_gro

        solv_top = complex_sys.top_file
        if not os.path.exists(solv_top):
            solv_top = os.path.join(solv_outdir, os.path.basename(solv_top))
        complex_sys.top_file = solv_top

        ndx_file = os.path.join(solv_outdir, "index.ndx")
        generate_custom_index_groups(solv_gro, ndx_file, resname)
        complex_sys.ndx = ndx_file

        chkpt_path = save_checkpoint(complex_sys, "solvate")
        logger.info(f"[Step 6 Complete] Solvated system in {solv_outdir}")

        return {
            "complex_sys": complex_sys,
            "solv_gro": solv_gro,
            "solv_top": complex_sys.top_file,
            "ndx_file": ndx_file,
            "checkpoint_path": chkpt_path,
        }


def step7_energy_minimization(
    nt: int = 1,
    gpu: bool = False,
    em_outdir: str = "em",
    em_steps: int = 10000,
    emtol: float = 10.0,
    emstep: float = 0.01,
    maxwarn: int = 5,
    checkpoint_in: Optional[Any] = None,
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 7: Energy minimisation (two-step).
    """
    _ensure_gromacs_py()

    with working_directory(work_dir):
        if isinstance(checkpoint_in, gmx.GmxSys):
            complex_sys = checkpoint_in
        elif isinstance(checkpoint_in, str) and os.path.exists(checkpoint_in):
            with open(checkpoint_in, "rb") as fh:
                complex_sys = pickle.load(fh)
        else:
            complex_sys, _ = load_latest_checkpoint("solvate")

        relocate_checkpoint_paths(complex_sys, ".")

        complex_sys.nt = nt
        complex_sys.ntmpi = 1
        complex_sys.gpu_id = "0" if gpu else None

        logger.info(f"[Step 7] Running two-step energy minimisation with {nt} CPU threads...")
        complex_sys.em_2_steps(
            out_folder=em_outdir,
            no_constr_nsteps=em_steps,
            constr_nsteps=em_steps,
            posres="",
            create_box_flag=False,
            emtol=emtol,
            emstep=emstep,
            maxwarn=maxwarn,
        )

        # Fix relative paths after EM
        if not os.path.exists(complex_sys.coor_file):
            complex_sys.coor_file = os.path.join(em_outdir, os.path.basename(complex_sys.coor_file))
        if hasattr(complex_sys, "tpr") and complex_sys.tpr and not os.path.exists(complex_sys.tpr):
            complex_sys.tpr = os.path.join(em_outdir, os.path.basename(complex_sys.tpr))

        em_png = os.path.join(em_outdir, "em_energy.png")
        try:
            ener_1 = complex_sys.sys_history[-1].get_ener(selection_list=["Potential"])
            ener_2 = complex_sys.get_ener(selection_list=["Potential"])
            ener_1["label"] = "no bond constraints"
            ener_2["label"] = "bond constraints"
            ener_2["Time (ps)"] = ener_2["Time (ps)"].values + ener_1["Time (ps)"].max()
            ener_all = pd.concat([ener_1, ener_2])

            fig, ax = plt.subplots(figsize=(10, 4))
            for label, grp in ener_all.groupby("label"):
                ax.plot(grp["Time (ps)"], grp["Potential"], label=label)
            ax.set_xlabel("Step")
            ax.set_ylabel("Potential energy (kJ mol⁻¹)")
            ax.set_title("Energy Minimisation")
            ax.legend()
            ax.grid(alpha=0.4)
            plt.tight_layout()
            fig.savefig(em_png, dpi=150)
            plt.close(fig)
        except Exception as exc:
            logger.warning(f"Could not plot EM energy profile: {exc}")

        chkpt_path = save_checkpoint(complex_sys, "em")
        logger.info(f"[Step 7 Complete] Energy minimisation complete: {em_outdir}")

        return {
            "complex_sys": complex_sys,
            "em_dir": em_outdir,
            "em_png": em_png,
            "checkpoint_path": chkpt_path,
        }


def step8_hpc_equi_prod(
    n_cores: Optional[int] = None,
    gpu_id: str = "0",
    ligand_resname: str = "UNL",
    tc_grps_cmplx: str = "Protein",
    prod_time_ns: float = 100.0,
    equi_outdir: str = "equi",
    prod_outdir: str = "prod",
    maxwarn: int = 10,
    checkpoint_in: Optional[Any] = None,
    work_dir: str = ".",
) -> Dict[str, Any]:
    """
    Step 8: Equilibration and production MD with post-processing.
    """
    _ensure_gromacs_py()

    with working_directory(work_dir):
        if n_cores is None:
            n_cores = int(
                os.environ.get("SLURM_CPUS_PER_TASK")
                or os.environ.get("PBS_NCPUS")
                or os.environ.get("OMP_NUM_THREADS")
                or multiprocessing.cpu_count()
            )

        if isinstance(checkpoint_in, gmx.GmxSys):
            complex_sys = checkpoint_in
        elif isinstance(checkpoint_in, str) and os.path.exists(checkpoint_in):
            with open(checkpoint_in, "rb") as fh:
                complex_sys = pickle.load(fh)
        else:
            complex_sys, _ = load_latest_checkpoint("em")

        relocate_checkpoint_paths(complex_sys, ".")

        index_file = getattr(complex_sys, "ndx", None)
        if not index_file or not os.path.exists(index_file):
            index_file = "complex_solv/index.ndx"
        if os.path.exists(index_file):
            complex_sys.ndx = os.path.abspath(index_file)

        complex_sys.nt = n_cores
        complex_sys.ntmpi = 1
        complex_sys.gpu_id = str(gpu_id) if gpu_id is not None else None

        dt = 0.002
        dt_ha = 0.001
        ha_step = int(1000 * 0.5 / dt_ha)
        ca_step = int(1000 * 1.0 / dt)
        ca_low_step = int(1000 * 4.0 / dt)

        prod_steps = int(1000 * prod_time_ns / dt)
        tc_grps = f"{tc_grps_cmplx} Water_and_ions"

        logger.info(f"[Step 8] Starting equilibration stage ({n_cores} cores, GPU: {gpu_id})...")
        complex_sys.em_equi_three_step_iter_error(
            out_folder=equi_outdir,
            no_constr_nsteps=5000,
            constr_nsteps=5000,
            nsteps_HA=ha_step,
            nsteps_CA=ca_step,
            nsteps_CA_LOW=ca_low_step,
            dt=dt,
            dt_HA=dt_ha,
            tc_grps=tc_grps,
            tau_t="0.1 0.1",
            ref_t="310 310",
            vsite='none',
            maxwarn=maxwarn,
            iter_num=1,
        )

        equi_chkpt = save_checkpoint(complex_sys, "equi")

        logger.info(f"[Step 8] Starting production MD ({prod_time_ns} ns)...")
        complex_sys.production(
            out_folder=prod_outdir,
            nsteps=prod_steps,
            tc_grps=tc_grps,
            tau_t="0.1 0.1",
            ref_t="310 310",
            dt=dt,
            vsite='none',
            maxwarn=1,
            nstlist=200,
        )

        prod_chkpt = save_checkpoint(complex_sys, "prod")

        logger.info(f"[Step 8] Post-processing trajectory with gromacs_py...")
        target_group = f"Protein_{ligand_resname.strip().upper()}"

        # 1. Center trajectory and wrap molecules in box
        try:
            logger.info(f"[Step 8] [PP1] Centering trajectory on protein and wrapping molecules...")
            complex_sys.center_mol_box(traj=True)
        except Exception as exc:
            logger.warning(f"[Step 8] gromacs_py center_mol_box failed: {exc}")

        # 2. Fit rotational & translational motion and extract dry complex
        try:
            logger.info(f"[Step 8] [PP2] Fitting on Backbone and extracting {target_group}...")
            complex_sys.convert_trj(
                select=f"Backbone\n{target_group}",
                fit="rot+trans",
                pbc="none",
            )
        except Exception as exc:
            logger.warning(f"[Step 8] gromacs_py convert_trj (Backbone fit) failed: {exc}. Trying fallback fit on Protein...")
            try:
                complex_sys.convert_trj(
                    select=f"Protein\n{target_group}",
                    fit="rot+trans",
                    pbc="none",
                )
            except Exception as exc2:
                logger.warning(f"[Step 8] gromacs_py convert_trj fallback failed: {exc2}")

        # 3. Resolve and generate standard dry trajectory (traj_dry.xtc) and dry TPR (traj_dry.tpr)
        prod_tpr = getattr(complex_sys, "tpr", None)
        if not prod_tpr or not os.path.exists(prod_tpr):
            candidates = sorted(glob.glob(os.path.join(prod_outdir, "*.tpr")))
            if candidates:
                prod_tpr = candidates[-1]

        fitted_xtc = getattr(complex_sys, "xtc", None)
        if not fitted_xtc or not os.path.exists(fitted_xtc):
            candidates = (
                sorted(glob.glob(os.path.join(prod_outdir, "*_compact_compact.xtc")))
                or sorted(glob.glob(os.path.join(prod_outdir, "*_compact.xtc")))
                or sorted(glob.glob(os.path.join(prod_outdir, "*_fit.xtc")))
                or sorted(glob.glob(os.path.join(prod_outdir, "*.xtc")))
            )
            if candidates:
                fitted_xtc = candidates[-1]

        dry_xtc = os.path.join(prod_outdir, "traj_dry.xtc")
        dry_tpr = os.path.join(prod_outdir, "traj_dry.tpr")

        if prod_tpr and os.path.exists(index_file):
            candidate_groups = [target_group, "Dry", f"Protein_{ligand_resname}", "Protein_Other", "Protein_LIG"]
            
            # Extract dry TPR
            for grp in candidate_groups:
                res_tpr = run_command(
                    ["gmx", "convert-tpr", "-s", prod_tpr, "-n", index_file, "-o", dry_tpr],
                    input_str=f"{grp}\n",
                    check=False,
                )
                if res_tpr.returncode == 0:
                    logger.info(f"[Step 8] Created dry TPR ({dry_tpr}) using group '{grp}'")
                    break

            # Extract / copy dry XTC
            if fitted_xtc and os.path.exists(fitted_xtc):
                for grp in candidate_groups:
                    res_xtc = run_command(
                        ["gmx", "trjconv", "-s", prod_tpr, "-f", fitted_xtc, "-n", index_file, "-o", dry_xtc],
                        input_str=f"{grp}\n",
                        check=False,
                    )
                    if res_xtc.returncode == 0:
                        logger.info(f"[Step 8] Created dry XTC ({dry_xtc}) using group '{grp}'")
                        break

        logger.info(f"[Step 8 Complete] Production finished. Output in {prod_outdir}")
        return {
            "complex_sys": complex_sys,
            "equi_chkpt": equi_chkpt,
            "prod_chkpt": prod_chkpt,
            "dry_xtc": dry_xtc,
            "dry_tpr": dry_tpr,
        }


# =============================================================================
# CLI Entrypoint functions
# =============================================================================

def step1_cli():
    parser = argparse.ArgumentParser(description="Step 1: Extract and prepare ligand from PDBQT.")
    parser.add_argument("--input-pdbqt", "--input_pdbqt", type=str, required=True, help="Input PDBQT path")
    parser.add_argument("--ligand-name", "--ligand_name", type=str, default="UNKNOWN", help="Ligand name")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory")
    parser.add_argument("--ph", type=float, default=None, help="pH for protonation")
    parser.add_argument("--full-protonate", "--full_protonate", action="store_true", help="Fully protonate ligand")
    args = parser.parse_args()

    step1_extract_and_prepare_ligand(
        input_pdbqt=args.input_pdbqt,
        ligand_name=args.ligand_name,
        residue_name=args.residue_name,
        ph=args.ph,
        full_protonate=args.full_protonate,
        ligand_dir=args.ligand_dir,
    )


def step2_cli():
    parser = argparse.ArgumentParser(description="Step 2: Run ACPYPE topology generation.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory")
    parser.add_argument("--net-charge", "--net_charge", type=int, default=0, help="Net charge")
    parser.add_argument("--charge-method", "--charge_method", type=str, default="bcc", help="Charge method")
    args = parser.parse_args()

    step2_run_acpype(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        net_charge=args.net_charge,
        charge_method=args.charge_method,
    )


def step3_cli():
    parser = argparse.ArgumentParser(description="Step 3: In-vacuo ligand simulation test.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory")
    parser.add_argument("--invacuo-dir", "--invacuo_dir", type=str, default="invacuo", help="Output directory")
    parser.add_argument("--ntomp", type=int, default=None, help="Number of OpenMP threads")
    args = parser.parse_args()

    step3_ligand_invacuo(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        invacuo_dir=args.invacuo_dir,
        ntomp=args.ntomp,
    )


def step4_cli():
    parser = argparse.ArgumentParser(description="Step 4: Prepare protein topology and combine with ligand.")
    parser.add_argument("--protein-pdb", "--protein_pdb", type=str, required=True, help="Protein PDB path")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory")
    parser.add_argument("--prepare-top", "--prepare_top", action="store_true", help="Run pdb2gmx to prepare topology")
    parser.add_argument("--ff", type=str, default="amber99sb-ildn", help="Force field")
    args = parser.parse_args()

    step4_prepare_protein(
        protein_pdb=args.protein_pdb,
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        prepare_top=args.prepare_top,
        ff=args.ff,
    )


def step5_cli():
    parser = argparse.ArgumentParser(description="Step 5: Merge protein and ligand topologies.")
    parser.add_argument("--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("--ligand-dir", "--ligand_dir", type=str, default="ligand", help="Ligand directory")
    parser.add_argument("--protein-top", "--protein_top", type=str, default=None, help="Protein topology path")
    parser.add_argument("--protein-pdb", "--protein_pdb", type=str, default=None, help="Protein PDB path")
    args = parser.parse_args()

    step5_merge_topologies(
        residue_name=args.residue_name,
        ligand_dir=args.ligand_dir,
        protein_top=args.protein_top,
        protein_pdb=args.protein_pdb,
    )


def step6_cli():
    parser = argparse.ArgumentParser(description="Step 6: Solvate, add ions, and create custom index groups.")
    parser.add_argument("-r", "--residue-name", "--residue_name", type=str, default="UNL", help="3-letter residue name")
    parser.add_argument("-c", "--ion-conc", "--ion_conc", type=float, default=0.15, help="Ion concentration (mol/L)")
    parser.add_argument("-b", "--box-dist", "--box_dist", type=float, default=1.5, help="Distance from box edge (nm)")
    parser.add_argument("--box-type", "--box_type", type=str, default="dodecahedron", help="Box type")
    args = parser.parse_args()

    step6_solvate_ions_index(
        residue_name=args.residue_name,
        ion_conc=args.ion_conc,
        box_dist=args.box_dist,
        box_type=args.box_type,
    )


def step7_cli():
    parser = argparse.ArgumentParser(description="Step 7: Energy minimisation (two-step).")
    parser.add_argument("--nt", type=int, default=1, help="Number of CPU threads")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for energy minimisation")
    args = parser.parse_args()

    step7_energy_minimization(
        nt=args.nt,
        gpu=args.gpu,
    )


def step8_cli():
    parser = argparse.ArgumentParser(description="Step 8: Equilibration and Production MD.")
    parser.add_argument("--n-cores", "--n_cores", type=int, default=None, help="Number of CPU cores")
    parser.add_argument("--gpu-id", "--gpu_id", type=str, default="0", help="GPU ID")
    parser.add_argument("--prod-time", "--prod_time", type=float, default=100.0, help="Production time (ns)")
    parser.add_argument("--tc-grps", "--tc_grps", type=str, default="Protein", help="TC groups")
    parser.add_argument("--ligand-resname", "--ligand_resname", type=str, default="UNL", help="Ligand residue name")
    args = parser.parse_args()

    step8_hpc_equi_prod(
        n_cores=args.n_cores,
        gpu_id=args.gpu_id,
        prod_time_ns=args.prod_time,
        tc_grps_cmplx=args.tc_grps,
        ligand_resname=args.ligand_resname,
    )
