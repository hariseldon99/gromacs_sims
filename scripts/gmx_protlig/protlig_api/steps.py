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

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gromacs_py import gmx

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
)

logger = logging.getLogger("gmx_protlig.steps")


def step1_extract_and_prepare_ligand(
    input_pdbqt,
    ligand_name="UnknownLigand",
    residue_name="UNL",
    ph=None,
    full_protonate=False,
    ligand_dir="ligand",
    work_dir=".",
):
    """
    Step 1: Extract best pose from PDBQT, convert to PDB, and protonate.
    """
    with working_directory(work_dir):
        residue_name = residue_name.upper()
        if len(residue_name) != 3:
            raise ValueError(f"Residue name must be exactly 3 letters: {residue_name}")

        abs_pdbqt = os.path.abspath(input_pdbqt)
        if not os.path.exists(abs_pdbqt):
            raise FileNotFoundError(f"Input PDBQT not found: {abs_pdbqt}")

        os.makedirs(ligand_dir, exist_ok=True)
        logger.info(f"[Step 1] Reading docked poses for {ligand_name} from {abs_pdbqt}...")

        with open(abs_pdbqt, "r") as fh:
            raw = fh.read()

        blocks = {}
        current_model = None
        current_lines = []

        for line in raw.splitlines(keepends=True):
            if line.startswith("MODEL "):
                current_model = int(line.split()[1])
                current_lines = [line]
            elif line.startswith("ENDMDL") and current_model is not None:
                current_lines.append(line)
                blocks[current_model] = current_lines
                current_model = None
                current_lines = []
            elif current_model is not None:
                current_lines.append(line)

        if not blocks:
            # Single model PDBQT without MODEL tags
            blocks[1] = raw.splitlines(keepends=True)

        scores = {}
        for mid, lines in blocks.items():
            for l in lines:
                if "VINA RESULT:" in l:
                    scores[mid] = float(l.split()[3])
                    break

        if scores:
            best_model = min(scores, key=scores.get)
            best_score = scores[best_model]
        else:
            best_model = 1
            best_score = 0.0

        best_pdbqt = os.path.join(ligand_dir, f"{ligand_name.lower()}_pose1.pdbqt")
        with open(best_pdbqt, "w") as fh:
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

        with open(prot_pdb, "r") as fh:
            content = fh.read()
        content = content.replace(" UNL ", f" {residue_name} ")
        with open(prot_pdb, "w") as fh:
            fh.write(content)

        logger.info(f"[Step 1 Complete] Prepared ligand PDB: {prot_pdb}")
        return {
            "ligand_dir": ligand_dir,
            "best_pdbqt": best_pdbqt,
            "raw_pdb": raw_pdb,
            "protonated_pdb": prot_pdb,
            "best_model": best_model,
            "best_score": best_score,
        }


def step2_run_acpype(
    residue_name="ERG",
    ligand_dir="ligand",
    net_charge=0,
    charge_method="bcc",
    work_dir=".",
):
    """
    Step 2: Run acpype on protonated ligand PDB.
    """
    with working_directory(work_dir):
        resname = residue_name.upper()
        ligand_pdb = os.path.join(ligand_dir, f"{resname}.pdb")
        if not os.path.exists(ligand_pdb):
            raise FileNotFoundError(f"Ligand PDB not found: {ligand_pdb}. Run Step 1 first.")

        logger.info(f"[Step 2] Running acpype on {ligand_pdb}...")

        # Run acpype directly inside ligand_dir so output stays contained
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
    residue_name="ERG",
    ligand_dir="ligand",
    invacuo_dir="invacuo",
    ntomp=None,
    work_dir=".",
):
    """
    Step 3: In-vacuo ligand simulation (energy minimisation + short MD).
    """
    with working_directory(work_dir):
        resname = residue_name.upper()
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

        with open(os.path.join(invacuo_dir, "em_vac.mdp"), "w") as fh:
            fh.write(em_mdp_content)
        with open(os.path.join(invacuo_dir, "md_vac.mdp"), "w") as fh:
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
    protein_pdb,
    residue_name="UNL",
    ligand_dir="ligand",
    prepare_top=False,
    ff="amber99sb-ildn",
    protein_outdir="protein",
    complex_outdir="complex_raw",
    work_dir=".",
):
    """
    Step 4: Prepare protein topology and combine protein + ligand coordinates.
    """
    with working_directory(work_dir):
        resname = residue_name.upper()
        abs_protein_pdb = os.path.abspath(protein_pdb)
        protein_name = os.path.basename(abs_protein_pdb).split(".")[0]

        ligand_gro = os.path.join(ligand_dir, f"{resname}.acpype", f"{resname}_GMX.gro")
        if not os.path.exists(abs_protein_pdb):
            raise FileNotFoundError(f"Protein PDB not found: {abs_protein_pdb}")
        if not os.path.exists(ligand_gro):
            raise FileNotFoundError(f"Ligand GRO missing: {ligand_gro}")

        os.makedirs(complex_outdir, exist_ok=True)
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

        prot_gro_path = os.path.join(protein_outdir, f"{protein_name}_pdb2gmx.gro")
        run_command(["gmx", "editconf", "-f", protein_sys.coor_file, "-o", prot_gro_path])

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
    residue_name="ERG",
    ligand_dir="ligand",
    protein_top="protein/topol.top",
    protein_pdb=None,
    complex_dir="complex",
    complex_gro_src="complex_raw/complex_raw.gro",
    work_dir=".",
):
    """
    Step 5: Merge protein and ligand topologies.
    """
    with working_directory(work_dir):
        resname = residue_name.upper()
        acpype_dir = os.path.join(ligand_dir, f"{resname}.acpype")

        if not os.path.exists(acpype_dir) or not os.path.exists(complex_gro_src):
            raise FileNotFoundError("Required inputs for topology merge not found.")

        os.makedirs(complex_dir, exist_ok=True)

        prot_top_path = protein_top if os.path.exists(protein_top) else None
        if prot_top_path is None:
            for root, _, files in os.walk(os.path.dirname(protein_top) or "protein"):
                for f in files:
                    if f.endswith(".top"):
                        prot_top_path = os.path.join(root, f)
                        break

        if not prot_top_path or not os.path.exists(prot_top_path):
            raise FileNotFoundError(f"Cannot locate protein topology: {protein_top}")

        lig_top = os.path.join(acpype_dir, f"{resname}_GMX.top")
        lig_itp = os.path.join(acpype_dir, f"{resname}_GMX.itp")

        with open(prot_top_path, "r") as fh:
            prot_text = fh.read()
        with open(lig_top, "r") as fh:
            lig_top_text = fh.read()
        with open(lig_itp, "r") as fh:
            lig_itp_text = fh.read()

        lig_atomtypes = extract_atomtypes_block(lig_top_text) or extract_atomtypes_block(lig_itp_text)
        lig_nonbonded = extract_nonbonded_params_block(lig_top_text) or extract_nonbonded_params_block(lig_itp_text)

        dst_itp = os.path.join(complex_dir, f"{resname}_GMX.itp")
        dst_posre = os.path.join(complex_dir, f"posre_{resname}.itp")

        clean_itp = strip_global_directives(lig_itp_text)
        with open(dst_itp, "w") as fh:
            fh.write(clean_itp)

        src_posre = os.path.join(acpype_dir, f"posre_{resname}.itp")
        if os.path.exists(src_posre):
            shutil.copy2(src_posre, dst_posre)
        else:
            with open(dst_posre, "w") as fh:
                fh.write(f"; Position restraints for {resname}\n[ position_restraints ]\n")

        for f in os.listdir(acpype_dir):
            if f.endswith(".itp") and f not in (f"{resname}_GMX.itp", f"posre_{resname}.itp"):
                shutil.copy2(os.path.join(acpype_dir, f), os.path.join(complex_dir, f))

        prot_text_rw = rewrite_include_paths(prot_text, prot_top_path, complex_dir)

        import re
        ff_include_re = re.compile(r'(#include\s+"[^"]*(?:forcefield|ff)[^"]*\.itp"[^\n]*\n)', re.IGNORECASE)
        m = ff_include_re.search(prot_text_rw)
        if m:
            insert_after = m.end()
            block = f"\n; ── Ligand {resname} atomtypes ──────────────────────────\n" + lig_atomtypes
            if lig_nonbonded:
                block += "\n" + lig_nonbonded
            block += "; ─────────────────────────────────────────────────────────────\n\n"
            prot_text_rw = prot_text_rw[:insert_after] + block + prot_text_rw[insert_after:]

        water_include_re = re.compile(r'(#include\s+"[^"]*(?:tip3p|tip4p|spc|water)[^"]*\.itp"[^\n]*\n)', re.IGNORECASE)
        lig_inc = f'#include "{resname}_GMX.itp"\n'
        m2 = water_include_re.search(prot_text_rw)
        if m2:
            prot_text_rw = prot_text_rw[:m2.start()] + f"; ── Ligand {resname} molecule type ───\n" + lig_inc + "\n" + prot_text_rw[m2.start():]

        mol_re = re.compile(r'(Protein\S*\s+\d+\s*\n)', re.IGNORECASE)
        matches = list(mol_re.finditer(prot_text_rw))
        if matches:
            last_match = matches[-1]
            prot_text_rw = prot_text_rw[:last_match.end()] + f"{resname:3s}            1\n" + prot_text_rw[last_match.end():]

        complex_top = os.path.join(complex_dir, "complex.top")
        with open(complex_top, "w") as fh:
            fh.write(prot_text_rw)

        complex_gro = os.path.join(complex_dir, "complex.gro")
        shutil.copy2(complex_gro_src, complex_gro)

        protein_name = os.path.basename(protein_pdb).split(".")[0] if protein_pdb else "protein"
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
    residue_name="UNL",
    ion_conc=0.15,
    box_dist=1.5,
    box_type="dodecahedron",
    solv_outdir="complex_solv",
    maxwarn=5,
    checkpoint_in=None,
    work_dir=".",
):
    """
    Step 6: Box creation, solvation, ion addition, and index group generation.
    """
    with working_directory(work_dir):
        resname = residue_name.upper()

        if isinstance(checkpoint_in, gmx.GmxSys):
            complex_sys = checkpoint_in
        elif isinstance(checkpoint_in, str) and os.path.exists(checkpoint_in):
            with open(checkpoint_in, "rb") as fh:
                complex_sys = pickle.load(fh)
        else:
            complex_sys, chkpt_file = load_latest_checkpoint("merge_topologies")

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

        solv_gro = complex_sys.coor_file
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
    nt=1,
    gpu=False,
    em_outdir="em",
    em_steps=10000,
    emtol=10.0,
    emstep=0.01,
    maxwarn=5,
    checkpoint_in=None,
    work_dir=".",
):
    """
    Step 7: Energy minimisation (two-step).
    """
    with working_directory(work_dir):
        if isinstance(checkpoint_in, gmx.GmxSys):
            complex_sys = checkpoint_in
        elif isinstance(checkpoint_in, str) and os.path.exists(checkpoint_in):
            with open(checkpoint_in, "rb") as fh:
                complex_sys = pickle.load(fh)
        else:
            complex_sys, _ = load_latest_checkpoint("solvate")

        complex_sys.nt = nt
        complex_sys.ntmpi = 1
        complex_sys.gpu_id = "0" if gpu else None

        logger.info(f"[Step 7] Running two-step energy minimisation with {nt} CPU threads...")
        complex_sys.em_2_steps(
            out_folder=em_outdir,
            no_constr_nsteps=em_steps,
            constr_nsteps=em_steps,
            posres=None,
            create_box_flag=False,
            emtol=emtol,
            emstep=emstep,
            maxwarn=maxwarn,
        )

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
    n_cores=None,
    gpu_id="0",
    tc_grps_cmplx="Protein",
    prod_time_ns=100.0,
    equi_outdir="equi",
    prod_outdir="prod",
    maxwarn=10,
    checkpoint_in=None,
    work_dir=".",
):
    """
    Step 8: Equilibration and production MD with post-processing.
    """
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

        index_file = "complex_solv/index.ndx"
        if os.path.exists(index_file):
            complex_sys.ndx = index_file

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

        logger.info(f"[Step 8] Starting equilibration stage...")
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
            vsite=None,
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
            vsite=None,
            maxwarn=1,
            nstlist=200,
        )

        prod_chkpt = save_checkpoint(complex_sys, "prod")

        logger.info(f"[Step 8] Post-processing trajectory...")
        complex_sys.center_mol_box(traj=True)
        complex_sys.convert_trj(
            select=f"{tc_grps_cmplx}\nSystem",
            fit="rot+trans",
            pbc="none",
        )

        prod_tpr = getattr(complex_sys, "tpr", None)
        if not prod_tpr:
            candidates = sorted(glob.glob(os.path.join(prod_outdir, "*.tpr")))
            if candidates:
                prod_tpr = candidates[-1]

        fitted_xtc = getattr(complex_sys, "xtc", None)
        if not fitted_xtc:
            candidates = sorted(glob.glob(os.path.join(prod_outdir, "*_fit.xtc"))) or sorted(glob.glob(os.path.join(prod_outdir, "*.xtc")))
            if candidates:
                fitted_xtc = candidates[-1]

        dry_xtc = os.path.join(prod_outdir, "traj_dry.xtc")
        dry_tpr = os.path.join(prod_outdir, "traj_dry.tpr")

        if prod_tpr and fitted_xtc and os.path.exists(index_file):
            run_command(
                ["gmx", "trjconv", "-s", prod_tpr, "-f", fitted_xtc, "-n", index_file, "-o", dry_xtc],
                input_str="Dry\n",
                check=False
            )
            run_command(
                ["gmx", "convert-tpr", "-s", prod_tpr, "-n", index_file, "-o", dry_tpr],
                input_str="Dry\n",
                check=False
            )

        logger.info(f"[Step 8 Complete] Production finished. Output in {prod_outdir}")
        return {
            "complex_sys": complex_sys,
            "equi_chkpt": equi_chkpt,
            "prod_chkpt": prod_chkpt,
            "dry_xtc": dry_xtc,
            "dry_tpr": dry_tpr,
        }
