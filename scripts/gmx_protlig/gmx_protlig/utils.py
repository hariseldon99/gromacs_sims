"""
Utility functions and context managers for gmx_protlig workflow.
"""

import os
import re
import glob
import pickle
import datetime
import subprocess
import contextlib
import logging
from typing import Dict, List, Tuple, Any, Optional

logger = logging.getLogger("gmx_protlig.utils")


@contextlib.contextmanager
def working_directory(path):
    """
    Context manager to temporarily change working directory.
    Creates directory if it does not exist.
    """
    prev_cwd = os.getcwd()
    target_path = os.path.abspath(path)
    os.makedirs(target_path, exist_ok=True)
    os.chdir(target_path)
    try:
        yield target_path
    finally:
        os.chdir(prev_cwd)


def run_command(
    cmd: Any,
    cwd: Optional[str] = None,
    check: bool = True,
    capture_output: bool = True,
    input_str: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """
    Run a shell/system command cleanly via subprocess with rich error reporting.
    """
    cmd_display = " ".join(cmd) if isinstance(cmd, list) else str(cmd)
    logger.info(f"Running command: {cmd_display}")
    try:
        res = subprocess.run(
            cmd,
            cwd=cwd,
            check=False,
            capture_output=capture_output,
            text=True,
            input=input_str,
        )
        if check and res.returncode != 0:
            err_msg = (
                f"Command failed with exit code {res.returncode}: {cmd_display}\n"
                f"STDOUT:\n{res.stdout}\n"
                f"STDERR:\n{res.stderr}"
            )
            logger.error(err_msg)
            raise subprocess.CalledProcessError(
                res.returncode, cmd, output=res.stdout, stderr=res.stderr
            )
        return res
    except Exception as e:
        if not isinstance(e, subprocess.CalledProcessError):
            logger.error(f"Error executing command '{cmd_display}': {e}")
        raise


def read_gro(path: str) -> Tuple[str, List[str], str]:
    """Parse a GRO file into (title, atom_lines, box_line)."""
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file: {path} (fewer than 3 lines)")
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms] if len(lines) > 2 + n_atoms else lines[-1]
    return title, atom_lines, box_line


def write_combined_gro(
    prot_gro: str,
    lig_gro: str,
    out_gro: str,
    title: str = "Complex",
) -> int:
    """Write a combined GRO of protein + ligand with atom renumbering."""
    _, prot_atoms, prot_box = read_gro(prot_gro)
    _, lig_atoms, _ = read_gro(lig_gro)

    total = len(prot_atoms) + len(lig_atoms)

    combined_atoms = []
    atom_idx = 1
    for line in prot_atoms + lig_atoms:
        new_line = line[:15] + f"{atom_idx % 100000:5d}" + line[20:]
        combined_atoms.append(new_line)
        atom_idx += 1

    out_dir = os.path.dirname(os.path.abspath(out_gro))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_gro, "w", encoding="utf-8") as fh:
        fh.write(f"{title}\n")
        fh.write(f"{total:5d}\n")
        fh.writelines(combined_atoms)
        fh.write(prot_box)

    return total


SECTION_RE = re.compile(r"^\s*\[\s*(\S+)\s*\]")


def extract_atomtypes_block(top_text: str) -> str:
    """Extract the full [ atomtypes ] block text from topology/ITP text."""
    lines = top_text.splitlines(keepends=True)
    in_atomtypes = False
    block_lines = []
    for line in lines:
        m = SECTION_RE.match(line)
        if m:
            if m.group(1).lower() == "atomtypes":
                in_atomtypes = True
                block_lines = [line]
            elif in_atomtypes:
                break
        elif in_atomtypes:
            block_lines.append(line)
    return "".join(block_lines)


def extract_nonbonded_params_block(top_text: str) -> str:
    """Extract [ nonbond_params ] if present."""
    lines = top_text.splitlines(keepends=True)
    in_block = False
    block_lines = []
    for line in lines:
        m = SECTION_RE.match(line)
        if m:
            if m.group(1).lower() == "nonbond_params":
                in_block = True
                block_lines = [line]
            elif in_block:
                break
        elif in_block:
            block_lines.append(line)
    return "".join(block_lines)


def strip_global_directives(itp_text: str) -> str:
    """Remove [ atomtypes ] and [ nonbond_params ] from an ITP file."""
    lines = itp_text.splitlines(keepends=True)
    out_lines = []
    skip_section = False
    for line in lines:
        m = SECTION_RE.match(line)
        if m:
            sec_name = m.group(1).lower()
            if sec_name in ("atomtypes", "nonbond_params"):
                skip_section = True
            else:
                skip_section = False

        if not skip_section:
            out_lines.append(line)
    return "".join(out_lines)


def rewrite_include_paths(top_text: str, from_path: str, to_dir: str) -> str:
    """
    Rewrite relative #include paths while ignoring system forcefields.
    Handles both file paths and directory paths for from_path.
    """
    base_from = from_path if os.path.isdir(from_path) else os.path.dirname(from_path)
    if not base_from:
        base_from = "."

    def fix_path(match):
        orig_path = match.group(1)
        if os.path.isabs(orig_path) or ".ff/" in orig_path or orig_path.endswith(".ff"):
            return match.group(0)
        abs_path = os.path.normpath(os.path.join(base_from, orig_path))
        try:
            rel = os.path.relpath(abs_path, to_dir)
        except ValueError:
            rel = abs_path
        return f'#include "{rel}"'

    return re.sub(r'#include\s+"([^"]+)"', fix_path, top_text)


def generate_custom_index_groups(gro_file: str, ndx_file: str, res_name: str):
    """
    Parses .gro file and appends required custom groups to index.ndx:
      - Protein_<RESNAME>
      - Dry
      - Water_and_ions
    """
    if not os.path.exists(gro_file):
        raise FileNotFoundError(f"Cannot generate index: GRO file not found at '{gro_file}'")

    ndx_dir = os.path.dirname(os.path.abspath(ndx_file))
    if ndx_dir:
        os.makedirs(ndx_dir, exist_ok=True)

    # 1. Generate default GROMACS index file first
    run_command(["gmx", "make_ndx", "-f", gro_file, "-o", ndx_file], input_str="q\n")

    # 2. Comprehensive solvent and ion residue recognition
    water_ion_resnames = {
        "SOL", "WAT", "HOH", "TIP3", "TIP4", "SPC", "SPCE", "TIP",
        "NA", "CL", "K", "MG", "CA", "ZN", "SOD", "CLA", "POT", "CAL",
        "MG2", "ZN2", "NA+", "CL-", "K+", "MG2+", "CA2+", "ZN2+",
        "ION", "IB+", "CL-", "CS+", "LI+", "RB+"
    }

    protein_atoms = []
    ligand_atoms = []
    water_ion_atoms = []

    with open(gro_file, "r", encoding="utf-8") as fh:
        lines = fh.readlines()

    atom_lines = lines[2:-1]

    target_res = res_name.strip().upper()
    for idx, line in enumerate(atom_lines, start=1):
        if len(line) < 10:
            continue
        res_name_in_gro = line[5:10].strip().upper()

        if res_name_in_gro == target_res:
            ligand_atoms.append(idx)
        elif res_name_in_gro in water_ion_resnames:
            water_ion_atoms.append(idx)
        else:
            protein_atoms.append(idx)

    protein_ligand_atoms = protein_atoms + ligand_atoms

    def format_group(group_name: str, atom_indices: List[int]) -> str:
        out = [f"\n[ {group_name} ]\n"]
        for i in range(0, len(atom_indices), 15):
            chunk = atom_indices[i : i + 15]
            out.append(" ".join(f"{atom:6d}" for atom in chunk) + "\n")
        return "".join(out)

    with open(ndx_file, "a", encoding="utf-8") as fh:
        fh.write(format_group(f"Protein_{target_res}", protein_ligand_atoms))
        fh.write(format_group("Dry", protein_ligand_atoms))
        fh.write(format_group("Water_and_ions", water_ion_atoms))

    logger.info(
        f"Generated custom index groups in {ndx_file}: "
        f"Protein={len(protein_atoms)}, Ligand={len(ligand_atoms)}, "
        f"Water_and_ions={len(water_ion_atoms)}"
    )


def save_checkpoint(obj: Any, prefix: str, search_dir: str = ".") -> str:
    """
    Save checkpoint file with date and timestamp.
    """
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    chkpt_name = f"checkpoint.{prefix}_{now}.pycpt"
    chkpt_path = os.path.join(search_dir, chkpt_name)
    with open(chkpt_path, "wb") as fh:
        pickle.dump(obj, fh)
    logger.info(f"Checkpoint saved: {chkpt_path}")
    return chkpt_path


def load_latest_checkpoint(prefix: str, search_dir: str = ".") -> Tuple[Any, str]:
    """
    Load latest checkpoint matching pattern, sorted chronologically.
    """
    pattern = os.path.join(search_dir, f"checkpoint.{prefix}_*.pycpt")
    chkpts = sorted(glob.glob(pattern))
    if not chkpts:
        raise FileNotFoundError(f"No checkpoint found matching {pattern}")
    latest = chkpts[-1]
    with open(latest, "rb") as fh:
        obj = pickle.load(fh)
    logger.info(f"Loaded checkpoint from: {latest}")
    return obj, latest


def relocate_checkpoint_paths(complex_sys: Any, new_base_dir: str = ".") -> Any:
    """
    Inspects and updates file paths stored inside a GmxSys object when the simulation
    directory has been moved across filesystems, mounts, or containers.
    """
    abs_base = os.path.abspath(new_base_dir)
    attrs = list(getattr(complex_sys, "__slots__", None) or vars(complex_sys).keys())

    for attr in attrs:
        val = getattr(complex_sys, attr, None)
        if not val or not isinstance(val, str):
            continue

        # If it looks like a path and the original path does not exist
        if ("/" in val or "\\" in val) and not os.path.exists(val):
            basename = os.path.basename(val)
            # Try searching under new_base_dir
            for root, _, files in os.walk(abs_base):
                if basename in files:
                    new_path = os.path.join(root, basename)
                    setattr(complex_sys, attr, new_path)
                    logger.info(f"Relocated {attr}: '{val}' -> '{new_path}'")
                    break

    return complex_sys
