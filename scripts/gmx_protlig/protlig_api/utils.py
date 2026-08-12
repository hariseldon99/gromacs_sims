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

logger = logging.getLogger("gmx_protlig")


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


def run_command(cmd, cwd=None, check=True, capture_output=True, input_str=None):
    """
    Run a shell/system command cleanly via subprocess.
    """
    logger.info(f"Running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    res = subprocess.run(
        cmd,
        cwd=cwd,
        check=check,
        capture_output=capture_output,
        text=True,
        input=input_str
    )
    return res


def read_gro(path):
    """Parse a GRO file into (title, atom_lines, box_line)."""
    with open(path, "r") as fh:
        lines = fh.readlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file: {path} (fewer than 3 lines)")
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms]
    return title, atom_lines, box_line


def write_combined_gro(prot_gro, lig_gro, out_gro, title="Complex"):
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

    os.makedirs(os.path.dirname(os.path.abspath(out_gro)), exist_ok=True)
    with open(out_gro, "w") as fh:
        fh.write(f"{title}\n")
        fh.write(f"{total:5d}\n")
        fh.writelines(combined_atoms)
        fh.write(prot_box)

    return total


SECTION_RE = re.compile(r"^\s*\[\s*(\S+)\s*\]")


def extract_atomtypes_block(top_text):
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


def extract_nonbonded_params_block(top_text):
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


def strip_global_directives(itp_text):
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


def rewrite_include_paths(top_text, from_dir, to_dir):
    """Rewrite relative #include paths while ignoring system forcefields."""
    def fix_path(match):
        orig_path = match.group(1)
        if os.path.isabs(orig_path) or ".ff/" in orig_path or orig_path.endswith(".ff"):
            return match.group(0)
        abs_path = os.path.normpath(os.path.join(os.path.dirname(from_dir), orig_path))
        try:
            rel = os.path.relpath(abs_path, to_dir)
        except ValueError:
            rel = abs_path
        return f'#include "{rel}"'

    return re.sub(r'#include\s+"([^"]+)"', fix_path, top_text)


def generate_custom_index_groups(gro_file, ndx_file, res_name):
    """
    Parses .gro file and appends required custom groups to index.ndx:
      - Protein_<RESNAME>
      - Dry
      - Water_and_ions
    """
    # 1. Generate default GROMACS index file first
    run_command(["gmx", "make_ndx", "-f", gro_file, "-o", ndx_file], input_str="q\n")

    # 2. Parse GRO file structure line-by-line
    water_ion_resnames = {"SOL", "NA", "CL", "NA+", "CL-", "WAT", "HOH", "ION"}
    
    protein_atoms = []
    ligand_atoms = []
    water_ion_atoms = []

    with open(gro_file, "r") as fh:
        lines = fh.readlines()

    atom_lines = lines[2:-1]

    for idx, line in enumerate(atom_lines, start=1):
        if len(line) < 10:
            continue
        res_name_in_gro = line[5:10].strip().upper()

        if res_name_in_gro == res_name.upper():
            ligand_atoms.append(idx)
        elif res_name_in_gro in water_ion_resnames:
            water_ion_atoms.append(idx)
        else:
            protein_atoms.append(idx)

    protein_ligand_atoms = protein_atoms + ligand_atoms

    def format_group(group_name, atom_indices):
        out = [f"\n[ {group_name} ]\n"]
        for i in range(0, len(atom_indices), 15):
            chunk = atom_indices[i : i + 15]
            out.append(" ".join(f"{atom:6d}" for atom in chunk) + "\n")
        return "".join(out)

    with open(ndx_file, "a") as fh:
        fh.write(format_group(f"Protein_{res_name.upper()}", protein_ligand_atoms))
        fh.write(format_group("Dry", protein_ligand_atoms))
        fh.write(format_group("Water_and_ions", water_ion_atoms))


def save_checkpoint(obj, prefix, search_dir="."):
    """Save checkpoint file with date stamp."""
    chkpt_name = f"checkpoint.{prefix}_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
    chkpt_path = os.path.join(search_dir, chkpt_name)
    with open(chkpt_path, "wb") as fh:
        pickle.dump(obj, fh)
    return chkpt_path


def load_latest_checkpoint(prefix, search_dir="."):
    """Load latest checkpoint matching pattern."""
    pattern = os.path.join(search_dir, f"checkpoint.{prefix}_*.pycpt")
    chkpts = sorted(glob.glob(pattern))
    if not chkpts:
        raise FileNotFoundError(f"No checkpoint found matching {pattern}")
    latest = chkpts[-1]
    with open(latest, "rb") as fh:
        obj = pickle.load(fh)
    return obj, latest
