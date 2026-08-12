#!/usr/bin/env python3
"""
Step 5 — Merge the protein topology (from pdb2gmx) with the ligand topology
(from acpype) into a single complex topology.
"""

import os
import sys
import re
import argparse
import shutil
import pickle
import datetime

from gromacs_py import gmx

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Merge protein and ligand topologies."
)
parser.add_argument(
    "--residue-name",
    type=str,
    default="ERG",
    help="3-letter GROMACS residue name for the ligand (default: ERG)",
)
parser.add_argument(
    "--ligand-dir",
    type=str,
    default="ligand",
    help="Ligand directory (default: ligand)",
)

parser.add_argument(
    "--protein-top",
    type=str,
    default="protein/topol.top",
    help="Path to protein topology (default: protein/topol.top)",
)

parser.add_argument(
    "--protein-pdb",
    type=str,
    help="Path to protein PDB file",
)
args = parser.parse_args()

RESIDUE_NAME = args.residue_name.upper()
LIGAND_DIR = args.ligand_dir
PROTEIN_TOP = args.protein_top

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROTEIN_TOPDIR  = os.path.dirname(PROTEIN_TOP)
ACPYPE_DIR      = os.path.join(LIGAND_DIR, f"{RESIDUE_NAME}.acpype")
COMPLEX_GRO_SRC = "complex_raw/complex_raw.gro"
COMPLEX_DIR     = "complex"

for p in (ACPYPE_DIR, COMPLEX_GRO_SRC):
    if not os.path.exists(p):
        print(f"[ERROR] Required path not found: {p}")
        sys.exit(1)

os.makedirs(COMPLEX_DIR, exist_ok=True)

# Locate protein topology
prot_top = PROTEIN_TOP if os.path.exists(PROTEIN_TOP) else None
if prot_top is None:
    for root, _, files in os.walk(PROTEIN_TOPDIR):
        for f in files:
            if f.endswith(".top"):
                prot_top = os.path.join(root, f)
                break
        if prot_top:
            break
if prot_top is None:
    print(f"[ERROR] Cannot find protein topology in {PROTEIN_TOPDIR}/")
    sys.exit(1)
print(f"[INFO] Using protein topology: {prot_top}")

lig_top = os.path.join(ACPYPE_DIR, f"{RESIDUE_NAME}_GMX.top")
lig_itp = os.path.join(ACPYPE_DIR, f"{RESIDUE_NAME}_GMX.itp")

for p in (lig_top, lig_itp):
    if not os.path.exists(p):
        print(f"[ERROR] Required file not found: {p}")
        sys.exit(1)

# ---------------------------------------------------------------------------
# Helpers: Section parsing and cleaning
# ---------------------------------------------------------------------------
SECTION_RE = re.compile(r"^\s*\[\s*(\S+)\s*\]")


def extract_atomtypes_block(top_text):
    """Extract the full [ atomtypes ] block text from a topology/ITP string."""
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
    """
    Remove [ atomtypes ] and [ nonbond_params ] from an ITP file so it only
    contains molecule-level directives (preventing GROMACS syntax errors).
    """
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


# ---------------------------------------------------------------------------
# 1. Read topology files
# ---------------------------------------------------------------------------
with open(prot_top) as fh:
    prot_text = fh.read()

with open(lig_top) as fh:
    lig_top_text = fh.read()

with open(lig_itp) as fh:
    lig_itp_text = fh.read()

# ---------------------------------------------------------------------------
# 2. Extract [ atomtypes ] from acpype (checking both TOP and ITP)
# ---------------------------------------------------------------------------
lig_atomtypes = extract_atomtypes_block(lig_top_text) or extract_atomtypes_block(lig_itp_text)
lig_nonbonded = extract_nonbonded_params_block(lig_top_text) or extract_nonbonded_params_block(lig_itp_text)

if not lig_atomtypes:
    print("[WARNING] Could not extract [ atomtypes ] from acpype topology.")
    print("          Proceeding without it — you may need to merge manually.")

# ---------------------------------------------------------------------------
# 3. Write cleaned ligand ITP into complex/ directory
# ---------------------------------------------------------------------------
dst_itp   = os.path.join(COMPLEX_DIR, f"{RESIDUE_NAME}_GMX.itp")
dst_posre = os.path.join(COMPLEX_DIR, f"posre_{RESIDUE_NAME}.itp")

# Clean the ITP by removing global directives [ atomtypes ]
clean_itp_text = strip_global_directives(lig_itp_text)
with open(dst_itp, "w") as fh:
    fh.write(clean_itp_text)
print(f"[OK] Cleaned and copied: {lig_itp} → {dst_itp}")

src_posre = os.path.join(ACPYPE_DIR, f"posre_{RESIDUE_NAME}.itp")
if os.path.exists(src_posre):
    shutil.copy2(src_posre, dst_posre)
    print(f"[OK] Copied: {src_posre} → {dst_posre}")
else:
    with open(dst_posre, "w") as fh:
        fh.write(f"; Position restraints for {RESIDUE_NAME} (auto-generated)\n")
        fh.write("[ position_restraints ]\n")
        fh.write("; atom  type      fx      fy      fz\n")
    print(f"[INFO] posre_{RESIDUE_NAME}.itp not found; placeholder created.")

# Copy extra ITP files
for f in os.listdir(ACPYPE_DIR):
    if f.endswith(".itp") and f not in (f"{RESIDUE_NAME}_GMX.itp", f"posre_{RESIDUE_NAME}.itp"):
        src = os.path.join(ACPYPE_DIR, f)
        dst = os.path.join(COMPLEX_DIR, f)
        shutil.copy2(src, dst)
        print(f"[OK] Copied extra ITP: {f}")

# ---------------------------------------------------------------------------
# 5. Build merged topology
# ---------------------------------------------------------------------------

def rewrite_include_paths(top_text, from_dir, to_dir):
    """Rewrite relative #include paths while ignoring system forcefields (.ff)."""
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


prot_text_rw = rewrite_include_paths(prot_text, prot_top, COMPLEX_DIR)

# Insert [ atomtypes ] after forcefield include
ff_include_re = re.compile(
    r'(#include\s+"[^"]*(?:forcefield|ff)[^"]*\.itp"[^\n]*\n)',
    re.IGNORECASE
)
m = ff_include_re.search(prot_text_rw)
if m:
    insert_after = m.end()
    atomtypes_block = (
        f"\n; ── Ligand {RESIDUE_NAME} atomtypes (from acpype) ──────────────────────────\n"
        + lig_atomtypes
    )
    if lig_nonbonded:
        atomtypes_block += "\n" + lig_nonbonded
    atomtypes_block += "; ─────────────────────────────────────────────────────────────\n\n"
    prot_text_rw = (
        prot_text_rw[:insert_after]
        + atomtypes_block
        + prot_text_rw[insert_after:]
    )
    print("[OK] Inserted [ atomtypes ] after force-field include.")
else:
    prot_text_rw = (
        f"; ── Ligand {RESIDUE_NAME} atomtypes (from acpype) ──────────────────────────\n"
        + lig_atomtypes
        + "\n"
        + prot_text_rw
    )

# Insert ligand ITP include before water topology
water_include_re = re.compile(
    r'(#include\s+"[^"]*(?:tip3p|tip4p|spc|water)[^"]*\.itp"[^\n]*\n)',
    re.IGNORECASE
)
lig_include_str = f'#include "{RESIDUE_NAME}_GMX.itp"\n'
m2 = water_include_re.search(prot_text_rw)
if m2:
    prot_text_rw = (
        prot_text_rw[:m2.start()]
        + f"; ── Ligand {RESIDUE_NAME} molecule type ───────────────────────────────────\n"
        + lig_include_str
        + "\n"
        + prot_text_rw[m2.start():]
    )
    print(f"[OK] Inserted {RESIDUE_NAME}_GMX.itp include before water topology.")

# Add ligand to [ molecules ]
mol_re = re.compile(r'(Protein\S*\s+\d+\s*\n)', re.IGNORECASE)
matches = list(mol_re.finditer(prot_text_rw))
if matches:
    last_match = matches[-1]
    prot_text_rw = (
        prot_text_rw[:last_match.end()]
        + f"{RESIDUE_NAME:3s}            1\n"
        + prot_text_rw[last_match.end():]
    )
    print(f"[OK] Added '{RESIDUE_NAME} 1' to [ molecules ].")

# ---------------------------------------------------------------------------
# Write Output and Checkpoint
# ---------------------------------------------------------------------------
complex_top = os.path.join(COMPLEX_DIR, "complex.top")
with open(complex_top, "w") as fh:
    fh.write(prot_text_rw)
print(f"\n[OK] Merged topology written: {complex_top}")

complex_gro = os.path.join(COMPLEX_DIR, "complex.gro")
shutil.copy2(COMPLEX_GRO_SRC, complex_gro)
print(f"[OK] Complex GRO copied: {complex_gro}")

protein_name = os.path.basename(args.protein_pdb).split(".")[0] if args.protein_pdb else "protein"
complex_name = f"complex_{protein_name.lower()}_{RESIDUE_NAME.lower()}"
complex_sys = gmx.GmxSys(name=complex_name, coor_file=complex_gro)
complex_sys.top_file = complex_top

chkpt = f"checkpoint.merge_topologies_{datetime.date.today().strftime('%Y%m%d')}.pycpt"
with open(chkpt, "wb") as fh:
    pickle.dump(complex_sys, fh)
print(f"\n[Checkpoint] Saved: {chkpt}")
