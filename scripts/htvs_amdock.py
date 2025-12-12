#!/usr/bin/env python3
"""
This script integrates the comprehensive AMDock workflow for High-Throughput Virtual Screening (HTVS), including dynamic binding site detection using **AutoLigand in Automatic Mode**, multi-processing control, and optimized PDB coordinate parsing designed to handle the fixed-width fields characteristic of PDB/FILL files, as seen in the source logs.

To address your query regarding using Bioconda for better PDB parsing, while the sources confirm AMDock and its dependencies (`pdb2pqr`, `openbabel`, `pymol`) are typically installed via **Conda**, the actual coordinate parsing of the generated `FILL` files is best handled using standard Python libraries optimized for array operations (`numpy`) to accurately read the fixed columns shown in the logs.
"""
### AMDock HTVS Replication Script with Dynamic Site Detection

import subprocess
import os
import sys
import json
import multiprocessing
import glob
import re
import shutil
import platform
import openbabel as ob
from openbabel import pybel  # common import name for OpenBabel's Python wrapper
from functools import partial
from tabulate import tabulate
import json
from collections import defaultdict


# Attempt to import numpy for robust numerical operations (recommended for Conda environment)
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    print("Warning: NumPy not found. Coordinate calculations will use slower Python built-ins.")
    NUMPY_AVAILABLE = False


# --- Configuration Constants (Derived from AMDock Logs) ---
FORCE_FIELD = str(os.getenv("FORCE_FIELD", "AMBER")) # Force field for PDB2PQR

# Protein preparation parameters (PDB2PQR/PROPKA)
# pH applied during PROPKA/PDB2PQR steps
PH_VALUE = float(os.getenv("PH_VALUE", 7.40)) # Allow override via env var
# AutoGrid4 Parameters
RECEPTOR_TYPES = "OA SA A C N NA HD" # Receptor atom types
LIGAND_TYPES = "OA A C N NA" # Ligand atom types
GRID_SPACING = 1.000 # Grid spacing (required 1 Å for AutoLigand)
SMOOTH_RADIUS = 0.5 # Potentials smoothing radius

# AD4 distance-dependent dielectric
DIELECTRIC_CONSTANT = float(os.getenv("DIELECTRIC_CONSTANT", -0.1465)) # Allow override via env var

# AutoLigand Parameters
AUTOLIGAND_FILL_POINTS = int(os.getenv("AUTOLIGAND_FILL_POINTS", 90)) # Fill points for AutoLigand (default 90, can be overridden via env var)

# Vina Box Size (Fixed for AutoLigand results)
VINA_BOX_SIZE_TUPLE = (30, 30, 30) # Fixed size used by AMDock for AutoLigand results

# AutoDock Vina Parallelization Control
TOTAL_HTVS_CORES = int(os.getenv('OMP_NUM_THREADS', multiprocessing.cpu_count()))  # Default to all available cores

VINA_CPU_PER_JOB = int(os.getenv('VINA_CPU_PER_JOB', multiprocessing.cpu_count())) # Cores allocated to each individual Vina run. Set low for high HTVS throughput. Otherwise default is all cores


# -----------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------

import re

# Define charge & radius table for your metals
DEFAULT_METAL_PARAMS = {
    "FE":  {"charge":  0.0, "radius": 1.40},   # Fe2+  (e.g., non-heme Fe(II))
    "MG":  {"charge": 0.0, "radius": 1.30},
    "CA":  {"charge": 0.0, "radius": 1.98},
    "MN":  {"charge": 0.0, "radius": 1.39},
}
#override these with environment variable to a json path if provided
METAL_PARAMS = json.loads(os.getenv("METAL_PARAMS_JSON", json.dumps(DEFAULT_METAL_PARAMS)))

def extract_metals_from_pdb(pdb_file):
    metals = []
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname in METAL_PARAMS:
                metals.append(line.rstrip("\n"))
    return metals

def convert_pdb_line_to_pqr(pdb_line):
    """
    Write a strict PDB-style PQR line with correct fixed columns:
    - residue name right-justified in cols 18-20
    - chain in col 22
    - resseq in cols 23-26
    - coords in cols 31-54
    - occupancy (charge) cols 55-60, bfactor (radius) cols 61-66
    - element in cols 77-78
    Returns a single line ending with '\n'.
    """
    import re

    ln = pdb_line.rstrip("\n")
    toks = ln.split()

    record = (ln[0:6].strip() or (toks[0] if toks else "HETATM")).upper()

    # serial may be fused with letters (e.g. "1396MG")
    raw_serial = ln[6:11] if len(ln) >= 11 else (toks[1] if len(toks) > 1 else "0")
    m = re.match(r"^\s*(\d+)([A-Za-z]{1,2})?\s*$", raw_serial)
    if m:
        serial = int(m.group(1))
        fused_element = (m.group(2) or "").upper()
    else:
        digits = re.sub(r"\D", "", raw_serial or "")
        serial = int(digits) if digits else (int(toks[1]) if len(toks) > 1 and toks[1].isdigit() else 0)
        fused_element = ""

    # atom name and resname (use fixed columns when available)
    atom_name = ln[12:16].strip() if len(ln) >= 16 and ln[12:16].strip() else (toks[2] if len(toks) > 2 else "")
    res_name = ln[17:20].strip() if len(ln) >= 20 and ln[17:20].strip() else (toks[3] if len(toks) > 3 else "")

    # chain id (col 22) and residue number (cols 23-26)
    chain = ln[21].strip() if len(ln) > 21 and ln[21].strip() else " "
    resno_field = ln[22:26] if len(ln) >= 26 else (toks[5] if len(toks) > 5 else "")
    resno_digits = re.sub(r"\D", "", resno_field or "")
    res_no = int(resno_digits) if resno_digits else (int(toks[5]) if len(toks) > 5 and toks[5].isdigit() else 0)

    # coordinates (fixed cols) with token fallback
    try:
        x = float(ln[30:38].strip()); y = float(ln[38:46].strip()); z = float(ln[46:54].strip())
    except Exception:
        x = y = z = 0.0
        for i in range(len(toks)-2):
            try:
                tx = float(toks[i]); ty = float(toks[i+1]); tz = float(toks[i+2])
                if all(-999.0 < v < 999.0 for v in (tx, ty, tz)):
                    x, y, z = tx, ty, tz
                    break
            except Exception:
                continue

    # charge/radius tail (PQR style) if present
    charge = None; radius = None
    if len(ln) > 54:
        tail = ln[54:].strip().split()
        if len(tail) >= 2:
            try:
                charge = float(tail[0]); radius = float(tail[1])
            except Exception:
                charge = radius = None

    # determine element: explicit col 77-78 > fused token > atom_name prefix > resname
    element = ""
    if len(ln) >= 78 and ln[76:78].strip():
        element = ln[76:78].strip().capitalize()
    elif fused_element:
        element = fused_element.capitalize()
    else:
        m2 = re.match(r"([A-Za-z]{1,2})", atom_name or "")
        if m2:
            element = m2.group(1).capitalize()
        else:
            element = (res_name[:2].capitalize() if res_name else "").capitalize()
    element = (element or "").rjust(2)[:2]

    # fallback defaults
    if charge is None or radius is None:
        params = METAL_PARAMS.get(element.strip().upper(), METAL_PARAMS.get(res_name.upper(), {"charge": 1.00, "radius": 0.00}))
        charge = params.get("charge", 1.00)
        radius = params.get("radius", 0.00)

    # prepare fields with exact justification: atom name left-justified in 4 cols,
    # residue name RIGHT-justified in 3 cols (cols 18-20), altLoc is blank (col 17)
    atom_name_field = (atom_name or element.strip()).capitalize()
    res_name_field = (res_name or "").upper()[:3].rjust(3)

    # Build fixed-column PDB line exactly placing fields in correct columns.
    # Layout:
    # 1-6 record, 7-11 serial, 12 blank, 13-16 atom name (left-just), 17 altLoc (space),
    # 18-20 resName (right-just), 21 blank, 22 chain, 23-26 resSeq, 27-30 padding,
    # 31-38 x, 39-46 y, 47-54 z, 55-60 occupancy, 61-66 tempFactor, pad to 76, element 77-78
    core = (
        f"{record:<6}"                  # 1-6
        f"{serial:5d}"                  # 7-11
        f" "                            # 12
        f"{atom_name_field:<4}"         # 13-16
        f" "                            # 17 altLoc
        f"{res_name_field:>3}"          # 18-20 (right-justified)
        f" "                            # 21
        f"{chain:1s}"                   # 22
        f"{res_no:4d}"                  # 23-26
        f"   "                          # 27-29 (in total make up to col30)
        f"{x:8.3f}{y:8.3f}{z:8.3f}"     # 31-54
        f"{charge:6.2f}{radius:6.2f}"   # 55-66 (occupancy & bfactor)
    )

    # pad to column 76 then put element at 77-78
    desired_core_len = 76
    if len(core) > desired_core_len:
        core = core[:desired_core_len]
    core = core.ljust(desired_core_len)
    element_field = element.capitalize().rjust(2)[:2]
    return core + element_field + "\n"

def merge_pqr_with_metals(pqr_file, metal_pqr_lines, output_file):
    with open(pqr_file) as f:
        pqr_lines = f.readlines()

    # Strip END lines to append cleanly
    pqr_lines = [ln for ln in pqr_lines if not ln.startswith("END")]

    with open(output_file, "w") as out:
        for ln in pqr_lines:
            out.write(ln)
        for m in metal_pqr_lines:
            out.write(m + "\n")
        out.write("END\n")


def execute_command(command, step_name, log_file, cwd, verbose=False):
    """Utility function to execute shell commands and log output."""
    if verbose:
        print(f"[{os.path.basename(cwd)}] Executing: {command}")
    try:
        # We assume command-line tools are in PATH, consistent with AMDock methodology
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True, 
            cwd=cwd
        )
        if verbose:
            print(f"[{os.path.basename(cwd)}] {step_name} completed successfully.")
        with open(log_file, 'a') as f:
            f.write(f"\n======== {step_name} STDOUT ========\n")
            f.write(result.stdout)
            f.write(f"\n======== {step_name} STDERR ========\n")
            f.write(result.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{os.path.basename(cwd)}] ERROR executing {step_name}.")
        with open(log_file, 'a') as f:
            f.write(f"\n======== {step_name} FAILED STDOUT ========\n")
            f.write(e.stdout)
            f.write(f"\n======== {step_name} FAILED STDERR ========\n")
            f.write(e.stderr)
        return False
    except FileNotFoundError:
        print(f"[{os.path.basename(cwd)}] ERROR: Tool required for {step_name} not found. Ensure tools like pdb2pqr, AutoLigand.py, autogrid4, and vina are accessible in your PATH, preferably installed via Conda/Bioconda.")
        return False
    
def run_pdb2pqr_and_propka(protein_pdb, pqr_path, log_file, cwd='.', verbose=False):
    """
    Run pdb2pqr with PROPKA titration integration and also try to run propka separately
    to reproduce the verbose PROPKA printout seen in the AMDock logs.  Results (stdout/stderr)
    are appended to `log_file`.
    """
    protonated_pdb = pqr_path.replace('.pqr', '_h.pdb')
    pdb2pqr_cmd = (
        f"pdb2pqr --keep-chain --noopt --ff={FORCE_FIELD} --with-ph {PH_VALUE} " f"{protein_pdb} {pqr_path} --nodebump -k --protonate-all "
        f"--titration-state-method=propka --pH {PH_VALUE} -f {protonated_pdb}"
    )
    # run pdb2pqr and capture output
    if verbose:
        print(f"Running pdb2pqr: {pdb2pqr_cmd}")
    try:
        res = subprocess.run(pdb2pqr_cmd, shell=True, check=True, capture_output=True, text=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a') as f:
            f.write("\n======== PDB2PQR/PROPKA FAILED STDOUT ========\n")
            f.write(getattr(e, 'stdout', '') or '')
            f.write("\n======== PDB2PQR/PROPKA FAILED STDERR ========\n")
            f.write(getattr(e, 'stderr', '') or '')
        if verbose:
            print("pdb2pqr failed; see log.")
        return False
    else:
        with open(log_file, 'a') as f:
            f.write("\n======== PDB2PQR/PROPKA STDOUT ========\n")
            f.write(res.stdout or '')
            f.write("\n======== PDB2PQR/PROPKA STDERR ========\n")
            f.write(res.stderr or '')
    return True

def detect_atom_types_from_pdbqt(pdbqt_path):
    """
    Return ordered list of unique AD4 atom-type tokens found in a PDBQT receptor file.
    Heuristics (in order):
    1. trailing whitespace-separated token on ATOM/HETATM lines
    2. element symbol in columns 77-78 (PDB element column)
    3. element inferred from atom name (cols 12-16)
    4. single-letter fallback tokens (preserve if present)
    """
    tokens = []
    seen = set()
    if not os.path.exists(pdbqt_path):
        return tokens

    with open(pdbqt_path, "r") as fh:
        for ln in fh:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            ln_stripped = ln.rstrip("\n\r")
            # 1) trailing token
            parts = ln_stripped.split()
            token = None
            if parts:
                last = parts[-1].strip()
                # accept if contains letters (e.g. "OA", "NA", "Mg", "HD", "A")
                if re.search(r"[A-Za-z]", last):
                    token = last
            # 2) element in fixed columns 77-78 (1-based)
            if not token:
                if len(ln_stripped) >= 78:
                    elem_col = ln_stripped[76:78].strip()
                    if elem_col and re.match(r"^[A-Za-z]{1,2}$", elem_col):
                        token = elem_col
            # 3) atom name-based fallback (cols 13-16)
            if not token:
                if len(ln_stripped) >= 16:
                    atom_name = ln_stripped[12:16].strip()
                    if atom_name:
                        m = re.match(r"^([A-Za-z]{1,2})", atom_name)
                        if m:
                            token = m.group(1)
            # 4) final fallback: single-letter tokens from anywhere in last few columns
            if not token:
                tail = ln_stripped[-10:].strip() if len(ln_stripped) >= 10 else ln_stripped.strip()
                m = re.search(r"([A-Za-z])\b", tail)
                if m:
                    token = m.group(1)

            if token:
                # normalize case to preserve what AutoGrid/prepare_receptor4 sees - keep original case
                if token not in seen:
                    seen.add(token)
                    tokens.append(token)

    return tokens


def detect_receptor_types_from_pdbqt(pdbqt_path):
    """Return a space-separated string of unique atom types found in a PDBQT (simple heuristic)."""
    types = []
    if not os.path.exists(pdbqt_path):
        return None
    with open(pdbqt_path, 'r') as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                parts = ln.strip().split()
                if not parts:
                    continue
                # heuristic: last token often contains the AD4 atom type / element token
                last = parts[-1]
                # ignore coordinate/charge numeric tokens
                if last.isalpha():
                    types.append(last)
                else:
                    # fallback: try to extract element from columns (element symbol usually at end)
                    # simple attempt: take 1-2 char element from atom name field
                    atom_name = parts[2] if len(parts) > 2 else ""
                    elem = re.sub(r'[^A-Za-z]', '', atom_name).strip()[:2].upper()
                    if elem:
                        types.append(elem)
    unique = sorted(set(types), key=lambda x: (len(x), x))
    return " ".join(unique) if unique else None

def generate_gpf_content(target_pdbqt, grid_center, npts, receptor_types=None, ligand_types=None):
    """
    Generate GPF content. Accepts receptor_types/ligand_types as lists (preferred) or space-separated strings.
    Emits exactly one 'map base.TYPE.map' line per ligand type token (AutoGrid expects maps for ligand_types).
    """
    cx, cy, cz = grid_center
    nx, ny, nz = npts

    def norm(x, fallback):
        if not x:
            x = fallback
        if isinstance(x, str):
            toks = [t for t in re.split(r'\s+', x.strip()) if t]
        else:
            toks = list(x)
        # preserve order & uniqueness
        out = []
        seen = set()
        for t in toks:
            if t not in seen:
                seen.add(t)
                out.append(t)
        return out

    receptor_tokens = norm(receptor_types, RECEPTOR_TYPES)
    ligand_tokens = norm(ligand_types, LIGAND_TYPES)

    base = os.path.splitext(os.path.basename(target_pdbqt))[0]
    # IMPORTANT: AutoGrid expects one "map" line per ligand atom type (not receptor types)
    map_lines = "\n".join(f"map {base}.{t}.map" for t in ligand_tokens)

    lines = [
        f"npts {nx} {ny} {nz} # num.grid points in xyz",
        f"gridfld {base}.maps.fld # grid_data_file",
        f"spacing {GRID_SPACING} # spacing(A)",
        f"receptor_types {' '.join(receptor_tokens)} # receptor atom types",
        f"ligand_types {' '.join(ligand_tokens)} # ligand atom types",
        f"receptor {base}.pdbqt # macromolecule",
        f"gridcenter {cx} {cy} {cz} # xyz-coordinates",
        f"smooth {SMOOTH_RADIUS} # store minimum energy w/in rad(A)",
        map_lines,
        f"elecmap {base}.e.map",
        f"dsolvmap {base}.d.map",
        f"dielectric {DIELECTRIC_CONSTANT}",
        ""
    ]
    return "\n".join(lines)


# -----------------------------------------------------------
# DYNAMIC BINDING SITE PARSING (Addressing user query)
# -----------------------------------------------------------

def parse_pdb_coordinates(file_path):
    """Reads a PDB/FILL_outX.pdb file and extracts fixed-width coordinates."""
    coordinates = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # PDB coordinate fields are fixed-width: X (31-38), Y (39-46), Z (47-54)
                if line.startswith(('ATOM', 'HETATM', 'HET')): 
                    try:
                        # Extract coordinates using fixed indices (0-indexed Python slicing)
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates.append((x, y, z))
                    except (ValueError, IndexError):
                        # Skip lines where parsing fails (e.g., bad formatting or short lines)
                        continue
    except FileNotFoundError:
        print(f"Error: Fill file not found at {file_path}")
        return []
    return coordinates

def calculate_geometric_center(coordinates):
    """Calculates the geometric center (average) of a list of coordinates, using NumPy if available."""
    if not coordinates:
        return None
    
    if NUMPY_AVAILABLE:
        # Use NumPy for efficient, robust geometric calculation (better parsing via Conda tools)
        coords_array = np.array(coordinates)
        center = np.mean(coords_array, axis=0)
        return tuple(center)
    else:
        # Fallback to pure Python calculation
        N = len(coordinates)
        x_sum, y_sum, z_sum = 0, 0, 0
        for x, y, z in coordinates:
            x_sum += x
            y_sum += y
            z_sum += z
        return (x_sum / N, y_sum / N, z_sum / N)

def compute_grid_npts(min_xyz, max_xyz, spacing=1.0, padding=3.0, require_even=True):
    """
    Return (nx, ny, nz) for AutoGrid given receptor min/max coords.
    - min_xyz, max_xyz: iterables of 3 floats
    - spacing: grid spacing in Å (default 1.0)
    - padding: margin added on each side in Å (total +2*padding)
    - require_even: if True, force npts to be even (some autogrid builds require this)
    """
    import math
    nx = math.ceil(((max_xyz[0] - min_xyz[0]) + 2*padding) / spacing)
    ny = math.ceil(((max_xyz[1] - min_xyz[1]) + 2*padding) / spacing)
    nz = math.ceil(((max_xyz[2] - min_xyz[2]) + 2*padding) / spacing)

    def adjust(v):
        v = int(v)
        if v < 2:
            v = 2
        if require_even and (v % 2 != 0):
            v += 1
        return v

    return (adjust(nx), adjust(ny), adjust(nz))

def run_autoligand_and_parse_sites(base_dir, protein_name, log_file, fill_points, verbose=False):
    """
    Executes AutoLigand (Automatic Mode) and dynamically parses the output PDB files
    to define the Vina box centers.
    """
    
    protein_pdbqt_basename = f"{protein_name}_h.pdbqt"
    
    if verbose:
        print(f"[{os.path.basename(base_dir)}] --- Running AutoLigand (Automatic Mode) ---")
    
    # AutoLigand.py requires the map base name and the fill points count.
    map_base_name = protein_pdbqt_basename.replace('.pdbqt', '')
    autoligand_cmd = f"AutoLigand -r {map_base_name} -p {fill_points}"
    
    if not execute_command(autoligand_cmd, "AutoLigand Site Search (Automatic)", log_file, base_dir):
        return []

    if verbose:       
        print(f"[{os.path.basename(base_dir)}] AutoLigand command executed. Parsing generated FILL files...")
    
    discovered_sites = []
    
    # 1. Find all generated FILL files (FILL_#outX.pdb)
    fill_files = glob.glob(os.path.join(base_dir, f"FILL_{fill_points}out*.pdb"))
    
    # Sort files numerically by rank (out1, out2, etc.)
    try:
        fill_files.sort(key=lambda x: int(re.search(r'out(\d+)\.pdb$', x).group(1)))
    except AttributeError:
        # Handle cases where fill files don't follow the expected naming convention
        print(f"[{os.path.basename(base_dir)}] Warning: Could not sort FILL files numerically.")
        pass

    if not fill_files:
        print(f"[{os.path.basename(base_dir)}] Error: No FILL files found after running AutoLigand. Site determination failed.")
        return []

    # 2. Process each fill file to calculate the geometric center
    for file_path in fill_files:
        coordinates = parse_pdb_coordinates(file_path)
        center = calculate_geometric_center(coordinates)

        if center:
            # Box size is fixed to 30x30x30, as observed in AMDock output for AutoLigand mode
            site_info = {
                "center": center, 
                "size": VINA_BOX_SIZE_TUPLE 
            }
            discovered_sites.append(site_info)

    if verbose:
        print(f"[{os.path.basename(base_dir)}] AutoLigand identified {len(discovered_sites)} binding site(s).")

    # AMDock typically reports the top 10 sites found.
    return discovered_sites[:10]


# -----------------------------------------------------------
# CORE DOCKING FUNCTION FOR HTVS
# -----------------------------------------------------------

def run_amdock_pipeline(job_config, verbose=False):
    """
    Performs the full AMDock-replicated pipeline (prep, grid, AutoLigand site search, Vina docking).
    """
    # Unpack configuration
    base_dir = job_config['working_dir']
    protein_pdb = job_config['protein_pdb']
    ligand_pdb = job_config['ligand_pdb']
    protein_name = job_config['protein_name']
    ligand_name = job_config['ligand_name']
    keep_metals = job_config.get('keep_metals', False)

    # Define intermediate file names
    protein_h_pdb = os.path.join(base_dir, f"{protein_name}_h.pdb")
    protein_pdbqt = f"{protein_name}_h.pdbqt"
    ligand_pdbqt = f"{ligand_name}.pdbqt"
    gpf_file = "grid.gpf"
    log_file = os.path.join(base_dir, f"{protein_name}_{ligand_name}_pipeline.log")
    
    os.makedirs(base_dir, exist_ok=True)
    
    # --- 1. Prepare Target Protein (PDB2PQR + PROPKA) ---
    # Run pdb2pqr in the working_dir and attempt a standalone propka run to reproduce AMDock verbosity.
    pqr_path = os.path.join(base_dir, f"{protein_name}_h.pqr")
    if not run_pdb2pqr_and_propka(protein_pdb, pqr_path, log_file, cwd=base_dir, verbose=verbose):
        print(f"[{os.path.basename(base_dir)}] ERROR: pdb2pqr/propka step failed. See {log_file}")
        return False

    # If metals are to be kept, extract and merge them into the PQR

    if keep_metals:
        metals = extract_metals_from_pdb(protein_pdb)
        metal_pqr_lines = [convert_pdb_line_to_pqr(m) for m in metals]
        merge_pqr_with_metals(pqr_path, metal_pqr_lines, pqr_path)
        if verbose:
            print(f"[{os.path.basename(base_dir)}] Merged {len(metal_pqr_lines)} metal ions into PQR.")
    
    # Uses Gasteiger charges by default
    pdb_h_path = os.path.join(base_dir, f"{protein_name}_h.pdb")

    # ensure pdb2pqr created the PQR
    if not os.path.exists(pqr_path):
        print(f"ERROR: expected PQR not found: {pqr_path}")
        return False

    # Control the verbosity of OpenBabel
    ob_outlvl = 2 if verbose else 0
    ob.obErrorLog.SetOutputLevel(ob_outlvl) 
    if not verbose:
        pybel.ob.obErrorLog.StopLogging()

    # Convert PQR to PDB using pybel (OpenBabel)
    # Maybe something wrong here!!!!!
    mols = list(pybel.readfile("pqr", pqr_path))
    if not mols:
        raise RuntimeError("No molecules read from PQR via pybel")
    mols[0].write("pdb", pdb_h_path, overwrite=True)

    if verbose:
        print(f"[{os.path.basename(base_dir)}] Converted PQR -> PDB using pybel: {pdb_h_path}")

    # then use pdb_h_path for prepare_receptor4
    if keep_metals:
        #Extract the names of the metal atoms
        metal_atom_names = [line[12:16].strip() for line in metals]
        metal_atoms_string = ' '.join(metal_atom_names)
        receptor_prepare_cmd = f"prepare_receptor4 -r {os.path.basename(pdb_h_path)} -C -p {metal_atoms_string} -U nphs_lps_waters_nonstdres_deleteAltB -o {os.path.basename(protein_pdbqt)}"
    else:
        receptor_prepare_cmd = f"prepare_receptor4 -r {os.path.basename(pdb_h_path)} -o {os.path.basename(protein_pdbqt)}"
    if not execute_command(receptor_prepare_cmd, "Prepare_Receptor4", log_file, base_dir, verbose=verbose): return False

    # --- 3. Convert Ligand PDB to PDBQT (Prepare_Ligand4) ---
    ligand_prepare_cmd = f"prepare_ligand4 -l {ligand_pdb} -o {ligand_pdbqt}"
    if not execute_command(ligand_prepare_cmd, "Prepare_Ligand4", log_file, base_dir, verbose=verbose): return False

    protein_coordinates = parse_pdb_coordinates(protein_pdb)
    if not protein_coordinates:
        print(f"[{os.path.basename(base_dir)}] ERROR: Unable to parse protein coordinates for grid setup.")
        return False    
    grid_center_full_protein = calculate_geometric_center(protein_coordinates)
    grid_npts = compute_grid_npts(
        min_xyz=(min(c[0] for c in protein_coordinates),
                min(c[1] for c in protein_coordinates),
                min(c[2] for c in protein_coordinates)),
        max_xyz=(max(c[0] for c in protein_coordinates),
                max(c[1] for c in protein_coordinates),
                max(c[2] for c in protein_coordinates)),
        spacing=GRID_SPACING,
        padding=3.0,
        require_even=True
    )
    # --- 4. Generate AutoGrid4 GPF and Run AutoGrid4 (Prerequisite for AutoLigand) ---
    # detect receptor atom types from the pdbqt and use them in the GPF to avoid autogrid mismatch
    receptor_tokens = detect_atom_types_from_pdbqt(os.path.join(base_dir, protein_pdbqt))
    ligand_tokens = detect_atom_types_from_pdbqt(os.path.join(base_dir, ligand_pdbqt))
    gpf = generate_gpf_content(protein_pdbqt, grid_center_full_protein, grid_npts, receptor_types=receptor_tokens, ligand_types=ligand_tokens)
    with open(os.path.join(base_dir, gpf_file), "w") as fh:
        fh.write(gpf)
    
    autogrid_cmd = f"autogrid4 -p {gpf_file} -l grid.glg"
    if not execute_command(autogrid_cmd, "AutoGrid4", log_file, base_dir, verbose=verbose): return False
    
    # --- 5. Binding Site Detection using AutoLigand (Automatic Mode) ---
    sites_to_dock = run_autoligand_and_parse_sites(base_dir, protein_name, log_file, AUTOLIGAND_FILL_POINTS, verbose=verbose)

    if not sites_to_dock:
        print(f"[{os.path.basename(base_dir)}] ERROR: AutoLigand site detection failed.")
        return False
        
    # --- 6. Run AutoDock Vina for all identified binding sites with Multiprocessing ---
    successful_runs = 0
    
    # Vina uses its internal threading (--cpu).
    for i, site in enumerate(sites_to_dock):
        cx, cy, cz = site['center']
        sx, sy, sz = site['size']
        
        output_pdbqt = f"{ligand_name}_pose_site{i+1}.pdbqt"
        log_file_vina = f"vina_site{i+1}.log"
        
        vina_cmd = (
            f"vina --receptor {protein_pdbqt} "
            f"--ligand {ligand_pdbqt} "
            f"--center_x {cx} --center_y {cy} --center_z {cz} "
            f"--size_x {sx} --size_y {sy} --size_z {sz} "
            f"--cpu {VINA_CPU_PER_JOB} "
            f"--out {output_pdbqt} "
            f"> {log_file_vina} 2>&1"
        )

        if execute_command(vina_cmd, f"AutoDock Vina (Site {i+1})", log_file, base_dir, verbose=verbose):
            successful_runs += 1
            mol_out = list(pybel.readfile("pdbqt", os.path.join(base_dir, output_pdbqt)))
            mol_out[0].write("pdb", os.path.join(base_dir, output_pdbqt.replace('.pdbqt', '.pdb')), overwrite=True)
        else:
            print(f"[{os.path.basename(base_dir)}] Vina docking failed for site {i+1}.")
            
    if successful_runs > 0:
        if verbose:
            print(f"[{os.path.basename(base_dir)}] Docking completed successfully.")
            print(f"\n*** Docking for {protein_name} + {ligand_name} Complete (Total sites docked: {len(sites_to_dock)}). ***")
        return True
    else:
        print(f"\n*** Docking for {protein_name} + {ligand_name} FAILED ENTIRELY. ***")
        return False

# -----------------------------------------------------------
# HTVS EXECUTION (PARALLELIZATION)
# -----------------------------------------------------------

def run_htvs(jobs_to_run, verbose=False):
    """
    Sets up Python's multiprocessing to run multiple docking jobs in parallel (HTVS).
    """
    
    # Calculate max parallel HTVS processes based on total cores and Vina's internal use.
    max_processes = max(1, TOTAL_HTVS_CORES // VINA_CPU_PER_JOB)
    
    if verbose:
        print(f"\n--- Starting HTVS Parallel Execution ---")
        print(f"Vina internal cores per job: {VINA_CPU_PER_JOB}")
        print(f"Maximum simultaneous HTVS jobs: {max_processes}")
        print(f"Total jobs requested: {len(jobs_to_run)}")

    with multiprocessing.Pool(processes=max_processes) as pool:
        results = pool.map(partial(run_amdock_pipeline, verbose=verbose), jobs_to_run)
        
    successful_jobs = sum(results)
    if verbose:
        print(f"\n--- HTVS Finished ---")
    print(f"Summary: {successful_jobs}/{len(jobs_to_run)} pairs successfully processed.")

# -----------------------------------------------------------
# ANALYZE HTVS RESULTS
# -----------------------------------------------------------


def read_file(path):
    with open(path, 'r', errors='replace') as f:
        return f.read().splitlines()

def find_tool_versions(lines):
    info = {}
    for ln in lines:
        m = re.search(r'PDB2PQR v?([0-9.]+)', ln, re.IGNORECASE)
        if m: info['pdb2pqr'] = m.group(1)
        m = re.search(r'propka(?:\s|3)?(?:\.\d+)?\s*([0-9.\-]*)', ln, re.IGNORECASE)
        if m and 'propka' not in info:
            # coarse capture (propka header lines are verbose)
            if 'propka' in ln.lower() or ln.strip().lower().startswith('propka'):
                info['propka_header'] = ln.strip()
        m = re.search(r'AMDOCK.*Version\s*([0-9.]+)', ln, re.IGNORECASE)
        if m: info['amdock'] = m.group(1)
        m = re.search(r'AutoDock Vina', ln, re.IGNORECASE)
        if m: info['vina'] = info.get('vina', 'present')
    return info

def extract_pdb2pqr_stats(lines):
    stats = {}
    for ln in lines:
        m = re.search(r'Created biomolecule object with\s+(\d+)\s+residues\s+and\s+(\d+)\s+atoms', ln, re.IGNORECASE)
        if m:
            stats['residues'] = int(m.group(1))
            stats['atoms'] = int(m.group(2))
        m = re.search(r'Unable to find amino or nucleic acid definition for\s+(\S+)', ln)
        if m:
            stats.setdefault('unknown_components', []).append(m.group(1))
        m = re.search(r'WARNING:Missing atom (\S+) in residue (\S+)\s+([A-Z])\s+(\d+)', ln)
        if m:
            stats.setdefault('missing_atoms', []).append({'atom':m.group(1),'residue':m.group(2),'chain':m.group(3),'resnum':m.group(4)})
        m = re.search(r'Added atom (\S+) to residue (\S+) at coordinates\s*([-\d\.]+),\s*([-\d\.]+),\s*([-\d\.]+)', ln)
        if m:
            stats.setdefault('added_atoms',[]).append({'atom':m.group(1),'residue':m.group(2),'coords':(float(m.group(3)),float(m.group(4)),float(m.group(5)))})
    return stats

def extract_propka_pkas(lines):
    # find the "INFO:pKa summary:" block and capture residue pKa lines following it
    pkas = []
    start = None
    for i, ln in enumerate(lines):
        if re.search(r'INFO:pKa summary', ln, re.IGNORECASE) or re.search(r'SUMMARY OF THIS PREDICTION', ln):
            start = i
            break
    if start is None:
        # fallback: find the dashed header that precedes the pKa summary table
        for i, ln in enumerate(lines):
            if re.match(r'-{5,}', ln):
                start = i+1
                break
    if start is None:
        return pkas
    # scan forward for lines matching residue pKa entries
    for ln in lines[start:]:
        if not ln.strip():
            # blank line often ends the table
            if pkas:
                break
            else:
                continue
        # match lines like: "ASP  36 A   3.19     0 % ..."
        m = re.match(r'\s*([A-Z]+)\s+(\d+)\s+([A-Z])\s+([0-9]+\.[0-9]+)', ln)
        if m:
            resname, resnum, chain, pka = m.group(1), m.group(2), m.group(3), float(m.group(4))
            pkas.append({'residue': f"{resname} {resnum}{chain}", 'pKa': pka})
        # also handle lines like "BTN  N1 A   4.20"
        else:
            m2 = re.match(r'\s*([A-Z0-9\-\+]+)\s+([A-Z0-9]+)\s+([A-Z])\s+([0-9]+\.[0-9]+)', ln)
            if m2:
                # coarse capture
                pkas.append({'residue': f"{m2.group(1)} {m2.group(2)}{m2.group(3)}", 'pKa': float(m2.group(4))})
    return pkas

def extract_propka_metadata(lines):
    meta = {}
    for ln in lines:
        m = re.search(r'propka([^\n]*)', ln, re.IGNORECASE)
        if m and 'propka_banner' not in meta:
            meta['propka_banner'] = ln.strip()
        m2 = re.search(r'Applying pKa values at a pH of\s*([0-9.]+)', ln, re.IGNORECASE)
        if m2:
            meta['pH_applied'] = float(m2.group(1))
        m3 = re.search(r'The pI is\s+([0-9.]+)', ln)
        if m3:
            meta['pI'] = float(m3.group(1))
    return meta

def extract_autoligand_outputs(lines):
    # parse lines like "Output # 1  Total Volume =  234.0  Total Energy per Vol, EPV =  -0.2140784829"
    outputs = []
    for ln in lines:
        m = re.search(r'Output\s*#\s*(\d+).*Total Volume\s*=\s*([0-9.]+).*EPV\s*=\s*([-\d\.]+)', ln)
        if m:
            outputs.append({'rank': int(m.group(1)), 'volume': float(m.group(2)), 'EPV': float(m.group(3))})
    return outputs

def extract_vina_affinities_from_log(lines):
    # detect typical vina affinity lines or the generic "-----+" table used in some logs.
    # Return a list of numeric affinities (floats) so callers can compute stats or format tables.
    affinities = []
    for i, ln in enumerate(lines):
        # Vina lines like: "1    -7.3    0.000    ...", but many formats exist
        m = re.match(r'^\s*\d+\s+([-\d\.]+)\s+', ln)
        if m:
            try:
                affinities.append(float(m.group(1)))
            except ValueError:
                pass
        # check for custom markers: previous pipeline used "-----+" followed by a line with affinity
        if ln.startswith("-----+"):
            if i+1 < len(lines):
                parts = lines[i+1].split()
                if len(parts) >= 2:
                    try:
                        affinities.append(float(parts[1]))
                    except Exception:
                        pass
    # Return raw numeric affinities for downstream processing/formatting
    return affinities

def collect_warnings(lines):
    warns = []
    for ln in lines:
        if 'WARNING:' in ln or ln.strip().startswith('WARNING') or ln.strip().upper().startswith('WARN'):
            warns.append(ln.strip())
    return warns

def summarize_log(path):
    lines = read_file(path)
    summary = {}
    summary['file'] = os.path.abspath(path)
    summary['tools'] = find_tool_versions(lines)
    summary['pdb2pqr'] = extract_pdb2pqr_stats(lines)
    summary['propka_meta'] = extract_propka_metadata(lines)
    pkas = extract_propka_pkas(lines)
    summary['propka_pkas_count'] = len(pkas)
    if pkas:
        # compute extremes
        sorted_pkas = sorted(pkas, key=lambda x: x['pKa'])
        summary['propka_lowest'] = sorted_pkas[:5]
        summary['propka_highest'] = sorted_pkas[-5:][::-1]
    summary['autoligand_outputs'] = extract_autoligand_outputs(lines)
    summary['vina_affinities_in_log'] = extract_vina_affinities_from_log(lines)
    summary['warnings'] = collect_warnings(lines)
    # try to find directory and parse vina_site*.log files if present
    log_dir = os.path.dirname(path) or '.'
    vina_files = sorted([os.path.join(log_dir, f) for f in os.listdir(log_dir) if re.match(r'vina.*\.log$', f, re.IGNORECASE)])
    all_vina_affs = {}
    for vf in vina_files:
        try:
            affs = extract_vina_affinities_from_log(read_file(vf))
            if affs:
                all_vina_affs[os.path.basename(vf)] = affs
        except Exception:
            pass
    if all_vina_affs:
        summary['vina_files'] = all_vina_affs
    return summary

def pretty_print(summary):
    print("\n=== DOCKING LOG SUMMARY ===\n")
    print(f"Log file: {summary.get('file')}\n")
    tools = summary.get('tools', {})
    if tools:
        print("Detected tools / headers:")
        for k,v in tools.items():
            print(f"  - {k}: {v}")
        print()
    pdb2 = summary.get('pdb2pqr', {})
    if pdb2:
        print("PDB2PQR / structure info:")
        if 'residues' in pdb2:
            print(f"  Residues: {pdb2['residues']}, Atoms: {pdb2['atoms']}")
        if pdb2.get('unknown_components'):
            print(f"  Unknown residue/component definitions: {', '.join(pdb2['unknown_components'])}")
        if pdb2.get('missing_atoms'):
            print(f"  Missing atoms detected: {len(pdb2['missing_atoms'])} (see warnings)")
        if pdb2.get('added_atoms'):
            print(f"  Atoms added by PDB2PQR: {len(pdb2['added_atoms'])}")
        print()
    prop_meta = summary.get('propka_meta', {})
    if prop_meta:
        print("PROPKA / pKa application:")
        if 'pH_applied' in prop_meta:
            print(f"  pH applied: {prop_meta['pH_applied']}")
        if 'pI' in prop_meta:
            print(f"  Calculated pI (folded): {prop_meta['pI']}")
        if 'propka_banner' in prop_meta:
            print(f"  {prop_meta['propka_banner']}")
        print(f"  pKa entries parsed: {summary.get('propka_pkas_count',0)}")
        if summary.get('propka_lowest'):
            print("  Lowest pKa residues (up to 5):")
            for r in summary['propka_lowest']:
                print(f"    {r['residue']:>12s}  pKa={r['pKa']:.2f}")
        if summary.get('propka_highest'):
            print("  Highest pKa residues (up to 5):")
            for r in summary['propka_highest']:
                print(f"    {r['residue']:>12s}  pKa={r['pKa']:.2f}")
        print()
    al = summary.get('autoligand_outputs', [])
    if al:
        print("AutoLigand (Automatic) site search:")
        for out in al:
            print(f"  Rank {out['rank']}: Volume={out['volume']:.1f}, EPV={out['EPV']:.4f}")
        print()
    vina_log_affs = summary.get('vina_affinities_in_log', [])
    if vina_log_affs:
        print("Affinities parsed from main log (may be aggregated):")
        rows = [(f"pose_{i+1}", f"{aff:.3f}") for i, aff in enumerate(vina_log_affs)]
        if rows:
            print(tabulate(rows, headers=["Pose ID", "Affinity (kcal/mol)"]))
        else:
            print("  (no affinities parsed)")
        print()

    if summary.get('vina_files'):
        print("Affinities found in vina log files:")
        best_affinities = []
        for fname, affs in summary['vina_files'].items():
            short_name = os.path.basename(fname)
            if affs:
                best_aff = min(affs)
                best_affinities.append(best_aff)
            else:
                print(f"  {short_name}: no affinities parsed")

        if best_affinities:
            rows = [(f"{i+1}", f"{float(aff):2.3f}") for i, aff in enumerate(best_affinities)]
            # sort aggregated best affinities by numeric value (increasing)
            rows.sort(key=lambda r: float(r[1]))
            print("\nAggregated best affinities per pose:")
            print(tabulate(rows, headers=["Pose ID", "Best Affinity (kcal/mol)"]))
        print()


def analyze_htvs_results(jobs_to_run):
    """
    Analyzes the log files from HTVS jobs and summarizes key information.
    """
    for job in jobs_to_run:
        working_dir = job['working_dir']
        log_file = os.path.join(working_dir, f"{job['protein_name']}_{job['ligand_name']}_pipeline.log")
        if os.path.exists(log_file):
            summary = summarize_log(log_file)
            print(f"\n--- Summary for {job['protein_name']} + {job['ligand_name']} ---")
            pretty_print(summary)
            
            print(f"--- End of Summary for {job['protein_name']} + {job['ligand_name']} ---\n")
        else:
            print(f"Log file not found for job: {log_file}")

if __name__ == '__main__':
    # =========================================================================
    # PREPARATION NOTES:
    # 1. ENVIRONMENT: Ensure a Conda environment is activated and contains: 
    #    pdb2pqr, prepare_receptor4.py, prepare_ligand4.py, autogrid4, vina, and AutoLigand.py.
    #    (See source for dependency list). numpy is highly recommended.
    # 2. INPUT FILES: Ensure actual 1hsg_protein.pdb and mk1.pdb files exist in the specified 'input' directory.
    # =========================================================================

    BASE_PROJECT_DIR = os.path.abspath("HTVS_AMDock_Project")
    
    MOCK_INPUT_DIR = os.path.join(BASE_PROJECT_DIR, "input")
    # These must be replaced with real PDB content for PDB2PQR to run correctly
    MOCK_PROTEIN_PDB = os.path.join(MOCK_INPUT_DIR, "1hsg_protein.pdb")
    MOCK_LIGAND_MK1 = os.path.join(MOCK_INPUT_DIR, "mk1.pdb")
    
    os.makedirs(MOCK_INPUT_DIR, exist_ok=True)
    
    # Create mock PDB files (Minimal content provided, MUST BE REAL PDB DATA TO RUN TOOLS)
    if not os.path.exists(MOCK_PROTEIN_PDB):
        with open(MOCK_PROTEIN_PDB, 'w') as f:
            f.write("ATOM      1  N   PHE A   1      29.361  39.686   5.862  1.00 66.45           N \n")
            f.write("ATOM   1849 OXT PHE B  99      26.380  42.929   4.902  1.00 66.45           O \n")
            f.write("END\n")
            
    if not os.path.exists(MOCK_LIGAND_MK1):
        with open(MOCK_LIGAND_MK1, 'w') as f:
            f.write("HETATM    1  N   LIG A   1      10.000  20.000  30.000  1.00 50.00           N\n")
            f.write("END\n")

    # Define jobs for the HTVS queue
    jobs = [
        # Job 1: 1hsg-mk1 docking
        {
            "protein_name": "1hsg_protein",
            "ligand_name": "mk1_A",
            "protein_pdb": MOCK_PROTEIN_PDB,
            "ligand_pdb": MOCK_LIGAND_MK1,
            "working_dir": os.path.join(BASE_PROJECT_DIR, "job_1hsg_mk1_A"),
        },
        # Job 2: Demonstrating parallelization on a second ligand/run
        {
            "protein_name": "1hsg_protein",
            "ligand_name": "ligand_Y",
            "protein_pdb": MOCK_PROTEIN_PDB,
            "ligand_pdb": MOCK_LIGAND_MK1, 
            "working_dir": os.path.join(BASE_PROJECT_DIR, "job_1hsg_ligY"),
        },
    ]
    
    run_htvs(jobs)
    analyze_htvs_results(jobs)