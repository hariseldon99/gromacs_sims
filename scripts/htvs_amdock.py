#!/usr/bin/env python3
"""
This script integrates the comprehensive AMDock workflow for High-Throughput Virtual Screening (HTVS), including dynamic binding site detection using **AutoLigand in Automatic Mode**, multi-processing control, and optimized PDB coordinate parsing designed to handle the fixed-width fields characteristic of PDB/FILL files, as seen in the source logs.

To address your query regarding using Bioconda for better PDB parsing, while the sources confirm AMDock and its dependencies (`pdb2pqr`, `openbabel`, `pymol`) are typically installed via **Conda**, the actual coordinate parsing of the generated `FILL` files is best handled using standard Python libraries optimized for array operations (`numpy`) to accurately read the fixed columns shown in the logs.
"""
### AMDock HTVS Replication Script with Dynamic Site Detection

import subprocess
import os
import multiprocessing
import glob
import re
import shutil
import platform
import openbabel as ob
from openbabel import pybel  # common import name for OpenBabel's Python wrapper
from functools import partial

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
TOTAL_HTVS_CORES = int(os.getenv('OMP_NUM_THREADS', 1))  # Default to serial

VINA_CPU_PER_JOB = int(os.getenv('VINA_CPU_PER_JOB', multiprocessing.cpu_count())) # Cores allocated to each individual Vina run. Set low for high HTVS throughput. Otherwise default is all cores


# -----------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------

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

def generate_gpf_content(target_pdbqt, grid_center, npts):
    """Generates the AutoGrid GPF file content replicating AMDock's parameters."""
    cx, cy, cz = grid_center
    nx, ny, nz = npts
    
    # Grid parameters must match those required for AutoLigand (1 Å spacing, specific maps)
    gpf_content = f"""
npts {nx} {ny} {nz} # num.grid points in xyz
gridfld {target_pdbqt.replace('.pdbqt', '.maps.fld')} # grid_data_file
spacing {GRID_SPACING} # spacing(A)
receptor_types {RECEPTOR_TYPES} # receptor atom types
ligand_types {LIGAND_TYPES} # ligand atom types
receptor {target_pdbqt} # macromolecule
gridcenter {cx} {cy} {cz} # xyz-coordinates
smooth {SMOOTH_RADIUS} # store minimum energy w/in rad(A)
map {target_pdbqt.replace('.pdbqt', '.OA.map')} # atom-specific affinity map
map {target_pdbqt.replace('.pdbqt', '.A.map')} # atom-specific affinity map
map {target_pdbqt.replace('.pdbqt', '.C.map')} # atom-specific affinity map
map {target_pdbqt.replace('.pdbqt', '.N.map')} # atom-specific affinity map
map {target_pdbqt.replace('.pdbqt', '.NA.map')} # atom-specific affinity map
elecmap {target_pdbqt.replace('.pdbqt', '.e.map')} # electrostatic potential map
dsolvmap {target_pdbqt.replace('.pdbqt', '.d.map')} # desolvation potential map
dielectric {DIELECTRIC_CONSTANT} # AD4 distance-dependent dielectric
"""
    return gpf_content


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

    # Define intermediate file names
    protein_h_pdb = os.path.join(base_dir, f"{protein_name}_h.pdb")
    protein_pdbqt = f"{protein_name}_h.pdbqt"
    ligand_pdbqt = f"{ligand_name}.pdbqt"
    gpf_file = "grid.gpf"
    log_file = os.path.join(base_dir, f"{protein_name}_{ligand_name}_pipeline.log")
    
    os.makedirs(base_dir, exist_ok=True)
    
    # --- 1. Prepare Target Protein (PDB2PQR + PROPKA) ---
    # AMDock uses PDB2PQR v3.6.1 and PROPKA v3.5.1 to apply pKa values at pH 7.40.
    pdb2pqr_cmd = (
        f"pdb2pqr --ff={FORCE_FIELD} --with-ph {PH_VALUE} {protein_pdb} "
        f"{protein_h_pdb.replace('.pdb', '.pqr')} --clean --nodebump -k --protonate-all --titration-state-method=propka --pH {PH_VALUE}"
    )
    if not execute_command(pdb2pqr_cmd, "PDB2PQR/PROPKA", log_file, base_dir, verbose=verbose): return False

    # --- 2. Convert prepared Receptor PDB to PDBQT (Prepare_Receptor4) ---
    # Uses Gasteiger charges by default
    pqr_path = os.path.join(base_dir, f"{protein_name}_h.pqr")
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
    mols = list(pybel.readfile("pqr", pqr_path))
    if not mols:
        raise RuntimeError("No molecules read from PQR via pybel")
    out_path = os.path.join(base_dir, f"{protein_name}_h.pdb")
    mols[0].write("pdb", out_path, overwrite=True)

    if verbose:
        print(f"[{os.path.basename(base_dir)}] Converted PQR -> PDB using pybel: {out_path}")

    # then use pdb_h_path for prepare_receptor4
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
    gpf_content = generate_gpf_content(protein_pdbqt, grid_center_full_protein, grid_npts)

    with open(os.path.join(base_dir, gpf_file), "w") as f:
        f.write(gpf_content)
    
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