#  Protein-Ligand MD Simulation — Step-by-Step Guide

**System:** Default CUEDC2 protein (AlphaFold model) + Default ERGOTAMINE (resname `ERG`) - can be readily adapted to any protein-ligand complex  
**Force field:** AMBER99SB-ILDN + AMBER (acpype) for ligand  
**Water model:** TIP3P  
**Production:** 100 ns at 310 K, 1 bar (NPT)  
**Singularity image (local):** `$HOME/images/gromacs_2024.5.1-GPU_20260718.sif`  
**Singularity image (HPC):** `$HOME/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif`

---

## Overview

| Step | Script | Where | Description |
|------|--------|-------|-------------|
| 1 | `01_extract_and_prepare_ligand.py` | Local | Extract best docked pose, protonate |
| 2 | `02_run_acpype.sh` | Local | Generate AMBER GROMACS topology |
| 3 | `03_ligand_invacuo.sh` | Local | In-vacuo test of ligand parameters |
| 4 | `04_prepare_protein.py` | Local | Protein topology + combine coordinates |
| 5 | `05_merge_topologies.py` | Local | Merge protein + ligand topologies |
| 6 | `06_solvate_ions_index.py` | Local | Solvate, add ions, create index |
| 7 | `07_energy_minimization.py` | Local | Energy minimisation + checkpoint |
| 8 | **Transfer to HPC** | — | Copy `simulations/` to cluster |
| 9 | `09_hpc_submit.sbatch` | HPC | Equilibration + production + post-processing |

---

## Before You Start

### Requirements on your local workstation
- Singularity (or Apptainer) installed.
- The GROMACS Singularity image at `$HOME/images/gromacs_2024.5.1-GPU_20260718.sif`.
- All input files are in the **parent directory** (`Docking Files/`):
  - `CUEDC2_corrected.pdb`
  - `CHEMBL442_ERGOTAMINE_out.pdbqt`

### SSH to remote

```bash
ssh USERNAME@REMOTE.IP.ADDRESS
```

### Open a terminal and navigate to the `simulations/` directory

```bash
cd "/home/<you>/workspace/Docking Files/simulations"
```

All subsequent commands are run **from this directory**.

---
### Start an interactive singularity session

Load up the singularity shell as an interactive slurm job:

```bash
srun --partition=GPU --gres=mps:50 --ntasks=1 --cpus-per-task=16 --pty singularity shell --nv -B "${PWD}:/host_pwd" --pwd /host_pwd /path/to/singularity/image/XXX.sif
```

* This starts an interactive singularity shell with 50% of the gpu and 16 processes. 
* Adjust as needed, especially the path of the sif image.

## Step 1 — Extract and Protonate the Ligand

**Script:** `01_extract_and_prepare_ligand.py`

This script accepts command-line arguments for flexible use with any ligand.

### Default (ERGOTAMINE)

```bash
bash run_local.sh 01_extract_and_prepare_ligand.py
```

### With custom ligand

```bash
bash run_local.sh 01_extract_and_prepare_ligand.py \
    --ligand-name DIGITOXIN \
    --residue-name DIG \
    --input-pdbqt ../CHEMBL254219_DIGITOXIN_out.pdbqt \
    --ph 7.4 \
    --ligand-dir ligand
```

**Arguments:**
- `--ligand-name` — Name of the ligand for output filenames (default: `ERGOTAMINE`)
- `--residue-name` — 3-letter GROMACS residue name, e.g. `ERG`, `DIG`, `RIF` (default: `ERG`)
- `--input-pdbqt` — Path to the input PDBQT file (default: `../CHEMBL442_ERGOTAMINE_out.pdbqt`)
- `--ph` — Protonation pH (default: `7.4`)
- `--ligand-dir` — Output directory (default: `ligand`)

**What this does:**
- Reads all docked poses from the PDBQT file.
- Finds the best pose (most-negative Vina score) and reports all poses.
- Extracts the best pose to `ligand_dir/<ligand_name>_pose1.pdbqt`.
- Converts it to PDB format using OpenBabel.
- Protonates at the specified pH and saves as `ligand_dir/<RESIDUE_NAME>.pdb`.

**Expected output (for ERGOTAMINE):**
```
ligand/
  ergotamine_pose1.pdbqt
  ergotamine_raw.pdb
  ERG.pdb                 ← protonated; residue name = ERG
```

**Expected output (for DIGITOXIN):**
```
ligand/
  digitoxin_pose1.pdbqt
  digitoxin_raw.pdb
  DIG.pdb                 ← protonated; residue name = DIG
```

---

## Step 2 — Generate AMBER Topology with acpype

**Script:** `02_run_acpype.sh`

```bash
bash run_local.sh 02_run_acpype.sh --residue-name ERG
```

or for another ligand:

```bash
bash run_local.sh 02_run_acpype.sh --residue-name DIG
```

**Arguments:**
- `--residue-name` — The 3-letter residue name matching Step 1 output (default: `ERG`)

**What this does:**
- Runs `acpype` on `ligand/<RESIDUE_NAME>.pdb` with AM1-BCC charges and AMBER force field.
- Generates GROMACS-compatible topology files in `ligand/<RESIDUE_NAME>.acpype/`.

**Expected output (for ERG, inside `ligand/ERG.acpype/`):**
```
ERG_GMX.gro       ← ligand coordinates
ERG_GMX.top       ← complete AMBER topology (includes atomtypes)
ERG_GMX.itp       ← molecule-type section (for #include)
posre_ERG.itp     ← position restraints for ligand
em.mdp, nvt.mdp   ← template MDP files
```

**Expected output (for DIG, inside `ligand/DIG.acpype/`):**
```
DIG_GMX.gro
DIG_GMX.top
DIG_GMX.itp
posre_DIG.itp
```

> **If acpype crashes:** Check the error message. Common causes:
> - Charged molecule — verify `ERG.pdb` has the correct number of H atoms.
> - OpenBabel protonation artifacts — try protonating with a different tool
>   (e.g. PyMOL → Builder → Add H, then save as PDB).

---

## Step 3 — In-Vacuo Ligand Test (Recommended)

**Script:** `03_ligand_invacuo.sh`

```bash
bash run_local.sh 03_ligand_invacuo.sh
```

**What this does:**
- Runs a quick energy minimisation of the ligand in vacuum.
- Runs 100 ps of vacuum MD.
- Uses GPU for non-bonded calculations (`-nb gpu`), 1 MPI thread, all available OpenMP threads.

**Expected output (`invacuo/`):**
```
em_vac.gro, em_vac.log
md_vac.trr, md_vac.log
```

**What to check:**
- Open `invacuo/md_vac.trr` in VMD: the molecule should move freely but not explode.
- Check `invacuo/md_vac.log` for "GROMACS reminds you" at the end (= successful run).

> **If the vacuum MD crashes:** The acpype topology has errors. Do NOT proceed to the
> protein-ligand simulation until this is resolved.

---

## Step 4 — Prepare Protein Topology and Combine Coordinates

**Script:** `04_prepare_protein.py`

```bash
bash run_local.sh 04_prepare_protein.py
```

**What this does:**
- Runs `pdb2gmx` on `CUEDC2_corrected.pdb` using AMBER99SB-ILDN.
- Uses standard (charged) N and C termini: NH₃⁺ and COO⁻.
- Existing hydrogen atoms from PDBFixer are ignored (`ignh`) and re-added
  by `pdb2gmx` to ensure correct AMBER naming.
- Combines the protein coordinates (GRO) with the acpype ligand coordinates
  into a single `complex_raw/complex_raw.gro`.

**Expected output:**
```
protein/          ← pdb2gmx output (topology, GRO, ITP files)
complex_raw/
  complex_raw.gro ← protein + ligand combined
checkpoint.prepare_protein_YYYYMMDD.pycpt
```

---

## Step 5 — Merge Protein and Ligand Topologies

**Script:** `05_merge_topologies.py`

```bash
bash run_local.sh 05_merge_topologies.py
```

**What this does:**
- Reads the protein topology from `protein/topol.top`.
- Reads the acpype topology from `ligand/ERG.acpype/ERG_GMX.top`.
- Inserts the ligand `[ atomtypes ]` AFTER the force-field `#include` and
  BEFORE any `[ moleculetype ]` directive (required by GROMACS).
- Inserts `#include "ERG_GMX.itp"` before the water topology.
- Adds `ERG 1` to the `[ molecules ]` section after the protein.
- Fixes all `#include` paths to be relative to `complex/`.
- Copies all ITP files to `complex/`.
- Creates a `gromacs_py` `GmxSys` object pointing to the merged files.

**Expected output:**
```
complex/
  complex.top      ← merged topology
  complex.gro      ← copy of complex_raw.gro
  ERG_GMX.itp
  posre_ERG.itp
  (other ITP files from acpype)
checkpoint.merge_topologies_YYYYMMDD.pycpt
```

### IMPORTANT — Manual Verification Required

After this step, **open `complex/complex.top` in a text editor** and verify:

1. Near the top, after `#include "amber99sb-ildn.ff/forcefield.itp"`, you see
   a `[ atomtypes ]` block with ERG-specific atom types (e.g. `ca`, `c3`, etc.).

2. Before the water `#include`, you see:
   ```
   #include "ERG_GMX.itp"
   ```

3. In the `[ molecules ]` section at the bottom, you see **both**:
   ```
   Protein_chain_A   1
   ERG               1
   ```
   (The protein entry name may differ — that is fine.)

If anything looks wrong, fix it by hand before proceeding. The [GROMACS
topology manual](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html)
is your reference.

---

## Step 6 — Solvate, Add Ions, Create Index

**Script:** `06_solvate_ions_index.py`

```bash
bash run_local.sh 06_solvate_ions_index.py
```

**What this does:**
- Creates a **rhombic dodecahedral** simulation box (1.5 nm minimum distance
  from protein to edge).
- Solvates with TIP3P water.
- Adds Na⁺/Cl⁻ ions to neutralise the system and reach 0.15 mol/L.
- Runs `gmx make_ndx` to create an index file with two custom groups:
  - `Protein_ERG` — for temperature coupling (protein + ligand treated together)
  - `Dry` — for post-processing (same group; removes water and ions from trajectory)

**Expected output:**
```
complex_solv/
  cuedc2_erg_water_ion.gro
  cuedc2_erg_water_ion.top
  index.ndx
checkpoint.solvate_YYYYMMDD.pycpt
```

> **If the index groups are not created automatically**, run `make_ndx` manually:
> ```bash
> bash run_local.sh  # then inside singularity:
> gmx make_ndx -f complex_solv/cuedc2_erg_water_ion.gro -o complex_solv/index.ndx
> # At the prompt, type:
> #   1 | <ERG_group_number>     (merge Protein with ERG)
> #   name <N> Protein_ERG
> #   1 | <ERG_group_number>
> #   name <N+1> Dry
> #   q
> ```
> Replace `<ERG_group_number>` with the group number of ERG shown in the output.

---

## Step 7 — Energy Minimisation (Local)

**Script:** `07_energy_minimization.py`

```bash
bash run_local.sh 07_energy_minimization.py
```

**What this does:**
- Two-stage steepest-descent minimisation:
  1. Without bond constraints (removes large forces first)
  2. With bond constraints (final fine minimisation)
- Convergence criterion: maximum force < 10 kJ mol⁻¹ nm⁻¹.
- Saves energy profile plot: `em/em_energy.png`.
- Saves a pickle checkpoint for resumption on HPC.

**Expected output:**
```
em/
  em_energy.png   ← energy vs. step plot
  *.gro, *.edr, *.log (gromacs_py-generated)
checkpoint.em_YYYYMMDD.pycpt   ← TRANSFER THIS TO HPC
```

**Check the energy plot:** The potential energy should decrease monotonically
and reach a plateau. A sudden spike or failure to converge indicates a problem
with the topology or starting coordinates.

---

## Step 8 — Transfer to HPC

Copy the **entire `simulations/` directory** to the HPC cluster:

```bash
# From your local machine:
rsync -avz --progress \
    "/path/to/Docking Files/simulations/" \
    <user>@<hpc_hostname>:/path/to/simulations/
```

Or use `scp`:
```bash
scp -r "/path/to/Docking Files/simulations/" \
    <user>@<hpc_hostname>:/path/to/simulations/
```

On the HPC, verify the Singularity image is available:
```bash
ls $HOME/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif
```

---

## Step 9 — Equilibration + Production (HPC)

### Navigate to `simulations/` on the HPC

```bash
cd /path/to/simulations
```

### Submit the SLURM job

```bash
sbatch 09_hpc_submit.sbatch
```

**Monitor the job:**
```bash
squeue -u $USER               # check job status
tail -f hpc_run_*.log         # follow the log file in real time
```

**What happens inside:**
1. Loads `checkpoint.em_YYYYMMDD.pycpt`.
2. **Equilibration** (three stages, ~5.5 ns total):
   - Stage 1 (0.5 ns): heavy-atom position restraints, dt = 0.001 ps
   - Stage 2 (1.0 ns): Cα position restraints, dt = 0.002 ps
   - Stage 3 (4.0 ns): weak Cα restraints, dt = 0.002 ps
3. **Production MD** (100 ns, dt = 0.002 ps, 50 million steps).
4. **Post-processing:**
   - Trajectory concatenation
   - Make molecules whole (PBC repair)
   - Center protein+ligand in box
   - Remove PBC jumps
   - Extract dry trajectory (protein + ERG only) → `prod/traj_dry.xtc`
   - Extract dry TPR → `prod/traj_dry.tpr`

**Checkpoints saved:**
```
checkpoint.equi_YYYYMMDD.pycpt
checkpoint.prod_YYYYMMDD.pycpt
```

---

## Step 10 — Inspect Results

### Check simulation health

```bash
# Inside Singularity — check equilibration energy
singularity exec \
    -B "${PWD}:/host_pwd" --pwd /host_pwd \
    $HOME/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif \
    gmx energy -f equi/<last_stage>.edr

# Check production energy
singularity exec \
    -B "${PWD}:/host_pwd" --pwd /host_pwd \
    $HOME/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif \
    gmx energy -f prod/md.edr
# Select: Temperature Pressure Volume Density (then 0 to quit)
```

### Visualise trajectory in VMD

Transfer `prod/traj_dry.xtc` and `prod/traj_dry.tpr` (or the final GRO) to your
local machine:

```bash
scp <user>@<hpc>:/path/to/simulations/prod/traj_dry.* .
```

Open in VMD:
1. File → New Molecule → select `traj_dry.tpr` (structure file)
2. File → Load Data Into Molecule → select `traj_dry.xtc` (trajectory)
3. Use the VMD trajectory player to inspect the simulation.

### Compute RMSD (optional)

```bash
singularity exec -B "${PWD}:/host_pwd" --pwd /host_pwd \
    $HOME/SIFDIR/gromacs/gromacs_2024.5.1-GPU_20260718.sif \
    gmx rms \
        -s prod/traj_dry.tpr \
        -f prod/traj_dry.xtc \
        -o prod/rmsd_backbone.xvg \
        -n complex_solv/index.ndx
# Select Backbone for fit, Backbone for RMSD
```

---

## Adapting for Other Ligands

All scripts now accept command-line arguments for ligand name, residue name, and input files.
To run the same protocol with a different ligand (DIGITOXIN, PARITAPREVIR, RIFAMYCIN, or MIDOSTAURIN):

### Quick example: Run full workflow for DIGITOXIN

```bash
# Step 1 — Extract and protonate
bash run_local.sh 01_extract_and_prepare_ligand.py \
    --ligand-name DIGITOXIN \
    --residue-name DIG \
    --input-pdbqt ../CHEMBL254219_DIGITOXIN_out.pdbqt

# Step 2 — Generate topology
bash run_local.sh 02_run_acpype.sh --residue-name DIG

# Step 3 — In-vacuo test
bash run_local.sh 03_ligand_invacuo.sh --residue-name DIG

# Step 4 onwards — use same residue name in arguments/topologies
bash run_local.sh 04_prepare_protein.py --residue-name DIG
bash run_local.sh 05_merge_topologies.py --residue-name DIG
bash run_local.sh 06_solvate_ions_index.py --residue-name DIG
bash run_local.sh 07_energy_minimization.py --residue-name DIG
```

### Ligand residue name reference

Choose a **unique 3-letter code** for each ligand (GROMACS requirement):

| Ligand | Residue Code | PDBQT File | Notes |
|--------|--------------|-----------|-------|
| ERGOTAMINE | `ERG` | `CHEMBL442_ERGOTAMINE_out.pdbqt` | Default |
| DIGITOXIN | `DIG` | `CHEMBL254219_DIGITOXIN_out.pdbqt` | |
| PARITAPREVIR | `PAR` | `CHEMBL3391662_PARITAPREVIR_out.pdbqt` | |
| RIFAMYCIN | `RIF` | `CHEMBL437765_RIFAMYCIN_out.pdbqt` | |
| MIDOSTAURIN | `MID` | `CHEMBL608533_MIDOSTAURIN_out.pdbqt` | |

### Full script parameter reference

**01_extract_and_prepare_ligand.py:**
```bash
--ligand-name NAME          # Name of the ligand (default: ERGOTAMINE)
--residue-name CODE         # 3-letter GROMACS residue code (default: ERG)
--input-pdbqt PATH          # Path to PDBQT file
--ph VALUE                  # Protonation pH (default: 7.4)
--ligand-dir DIR            # Output directory (default: ligand)
```

**02_run_acpype.sh:**
```bash
--residue-name CODE         # 3-letter residue code matching Step 1 (default: ERG)
--ligand-dir DIR            # Ligand directory (default: ligand)
```

**03_ligand_invacuo.sh:**
```bash
--residue-name CODE         # 3-letter residue code (default: ERG)
--ligand-dir DIR            # Ligand directory (default: ligand)
--invacuo-dir DIR           # Output directory (default: invacuo)
```

**04_prepare_protein.py, 05_merge_topologies.py, 06_solvate_ions_index.py, 07_energy_minimization.py:**
```bash
--residue-name CODE         # 3-letter residue code (default: ERG)
```

**08_hpc_equi_prod.py:**
```bash
--residue-name CODE         # 3-letter residue code (default: ERG)
```

---

## Troubleshooting

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `acpype` crashes with "Unusual charge" | Wrong formal charge | Verify protonation state; try `-n 0` or `-n 1` for neutral/cation |
| `pdb2gmx` fails with "Atom not found" | Non-standard residue / wrong PDB | Use `ignh=True`; inspect PDB for HET atoms |
| "Fatal error: Invalid order for directive atomtypes" | Ligand atomtypes missing at top | Check topology merge (Step 5) |
| EM does not converge | Clash in starting structure | Try more EM steps; check acpype GRO |
| Production MD crashes immediately | Index group mismatch | Verify `Protein_ERG` group in `index.ndx` matches `TC_GRPS` |
| `traj_dry.xtc` not created | `Dry` index group missing | Run `make_ndx` manually; create `Dry = Protein_ERG` |

---

## File Structure Reference

After all steps complete locally, `simulations/` should look like:

```
simulations/
  run_local.sh
  01_extract_and_prepare_ligand.py
  02_run_acpype.sh
  03_ligand_invacuo.sh
  04_prepare_protein.py
  05_merge_topologies.py
  06_solvate_ions_index.py
  07_energy_minimization.py
  08_hpc_equi_prod.py
  09_hpc_submit.sbatch
  INSTRUCTIONS.md          ← this file
  ligand/
    ERG.pdb
    ERG.acpype/
  invacuo/                 ← in-vacuo test results
  protein/                 ← pdb2gmx output
  complex_raw/             ← combined GRO (before merge)
  complex/                 ← merged topology + GRO
  complex_solv/            ← solvated system + index
  em/                      ← energy minimisation
  checkpoint.em_YYYYMMDD.pycpt    ← transfer this to HPC
```

After HPC runs:

```
simulations/
  equi/              ← equilibration output
  prod/
    md.gro, md.xtc, md.edr, md.log
    traj_dry.xtc     ← dry production trajectory
    traj_dry.tpr     ← dry topology for VMD
  checkpoint.equi_YYYYMMDD.pycpt
  checkpoint.prod_YYYYMMDD.pycpt
```
