# AMBER Topology Build for PDL — Pd(II) N₄ Square-Planar Complex

## 🚀 Palladium-Ligand Topology & Simulation Workflow
This guide details the high-fidelity protocol for building a classical bonded model of a Palladium-ligand complex for GROMACS. This workflow uses Gaussian 16, Multiwfn, and acpype to create a custom-parameterized force field.

## Ligand Summary

| Property | Value |
|---|---|
| Residue name | PDL |
| Metal | Pd(II), square-planar d⁸, diamagnetic |
| Coordination | 4× N at 2.01–2.03 Å — bis(acylhydrazone) motif |
| Coordinating atoms | Two imine N (N8, N9) + two hydrazone N (N20, N21) |
| Total charge | 0 (Pd²⁺ balanced by two anionic hydrazone arms) |
| Spin multiplicity | 1 (singlet) |
| Heavy atoms after prep | 41 (47 original − 6 AutoDock dummy atoms) |
| Total atoms with H | 65 (41 heavy + 24 H: 4×CH₃ + 8 aryl CH + 2×NH₂) |

## Software Versions

| Tool | Version | License | Download |
|---|---|---|---|
| acpype|  2022.6.6 or later | FOSS | [https://alanwilter.github.io/acpype/](https://alanwilter.github.io/acpype/)
| Gaussian 16 | Rev A.03 or later| Proprietary |[https://gaussian.com/gaussian16/](https://gaussian.com/gaussian16/) |
| GROMACS | 2023 or later | FOSS | [https://anaconda.org/channels/conda-forge/packages/gromacs/overview](https://anaconda.org/channels/conda-forge/packages/gromacs/overview)
| Python | ≥ 3.8 | FOSS| [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge) |
| Multiwfn |2026.4.10 or later | Free | [http://sobereva.com/multiwfn/](http://sobereva.com/multiwfn/)|

---

## Prelude: Stage 0 — PDB Preparation

### 0a. Strip AutoDock Dummy Atoms

```bash
python3 00_prep_ligand.py best_docked_ligandonly.pdb   # → PDL_prep.pdb (41 atoms)
```

The `AutoDockZn` docking output contains six non-physical atoms that must be removed. The script identifies them by their sub-Ångström inter-element distances and removes them, then rebuilds element columns, unique AMBER-compatible atom names, and CONECT records.

| Serials | Element | Distance to parent | Interpretation |
|---|---|---|---|
| 28, 29 | `BH` | 0.86 Å to N8, N9 | AutoDockZn pseudo-boron encoding N-donor directionality |
| 37, 38 | `H` | 0.86 Å to N36 | AutoDock pseudo-H on terminal NH₂ |
| 46, 47 | `H` | 0.86 Å to N45 | AutoDock pseudo-H on terminal NH₂ |

The stripping criterion is: element is non-standard (`BH`, `DU`, `EP`, etc.) **or** element is H with an inter-element distance below 0.95 Å. Real X–H covalent bonds are always ≥ 0.96 Å. **Do not run `00_prep_ligand.py` on a file that already has real H atoms added** — it will strip N–H bonds at 1.03 Å because they fall below the 1.05 Å outer threshold.

Output: `PDL_prep.pdb`, 41 atoms, Pd at serial 1.

### 0b. Add Hydrogen Atoms

Automated protonation tools all fail for this system:

- `reduce` operates from a compiled residue dictionary and adds zero H atoms to non-standard `HETATM` residues, reporting `Added 0 hydrogens (0 hets)`.
- `obabel --gen3d` and `RDKit EmbedMolecule` regenerate all 3D coordinates from scratch using distance geometry, destroying the docked pose entirely.

**Protonation must be performed in Avogadro** (or equivalent visualiser). After adding H atoms, inspect the structure carefully against the known chemistry before proceeding.

The 24 correct H atoms for this ligand are:

| Group | Count | Atoms | Notes |
|---|---|---|---|
| 4 × CH₃ (methyl on N-amino arms) | 12 H | C24–C27 | sp3, 3H each |
| Aromatic ring CH (ring 1: C29, C30, C32, C33) | 4 H | para-aminophenyl | 1H per unsubstituted aromatic C |
| Aromatic ring CH (ring 2: C36, C37, C39, C40) | 4 H | para-aminophenyl | 1H per unsubstituted aromatic C |
| 2 × NH₂ (terminal aniline) | 4 H | N34, N41 | sp3, 2H each |

The following positions must carry **no hydrogen** in the neutral deprotonated complex:

- **Imine carbons** C2, C3 (C=N at 1.274 Å), **hydrazone carbons** C6, C7 (C=N at 1.298 Å), **carbonyl carbons** C4, C5, C22, C23 (C=O at 1.21–1.25 Å) — all sp2, zero H.
- **Imine nitrogen** N8, N9 — these coordinate Pd and carry no H.
- **Hydrazone alpha-nitrogen** N10, N11 (the N of the N–N= unit) — no H.
- **Coordinating hydrazone nitrogen** N20, N21 — these are the *deprotonated* donor atoms that give the ligand its –2 charge; they must carry no H.
- **Carbonyl oxygens** O16–O19 — C=O, no H.
- **Aromatic substituted carbons** C28, C31 (ring 1) and C35, C38 (ring 2) each carry three heavy-atom bonds and must have zero H.

> **Critical chemistry note.** The complex is neutral because the two hydrazone arms are *deprotonated*: each anionic N–H bond is lost upon coordination, contributing a formal charge of −1 per arm to balance Pd²⁺. Protonating N20 or N21 changes the total charge to +2 and completely alters the RESP charge distribution. The anomalous −4.17 B-factor value on N11 in the original AutoDock PDB is a scoring artefact and carries no chemical meaning.

Output: `PDL_prep_H.pdb`, 65 atoms (41 heavy + 24 H).

### 0c. Structure Correction Script

If the protonation step introduced any of the forbidden H atoms listed above, run the structure correction script:

```bash
# Requires PDL_prep_H_corrected.pdb in the working directory
bash fix_structure_v2.sh
```

This script bypasses `00_prep_ligand.py` (which must not be called on an H-added file) and directly constructs `PDL_final_v2.pdb` by splitting the 65-atom structure into ZN (metal) and PDL (ligand) residues with correct PDB column alignment, then runs antechamber and patches carbonyl O types.

#### User Note
Avogadro just put an extra H on N10, so I manually got rid of it in `PyMol` and saved as `PDL_prep_H_corrected.pdb` in the working directory. 


### 📋 0d. Prerequisites & Installation## Multiwfn Installation

   1. Download: Get the Linux binary from [sobereva.com](http://sobereva.com/multiwfn/).
   2. Setup: Extract the folder and add export Path_to_Multiwfn/Multiwfn to your .bashrc.
   3. Memory Management (CRITICAL):
   To prevent "Segmentation Faults" during ESP calculation for large systems or metals, always increase the stack size limit in your terminal before launching:
  ```bash
   ulimit -s unlimited
  ```
   
   

## 🛠 Step 1: Geometry Optimization (Gaussian + SLURM)

   1. Load full exact protonated pdb in avogadro2 and Click Input --> Gaussian. 
   2. Keywords: Use B3LYP/def2SVP with an Effective Core Potential (ECP) for Pd.
 
 ```bash
 %Chk=PDL_small_opt.chk
 %CPU=
 %Mem=
 # B3LYP/def2SVP Integral=(Grid=UltraFine) NoSymm
 Opt=(MaxCycles=500,MaxStep=10,CalcFC)
 SCF=(MaxCycles=256,XQC,NoVarAcc)

 PDL small model B3LYP/def2SVP optimisation

 0  1
 [Atom Co-ordinates. Can export from avogadro and edit the above config ...]
 ...
  ```
 
   3. SLURM Batch Submission: Use a .sh script to request enough memory for the Hessian calculation.
```bash
#!/bin/bash
#SBATCH --job-name=g16_smopt
#SBATCH --output=logs/g16_sm_%j.out
#SBATCH --error=logs/g16_sm_%j.err
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --partition=CPU
#SBATCH --qos=elevated

mkdir -p logs
echo "Starting Gaussian jobs"
# 2. Gaussian Environment (Source without strict mode)
source /etc/g16setup

# 3. Scratch Directory Setup
export GAUSS_SCRDIR=$(pwd)/scratch_g16
mkdir -p "$GAUSS_SCRDIR"

# 4. Correct Filenames (Ensuring they match your 'ls')
SM_OPT="PDL_small_opt.gjf"
CPUS=$(taskset -cp $$ | awk -F':' '{print $2}')

echo "[$(date)] Starting Gaussian Job"
echo "Using File: $SM_OPT"

# 5. Patch Files in-place
sed -i -e "/^.*%cpu/I s/=.*$/=$CPUS/" "$SM_OPT"
sed -i "s|%Mem=.*|%Mem=32GB|" "$SM_OPT"
echo "Patched $SM_OPT"


echo "[$(date)] Small model opt restart — MaxStep=5, MaxCycles=500"
echo "  Checkpoint: PDL_small_opt.chk"
echo "  GAUSS_SCRDIR: $GAUSS_SCRDIR"

# ── Sanity checks ────────────────────────────────────────────────────────────
if [[ ! -f PDL_small_opt.chk ]]; then
    echo "[WARNING] PDL_small_opt.chk not found."
    echo "        This job must continue from the existing checkpoint."
fi
if [[ ! -f PDL_small_opt.gjf ]]; then
    echo "[WARNING] PDL_small_opt.gjf not found."
    exit 1
fi

# ── Run optimisation ─────────────────────────────────────────────────────────
echo "[$(date)] Running g16 < PDL_small_opt.gjf ..."
g16 < PDL_small_opt.gjf > logs/PDL_small_opt2.log 2>&1

# ── Check result ─────────────────────────────────────────────────────────────
echo ""
if grep -q "Normal termination" logs/PDL_small_opt2.log; then
    echo "[$(date)] [OK] Optimisation converged — Normal termination."
    FINAL_STEP=$(grep "Step number" logs/PDL_small_opt2.log | tail -1)
    FINAL_E=$(grep "SCF Done" logs/PDL_small_opt2.log | tail -1 | awk '{print $5}')
    echo "  $FINAL_STEP"
    echo "  Final SCF energy: $FINAL_E Ha"
else
    LAST_STEP=$(grep "Step number" logs/PDL_small_opt2.log | tail -1)
    LAST_CONV=$(grep -A5 "Item.*Value.*Threshold" logs/PDL_small_opt2.log | tail -8)
    echo "[$(date)] [WARN] Did not reach Normal termination."
    echo "  Last: $LAST_STEP"
    echo "  Last convergence table:"
    echo "$LAST_CONV"
    echo ""
    echo "  If oscillating: reduce MaxStep further (try MaxStep=3 in .gjf)"
    echo "  If still descending: increase walltime and resubmit from checkpoint"
fi

rm -rf "$GAUSS_SCRDIR"
echo "[$(date)] Done."
```


## 🧪 Step 2: RESP Charge Calculation (Multiwfn)

   1. Convert gaussian chk file to formatted chk
 
  ```bash
   formchk PDL_small_opt.chk PDL_small_opt.fchk
  ```
   
   2. Fitting:
      * Launch: Multiwfn ligand_resp.fchk (after ulimit -s unlimited).
      * Path: 7 → 18 → 1.
      * Pd Radius: When prompted, input 1.63 Å.
   3. Output: Save the results as PDL_small_opt.chg.


## 🏗 Step 3: Topology & Coordinate Automation
## 3.1. Generate the Organic Scaffold
Remove the Pd atom from your optimized PDB (manually) and run acpype:

```bash
acpype -i organic_only.pdb -c bcc -n 0 -f
```

## 3.2. Automated Topology Grafting
Use this script to inject the Palladium, shift all 64 organic atom indices by +1, and apply RESP charges.

```python
import re

def sync_topol_and_gro(chg_file, original_itp):
    # 1. Load data from .chg
    with open(chg_file, 'r') as f:
        data = [line.split() for line in f if line.strip()]
    
    # 2. Extract organic framework from original ITP
    with open(original_itp, 'r') as f:
        itp_content = f.read()
    
    # Extract original atom types and interactions
    atom_types = re.findall(r'\[ atoms \](.*?)\n\[', itp_content, re.S)[0].strip().split('\n')
    organic_types = [line.split()[1] for line in atom_types if line.strip() and not line.strip().startswith(';')]
    
    # 3. Create the NEW .gro file (nm units)
    with open("PDL_final.gro", "w") as f:
        f.write("PDL optimized complex\n65\n")
        for i, (parts) in enumerate(data):
            symbol = "PD" if i == 0 else f"{parts[0]}{i+1}" # Unique names: PD, C2, C3...
            x, y, z = float(parts[1])/10, float(parts[2])/10, float(parts[3])/10
            f.write(f"{1:5d}{'PDL':5s}{symbol:>5s}{i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(f"{10.0:10.5f}{10.0:10.5f}{10.0:10.5f}\n")

    # 4. Create the NEW .itp file
    with open("PDL_final.itp", "w") as f:
        f.write("[ atomtypes ]\n")
        # Add all types used in the file
        f.write("Pd   46   106.420   0.000   A   0.24500   0.21338\n")
        # Extract existing types from original ITP to be safe
        types_section = re.findall(r'\[ atomtypes \](.*?)\n\[', itp_content, re.S)[0]
        f.write(types_section.strip() + "\n\n")

        f.write("[ moleculetype ]\nPDL   3\n\n[ atoms ]\n")
        for i, parts in enumerate(data):
            symbol = "PD" if i == 0 else f"{parts[0]}{i+1}"
            charge = parts[4]
            at_type = "Pd" if i == 0 else organic_types[i-1]
            mass = 106.42 if i == 0 else (12.01 if 'c' in at_type else (14.01 if 'n' in at_type else 1.008))
            f.write(f"{i+1:6d} {at_type:5s} 1  PDL  {symbol:5s} {i+1:6d} {float(charge):10.6f} {mass:8.3f}\n")
        
        # Shift all interaction indices by +1
        f.write("\n[ bonds ]\n1 8 1 0.20131 250000\n1 9 1 0.20146 250000\n1 20 1 0.20297 250000\n1 21 1 0.20299 250000\n")
        
        # Simplified logic for remaining sections
        for section in ["bonds", "pairs", "angles", "dihedrals"]:
            if section != "bonds": f.write(f"\n[ {section} ]\n")
            lines = re.findall(rf'\[ {section} \](.*?)\n\[', itp_content + "\n[", re.S)[0].strip().split('\n')
            n = 4 if section == "dihedrals" else (3 if section == "angles" else 2)
            for line in lines:
                if not line.strip() or line.strip().startswith(';'): continue
                p = line.split()
                shifted = [str(int(x) + 1) if j < n else x for j, x in enumerate(p)]
                f.write("  " + "  ".join(shifted) + "\n")

    print("✅ Created PDL_final.gro and PDL_final.itp. They are now perfectly synced.")

sync_topol_and_gro("PDL_small_opt.chg", "PDL_H_optimized_nopd.acpype/PDL_H_optimized_nopd_GMX.itp")

```


## 🔗 Step 4: Position Restraints (posre_PDL_GMX.itp)
To ensure the ligand (including the metal) can be equilibrated without flying away, generate a position restraint file specifically for the PDL residue.

   1. Create the file:
   ```bash
   gmx genrestr -f ligand_final.gro -o posre_PDL_GMX.itp -fc 1000 1000 1000
   ``` 
 
   2. Integration: Ensure your PDL_final.itp includes the conditional block at the very end:
   3. 
   ```bash
   ; Ligand position restraints
   #ifdef POSRES_LIG
   #include "posre_PDL_GMX.itp"
   #endif
   ```
   

## 🚀 Step 6: Launching the Simulation of ligand in-vacuo

   1. Grompp:
   ```bash 
   gmx grompp -f em.mdp -c boxed.gro -p topol.top -o em.tpr -v
   ```
   
   2. GPU Run:
   ```bash
   gmx mdrun -v -deffnm em -nt 12 -nb gpu -pme gpu -pmefft gpu
   ```

## TODO

Try using FOSS QM software such as `ORCA` or `GAMESS` instead of `Gaussian16` for the optimization.