# Transition Metal Complex Parametrization Pipeline


This guide outlines the successful workflow for parametrizing a small transition metal complex (**Nickel tetracarbonyl, [Ni(CO)₄]**) using a quantum mechanics (QM) driven force field pipeline. The pipeline transitions from **ORCA** (QM) to **easyPARM** (Seminario-based force constants) to **AMBER/LEaP** (topology building), **ParmEd** (conversion to GROMACS), and finally to **GROMACS Production Runs**.

---

## Workflow Overview

[1. ORCA QM Run] ──> [2. easyPARM] ──> [3. AMBER tleap] ──> [4. ParmEd] ──> [5. GROMACS Ready]

(Opt + Freq + (Seminario Bond/ (Builds Amber (Converts to (Topology &
CHELPG Charges) Angle Constants) Topology Files) .top & .gro) Coordinates)


---

## Step 1: Quantum Mechanics (ORCA)

We used a multi-step ORCA job (`$new_job`) to first optimize the geometry using a pure GGA functional (`PBE`) for speed, followed by a frequency and electrostatic potential charge (`CHELPG`) calculation on the optimized structure.

### Input File (`ni_co4.inp`)
```bash
# STEP 1: GEOMETRY OPTIMIZATION
! PBE D3BJ def2-SVP def2/J opt
%pal nprocs 12 end
%maxcore 2000
* xyz 0 1
Ni       0.00000000      0.00000000      0.00000000
C        1.05000000      1.05000000      1.05000000
O        1.75000000      1.75000000      1.75000000
C       -1.05000000     -1.05000000      1.05000000
O       -1.75000000     -1.75000000      1.75000000
C        1.05000000     -1.05000000     -1.05000000
O        1.75000000     -1.75000000     -1.05000000
C       -1.05000000      1.05000000     -1.05000000
O       -1.75000000      1.75000000     -1.05000000
*
```

### Execution Command
```bash
orca ni_co4.inp > ni_co4.out
```
*Critical Outputs Generated:* `ni_co4.xyz` (optimized structure), 

### Input file (`ni_co4_freq.inp`)
```bash 
# STEP 2: FREQUENCY CALCULATION
! PBE D3BJ def2-SVP def2/J freq
%pal nprocs 12 end
%scf MaxIter 1000 end
! CHELPG
%maxcore 2000
* xyzfile 0 1 ni_co4.xyz

```
### Execution Command
```bash
orca ni_co4_freq.inp > ni_co4_freq.out
```

*Critical Outputs Generated:* `ni_co4_freq.out` (charges), and `ni_co4_freq.hess` (vibrational Hessian matrix).

---

## Step 2: Parameter Generation (easyPARM)

We executed easyPARM in **Interactive Mode** (Option 1: *Generate molecular complex parameters*) to map GAFF2 atom types and generate custom metal force field parameters via the Seminario method.

### Inputs Provided to easyPARM Prompts:
*   **Structure File:** `ni_co4.xyz`
*   **Charge Method:** `CHELPG (ORCA)`
*   **Charge Output Source:** `ni_co4_freq.out`
*   **QM Frequency Format:** `Orca Output`
*   **Hessian Source File:** `ni_co4_freq.hess`
*   **Atom Type Strategy:** `GAFF2`

*Outputs Generated:* `COMPLEX.frcmod` (custom bond/angle constants), `COMPLEX.lib` (off-standard residue definitions), and `tleap.input` (automated build script template).

---

## Step 3: Topology Generation (AMBER LEaP)

Before running LEaP, the placeholder string `Your_Whole_System.pdb` inside the automatically generated script was updated to match our actual molecule system coordinates (`COMPLEX.pdb` or `NIC.pdb`).

### Script Correction Commands
```bash
sed -i 's/Your_Whole_System.pdb/COMPLEX.pdb/g' tleap.input
sed -i 's/Your_Whole_System.prmtop/NIC.prmtop/g' tleap.input
sed -i 's/Your_Whole_System.inpcrd/NIC.inpcrd/g' tleap.input
```

### Execution Command
```bash
tleap -f tleap.input
```
*Outputs Generated:* `NIC.prmtop` (Amber Topology) and `NIC.inpcrd` (Amber Coordinates).

---

## Step 4: GROMACS Conversion (ParmEd)

We bypassed the interactive ParmEd CLI by invoking its underlying API layer via a quick Python execution hook. This cleanly translated the metal coordination bonding terms into GROMACS native architectures.

### Conversion Command
```bash
python -c 'import parmed as pmd; parm = pmd.load_file("NIC.prmtop", "NIC.inpcrd"); parm.save("NIC.top"); parm.save("NIC.gro")'
```
*Outputs Generated:* `NIC.top` (GROMACS System Topology) and `NIC.gro` (GROMACS Structure Coordinates).

---

## Step 5: GROMACS MD Parameter Settings (MDP)

Save these parameter blocks to your directory to control the GROMACS execution engine framework.

### 1. Energy Minimization Input (`em.mdp`)
```ini
; Parameters describing what to do, when to stop and what to save
integrator      = steep     ; Algorithm (steep = steepest descent minimization)
nsteps          = 500       ; Maximum number of (minimization) steps to perform
nstxout         = 10

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1             ; Frequency to update the neighbour list and long range forces
cutoff-scheme   = Verlet
rlist           = 1.2           ; Cut-off for making neighbour list (short range forces)
coulombtype     = PME           ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2           ; long range electrostatic cut-off
vdw-type        = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2           ; long range Van der Waals cut-off
pbc             = xyz           ; Periodic Boundary Conditions
DispCorr        = no
```

### 2. Production MD Input (`md.mdp`)
```ini
integrator               = md
nsteps                   = 10000
nstxout                  = 10
cutoff-scheme            = verlet
coulombtype              = PME
constraints              = h-bonds
vdwtype                  = cutoff
vdw-modifier             = force-switch
rlist                    = 1.0
rvdw                     = 1.0
rvdw-switch              = 0.9
rcoulomb                 = 1.1
DispCorr                 = EnerPres
lincs-iter               = 2
fourierspacing           = 0.25
gen-vel                  = yes
```

---

## Step 6: GROMACS Production Pipeline

With the native GROMACS configurations ready, use the following production workflow to define the box shape, minimize the system, and execute the Molecular Dynamics run utilizing hardware acceleration.

### 1. Box and Centering
```bash
echo 0 | gmx editconf -f NIC.gro -bt octahedron -d 1.5 -c -princ -o boxed.gro
```

### 2. Energy Minimization (EM)
```bash
gmx grompp -f em.mdp -c boxed.gro -p NIC.top -o em.tpr -v
# Run EM using 6 available threads
gmx mdrun -nt 6 -v -deffnm em
```

### 3. Production MD (GPU PME Offload)
```bash
gmx grompp -f md.mdp -c em.gro -p NIC.top -o md.tpr -r em.gro
# -nt 6: Use 6 CPU threads
# -nb gpu: Non-bonded interactions on GPU
# -pme gpu: PME calculations on GPU
# -pmefft gpu: PME FFT on GPU (if supported by your GMX version)
gmx mdrun -nt 6 -nb gpu -pme gpu -pmefft gpu -v -deffnm md
```

---

## Step 7: Visualization & Troubleshooting (VMD)

Because GROMACS `.gro` files lack explicit bond information, VMD will initially hide the long metal coordination bonds (`Ni-C`, ~1.78 Å) since they fall outside default organic distance thresholds. The underlying physics inside your simulation topology (`NIC.top`) remain intact.

### Manual Bonding in VMD Tk Console
To visually draw the coordination bonds inside VMD, open **Extensions ➔ Tk Console** and map the 0-based atom indices:
```tcl
topo addbond 0 1   ;# Ni1 to C1
topo addbond 0 3   ;# Ni1 to C2
topo addbond 0 5   ;# Ni1 to C3
topo addbond 0 7   ;# Ni1 to C4
```

### Resetting or Clearing Custom Bonds
If you ever need to purge these manual coordinates to rewrite them for a solvated system box, use the cleanup utilities inside the console window:
```tcl
# Clear all current bonds
topo clearbonds

# Recalculate VMD's standard organic guesses
vmd_calculate_bonds 0
```
*Click **Apply** in the Graphics ➔ Representations screen to force visual viewport rendering updates.*