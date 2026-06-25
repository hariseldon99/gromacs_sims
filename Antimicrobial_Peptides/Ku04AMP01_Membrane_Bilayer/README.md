# Simulation of KU04AMP01 Peptide with Inner Cell Membrane of Gram-Negative Bacteria

## Overview

This README documents the workflow for simulating the antimicrobial peptide **KU04AMP01**
with the inner membrane (**G-IM**) of Gram-negative bacteria using GROMACS with the
CHARMM36 force field.

**KU04AMP01** is a 36-residue linear cationic AMP (structure: [KU04AMP01.pdb](KU04AMP01.pdb)):

```
D-N-K-A-K-S-K-K-R-D-K-E-K-P-S-S-G-R-P-G-Q-T-N-S-V-P-N-A-A-I-Q-V-Y-K-E-D
```

Net charge ≈ **+4** at physiological pH (7 Lys + 2 Arg − 3 Asp − 2 Glu, plus N/C termini).
All residues are standard amino acids, so topology generation uses `pdb2gmx` directly —
no CGenFF or Ligand Reader is needed.

---

## 1. G-IM Membrane Topology

The equilibrated G-IM model is stored in `GIM_topol/`:

- [`GIM_topol/last.pdb`](GIM_topol/last.pdb) — final snapshot of the CHARMM step7 production run (use this as the starting membrane)
- `GIM_topol/G-IM.tar.gz` — full simulation archive (CHARMM format: `.inp`, `.psf`, `.crd`, `.dcd`)
- `GIM_topol/tarball/` — extracted archive contents

Topology files can also be downloaded from the [CHARMM-GUI Biomembrane Library](https://charmm-gui.org/?doc=archive&lib=biomembrane):

- PDB of last snapshot: <https://charmm-gui.org/archive/biomembrane/G-IM/last.pdb>
- Full archive: <https://charmm-gui.org/archive/biomembrane/G-IM/G-IM.tar.gz>

### 1.1 G-IM System Properties

Size (Å³): 79.79 × 79.79 × 100  
\# water / Na⁺ / Cl⁻: 11550 / 76 / 30  
Temperature (K): 310.15  
**Lipid composition** (symmetric bilayer; counts are per leaflet):

| Lipid Name | Lipid Head/Tail | # Lipids per Leaflet | APL (Å²) | Charge (e) |
|---|---|---|---|---|
| PMPE | PE (16:0/cy17:0) | 46 | 62.25 | 0 |
| POPE | PE (16:0/18:1 9Z) | 13 | 62.18 | 0 |
| QMPE | PE (15:0/cy17:0) | 12 | 62.19 | 0 |
| OYPE | PE (18:1 9Z/16:1 9Z) | 8 | 62.35 | 0 |
| PMPG | PG (16:0/cy17:0) | 10 | 64.22 | −1 |
| PYPG | PG (16:0/16:1 9Z) | 9 | 64.04 | −1 |
| PVCL2 | CL (1′-[16:0/18:1 11Z],3′-[16:0/18:1 11Z]) | 2 | 127.43 | −2 |
| **TOTAL** | | **100** | | |

The cy17:0 (cyclopropane) chains in PMPE, QMPE, and PMPG are characteristic of
*E. coli* inner membrane. Total net lipid charge per leaflet:
10 × (−1) + 2 × (−2) = −14 e (both leaflets: −28 e), partially neutralised
by the pre-existing 76 Na⁺ / 30 Cl⁻.

### 1.2 G-IM CHARMM Archive Stages

The archive contains CHARMM-format files only; there are no GROMACS input files. For
reference:

| Stage | Key files | What it represents |
|---|---|---|
| step3 | `step3_packing.pdb`, `step3_packing.inp` | Lipid packing |
| step4 | `step4_lipid.psf`, `step4.2_waterbox.*`, `step4.3_ion.*` | Bilayer assembly, solvation, ion placement |
| step5 | `step5_assembly.pdb/.psf/.crd` | Fully assembled system before MD |
| step6.1–6.6 | `step6.x_equilibration.inp` | Six-stage NVT/NPT equilibration in CHARMM |
| step7 | `step7_production.inp`, `step7_1000.dcd` | Production MD, 250,000 steps at 310.15 K |
| **last snapshot** | [`GIM_topol/last.pdb`](GIM_topol/last.pdb) | **Final equilibrated frame — start here** |

---

## 2. KU04AMP01 Peptide Topology

KU04AMP01 consists entirely of standard amino acid residues. Its GROMACS topology is
generated directly with `pdb2gmx` — no Ligand Reader or CGenFF is needed.

### 2.1 Generate topology with `pdb2gmx`

```bash
gmx pdb2gmx \
  -f KU04AMP01.pdb \
  -o ku04amp01.gro \
  -p ku04amp01.top \
  -ff charmm36-jul2022 \
  -water tip3p \
  -ignh
```

| Output file | Contents |
|---|---|
| `ku04amp01.gro` | Peptide coordinates in GROMACS GRO format |
| `ku04amp01.top` | Standalone topology with `[ moleculetype ]` |
| `posre.itp` | Position restraint file |

The `-ignh` flag discards the PDB hydrogens and lets `pdb2gmx` re-add them with
CHARMM36 naming conventions. Verify the charge:

```bash
grep "qtot" ku04amp01.top | tail -1
# Expected: qtot ≈ +4
```

### 2.2 Extract the peptide ITP

For merging with the membrane topology, extract the `[ moleculetype ]` section from
`ku04amp01.top` into a standalone `KU04AMP01.itp`:

```bash
python3 - <<'EOF'
import re
with open("ku04amp01.top") as f:
    content = f.read()
match = re.search(r'(\[ moleculetype \].*?)(?=\[ system \])', content, re.DOTALL)
if match:
    with open("KU04AMP01.itp", "w") as out:
        out.write(match.group(1))
    print("KU04AMP01.itp written")
else:
    print("ERROR: could not find [ moleculetype ] section")
EOF
```

---

## 3. Build G-IM in CHARMM-GUI Bilayer Builder with GROMACS Output

The preferred approach is to rebuild the G-IM bilayer fresh in CHARMM-GUI Bilayer Builder
requesting GROMACS output. This correctly handles the CHARMM→GROMACS conversion for the
cyclopropane-chain lipids and generates the staged equilibration MDP files.

> **Note on uploading `last.pdb` directly:** Using the "Add Molecules to Bilayer System"
> sub-option and uploading `GIM_topol/last.pdb` would cause CHARMM-GUI to run CGenFF on
> every residue in the PDB. This fails for the non-standard lipid types. Build a fresh
> bilayer instead.

Go to [CHARMM-GUI Membrane Builder → Bilayer Builder](https://charmm-gui.org/?doc=input/membrane.bilayer)
and select **"Build New Bilayer"**.

**Membrane composition** (symmetric; enter same counts for both leaflets):

| Lipid | Count per leaflet |
|---|---|
| PMPE | 46 |
| POPE | 13 |
| QMPE | 12 |
| OYPE | 8 |
| PMPG | 10 |
| PYPG | 9 |
| PVCL2 | 2 |

**System settings:**
- Box size: 79.79 × 79.79 Å lateral (X/Y); Z set automatically (~100 Å)
- Water model: TIP3P
- Temperature: 310.15 K
- Ion concentration: 0.15 M NaCl (no Ca²⁺ needed; G-IM contains no LPS)

**Force field / output page:**
- Force field: **CHARMM36**
- Output format: **GROMACS**

Download the package. Key GROMACS files in `gromacs/`:

| File | Contents |
|---|---|
| `gromacs/step5_input.gro` | Assembled coordinates — **use this, not the PDB** |
| `gromacs/topol.top` | Master topology |
| `gromacs/step6.x_equilibration.mdp` | Six staged equilibration MDP files |
| `gromacs/step7_production.mdp` | Production MDP |

---

## 4. Insert KU04AMP01 into the GROMACS Membrane System

Copy the contents of `gromacs/` to a new `simulation/` directory and work from there.

### Step 4a — Use `step5_input.gro` as the membrane GRO

CHARMM-GUI writes both `step5_input.pdb` and `step5_input.gro`. The PDB has no `CRYST1`
record, so converting it with `gmx editconf` produces a GRO with zero box vectors, causing a
fatal `Empty diagonal for a 3-dimensional periodic box` error in `grompp`. Use the GRO directly:

```bash
cp step5_input.gro membrane.gro
```

Verify the box line (last line) is non-zero:

```bash
tail -1 membrane.gro
# Expected e.g.:   7.97900   7.97900  10.00000   0.00000 ...
```

### Step 4b — Place KU04AMP01 above a leaflet

Determine the z-coordinate of the upper phosphate plane:

```bash
# Phosphate atoms in PE/PG lipids are named P
grep " P " membrane.gro | awk '{print $NF}' | sort -n | tail -5
```

Translate the peptide so its N-terminus is ~15–20 Å (1.5–2.0 nm) above that plane:

```bash
# Example: if the topmost phosphate is at z = 7.5 nm, offset by 1.5 nm → z = 9.0
gmx editconf -f ku04amp01.gro \
    -o amp01_placed.gro -translate 0 0 <z_offset_in_nm>
```

### Step 4c — Merge coordinates

```bash
python3 - <<'EOF'
with open("membrane.gro") as f:
    lines = f.readlines()
natoms_mem = int(lines[1].strip())
mem_atoms = lines[2:2+natoms_mem]
box = lines[2+natoms_mem]

with open("amp01_placed.gro") as f:
    plines = f.readlines()
natoms_amp = int(plines[1].strip())
amp_atoms = plines[2:2+natoms_amp]

with open("system.gro", "w") as out:
    out.write(lines[0])
    out.write(f"{natoms_mem + natoms_amp}\n")
    out.writelines(mem_atoms)
    out.writelines(amp_atoms)
    out.write(box)
print(f"system.gro: {natoms_mem + natoms_amp} atoms")
EOF
```

The box vector is copied from `membrane.gro`, ensuring correct periodic boundary conditions.

### Step 4d — Update the topology

CHARMM36 atomtypes for standard amino acids are already included in the membrane
`forcefield.itp`, so **no additional atomtype file is needed**.

Copy the peptide ITP and position restraint file into the simulation directory:

```bash
cp ../KU04AMP01.itp .
cp posre.itp .
```

Edit `topol.top`:

1. After the last lipid/water `#include` line, add:
   ```
   #include "KU04AMP01.itp"
   ```
2. At the end of the `[ molecules ]` block, add:
   ```
   KU04AMP01         1
   ```

---

## 5. Verify, Neutralize, and Build Index Groups

### Neutralize the system

KU04AMP01 adds approximately +4 to the system net charge. PME electrostatics requires
a neutral system. Add 4 Cl⁻ counter-ions by replacing water molecules.

First produce a TPR for `genion`:

```bash
gmx grompp -f step6.1_equilibration.mdp -c system.gro -r system.gro \
    -p topol.top -n index.ndx -o genion.tpr -maxwarn 5
```

Add neutralising ions. Do **not** pass `-n index.ndx` to `genion` — select `TIP3` when prompted:

```bash
gmx genion -s genion.tpr -o system.gro -p topol.top \
    -pname SOD -nname CLA -neutral
# At the prompt: select the group corresponding to TIP3 (water)
```

`topol.top` is updated automatically (4 fewer TIP3, 4 more CLA). Regenerate the index file
since atom numbers have changed.

### Generate the index file

The CHARMM-GUI equilibration MDP files reference two custom groups — `MEMB` (all lipid and
peptide atoms) and `SOLV` (water + ions) — not present in the default `make_ndx` output.

```bash
python3 make_index.py
# MEMB residues: PMPE, POPE, QMPE, OYPE, PMPG, PYPG, PVCL2, KU04AMP01
# SOLV residues: TIP3, SOD, CLA
```

> Unlike G-OM, G-IM has no LPS glycan residues (ECLIP, ARHM, AGAL, etc.), so the index
> script is simpler and no special glycan residue names need to be listed.

### Sanity check

`grompp` should now complete without a net-charge warning:

```bash
gmx grompp -f step6.1_equilibration.mdp -c system.gro -r system.gro \
    -p topol.top -n index.ndx -o test.tpr -maxwarn 5
```

All three extra flags are required:
- `-r system.gro` — reference coordinates for position restraints (all equilibration MDPs use `-DPOSRES`)
- `-n index.ndx` — custom `MEMB`/`SOLV` groups referenced by `tc_grps` and `comm_grps`

---

## 6. Energy Minimization

Run EM locally before transferring to HPC:

```bash
gmx grompp \
    -f step6.0_minimization.mdp \
    -c system.gro -r system.gro \
    -p topol.top \
    -n index.ndx \
    -o em.tpr -maxwarn 5
gmx mdrun -v -deffnm em -nt 12
```

Check convergence of potential energy:

```bash
echo "Potential 0" | gmx energy -f em.edr -o em_potential.xvg
```

```python
import matplotlib.pyplot as plt
x, y = [], []
with open('em_potential.xvg') as f:
    for line in f:
        if line.startswith(('#','@')): continue
        cols = line.split()
        x.append(float(cols[0])); y.append(float(cols[1]))
plt.plot(x, y)
plt.xlabel('Step'); plt.ylabel('Potential Energy (kJ/mol)')
plt.title('EM convergence')
plt.tight_layout(); plt.show()
```

Prepare a visualizable starting structure by removing PBC artifacts:

```bash
gmx trjconv -s em.tpr -f em.gro -o em_whole.gro \
    -pbc mol -center -n index.ndx
```

Then transfer `em.gro`, `em.edr`, `em.log` to HPC and start equilibration from there.

---

## 7. Equilibration and Production

### Equilibration protocol

The six CHARMM-GUI equilibration stages total ~1.875 ns. Steps 6.1–6.2 (NVT, 125 ps each,
dt = 1 fs) use V-rescale temperature coupling at 310.15 K applied separately to `MEMB` and
`SOLV`; velocities are generated at the start of step 6.1. Steps 6.3–6.6 switch to NPT
(dt = 2 fs, 500 ps each) with C-rescale semi-isotropic pressure coupling (τ_p = 5 ps,
P = 1 bar). Lipid head-group restraints decay progressively:
1000/1000 → 400/400 → 400/200 → 200/200 → 40/100 → 0/0 kJ mol⁻¹ nm⁻² / kJ mol⁻¹ rad⁻².

**KU04AMP01 is unrestrained throughout all equilibration stages.** Its ITP contains no
`#ifdef POSRES` block. It is placed above a leaflet and is free to approach the membrane
naturally.

> **Note on `-r` (reference coordinates):** All equilibration stages use `system.gro` as
> the fixed `-r` reference, not the output of the previous stage. This ensures position
> restraints always pull atoms back toward the original well-placed starting structure and
> do not drift across stages.

```bash
nprocs=48

# Staged equilibration (steps 6.1 through 6.6)
for step in 6.1 6.2 6.3 6.4 6.5 6.6; do
    prev=$(echo $step | awk -F. '{if ($2==1) print "em"; else printf "step%s.%s", $1, $2-1}')
    echo ">>> Starting equilibration step${step} (prev: ${prev})"
    gmx grompp \
        -f step${step}_equilibration.mdp \
        -c ${prev}.gro -r system.gro \
        -p topol.top \
        -n index.ndx \
        -o step${step}.tpr -maxwarn 5
    gmx mdrun -v -deffnm step${step} -nb gpu -pme gpu -bonded gpu -nt ${nprocs} -gpu_id 0
    echo ">>> Completed step${step}. Checkpoint: step${step}.cpt"
done
echo ">>> All equilibration stages complete."
```

### Production MD

Production is split into chunks to fit HPC queue time limits. Each chunk restarts from the
checkpoint of the previous one. The total simulation time is:

$$t_{\text{total}} = \text{set\_cntmax} \times \text{nsteps} \times \text{dt}$$

For example, `nsteps = 5000000` (2 fs timestep → 10 ns/chunk) with `set_cntmax = 10`
gives 100 ns total. Adjust `nsteps` in `step7_production.mdp` and `set_cntmax` to reach
your target. To resume after interruption, set `cnt_start` to the last completed chunk + 1.

```bash
cnt_start=1
set_cntmax=10
for cnt in $(seq ${cnt_start} ${set_cntmax}); do
    pcnt=$((cnt - 1))
    istep=step7_${cnt}
    echo ">>> Starting production chunk ${cnt}/${set_cntmax}: ${istep}"
    if [ ${cnt} -eq 1 ]; then
        pstep=step6.6
        gmx grompp \
            -f step7_production.mdp \
            -c ${pstep}.gro \
            -p topol.top \
            -n index.ndx \
            -o ${istep}.tpr -maxwarn 2
    else
        pstep=step7_${pcnt}
        gmx grompp \
            -f step7_production.mdp \
            -c ${pstep}.gro -t ${pstep}.cpt \
            -p topol.top \
            -n index.ndx \
            -o ${istep}.tpr -maxwarn 2
    fi
    gmx mdrun -v -deffnm ${istep} -nb gpu -pme gpu -bonded gpu -nt ${nprocs} -gpu_id 0
    echo ">>> Completed production chunk ${cnt}/${set_cntmax}. Checkpoint: ${istep}.cpt"
done
echo ">>> All production chunks complete."
```

On an HPC cluster, replace the `gmx mdrun` calls with the appropriate scheduler script
(PBS/Slurm). Use NPT ensemble with semi-isotropic pressure coupling
(`Pcoupltype = semiisotropic`) for a membrane system.

### Simulation timescales

| Goal | Time needed |
|---|---|
| Initial membrane adsorption | 10–100 ns |
| Partial insertion into headgroup region | 100–500 ns |
| Pore formation / full translocation | 1–10 µs (consider enhanced sampling) |

**Practical recommendation: 100–500 ns** for a first test. At 100 ns you should see whether
KU04AMP01 adsorbs and begins interacting with the PG/PE head groups; at 500 ns deeper
penetration is possible.

If translocation within a feasible compute budget is needed, consider:
- **Steered MD (SMD)** — pull the peptide through the membrane along Z with a constant force
- **Umbrella sampling** — compute the translocation PMF along Z
- **Metadynamics or REST2** — accelerate sampling without a fixed reaction coordinate

### Runtime estimate

| Component | Count | ~Atoms |
|---|---|---|
| Lipids (200 lipids × ~50 atoms) | 200 | ~10,000 |
| Water (11,550 molecules) | — | ~34,650 |
| Ions (Na⁺/Cl⁻) | ~110 | ~110 |
| KU04AMP01 | 1 | ~560 |
| **Total** | | **~45,000 atoms** |

G-IM is smaller and simpler than G-OM (no LPS glycan topology). Expected performance on a
V100 + 48 CPUs: **~150–250 ns/day**.

| Target simulation | Estimated wall time |
|---|---|
| 100 ns | ~10–16 hours |
| 500 ns | ~2–3 days |

Run a short benchmark first:

```bash
gmx mdrun -nsteps 5000 -resetstep 4000 -v
# Prints "Performance: X ns/day"
```

Use `-nb gpu -pme gpu -bonded gpu` in `mdrun` to offload as much as possible to the GPU.

---

## 8. Post-Simulation Analysis

### 8.1 Insertion depth

Track the z-coordinate of the KU04AMP01 centre of mass relative to the membrane midplane:

```bash
gmx traj -f md.xtc -s md.tpr -n index.ndx -com -oz -o amp01_z.xvg
```

### 8.2 Interaction energies

To compute per-residue (or per-group) interaction energies between KU04AMP01 and lipid types:

1. Build index groups for peptide residues and lipid species:
   ```bash
   gmx make_ndx -f md.tpr -o interaction_groups.ndx
   ```

2. Prepare a rerun MDP with `energygrps` set to the desired groups:
   ```ini
   energygrps = KU04AMP01 PMPE POPE PMPG PYPG PVCL2
   ```

3. Recompute energies on the saved trajectory:
   ```bash
   gmx grompp -f rerun_energy.mdp -c md.gro -p topol.top \
       -n interaction_groups.ndx -t md.cpt -o rerun_energy.tpr
   gmx mdrun -s rerun_energy.tpr -rerun md.xtc -deffnm rerun_energy
   gmx energy -f rerun_energy.edr -o amp01_membrane_interactions.xvg
   ```

   The descriptive nonbonded interaction energy between groups A and B is:
   ```
   E_interaction = Coul-SR(A, B) + LJ-SR(A, B)
   ```

4. The rerun reports short-range Coulomb and LJ terms only. It does not include entropy,
   conformational reorganisation, or membrane deformation. Use it to rank which residues
   contact the membrane most strongly, paired with H-bond and contact analysis.
   MM/PBSA can be applied as a follow-up for a more binding-like free-energy decomposition.

### 8.3 Hydrogen bonds

```bash
gmx hbond -f md.xtc -s md.tpr -n interaction_groups.ndx \
    -num hbond_amp01_membrane.xvg
```

### 8.4 Tilt angle

Monitor the angle between the KU04AMP01 helix axis and the membrane normal (Z-axis) over
time to quantify orientation and insertion state.

