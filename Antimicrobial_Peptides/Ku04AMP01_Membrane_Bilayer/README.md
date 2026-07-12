# KU04AMP01 + G-IM Membrane MD Workflow (CHARMM-GUI auto-build + peptide height correction)

## Overview

This workflow uses **CHARMM-GUI Membrane Builder** to automatically generate a full **peptide + membrane** GROMACS system (topology and coordinates consistent), then corrects peptide placement if it starts too deep in the bilayer, and finally runs minimization/equilibration/production.

This avoids topology mismatch issues from manually mixing peptide and membrane parameter stacks.

**Note:** This README is for the all-atom simulation. Usually membrane sims need to run for longer times than AA sims will allow, requiring coarsening. For that option, see [MARTINI_vs_EnhancedSampling.md](MARTINI_vs_EnhancedSampling.md)

---

## 1) Build the full system in CHARMM-GUI (peptide already included)

1. Open **CHARMM-GUI Membrane Builder → Bilayer Builder**.
2. Build the G-IM membrane composition:
   - PMPE: 46
   - POPE: 13
   - QMPE: 12
   - OYPE: 8
   - PMPG: 10
   - PYPG: 9
   - PVCL2: 2  
   (counts per leaflet, symmetric bilayer)
3. Add **KU04AMP01** in the builder (upload peptide PDB).
4. Force field/output:
   - Force field: **CHARMM36**
   - Output format: **GROMACS**
5. Download and extract the package.

Use these generated files as base:
- `gromacs/step5_input.gro`
- `gromacs/topol.top`
- `gromacs/toppar/*`
- `gromacs/step6.*_equilibration.mdp`
- `gromacs/step7_production.mdp`

---

## 2) Prepare simulation folder

```bash
mkdir -p simulation
cp -r charmm-gui-*/gromacs/* simulation/
cd simulation
```

---

## 3) Check whether peptide is buried (quick sanity check)

If peptide starts inside bilayer, reposition before energy minimization.

### 3.1 Find lipid phosphate Z planes (headgroups)

```bash
grep -E 'PMPE|POPE|QMPE|OYPE|PMPG|PYPG|PVCL2' step5_input.gro | grep ' P ' | awk '{print $6}' | sort -n | head
grep -E 'PMPE|POPE|QMPE|OYPE|PMPG|PYPG|PVCL2' step5_input.gro | grep ' P ' | awk '{print $6}' | sort -n | tail
```

Typical result:
- lower leaflet phosphate plane: \(z \approx 3.5\)–\(3.7\) nm
- upper leaflet phosphate plane: \(z \approx 7.3\)–\(7.5\) nm

### 3.2 Inspect peptide Z range

If peptide is first molecule in `step5_input.gro`, inspect its Z span (example range shown below; adapt if needed):

```bash
awk 'NR>=3 && NR<=571 {print $6}' step5_input.gro | sort -n | head
awk 'NR>=3 && NR<=571 {print $6}' step5_input.gro | sort -n | tail
```

If peptide spans both leaflet planes, it is embedded and should be moved upward.

---

## 4) Reposition peptide upward (recommended)

### 4.1 Build index groups and extract peptide + non-peptide

```bash
gmx make_ndx -f step5_input.gro -o split.ndx
```

In the prompt:
- Use existing `Protein` group for peptide.
- Create non-protein group:
  - `! "Protein"`
- `q`

Extract peptide:

```bash
gmx trjconv -f step5_input.gro -s step5_input.gro -n split.ndx -o peptide_only.gro
```

Select: `Protein`

Extract rest of system:

```bash
gmx trjconv -f step5_input.gro -s step5_input.gro -n split.ndx -o rest_only.gro
```

Select: non-protein group (`!Protein`)

### 4.2 Choose target height and translate peptide

Set target peptide COM above upper phosphate plane:

$$
z_{\text{target}} = z_{\text{P,upper}} + (1.5\text{ to }2.0)\,\text{nm}
$$

Optional: get peptide COM before move

```bash
gmx traj -f peptide_only.gro -s peptide_only.gro -com -ox pep_com.xvg
```

Translate peptide by chosen \(\Delta z\): I chose 2 nm. Also, need to rotate the peptide to make sure it is outside the membrane completely.

```bash
# Rotate 90° around Y-axis (flattens peptide to XY plane)
gmx editconf -f peptide_only.gro -o peptide_rotated.gro \
  -rotate 0 90 0

# Translate upward by 1.5 nm (to clear upper leaflet)
gmx editconf -f peptide_rotated.gro -o peptide_positioned.gro \
  -translate 0 0 2.0
```

**Note:** Charmm gui yielded a peptide that was slightly inside the membrane. The above prescription rotated it out and made it lie roughly flat along the phosphate plane. If we want spontaneous insertion into the membrane, then an amphipathic peptide should be oriented parallel. Alternatively, a hydrophobic peptide should be perpendicular, with the N-terminal facing the membrane (downward), since the Ku04AMP01 has a chemically active N-terminal. For that, run the following:

```python
import argparse
import numpy as np
import MDAnalysis as mda

def main():
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(
        description="Rotate a peptide around its Center of Mass (COM) using MDAnalysis."
    )
    
    # Required positional arguments for file paths
    parser.add_argument("input_file", help="Path to the input file (e.g., input.gro, input.pdb)")
    parser.add_argument("output_file", help="Path to save the rotated output file (e.g., output.gro)")
    
    # Optional arguments to easily tweak the rotation
    parser.add_argument("-a", "--angle", type=float, default=90.0, help="Rotation angle in degrees (default: 90.0)")
    parser.add_argument("-x", "--axis", choices=['x', 'y', 'z'], default='y', help="Axis to rotate around (default: y)")
    
    args = parser.parse_args()

    # 1. Load the universe
    u = mda.Universe(args.input_file)
    
    # 2. Select the peptide atoms
    peptide = u.select_atoms("protein")
    if len(peptide) == 0:
        raise ValueError("No protein atoms found in the selection! Check your input file formatting.")

    # 3. Calculate Center of Mass
    peptide_com = peptide.center_of_mass()
    
    # 4. Map the chosen axis letter to a directional vector
    axis_mapping = {
        'x': [1, 0, 0],
        'y': [0, 1, 0],
        'z': [0, 0, 1]
    }
    direction_vector = axis_mapping[args.axis]

    # 5. Perform the in-place rotation (Fixed: changed direction= to axis=)
    peptide.rotateby(angle=args.angle, axis=direction_vector, point=peptide_com)
    
    # 6. Save the new coordinates
    peptide.write(args.output_file)
    print(f"Successfully rotated {args.input_file} by {args.angle}° around the {args.axis}-axis via its COM.")
    print(f"Output saved to: {args.output_file}")

if __name__ == "__main__":
    main()
```


If the N-terminal ends up pointing upwards instead of downwards after this command, rerun it using `-rotate 0 -90 0` or `-rotate 0 270 0`. To check this, load the PDB in `PyMol` and run in colsole

```bash

show cell
color gray, all
color red, resi 1 
```

Now the N-terminus is red.

### 4.3 Merge back into one system (peptide first, rest after)
Adjust the name of the peptide file as needed.

```bash
python3 - <<'PY'
def read_gro(path):
    with open(path) as f:
        lines = f.readlines()
    title = lines[0]
    nat = int(lines[1].strip())
    atoms = lines[2:2+nat]
    box = lines[2+nat].rstrip("\n")
    return title, nat, atoms, box

t1,n1,a1,b1 = read_gro("peptide_perp_perfect.gro")
t2,n2,a2,b2 = read_gro("rest_only_tall.gro")

with open("system.gro","w") as out:
    out.write("KU04AMP01 + G-IM (peptide height corrected)\n")
    out.write(f"{n1+n2}\n")
    out.writelines(a1)   # peptide first
    out.writelines(a2)   # rest
    out.write(b2 + "\n") # keep box from rest system
print(f"Wrote system.gro with {n1+n2} atoms")
PY
```

topol.top` generally does **not** need changes if no molecules were added/removed.

---

## 5) Preprocess, neutralize, and index groups

### 5.1 First grompp (for genion)

```bash
gmx grompp -f step6.1_equilibration.mdp -n index.ndx -c system.gro -r system.gro -p topol.top -o genion.tpr -maxwarn 5
```

If grompp crashes with unknown index group ids, then run the `make_index.py` script from section 5.3 below.

### 5.2 Neutralize

```bash
gmx genion -s genion.tpr -o system.gro -p topol.top -pname SOD -nname CLA -neutral
```

Select water group (`TIP3`) at prompt.

### 5.3 Build/update `index.ndx` groups required by CHARMM-GUI mdp

Need at minimum:
- `MEMB` = lipids + peptide
- `SOLV` = TIP3 + ions

If you have your helper script:

```bash
python3 make_index.py
```

Otherwise, create manually with `gmx make_ndx`.

---

## 6) Energy minimization

```bash
gmx grompp -f step6.0_minimization.mdp -c system.gro -r system.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 5
gmx mdrun -deffnm em -v
```

---

## 7) Equilibration (same protocol as original README)

Run steps `6.1` to `6.6` in sequence.  
Use:
- `-c` previous step output
- `-r system.gro` (fixed restraint reference)
- `-p topol.top`
- `-n index.ndx`

Example loop:

```bash
nprocs=48

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
    echo ">>> Completed step${step}"
done
```

---

## 8) Production MD (same as original workflow)

Run chunked production with checkpoint continuation:

```bash
cnt_start=1
set_cntmax=10
nprocs=48

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
    echo ">>> Completed production chunk ${cnt}/${set_cntmax}"
done
```

---

## 9) Sanity checks

### 9.1 Visual check

Before EM:
- open `system.gro` in VMD/PyMOL
- confirm peptide is above upper leaflet (not crossing membrane center)
- check no severe steric overlap with lipids/water

### 9.2 grompp check

```bash
gmx grompp -f step6.1_equilibration.mdp -c system.gro -r system.gro -p topol.top -n index.ndx -o test.tpr -maxwarn 5
```

---

## Notes / Best practices

- Keep topology/parameters strictly from one CHARMM-GUI package (no mixing stacks).
- If peptide was moved, always do minimization before equilibration.
- Do not change atom counts during reposition unless you also update topology.
- If you regenerate ions (`genion`), rebuild index groups because atom numbers change.
- For adsorption-first physics, peptide should start in water ~1.5–2.0 nm above phosphate plane.
- Because the peptide is 6.6 nm long, the Z-box must be expanded to 17.0 nm.
    1.  **Split** the system into `peptide_alone.gro` and `membrane_water.gro`.
    2.  **Use `gmx editconf`** to set the box to `8.08472 8.08472 17.0` for both, centering the membrane at $Z=5.5$ and the peptide at $Z=12.8$.
    3.  **Merge** the `.gro` files using the Python script provided in the original `charmm-gui` output to handle coordinate concatenation.
- Solvation and Topology Update
    1.  Fill the vacuum: `gmx solvate -cp merged_tall.gro -cs spc216.gro -o final_solvated.gro -p topol.top`.
    2.  **Update `topol.top`**: Manually update the `[ molecules ]` section, combining the original water count with new SOL molecules.
    3.  **Ensure `index.ndx`** is updated to include `SOL` in the water group (e.g., `water_residues = ['TIP3', 'SOLV', 'SOL']`).
    4.  May need to rotate the peptide again to ensure that N-terminals face the membrane

---