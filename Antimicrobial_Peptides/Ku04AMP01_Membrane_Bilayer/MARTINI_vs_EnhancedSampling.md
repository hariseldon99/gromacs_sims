# KU04AMP01–Gram-negative Inner Membrane Coarse-Grained MD Simulation

> **Martini 3 / GROMACS workflow for simulating the antimicrobial peptide KU04AMP01 interacting with a physiologically representative *Escherichia coli* Gram-negative inner membrane.**

---

## Table of Contents

- [KU04AMP01–Gram-negative Inner Membrane Coarse-Grained MD Simulation](#ku04amp01gram-negative-inner-membrane-coarse-grained-md-simulation)
  - [Table of Contents](#table-of-contents)
    - [Appendices](#appendices)
- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Simulation Strategy](#simulation-strategy)
- [Membrane Composition](#membrane-composition)
- [Lipid Substitutions](#lipid-substitutions)
  - [Motivation](#motivation)
  - [Design principles](#design-principles)
  - [Final substitutions](#final-substitutions)
- [Final Martini 3 Membrane Composition](#final-martini-3-membrane-composition)
- [Workflow](#workflow)
  - [1. Prepare the Peptide Structure](#1-prepare-the-peptide-structure)
    - [Objective](#objective)
    - [Remove hydrogen atoms](#remove-hydrogen-atoms)
    - [Orient the peptide](#orient-the-peptide)
    - [Expected output](#expected-output)
- [2. Generate the Martini 3 Peptide](#2-generate-the-martini-3-peptide)
    - [Objective](#objective-1)
    - [Expected output](#expected-output-1)
- [3. Position the Peptide](#3-position-the-peptide)
    - [Objective](#objective-2)
    - [Expected output](#expected-output-2)
- [4. Construct the Membrane](#4-construct-the-membrane)
    - [Objective](#objective-3)
  - [4.1 Prepare the topology includes](#41-prepare-the-topology-includes)
  - [4.2 Generate TOCL coordinates](#42-generate-tocl-coordinates)
    - [Expected output](#expected-output-3)
  - [4.3 Build the membrane](#43-build-the-membrane)
    - [Expected output](#expected-output-4)
- [5. Create the Index File](#5-create-the-index-file)
    - [Objective](#objective-4)
    - [Expected output](#expected-output-5)
- [6. Energy Minimization](#6-energy-minimization)
    - [Objective](#objective-5)
    - [Expected output](#expected-output-6)
- [7. Equilibration](#7-equilibration)
    - [Objective](#objective-6)
    - [Expected output](#expected-output-7)
- [8. Production Molecular Dynamics](#8-production-molecular-dynamics)
    - [Objective](#objective-7)
    - [Expected output](#expected-output-8)
  - [Simulation Parameters](#simulation-parameters)
  - [Project Notes](#project-notes)
    - [Elastic network](#elastic-network)
    - [Membrane composition](#membrane-composition-1)
    - [Deprecated topologies](#deprecated-topologies)
    - [Coordinate generation](#coordinate-generation)
    - [Simulation protocol](#simulation-protocol)
  - [Limitations](#limitations)
  - [References](#references)
  - [Citation](#citation)
  - [Acknowledgements](#acknowledgements)
- [Appendix A — Comparison of Enhanced Sampling Methods](#appendix-a--comparison-of-enhanced-sampling-methods)
- [Overview](#overview-1)
- [Scientific Objective](#scientific-objective)
- [Candidate Simulation Strategies](#candidate-simulation-strategies)
- [Atomistic Molecular Dynamics](#atomistic-molecular-dynamics)
- [Umbrella Sampling](#umbrella-sampling)
  - [Advantages](#advantages)
  - [Limitations](#limitations-1)
- [Steered Molecular Dynamics](#steered-molecular-dynamics)
  - [Advantages](#advantages-1)
  - [Limitations](#limitations-2)
- [Metadynamics](#metadynamics)
  - [Advantages](#advantages-2)
  - [Limitations](#limitations-3)
- [Martini 3 Coarse-Grained Molecular Dynamics](#martini-3-coarse-grained-molecular-dynamics)
- [Advantages for the Present Study](#advantages-for-the-present-study)
- [Why Unbiased Martini 3 Was Chosen](#why-unbiased-martini-3-was-chosen)
- [Future Directions](#future-directions)
- [Conclusions](#conclusions)
- [Appendix B — Lipid Substitutions and Cyclopropane Lipid Parameterization](#appendix-b--lipid-substitutions-and-cyclopropane-lipid-parameterization)
- [Overview](#overview-2)
- [Native Membrane Composition](#native-membrane-composition)
- [The Cyclopropane Problem](#the-cyclopropane-problem)
- [Why Custom Parameters Were Not Used](#why-custom-parameters-were-not-used)
- [Design Criteria](#design-criteria)
- [Candidate Replacement Strategies](#candidate-replacement-strategies)
- [Why Saturated Lipids Were Rejected](#why-saturated-lipids-were-rejected)
- [Adopted Lipid Substitutions](#adopted-lipid-substitutions)
- [Rationale for Individual Substitutions](#rationale-for-individual-substitutions)
  - [PMPE → DVPE](#pmpe--dvpe)
  - [QMPE → POPE](#qmpe--pope)
  - [OYPE → DOPE](#oype--dope)
  - [PMPG and PYPG → POPG](#pmpg-and-pypg--popg)
  - [PVCL2 → TOCL](#pvcl2--tocl)
- [Consequences of the Substitutions](#consequences-of-the-substitutions)
- [Future Improvements](#future-improvements)
- [Conclusions](#conclusions-1)
- [Appendix C — Repository Maintenance Notes and Deprecated Topologies](#appendix-c--repository-maintenance-notes-and-deprecated-topologies)
- [Overview](#overview-3)
- [Development Strategy](#development-strategy)
- [Custom Lipid Topologies](#custom-lipid-topologies)
- [COBY Compatibility Issues](#coby-compatibility-issues)
- [Manual TOCL Coordinate Generation](#manual-tocl-coordinate-generation)
- [Current Production Workflow](#current-production-workflow)
- [Deprecated Files](#deprecated-files)
- [Recommendations for Future Development](#recommendations-for-future-development)
- [Conclusions](#conclusions-2)

### Appendices

- [Appendix A — Comparison of Enhanced Sampling Methods](#appendix-a-comparison-of-enhanced-sampling-methods)
- [Appendix B — Lipid Substitutions and cy17:0 Parameterization](#appendix-b-lipid-substitutions-and-cy17:0-parameterization)
- [Appendix C — Historical Notes and Deprecated Topologies](appendix-c-historical-notes-and-deprecated-topologies)

---

# Overview

This repository contains all files required to reproduce coarse-grained molecular dynamics simulations of the 36-residue antimicrobial peptide **KU04AMP01** interacting with a physiologically representative **Gram-negative inner membrane (G-IM)** of *Escherichia coli* using the **Martini 3** coarse-grained force field.

The workflow includes

- coarse-graining of the peptide using **martinize2**,
- construction of a heterogeneous seven-lipid membrane using **COBY**,
- preparation of GROMACS topologies,
- equilibration,
- long unbiased production simulations.

The simulation protocol was developed as a computationally tractable alternative to multi-microsecond atomistic CHARMM36 simulations while preserving the essential membrane physics relevant to peptide adsorption and spontaneous membrane insertion.

---

# Scientific Background

Spontaneous translocation of antimicrobial peptides across bacterial membranes typically occurs on the **microsecond to millisecond** timescale.

Although all-atom simulations using CHARMM36 provide excellent structural accuracy, they are generally unable to reach these timescales without enhanced sampling methods.

The Martini 3 coarse-grained force field reduces the number of interaction sites by approximately a factor of four while allowing a larger integration timestep and a smoother effective energy landscape. In practice, this provides roughly two orders of magnitude greater sampling efficiency than atomistic simulations while retaining realistic membrane mechanics and peptide–lipid interactions [1].

For the present study this offers three important advantages:

- spontaneous insertion can be observed without imposing an external reaction coordinate;
- multiple insertion pathways can be sampled within practical computational cost;
- membrane deformation, peptide bending and lipid rearrangement emerge naturally during the simulation.

A detailed comparison with umbrella sampling, steered molecular dynamics and metadynamics is provided in **[Appendix A](#appendix-a-comparison-of-enhanced-sampling-methods).**

---

# Simulation Strategy

The overall workflow is

```text
Atomistic peptide
        │
        ▼
martinize2
        │
        ▼
Coarse-grained peptide
        │
        ▼
COBY membrane construction
        │
        ▼
Energy minimization
        │
        ▼
NPT equilibration
        │
        ▼
Unbiased production simulation
```

Only validated Martini 3 lipid parameters are used in the final production system.

Several custom lipid topologies explored during method development are retained only for documentation purposes and **must not** be used for production simulations. Their history is described in **[Appendix C](appendix-c-historical-notes-and-deprecated-topologies).**

---

# Membrane Composition

The target membrane is a heterogeneous model of the ***E. coli* Gram-negative inner membrane (G-IM)** containing seven lipid species.

| Native lipid | Count / leaflet |
|--------------|---------------:|
| PMPE | 46 |
| POPE | 13 |
| QMPE | 12 |
| OYPE | 8 |
| PMPG | 10 |
| PYPG | 9 |
| PVCL2 | 2 |

The bilayer is symmetric, giving twice these numbers in the complete simulation box.

---

# Lipid Substitutions

## Motivation

Three native membrane lipids contain **cyclopropane fatty acids (cy17:0)**.

At the time of writing, the official Martini 3 lipid parameter library does **not** include validated parameters for cyclopropane-containing phospholipids [2].

Consequently, direct coarse-grained representations of

- PMPE,
- QMPE,
- PMPG,

are unavailable.

The objective of the substitutions described below is therefore **not** to reproduce cyclopropane chemistry exactly, but to preserve the physical state of the membrane while remaining entirely within validated Martini 3 parameters.

A detailed discussion of this problem, together with alternative substitution strategies that were considered and rejected, is given in **[Appendix B](#appendix-b-lipid-substitutions-and-cy17:0-parameterization).**

---

## Design principles

The substitutions were selected according to four criteria.

1. Preserve the liquid-disordered phase at 310 K.

2. Preserve overall phospholipid headgroup composition.

3. Minimize perturbation of membrane charge.

4. Use only officially parameterized Martini 3 lipids.

Substitutions that introduced gel-phase lipids (for example DPPE or DPPG) were rejected because they artificially suppress membrane insertion by substantially increasing bilayer order (Appendix B).

---

## Final substitutions

| Native lipid | Martini 3 substitute | Rationale |
|--------------|---------------------|-----------|
| PMPE | DVPE | Vaccenic-acid analogue; preserves membrane fluidity |
| QMPE | POPE | Closest validated PE lipid available |
| OYPE | DOPE | Closest validated unsaturated PE lipid compatible with COBY workflow |
| PMPG | POPG | Used to avoid current COBY incompatibilities with PYPG |
| PYPG | POPG | Equivalent PG used throughout membrane construction |
| PVCL2 | TOCL | Validated Martini cardiolipin with manually generated coordinates |

The substitution of **POPG** for **PYPG** is an implementation workaround for current COBY limitations rather than a force-field limitation.

Likewise, **TOCL** coordinates are supplied externally because COBY does not currently construct this lipid correctly.

These workarounds affect coordinate generation only; all production simulations employ standard Martini 3 interaction parameters.

---

# Final Martini 3 Membrane Composition

The production membrane contains the following lipid composition per leaflet.

| Martini lipid | Count |
|---------------|------:|
| DVPE | 46 |
| POPG | 19 |
| POPE | 25 |
| DOPE | 8 |
| TOCL | 2 |

The membrane is constructed as a symmetric bilayer.

---
# Workflow

This section describes the complete procedure used to prepare, equilibrate, and simulate the coarse-grained KU04AMP01–membrane system.

The workflow consists of

1. Preparing the peptide structure.
2. Generating the Martini 3 coarse-grained peptide.
3. Constructing the membrane using COBY.
4. Assembling the complete simulation system.
5. Energy minimization.
6. Equilibration.
7. Production molecular dynamics.

---

## 1. Prepare the Peptide Structure

### Objective

Prepare a hydrogen-free peptide structure suitable for Martini coarse-graining and orient it with respect to the membrane.

### Remove hydrogen atoms

Remove hydrogen atoms and terminal OXT atoms before coarse-graining.

```bash
grep -v -E " H   | H[A-Z]| OXT" KU04AMP01.pdb > KU04AMP01_clean.pdb
```

> **Note**
>
> Hydrogen atoms and terminal OXT atoms are discarded during Martini mapping and therefore should be removed before running `martinize2`.

---

### Orient the peptide

Orient the peptide with respect to the membrane normal before coarse-graining. In this work, the orientation was obtained using the **PPM server** of the **Orientations of Proteins in Membranes (OPM)** database (Lomize *et al.*, 2012) [8].

Recommended settings:

- Upload `KU04AMP01_clean.pdb`
- Set **Topology (N-ter)** to **in**
- Leave all remaining settings at their defaults.

After downloading the oriented structure, overwrite

```text
KU04AMP01_clean.pdb
```

Remove the dummy membrane atoms inserted by OPM.

```bash
grep -v HETATM KU04AMP01_clean.pdb > temp.pdb
mv temp.pdb KU04AMP01_clean.pdb
```

> **Note**
>
> The peptide was oriented with its positively charged N-terminus facing the membrane to reproduce the initial adsorption geometry used in the preceding all-atom CHARMM36 simulations.

### Expected output

```text
KU04AMP01_clean.pdb
```

---

# 2. Generate the Martini 3 Peptide

### Objective

Generate a Martini 3 representation of the peptide using `martinize2`.

Run

```bash
martinize2 \
    -f KU04AMP01_clean.pdb \
    -dssp $MKDSSP \
    -nt \
    -cys auto \
    -p backbone \
    -sep \
    -ff martini3001 \
    -elastic \
    -x KU04AMP01_cg.pdb \
    -o KU04.top \
    -scfix \
    -ef 500.0 \
    -el 0.5 \
    -eu 0.9 \
    -ea 0 \
    -ep 0 \
    -maxwarn 5
```

The elastic-network parameters follow the ElNeDyn protocol of Periole *et al.* [6].

Replace the automatically generated force-field include

```cpp
#include "martini_v3.0.0.itp"
```

with the appropriate installation-specific Martini force-field path.

> **Tip**
>
> Rename the default molecule name (`molecule_0`) generated by `martinize2` to `KU04` to improve the readability of topology files and GROMACS output.

### Expected output

```text
KU04AMP01_cg.pdb
KU04.itp
KU04.top
```

---

# 3. Position the Peptide

### Objective

Move the peptide above the membrane before membrane construction.

COBY centers proteins within the membrane by default. For the present simulations the peptide was initially placed above the upper leaflet.

```bash
gmx editconf \
    -f KU04AMP01_cg.pdb \
    -o KU04AMP01_cg_shifted.pdb \
    -translate 0 0 2.4
```

> **Warning**
>
> Unless the peptide is translated after coarse-graining, the simulation will begin with the peptide embedded inside the membrane rather than adsorbed above its surface.

The translation assumes the default simulation box

```text
14 × 14 × 8 nm³
```

used in this simulation. This differs from the more standard 14 × 14 × 84 nm³ box in order to shorten the z-dimension. When filled with water, less water 
height leaves no room for the peptide to escape into a bulk-water buffer zone.

> **Note**
> 
>You may have to experiment with this z-shift a bit. Even 2.4 nm is too far for insertion. 

### Expected output

```text
KU04AMP01_cg_shifted.pdb
```

---

# 4. Construct the Membrane

### Objective

Construct the heterogeneous Gram-negative inner membrane using COBY.

---

## 4.1 Prepare the topology includes

Create

```text
top_for_COBY.itp
```

containing

```cpp
#include "martini3/martini_v3.0.0.itp"
#include "martini3/martini_v3.0.0_ffbonded_v2.itp"

#include "martini3/martini_v3.0.0_phospholipids_PE_v2.itp"
#include "martini3/martini_v3.0.0_phospholipids_PG_v2.itp"
#include "martini3/martini_v3.0.0_phospholipids_CL_v2.itp"

#include "martini3/martini_v3.0.0_solvents_v1.itp"
#include "martini3/martini_v3.0.0_ions_v1.itp"
#include "martini3/martini_v3.0.0_sterols_v1.itp"

#include "KU04.itp"
```

---

## 4.2 Generate TOCL coordinates

Generate the TOCL coordinate file supplied with this repository.

```bash
python simulation_cg/generate_tocl.py
```

> **Warning**
>
> Current COBY releases do not correctly generate TOCL coordinates. This step must be completed before constructing the membrane.

### Expected output

```text
tocl_single.gro
```

---

## 4.3 Build the membrane

Configure

```text
simulation_cg/memb_build.py
```

using the membrane composition described in the section
[Final Martini 3 Membrane Composition](#final-martini-3-membrane-composition).

Construct the complete system

```bash
python simulation_cg/memb_build.py
```

Note: The memb-build script looks like this:

```python
#!/usr/bin/env python
import COBY

sysname = "ku04_gim_multi"

COBY.COBY(
    box=[14, 14, 8],
    # Leaflet compositions remain unchanged
    membrane="lipid:DVPE:46 lipid:POPG:19 lipid:POPE:12 lipid:DOPE:8 lipid:TOCL:2 apl:0.61",
    
    # POINT TO THE PEPTIDE FILE HERE
    protein="file:KU04AMP01_cg.pdb moleculetypes:KU04 center_protein:False", 
    
    solvation="default",
    itp_input="file:top_for_COBY.itp",
    molecule_import="file:tocl_single.gro moleculetypes:TOCL",
    out_sys=sysname,
    out_top=sysname + ".top",
    out_log=sysname + ".log",
    sn=sysname,
)
```

> **Warning**
>
> If COBY spits out a warning that the protein is too high and jumps across the bounding box, then reduce the z-shift of KU04 a bit in the previous step. Alternatively, increase the z-height of the box.

### Expected output

```text
ku04_gim.gro
ku04_gim.top
ku04_gim.log
```


---

# 5. Create the Index File

### Objective

Create temperature-coupling groups for GROMACS. Replace `ku04_gim.gro` with whatever final gro file you got, if running multi-peptide sims.

```bash
gmx make_ndx \
    -f ku04_gim.gro \
    -o index.ndx
```

Use

```text
"DVPE" | "POPG" | "POPE" | "DOPE" | "TOCL"
name 22 Membrane

"W" | "Ion"
name 23 Solvent

q
```

The resulting coupling groups are

| Group | Contents |
|--------|----------|
| Protein | KU04 peptide |
| Membrane | All phospholipids |
| Solvent | Water and ions |

### Expected output

```text
index.ndx
```

---

# 6. Energy Minimization

### Objective

Update the default `minim.mdp`: update the top parameters and append the Martini 3 interaction blocks to look like this:

```bash
define                  = -DFLEXIBLE   ; Allows water/peptides to relax fully
integrator              = steep
emtol                   = 10.0         ; Tighter tolerance ensures a very stable starting structure
emstep                  = 0.01         ; Faster, robust relaxation for crowded boxes
nsteps                  = 50000

; MANDATORY MARTINI 3 NON-BONDED PARAMETERS
cutoff-scheme            = Verlet
nstlist                  = 20
ns-type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005

coulombtype              = cutoff
coulomb-modifier         = Potential-shift-verlet
rcoulomb                 = 1.1
epsilon_r                = 15          ; Martini 3 screening constant
vdwtype                  = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

constraints              = none
```

Remove steric clashes and relax the initial configuration.

Prepare the run

```bash
gmx grompp \
    -f minim.mdp \
    -c ku04_gim.gro \
    -r ku04_gim.gro \
    -p ku04_gim.top \
    -o em.tpr
```

Run

```bash
gmx mdrun \
    -deffnm em
```

### Expected output

```text
em.gro
em.edr
em.log
em.cpt
```

---

# 7. Equilibration

### Objective

Equilibrate the system under NPT conditions.

```bash
gmx grompp \
    -f npt.mdp \
    -c em.gro \
    -r em.gro \
    -p ku04_gim.top \
    -n index.ndx \
    -o npt.tpr
```

```bash
gmx mdrun \
    -deffnm npt
```

The supplied protocol performs approximately **50 ns** of restrained equilibration.

### Expected output

```text
npt.gro
npt.cpt
npt.edr
```

---

# 8. Production Molecular Dynamics

### Objective

Perform an unbiased Martini 3 production simulation.

Prepare the production run

```bash
gmx grompp \
    -f md.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p ku04_gim.top \
    -n index.ndx \
    -o md.tpr
```

Run

```bash
gmx mdrun \
    -deffnm md
```

For HPC environments, an example PBS submission script is provided:

```text
simulation_cg/run_eq_prod_martini.pbs
```

### Expected output

```text
md.xtc
md.tpr
md.cpt
md.edr
md.log
md.gro
```

---

## Simulation Parameters

The repository includes three GROMACS parameter files:

| File | Purpose |
|------|---------|
| `minim.mdp` | Energy minimization |
| `npt.mdp` | Equilibration under NPT conditions |
| `md.mdp` | Production molecular dynamics |

These files were optimized for Martini 3 membrane simulations following established Martini simulation practices, primarily those of de Jong *et al.* (2016) [5].

Users intending to modify the simulation protocol should consult that reference before altering integration parameters, pressure coupling, or thermostat settings.

---

## Project Notes

### Elastic network

The peptide structure is stabilized using the ElNeDyn elastic network generated automatically by `martinize2`.

The elastic-network parameters are

| Parameter | Value |
|-----------|------:|
| Force constant | 500 kJ mol⁻¹ nm⁻² |
| Lower cutoff | 0.5 nm |
| Upper cutoff | 0.9 nm |

These values follow the recommendations of Periole *et al.* (2009) [6].

---

### Membrane composition

The membrane composition reproduces the phospholipid distribution of the *Escherichia coli* Gram-negative inner membrane as closely as possible using currently available Martini 3 lipid parameters.

Several substitutions were required because validated Martini 3 parameters for cyclopropane-containing phospholipids are not presently available.

The complete scientific justification is provided in
[Appendix B](#appendix-b-lipid-substitutions-and-cy17:0-parameterization).

---

### Deprecated topologies

Early versions of this project explored custom Martini topologies for cyclopropane-containing lipids.

These topologies are retained only for historical reference and **must not** be used for production simulations.

See
[Appendix C](appendix-c-historical-notes-and-deprecated-topologies)
for details.

---

### Coordinate generation

Current versions of COBY do not correctly generate coordinates for Martini TOCL.

Accordingly, the repository provides a helper script

```text
simulation_cg/generate_tocl.py
```

which generates the required coordinate file prior to membrane construction.

Only coordinate generation is affected; all production simulations use the standard Martini 3 TOCL interaction parameters.

---

### Simulation protocol

The production simulations described here are entirely unbiased.

No external steering forces, umbrella restraints, metadynamics bias, or enhanced-sampling methods are applied.

Consequently, all observed peptide adsorption and membrane insertion events arise spontaneously from the Martini dynamics.

The rationale for this approach is discussed in
[Appendix A](#appendix-a-comparison-of-enhanced-sampling-methods).

---

## Limitations

The present workflow has the following limitations.

- Cyclopropane-containing phospholipids are represented by validated Martini 3 substitutes rather than explicit coarse-grained models.
- The membrane is symmetric and therefore does not reproduce leaflet asymmetry.
- Polarizable water is not employed.
- Long-timescale insertion events remain stochastic and may require multiple independent trajectories for robust statistical analysis.

These limitations reflect the current state of the Martini 3 ecosystem rather than restrictions imposed by the workflow itself.

---

## References

1. Souza, P. C. T., Alessandri, R., Barnoud, J., *et al.* Martini 3: A General Purpose Force Field for Coarse-Grained Molecular Dynamics. *Nature Methods* **18**, 382–388 (2021). https://doi.org/10.1038/s41592-021-01098-3

2. Martini Force Field. Martini 3 Lipid Parameter Library. https://cgmartini.nl/

3. Wassenaar, T. A., Ingólfsson, H. I., Böckmann, R. A., Tieleman, D. P. & Marrink, S. J. Computational Lipidomics with INSANE: A Versatile Tool for Generating Custom Membranes for Molecular Simulations. *Journal of Chemical Theory and Computation* **11**, 2144–2155 (2015).

4. Grime, J. M. A. & Madsen, J. J. COBY: A Python-Based Membrane Builder for Complex Lipid Bilayers. *(Use the exact publication corresponding to the COBY version used in this repository.)*

5. de Jong, D. H., Baoukina, S., Ingólfsson, H. I., & Marrink, S. J. Martini Straight: Boosting Performance Using a Shorter Cutoff and GPUs. *Computer Physics Communications* **199**, 1–7 (2016).

6. Periole, X., Cavalli, M., Marrink, S. J. & Ceruso, M. A. Combining an Elastic Network with a Coarse-Grained Molecular Force Field: Structure, Dynamics and Intermolecular Recognition. *Journal of Chemical Theory and Computation* **5**, 2531–2543 (2009).

7. Abraham, M. J., Murtola, T., Schulz, R., *et al.* GROMACS: High Performance Molecular Simulations through Multi-Level Parallelism from Laptops to Supercomputers. *SoftwareX* **1–2**, 19–25 (2015).

8. Lomize, M. A., Pogozheva, I. D., Joo, H., Mosberg, H. I. & Lomize, A. L. OPM Database and PPM Web Server: Resources for Positioning Proteins in Membranes. *Nucleic Acids Research* **40**, D370–D376 (2012).

---

## Citation

If this workflow contributes to published work, please cite the Martini 3 force field [1], GROMACS [7], and the original publications describing any additional software used in your simulations.

---

## Acknowledgements

This workflow integrates established tools from the Martini and GROMACS communities, together with project-specific scripts developed for constructing physiologically representative Gram-negative inner membranes.



# Appendix A — Comparison of Enhanced Sampling Methods


---

# Overview

The objective of this project is to investigate the spontaneous interaction of the antimicrobial peptide **KU04AMP01** with a physiologically representative **Gram-negative inner membrane (G-IM)** of *Escherichia coli*. The primary quantity of interest is the natural pathway by which the peptide adsorbs to, inserts into, and perturbs the membrane.

Several simulation strategies were evaluated before selecting the final protocol. This appendix summarizes those approaches and explains why **unbiased Martini 3 molecular dynamics** was ultimately chosen.

---

# Scientific Objective

The simulations were designed to address the following questions.

- Does KU04AMP01 spontaneously adsorb to the membrane?
- Does spontaneous membrane insertion occur?
- What structural rearrangements accompany insertion?
- Which lipid species participate in peptide binding?
- Does the peptide induce local membrane deformation?

These questions require the dynamics to evolve naturally, without imposing a predefined reaction coordinate.

---

# Candidate Simulation Strategies

The following approaches were considered.

| Method | Bias applied | Can observe spontaneous pathway? | Free-energy calculation | Selected |
|---------|--------------|----------------------------------|-------------------------|----------|
| Atomistic MD | None | Yes | No | ✗ |
| Umbrella sampling | Harmonic restraint | No | Yes | ✗ |
| Steered MD | External force | No | No | ✗ |
| Metadynamics | History-dependent bias | Partially | Yes | ✗ |
| Martini 3 CG MD | None | Yes | No | ✓ |

---

# Atomistic Molecular Dynamics

Conventional atomistic molecular dynamics provides the highest structural resolution and the most accurate description of peptide–lipid interactions.

However, spontaneous membrane insertion typically occurs on microsecond to millisecond timescales, whereas simulations of comparable atomistic systems are usually limited to hundreds of nanoseconds or a few microseconds.

For the present system, unbiased atomistic simulations were therefore unlikely to sample multiple spontaneous insertion events within practical computational cost.

Atomistic simulations remain valuable for

- validating coarse-grained results,
- analysing specific binding configurations,
- studying local hydrogen-bonding networks, and
- refining structural models obtained from coarse-grained simulations.

---

# Umbrella Sampling

Umbrella sampling computes the free-energy profile along a predefined reaction coordinate by restraining the system within overlapping sampling windows.

For membrane insertion studies, the reaction coordinate is typically the peptide centre-of-mass distance from the membrane.

## Advantages

- Quantitative free-energy profiles.
- Potential of mean force (PMF).
- Well-established methodology.
- Excellent convergence for simple reaction coordinates.

## Limitations

Umbrella sampling assumes that the important physics can be described by the chosen reaction coordinate.

For antimicrobial peptides this assumption is often inadequate because insertion involves

- peptide rotation,
- peptide bending,
- membrane deformation,
- lipid rearrangement,
- transient pore formation, and
- cooperative structural fluctuations.

These processes are not uniquely determined by the peptide's distance from the membrane.

Consequently, umbrella sampling is well suited for determining free-energy barriers **after** the insertion mechanism is known, but is less appropriate for discovering that mechanism.

---

# Steered Molecular Dynamics

Steered molecular dynamics (SMD) accelerates rare events by applying an external force to selected atoms or collective variables.

Typical applications include

- pulling peptides into membranes,
- extracting ligands,
- unfolding proteins, and
- estimating non-equilibrium work.

## Advantages

- Efficient generation of insertion trajectories.
- Computationally inexpensive.
- Useful for generating initial configurations.

## Limitations

The applied force changes the natural dynamics.

The observed pathway therefore depends on

- pulling direction,
- pulling velocity,
- spring constant, and
- choice of reaction coordinate.

Different parameter choices can produce different insertion mechanisms.

Because the objective of the present study is to observe spontaneous membrane insertion, externally driven trajectories were considered unsuitable.

---

# Metadynamics

Metadynamics accelerates sampling by constructing a history-dependent bias potential along one or more collective variables.

Unlike umbrella sampling, the bias evolves dynamically during the simulation.

## Advantages

- Efficient exploration of rugged free-energy landscapes.
- Can estimate multidimensional free-energy surfaces.
- Requires less manual setup than umbrella sampling.

## Limitations

The success of metadynamics depends critically on the choice of collective variables.

For peptide–membrane systems, no small set of collective variables fully captures

- peptide orientation,
- insertion depth,
- membrane deformation,
- lipid rearrangement, and
- conformational flexibility.

Poorly chosen collective variables may therefore accelerate an unphysical pathway while suppressing the one of genuine interest.

---

# Martini 3 Coarse-Grained Molecular Dynamics

Martini 3 represents approximately four heavy atoms by a single coarse-grained interaction site, substantially reducing the number of degrees of freedom while preserving essential structural and thermodynamic properties.

This reduction permits

- larger integration timesteps,
- smoother effective energy landscapes,
- faster diffusion,
- substantially longer simulations than atomistic models.

Importantly, no external bias is required.

The peptide is free to

- diffuse,
- rotate,
- adsorb,
- insert,
- deform the membrane,
- recruit lipids,

according to the force field alone.

---

# Advantages for the Present Study

The Martini 3 approach satisfies the principal requirements of this project.

- No reaction coordinate is imposed.
- No external force is applied.
- Multiple insertion pathways remain possible.
- Membrane deformation emerges naturally.
- Lipid rearrangements occur without constraint.
- Multi-microsecond trajectories are computationally practical.

Although the absolute kinetics of coarse-grained simulations should not be interpreted quantitatively, the method is well suited for identifying physically plausible insertion mechanisms.

---

# Why Unbiased Martini 3 Was Chosen

The principal objective of this work is **mechanistic**, rather than thermodynamic.

Specifically, the goal is to determine *how* KU04AMP01 interacts with the membrane rather than to compute the free-energy barrier associated with a predefined insertion pathway.

Unbiased Martini 3 simulations allow the peptide to explore configurational space without externally imposed forces or restraints. Consequently,

- adsorption,
- insertion,
- membrane deformation,
- lipid recruitment, and
- peptide conformational changes

all emerge naturally from the dynamics.

Once representative insertion pathways have been identified, enhanced-sampling methods such as umbrella sampling may be employed in future work to quantify the corresponding free-energy landscapes.

---

# Future Directions

The workflow presented in this repository provides a starting point for subsequent quantitative studies.

Possible extensions include

- umbrella sampling of insertion pathways identified here;
- metadynamics using collective variables informed by unbiased trajectories;
- atomistic backmapping of representative coarse-grained configurations;
- comparison with experimental measurements of membrane permeabilization;
- investigation of peptide oligomerization at higher peptide concentrations.

---

# Conclusions

The final simulation protocol employs **unbiased Martini 3 coarse-grained molecular dynamics** because it provides the best compromise between computational efficiency and physical realism for studying spontaneous peptide–membrane interactions.

Enhanced-sampling methods remain valuable tools for subsequent thermodynamic analyses but were not adopted as the primary simulation strategy because they require assumptions regarding the reaction coordinates governing membrane insertion.

---
# Appendix B — Lipid Substitutions and Cyclopropane Lipid Parameterization

[← Back to README](README.md)

---

# Overview

The experimentally determined phospholipid composition of the *Escherichia coli* Gram-negative inner membrane (G-IM) includes several phospholipids containing **cyclopropanated fatty acid chains**. These cyclopropane modifications play an important role in regulating membrane fluidity, permeability, and resistance to environmental stress.

At present, however, the official Martini 3 lipid library does not provide validated coarse-grained models for cyclopropane-containing phospholipids. Consequently, several substitutions were required to construct a membrane that remains both physically realistic and fully compatible with the published Martini 3 force field.

This appendix documents the scientific rationale behind those substitutions.

---

# Native Membrane Composition

The target membrane reproduces the experimentally determined phospholipid composition of the *E. coli* Gram-negative inner membrane.

| Native lipid | Description |
|--------------|-------------|
| PMPE | 16:0/cy17:0 phosphatidylethanolamine |
| POPE | 16:0/18:1 phosphatidylethanolamine |
| QMPE | 15:0/cy17:0 phosphatidylethanolamine |
| OYPE | 18:1/18:1Δ11 (vaccenoyl) phosphatidylethanolamine |
| PMPG | 16:0/cy17:0 phosphatidylglycerol |
| PYPG | 16:0/18:1Δ11 phosphatidylglycerol |
| PVCL2 | Mixed-chain cardiolipin containing palmitoyl and vaccenoyl chains |

Only **POPE** is directly available in the Martini 3 phospholipid library without modification.

---

# The Cyclopropane Problem

Cyclopropane fatty acids are produced by enzymatic modification of unsaturated fatty acids, replacing a carbon–carbon double bond with a three-membered cyclopropane ring.

Although this modification changes only a single bond, it influences

- chain flexibility,
- lipid packing,
- membrane fluidity,
- lateral diffusion,
- permeability.

Accurately reproducing these effects requires dedicated coarse-grained parameterization.

At present, the Martini 3 lipid library contains no validated cyclopropane phospholipids. Consequently,

- PMPE,
- QMPE,
- PMPG,

cannot currently be represented explicitly.

---

# Why Custom Parameters Were Not Used

During the early stages of this project, several custom Martini lipid topologies were developed by modifying existing phospholipid definitions.

Although these topologies were chemically plausible, they were ultimately abandoned for production simulations because

- they were not validated against atomistic simulations or experimental observables;
- small bead-type changes can significantly alter membrane mechanics;
- custom parameters reduce reproducibility;
- future Martini releases are expected to provide official cyclopropane lipid models.

Accordingly, the final production simulations use **only officially parameterized Martini 3 lipids**.

The historical development of these custom topologies is documented in
[Appendix C](appendix-c-historical-notes-and-deprecated-topologies).

---

# Design Criteria

Any replacement membrane was required to satisfy the following conditions.

1. Remain in the liquid-disordered phase at 310 K.
2. Preserve the phosphatidylethanolamine/phosphatidylglycerol ratio.
3. Preserve the overall membrane charge.
4. Use only validated Martini 3 lipid parameters.

No single substitution satisfies all four criteria perfectly. The adopted membrane therefore represents a compromise between chemical fidelity and force-field reliability.

---

# Candidate Replacement Strategies

Three broad approaches were considered.

| Strategy | Status | Reason |
|----------|--------|--------|
| DPPE/DPPG substitution | Rejected | Artificially rigid membrane |
| Custom cyclopropane topologies | Rejected | Not validated |
| Official Martini 3 substitutes | Adopted | Fully validated and reproducible |

---

# Why Saturated Lipids Were Rejected

Replacing cyclopropane-containing lipids with fully saturated phospholipids initially appears attractive because neither contains a carbon–carbon double bond.

In practice, however, this substitution produces a membrane substantially more ordered than the native *E. coli* inner membrane.

DPPE and DPPG possess relatively high gel-to-fluid transition temperatures and exhibit markedly reduced lateral mobility near physiological temperatures.

Their use would therefore

- decrease membrane fluidity,
- reduce lipid diffusion,
- suppress local membrane deformation,
- increase the energetic barrier to peptide insertion.

Because spontaneous membrane insertion is the principal phenomenon investigated in this work, introducing gel-phase lipids would systematically bias the simulation.

For this reason, saturated substitutes were rejected.

---

# Adopted Lipid Substitutions

The final membrane uses the following substitutions.

| Native lipid | Martini 3 substitute | Primary justification |
|--------------|---------------------|-----------------------|
| PMPE | DVPE | Preserves a fluid PE environment using validated parameters |
| QMPE | POPE | Closest validated PE while avoiding unnecessary unsaturation |
| OYPE | DOPE | Closest validated unsaturated PE analogue |
| PMPG | POPG | Validated PG substitute |
| PYPG | POPG | Equivalent PG used consistently throughout membrane construction |
| PVCL2 | TOCL | Official Martini cardiolipin |

These substitutions preserve

- membrane charge,
- phospholipid headgroup composition,
- liquid-disordered membrane behaviour,

while remaining entirely within the published Martini 3 force field.

---

# Rationale for Individual Substitutions

## PMPE → DVPE

PMPE contains one saturated palmitoyl chain and one cyclopropanated C17 chain.

The cyclopropane ring originates biosynthetically from a cis double bond and retains many of the packing characteristics of an unsaturated chain. Consequently, replacing the cyclopropane chain with a validated vaccenoyl chain provides a closer approximation to the native membrane than replacing it with a fully saturated chain.

DVPE therefore preserves a fluid phosphatidylethanolamine environment while remaining entirely within the official Martini 3 lipid library.

---

## QMPE → POPE

QMPE differs from PMPE in an important respect: its **sn1 chain is a saturated C15:0 fatty acid**, whereas PMPE contains a C16:0 chain.

Because QMPE contains only **one** cyclopropane chain, mapping it to DVPE would introduce **two** chemically significant changes simultaneously:

1. replacing the cyclopropane ring with a cis double bond, **and**
2. replacing the native C15 saturated chain with an unsaturated vaccenoyl chain.

This would increase the overall unsaturation of the molecule beyond that present in the native lipid.

POPE, although not an exact chemical analogue, introduces only a single unsaturated chain while preserving the general phosphatidylethanolamine character of the membrane. It therefore represents the more conservative substitution and avoids overestimating membrane fluidity.

---

## OYPE → DOPE

OYPE contains oleoyl and vaccenoyl chains.

Because Martini does not distinguish between positional isomers of cis monounsaturated chains at this level of coarse graining, DOPE provides the closest validated phosphatidylethanolamine analogue available in the standard Martini 3 library.

---

## PMPG and PYPG → POPG

Current versions of COBY do not directly support PYPG.

Using POPG throughout the phosphatidylglycerol fraction simplifies membrane construction while preserving the total concentration of anionic phospholipids.

This is an implementation decision rather than a limitation of the Martini force field.

---

## PVCL2 → TOCL

The experimental membrane contains a mixed-chain cardiolipin.

Martini 3 provides validated parameters for tetraoleoyl cardiolipin (TOCL), which was adopted for the present simulations.

Only coordinate generation required a project-specific workaround; all interaction parameters are the standard Martini 3 values.

---

# Consequences of the Substitutions

The adopted membrane is not chemically identical to the experimental lipid composition.

Nevertheless, it preserves the properties most relevant to the present study.

- Overall membrane charge.
- Approximate phospholipid class distribution.
- A predominantly liquid-disordered membrane.
- Exclusive use of validated Martini 3 parameters.
- Full reproducibility using publicly available software and force fields.

The principal missing feature is the explicit representation of cyclopropane fatty acids.

---

# Future Improvements

Once validated Martini 3 parameters for cyclopropane-containing phospholipids become available, the preferred strategy will be to replace the corresponding substitute lipids while leaving the remainder of the workflow unchanged.

Because the simulation protocol relies exclusively on standard Martini infrastructure, such an update should require only minimal modification to the membrane-construction stage.

---

# Conclusions

The membrane composition adopted in this repository is not intended to be a chemically exact coarse-grained representation of the *E. coli* Gram-negative inner membrane.

Instead, it represents the closest approximation currently achievable while

- preserving membrane fluidity,
- preserving phospholipid composition,
- avoiding unvalidated force-field modifications,
- remaining fully reproducible using published Martini 3 parameters.

This conservative approach prioritizes reproducibility and force-field reliability over speculative parameter development while providing a robust framework for future refinement as the Martini lipid library evolves.

---

# Appendix C — Repository Maintenance Notes and Deprecated Topologies

[← Back to README](README.md)

---

# Overview

During the development of this workflow, several approaches were investigated for representing the experimentally observed *Escherichia coli* Gram-negative inner membrane within the Martini 3 force field.

Some of these approaches required project-specific topology modifications or coordinate-generation workarounds. As the workflow evolved, these were replaced wherever possible by officially parameterized Martini 3 lipids.

This appendix documents those development decisions and identifies files that are retained solely for historical reference.

---

# Development Strategy

The objective throughout this project was to reproduce the experimentally determined membrane composition while remaining compatible with the published Martini 3 force field.

Three practical difficulties were encountered.

1. Cyclopropane-containing phospholipids are not currently available in the official Martini 3 lipid library.
2. COBY does not currently generate coordinates correctly for all lipids required by the target membrane.
3. Some custom lipid definitions proved incompatible with COBY's internal validation routines.

These issues were addressed iteratively during the development of the workflow.

---

# Custom Lipid Topologies

To reproduce the native membrane composition as closely as possible, custom Martini 3 lipid topology files were initially constructed for the missing lipid species.

These topologies were intended only to bridge gaps in the available Martini 3 lipid library and were **not** derived from the Martini developers.

No production simulations described in this repository employ these custom topologies.

---

# COBY Compatibility Issues

During membrane construction, COBY repeatedly terminated with exceptions reporting bead-type mismatches.

Investigation showed that these failures originated from inconsistencies between the custom lipid definitions and the bead types expected by COBY.

In particular, one phosphatidylethanolamine topology (QMPE) was found to contain incorrect bead assignments.

Because COBY validates lipid definitions against its internal Martini specifications, these inconsistencies prevented successful membrane construction.

Rather than attempting to maintain an independent set of modified lipid topologies, the workflow was redesigned to use only officially parameterized Martini 3 lipids.

This decision substantially improved reproducibility while eliminating compatibility problems.

---

# Manual TOCL Coordinate Generation

One exception remains in the final workflow.

Although Martini 3 provides validated interaction parameters for tetraoleoyl cardiolipin (TOCL), current versions of COBY do not correctly generate TOCL coordinates during membrane construction.

To avoid modifying the force field itself, TOCL coordinates are generated separately using

```text
simulation_cg/generate_tocl.py
```

The generated coordinates are subsequently incorporated into the membrane during system construction.

Only coordinate generation is affected.

All interaction parameters remain those of the official Martini 3 force field.

---

# Current Production Workflow

The production simulations described in this repository use

- officially parameterized Martini 3 phospholipids;
- officially parameterized Martini 3 cardiolipin (TOCL);
- project-specific coordinate generation only where required;
- no custom interaction parameters.

Consequently, all production simulations are reproducible using publicly available Martini 3 force-field files together with the scripts provided in this repository.

---

# Deprecated Files

Some files retained within the repository document earlier stages of workflow development.

These files are preserved for transparency and reproducibility only.

They should **not** be used when preparing new simulation systems.

If future versions of Martini 3 provide validated cyclopropane phospholipids or if COBY adds native support for the missing lipid species, these historical files may be removed entirely.

---

# Recommendations for Future Development

Future versions of this workflow should preferentially adopt official Martini lipid parameterizations as they become available.

In particular,

- validated cyclopropane phospholipids should replace the current substitute lipids;
- manual TOCL coordinate generation should be removed once supported by COBY;
- any remaining project-specific workarounds should be retired in favour of upstream implementations.

This policy minimizes maintenance effort while ensuring continued compatibility with future Martini releases.

---

# Conclusions

The current workflow deliberately prioritizes reproducibility over exact chemical fidelity.

Early attempts to reproduce the native membrane using project-specific Martini lipid definitions demonstrated the practical difficulties of maintaining custom force-field components alongside evolving software such as COBY.

The final workflow therefore adopts a conservative strategy:

- use official Martini 3 interaction parameters wherever available;
- introduce project-specific code only for coordinate generation;
- avoid maintaining independent force-field modifications.

This approach provides a robust, reproducible foundation while remaining straightforward to update as the Martini ecosystem evolves.

---