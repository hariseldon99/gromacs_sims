Good news: the docking has actually gone well on the whole. Let me walk through the evidence systematically.

## Summary of Docking Run
The main issue that arose was the palladium atom in the ligand, which autodock4Zn does not support by default. A new [AD4Zn.dat](AD4Zn.dat) file was created from the default one used by `AMDock` with a line added for Palladium (See comment in the dat file), then loaded manually using the `AMDock` gui.

Additionally, the protein had to be manually cleaned of water and malformatted comments with

```bash
$ grep -v "BIOVIA\|REMARK.*BIOVIA" 1LML.pdb > 1LML_clean.pdb
```

Default parameters yielded critically low sampling — only 6 GA runs

The docking log showed exactly 6 conformations across 2 clusters:
```
Number of conformations in this cluster = 5   (cluster 1)
Number of conformations in this cluster = 1   (cluster 2)
```


AutoDock4's default is `ga_run 10`; best practice for a metalloprotein active site is **25–50 runs**. With only 6 runs, the conformational space is severely undersampled — you cannot trust that the global minimum has been found. The number of runs can be set in the configuration tab of `AMDOCK`. This was done during the final docking run.


## What the log shows: excellent convergence

This is the best-quality run so far across all three sessions:

```
Cluster 1: -8.89 kcal/mol  | 21/25 runs  | ######################
Cluster 2: -7.89 kcal/mol  |  3/25 runs  | ###
Cluster 3: -6.89 kcal/mol  |  1/25 run   | #

Information entropy = 0.16   (vs 0.70 in previous run)
```

**21 out of 25 runs converge to the same pose**, with internal cluster RMSD < 0.55 Å. Information entropy of 0.16 is close to the theoretical minimum (0 = all runs identical). The AD4ZnPd.dat file is confirmed active at line 77. This is a trustworthy result.

---

## Improvements vs previous runs

| Parameter | Original | Improved | Pd ligand |
|---|---|---|---|
| GA runs | 6 | 10 | **25** |
| Box size (Å) | 18 | 18 | **27** |
| Clusters | 2 | 6 | **3** |
| Dominant cluster population | 83% (5/6, artefact) | 10% | **84% (21/25, real)** |
| Entropy | ~0.65 | 0.70 | **0.16** |
| Best ΔG (kcal/mol) | −6.54 | −6.28 | **−8.89** |

---

## The "distortion" is physically meaningful, not an artefact

Looking at atom 11 in the best pose:

```
ATOM  11  N00  NA  (20.901, 42.397, 17.746)  vdW = -4.17
```

The protein Zn²⁺ sits at (19.396, 43.866, 16.941). The distance is:

$$d = \sqrt{(1.50)^2 + (1.47)^2 + (0.81)^2} \approx 2.25 \text{ Å}$$

This NA-type nitrogen is **coordinating directly to the protein Zn²⁺**, at a distance consistent with a Zn–N coordinate bond (~2.1 Å). The large −4.17 kcal/mol vdW on atom 11 is **not a clash** — it is the TZ map attraction being correctly triggered by the AutoDockZn force field for an NA nitrogen approaching Zn. This is the intended mechanism of the AD4Zn potential.

The Pd centre in the best pose sits at (21.696, 39.902, 19.187), giving a Pd–Zn separation of ~5.1 Å. The physical picture is:

```
[Protein Zn²⁺] ←~2.25Å→ [NA nitrogen] ← (NHC core, rigid) → [Pd²⁺] ← flexible arms →
```

The visual "distortion" you see in the screenshot is the **two pendant aminophenyl arms** (the 4 flexible torsions) adopting different conformations as they fill the binding pocket independently. The Pd–NHC core itself is treated as a rigid ROOT block throughout the docking and is not deforming.

---

## Remaining issues

**1. Pd charge = 0.000 (significant)**

```
ATOM  1  Pd02  ...  vdW = -0.02  Elec = -0.00  q = +0.000
```

Pd²⁺ carries a formal charge of +2. The zero RESP charge means the entire electrostatic contribution of the metal is missing, which is why the electrostatic energy is slightly *repulsive* (+0.07 kcal/mol) rather than attractive. For the docking score this is a secondary effect since vdW dominates, but it means the binding energy is likely **underestimated** — the real interaction with the negatively charged active site residues (HIS264, GLU265) would be more favourable.

To fix: assign a RESP-derived partial charge to Pd in the PDBQT file. For Pd(II) in a square-planar NHC complex, a reasonable starting value from the literature is **+0.4 to +0.6 e** (the charge is partially donated back by the NHC carbons). Edit the PDBQT directly:

```
ATOM  1  Pd02  UNK  ...  +0.500  Pd
```

**2. Still only 25 GA runs**

For a 41-heavy-atom ligand with 4 torsions, 25 runs is sufficient given the excellent convergence seen. But for a publication, increasing to 50 runs is advisable to confirm the dominant cluster holds.

**3. Ligand efficiency is low (−0.22 kcal/mol/HA)**

A LE of −0.22 for 41 heavy atoms is modest. The −8.89 kcal/mol binding energy is largely driven by the large NHC+Pd scaffold filling the pocket, rather than high per-atom efficiency. This is typical for metal-chelating ligands but worth noting if optimising the scaffold.
