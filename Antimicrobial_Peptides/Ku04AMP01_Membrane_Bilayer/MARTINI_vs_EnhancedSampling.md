# MARTINI CG vs. Enhanced Sampling Methods for KU04AMP01 Translocation

## Background

Spontaneous peptide translocation across a lipid bilayer occurs on the µs–ms timescale,
well beyond what unbiased atomistic CHARMM36 MD can reach at 150–250 ns/day on a V100 GPU.
This document compares four approaches for accessing translocation within a feasible compute
budget, with specific attention to the G-IM system composition and KU04AMP01's sequence.

---

## The Core Trade-off

All four methods sacrifice something to compress time. They differ in *what* is sacrificed:
free-energy accuracy, mechanistic detail, reaction-coordinate assumptions, or physical
realism of the force field.

---

## Method-by-Method Analysis

### 1. Steered MD (SMD / Constant-Velocity Pulling)

Pull KU04AMP01 along the membrane normal (Z) at a constant velocity using a harmonic spring
(`pull = umbrella` in GROMACS MDP).

**Pros:**
- Simplest to set up — one extra `pull` block in the MDP
- Fast: a single pulling run takes 1–5 ns
- Immediately visualises a translocation pathway

**Cons for this system:**
- Pulling velocity is always orders of magnitude faster than biology. For a 36-residue
  flexible peptide with Pro/Gly break-points, the backbone does not have time to rearrange —
  you drag a misfolded structure through the membrane.
- Work measured by Jarzynski's equality converges very slowly for a large, flexible solute.
  Dozens to hundreds of trajectories are needed for a reliable ΔG; error bars remain large.
- Severely overestimates the free energy barrier because lateral peptide–lipid rearrangement
  is suppressed.
- Best used only for generating seed configurations for umbrella sampling windows, not as a
  standalone answer.

**Verdict:** Pathway-generation tool only. Do not use for energetics.

---

### 2. Umbrella Sampling (US)

Generate a series of configurations at evenly spaced Z positions (via SMD), restrain the
peptide at each with a harmonic bias, run independent simulations at each window, then
reconstruct the potential of mean force (PMF) with WHAM or MBAR.

**Pros:**
- Gold standard for translocation PMFs — widely published for AMPs
- Statistical error is controllable and quantifiable
- Works directly with CHARMM36, preserving all atomistic detail

**Cons for this system:**
- The 1D Z reaction coordinate is valid only if the peptide inserts stiffly and axially.
  KU04AMP01 has two Pro and two Gly break-points (residues 14–21, 26) and is likely to
  tilt, bend, or lie flat at intermediate depths. A 1D Z PMF mixes physically distinct
  configurations into the same window — a well-known artefact for flexible AMPs.
- Window count: ~40–60 windows at 0.2 nm spacing × 50–200 ns each =
  **2,000–12,000 ns of total compute time** → 8–80 GPU-days on a V100.
- Hysteresis between insertion and extraction directions is common for amphipathic peptides.
- Managing 40–60 parallel `pull` jobs and converging WHAM is operationally complex.

**Verdict:** Most rigorous free-energy method available, but the 1D reaction coordinate is a
genuine scientific concern for a flexible peptide, and the compute cost is the highest of all
four approaches.

---

### 3. Metadynamics (Well-Tempered / PLUMED)

Deposit Gaussian bias hills along one or more collective variables (CVs) — typically Z
distance and optionally a tilt or secondary-structure CV — to accelerate escape from
free-energy minima. Implemented via PLUMED, natively supported by GROMACS.

**Pros:**
- Does not require a pre-defined pathway; the system finds its own route.
- Well-tempered metadynamics (WTmetaD) converges to the true free energy surface.
- 2D CVs (Z + helix tilt, or Z + buried hydrophobic SASA) can capture conformational changes
  simultaneously with translocation.
- Funnel metadynamics or parallel-tempering WTmetaD (PT-WTmetaD) improve convergence
  further.

**Cons for this system:**
- CV choice is critical. A poor CV gives a PMF that is converged but physically meaningless.
  For a 36-residue disordered peptide, no single pair of CVs fully captures conformational
  space.
- Hill-filling convergence time is hard to predict a priori: 500 ns to 5 µs depending on
  barrier height and CV quality.
- PT-WTmetaD and funnel metadynamics add significant setup complexity.
- At CHARMM36 resolution, confident PMF convergence likely still requires multi-µs runs.

**Verdict:** Scientifically the most flexible of the three enhanced-sampling atomistic
methods (especially with 2D CVs), but CV design is non-trivial for a flexible peptide
and convergence monitoring is burdensome.

---

### 4. MARTINI 3 Coarse-Grained MD (Unbiased)

Reduce the ~45,000-atom system to ~11,000 CG beads (≈4:1 heavy-atom mapping). Combined with
a 20 fs timestep and a smoother energy landscape, the effective speedup is ~100×:

$$t_{\text{CG/day}} \approx 150\text{–}250 \frac{\text{ns}}{\text{day}} \times 100 \approx 15\text{–}25\ \mu\text{s/day}$$

At 15 µs/day, a 10 µs translocation trajectory takes < 1 day. This makes spontaneous,
unforced translocation genuinely accessible.

**Pros:**
- **No reaction coordinate needed** — the single most important advantage for a flexible,
  disordered peptide.
- Captures multiple translocation events per run, providing ensemble statistics and pathway
  diversity rather than a single forced trajectory.
- Lateral lipid rearrangement, peptide bending, and membrane thinning all happen naturally.
- Can be combined with backmapping to CHARMM36 for atomic-resolution analysis of key
  intermediates.

**Cons for this system:**
- **Cyclopropane-chain lipids (cy17:0) have no official MARTINI 3 parameters** in the
  M3-Lipid-Parameters release (Pedersen et al. 2025). PMPE, QMPE, and PMPG
  (68 of 100 lipids per leaflet) cannot be used directly. See workaround below.
- MARTINI tends to underestimate the free-energy penalty for charged residues crossing the
  membrane core, so translocation barriers may be lower than in reality. The Asp and Glu
  residues on KU04AMP01 are most affected.
- Backbone secondary structure must be constrained via an elastic network (ElNeDyn)
  [(Periole et al. 2009)](https://doi.org/10.1021/ct9002114); helix-coil transitions during
  insertion are artificially restricted.
- Free energies are qualitative — reporting ΔG‡ in kJ/mol requires extensive validation.

---

## The Cyclopropane Problem and Workaround

Three of the seven G-IM lipid species contain **cy17:0** (cyclopropane) acyl chains, for
which MARTINI 3 has no validated bead type. These account for 68 of the 100 lipids per
leaflet. cy17:0 (cis-9,10-methylene-hexadecanoic acid) is a 17-carbon fatty acid bearing a
cyclopropane ring at C9–C10 derived biosynthetically from cis-vaccenic acid (18:1Δ11)
[(Grogan & Cronan 1997)](https://doi.org/10.1128/mr.61.4.429-441.1997).

| G-IM Lipid | cy17:0 species | Replace with | MARTINI 3 name |
|---|---|---|---|
| PMPE (16:0/cy17:0) | cy17:0 | dipalmitoyl (16:0/16:0) | `DPPE` |
| QMPE (15:0/cy17:0) | cy17:0 | dipentadecanoyl → dipalmitoyl (approx) | `DPPE` |
| PMPG (16:0/cy17:0) | cy17:0 | dipalmitoyl (16:0/16:0) | `DPPG` |

> **Note on QMPE sn1 chain identity**: "15:0" (pentadecanoic acid) is fully saturated;
> however the custom QMPE topology assigns a D3/double-bond bead to the sn1 chain. If
> QMPE truly contains 15:0 on sn1, the D3 bead is incorrect. This ambiguity should be
> resolved before use.

**Availability in the Singularity image's M3-Lipid-Parameters (Pedersen et al. 2025):**

| Lipid | Available? | Notes |
|---|---|---|
| POPE | **Yes** | `martini_v3.0.0_phospholipids_PE_v2.itp` |
| PYPG | **Yes** | `martini_v3.0.0_phospholipids_PG_v2.itp` (16:0/16:1∆9-PG) |
| DPPE | **Yes** | `martini_v3.0.0_phospholipids_PE_v2.itp` — **gel phase at 37°C; do not use as cy17:0 substitute** |
| DPPG | **Yes** | `martini_v3.0.0_phospholipids_PG_v2.itp` — **near Tm at 37°C; avoid** |
| DVPE | **Yes** | di-C18:1Δ11 (cis-vaccenic) PE in PE ITP |
| DVPG | **Yes** | di-C18:1Δ11 (cis-vaccenic) PG in PG ITP |
| OYPE | **No** | Not in M3-Lipid-Parameters. **YOPE** (16:1∆9/18:1∆9-PE, sn1/sn2 swapped) is available |
| PVCL2 | **Likely** | Cardiolipin ITP (`martini_v3.0.0_phospholipids_CL_v2.itp`) exists; verify exact residue name |
| PMPE, QMPE, PMPG | **No** | Not in any official release; cy17:0 unparameterized |

---

### ⚠️ Phase Behaviour Warning: DPPE/DPPG Is the Wrong Substitute for a Translocation Study

The DPPE/DPPG substitution listed above is **not safe** for a peptide translocation study.
Replacing 68% of the membrane with gel-forming lipids will give a false-negative result —
the peptide becomes trapped not because of biology, but because the membrane is
artificially frozen.

**Why DPPE is in the gel phase at simulation temperature.**
DPPE (16:0/16:0-PE) has a main-phase transition temperature Tm ≈ 63°C; DPPG Tm ≈ 40°C.
At 310 K (37°C) these lipids are solidly in the **gel (Lβ) phase**, characterised by
near-all-*trans* packed chains that exclude lateral void space. Substituting 68 of 100
lipids per leaflet with gel-forming surrogates raises the effective bilayer Tm above
simulation temperature, making the membrane:

- ~10-fold less permeable to solutes than a fluid membrane
- Resistant to the lipid-chain displacement required for peptide insertion
- Unable to support the membrane thinning and pore nucleation that AMPs exploit

The result is the peptide adsorbs but cannot insert — a false-negative that matches
experiment in outcome but for the wrong reason.

**cy17:0 is physiologically fluid at 37°C.**
The cyclopropane ring acts as a partial packing disruptor — it reduces fluidity relative
to its biosynthetic precursor (18:1Δ11, cis-vaccenic acid), but keeps the membrane well
above the gel transition [(Poger & Mark 2015)](https://doi.org/10.1021/jp5092717).
The *E. coli* inner membrane is unambiguously in the **liquid-disordered (Ld) phase** at
37°C regardless of cyclopropane content.

**Recommended substitutions that preserve fluid-phase behaviour:**

cy17:0 is derived biosynthetically from 18:1Δ11 (cis-vaccenic acid) by CFA synthase
[(Grogan & Cronan 1997)](https://doi.org/10.1128/mr.61.4.429-441.1997), so vaccenic-acid
chains are the closest biochemically motivated and thermodynamically valid surrogates:

| G-IM Lipid | Preferred substitute | Rationale | MARTINI 3 name |
|---|---|---|---|
| PMPE (16:0/cy17:0) | 16:0/18:1Δ11-PE (if available) or POPE | Vaccenic precursor; fluid at 37°C ✓ | `PVPE`\* or `POPE` |
| QMPE (?/cy17:0) | POPE (16:0/18:1Δ9-PE) | Best available fluid PE; already in composition | `POPE` |
| PMPG (16:0/cy17:0) | 16:0/18:1Δ11-PG (if available) or PYPG | Vaccenic precursor; fluid ✓ | `PVPG`\* or `PYPG` |

\* `PVPE` and `PVPG` may not exist as named entries; check the ITP files.
`DVPE` (di-18:1Δ11-PE) and `DVPG` (di-18:1Δ11-PG) are confirmed present and are also
acceptable fluid-phase surrogates.

**Fallback.** If only symmetric-chain ITPs are convenient, use **POPE** for all PE
surrogates and **PYPG** for PG surrogates — both already present in the membrane,
unambiguously fluid at 37°C.

The qualitative ranking of substitution quality:

$$\text{cy17:0-native} > \text{vaccenic (18:1}\Delta11\text{)} > \text{oleoyl (18:1}\Delta9\text{)} \gg \text{dipalmitoyl (16:0/16:0 — gel phase)}$$

Re-parametrizing cy17:0 from scratch using the MARTINI 3 small-molecule workflow
(QM target data → iterative Boltzmann inversion, Pedersen et al. 2025) is publishable
work in itself and is not recommended unless cyclopropane specificity is a primary
scientific question.

---

## Comparison Summary

| Method | Compute cost (V100) | Free-energy accuracy | Reaction coord needed | Handles peptide flexibility | cy17:0 issue |
|---|---|---|---|---|---|
| Steered MD | Very low (1–5 ns) | Poor | Yes (Z) | Poor | None |
| Umbrella sampling | Very high (8–80 GPU-days) | High (if CV is right) | Yes (Z) | Poor for flexible peptides | None |
| Metadynamics | High (0.5–5 µs) | Moderate–high | Yes (2D helps) | Moderate | None |
| **MARTINI CG** | **Low (1–2 GPU-days / 10 µs)** | Qualitative | **No** | **Good** | Substitute DPPE/DPPG |

---

## Practical Recommendation

**For a first-pass translocation study within a 1–2 week compute window**, the MARTINI CG
approach (Stage 1 alone) is the most cost-effective choice. It does not require prior
knowledge of the translocation pathway, handles KU04AMP01's flexibility naturally, and
runs to completion within 3 GPU-days. Backmapping (Stage 2) adds < 1 GPU-day per
intermediate and is always worthwhile. Stage 3 (metadynamics PMF) is recommended only if
a peer-review-quality quantitative free-energy estimate is required.

**Initial Conditions** The all-atom CHARMM-ff simulation seems to indicate adsorption in a 50 ns time window.
Choose a frame after adsorption has clearly stabilized (e.g., 25–30 ns in your case).Convert to Martini, then energy-minimize, equilibrate and run for many microseconds.
One robust suggestion is to make multiple runs (each with different random initial velocities) and get the fraction of trajectories that do translocate.

# Coarse-Grained MD Setup: *E. coli* Inner Membrane with KU04 Peptide
This repository contains the workflow, coordinate construction steps, and custom topology topologies required to simulate a physiologically accurate, seven-lipid *Escherichia coli* inner membrane containing an embedded 36-residue KU04 peptide using the **Martini 3 Coarse-Grained Force Field**.
---
## 🛠️ Step-by-Step Execution Workflow
### Step 1: Coarse-Grain the KU04 Peptide

First, clean out the hydrogens in the peptide. They're lost while coarse-graining anyways

```bash
grep -v -E " H   | H[A-Z]| OXT" KU04AMP01.pdb > KU04AMP01_clean.pdb
```

Align the peptide so that it is properly aligned along the z axis. This alignment is crucial for the subsequent simulation. The N-terminus must be facing downward towards the phospholipids, as the opposite charges will attract, facilitating membrane insertion.

We can do this using the [OPM web server](https://opm.phar.umich.edu/ppm_server2_cgopm). There, upload the clean PDB and hit submit  after setting setting `Topology (N-ter)` to `in` *i.e.* down. You’ll be directed to a page where you’ll need to wait a few minutes. Eventually, you’ll be able to download the oriented protein structure, which you can safely overwrite into the `KU04AMP01_clean.pdb` locally. Remove the dummy membrane atoms that OPM added by 

```bash
grep -v HETATM KU04AMP01_clean.pdb > temp.pdb
mv temp.pdb KU04AMP01_clean.pdb
```

Next, run `martinize2` to coarse-grain the peptide.

```bash
martinize2 -f KU04AMP01_clean.pdb \
           -dssp $MKDSSP \
           -nt -cys auto -p backbone -sep \
           -ff martini3001 -elastic \
           -x KU04AMP01_cg.pdb -o KU04.top \
           -scfix -ef 500.0 -el 0.5 -eu 0.9 -ea 0 -ep 0 -maxwarn 5
```

This will create a bare topology file KU04.top that transcludes the CG-itp file. If you want, you can systematically change the default name (molecule_0) to the peptide name (KU04).

Change the transclusion of the martini itp file to 

```bash
#include <martini-ff-dir>/martini_v3.0.0.itp
```
Note that you need not include the full absolute path if the directory is in `$PWD` or the standard gromacs topology directory. For visualization, however, it is best to put the absolute path.

For visualization, install `martiniglass` outside the sif image (didn;t wanna bloat it) and folow the instructions [in their docs](https://martiniglass.readthedocs.io/). Note that you can supply `cg_bonds` with your martinize2 topology file and it'll work fine.

### Step 2: Build Membrane Coordinates with COBY

To prevent coordinate construction mismatches, configure a Python script (`memb_build.py`) to bypass COBY's outdated 11-bead internal library for `PYPG` and specify the explicit 12-bead structure layout alongside your project's exact target lipid ratios.

Create a file named `memb_build.py` with the following content:

```python
#!/usr/bin/env python
import COBY

sysname = "ku04amp01_gim"

# Custom structural templates to force correct bead layouts matching Martini 3 ITPs
custom_lipid_definitions = {
    "PYPG": {
        "headgroup": ["GL0", "PO4", "GL1", "GL2"],
        "tails": [["C1A", "C2A", "D3A", "C4A"], ["C1B", "C2B", "C3B", "C4B"]]
    },
    "PVCL2": {
        "headgroup": ["GL0", "PO41", "GL11", "GL12", "PO42", "GL21", "GL22"],
        "tails": [
            ["C1A", "C2A", "C3A", "C4A"], ["C1B", "C2B", "D3B", "C4B"],
            ["C1C", "C2C", "C3C", "C4C"], ["C1D", "C2D", "D3D", "C4D"]
        ]
    }
}

COBY.COBY(
    box =, 
    custom_lipids = custom_lipid_definitions,
    membrane = "lipid:PMPE:46:charge:lib lipid:POPE:13 lipid:QMPE:12:charge:lib lipid:OYPE:8:charge:lib lipid:PMPG:10:charge:lib lipid:PYPG:9:charge:lib lipid:PVCL2:2:charge:lib apl:0.5",
    protein = "file:KU04AMP01_cg.pdb moleculetypes:KU04",
    solvation = "default",
    itp_input = "file:top_for_COBY.itp",
    out_sys = sysname,
    out_top = sysname + ".top",
    out_log = sysname + ".log",
    sn = sysname,
)
```

Execute the builder:

```bash
python3 memb_build.py
```

### Step 3: Append Missing Topologies to `top_for_COBY.itp`

**Warning: Do not do this! See notes below.**

COBY is strictly a coordinate generator and does not produce physics data. To clear GROMACS `No such moleculetype` exceptions, append the custom topologies generated below directly to your `top_for_COBY.itp` file.

```bash
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CUSTOM BENCHMARK E. COLI INNER MEMBRANE PARAMETERS (MARTINI 3)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

[ moleculetype ]
PMPE              3
[ atoms ]
  1  Q1   1     PMPE    NH3  1    0      74
  2  Q5   1     PMPE    PO4  2    0      74
  3  P1   1     PMPE    GL1  3    0      74
  4  P1   1     PMPE    GL2  4    0      74
  5  C1   1     PMPE    C1A  5    0      74
  6  C1   1     PMPE    C2A  6    0      74
  7  C1   1     PMPE    C3A  7    0      74
  8  C1   1     PMPE    C4A  8    0      74
  9  C1   1     PMPE    C1B  9    0      74
 10  C1   1     PMPE    C2B  10   0      74
 11  C1   1     PMPE    C3B  11   0      74
[ bonds ]
  1  2   1   0.47  12500
  2  3   1   0.47  12500
  3  4   1   0.37  25000
  3  5   1   0.47  12500
  5  6   1   0.47   7500
  6  7   1   0.47   7500
  7  8   1   0.47   7500
  4  9   1   0.47  12500
  9 10   1   0.47   7500
 10 11   1   0.47   7500
[ angles ]
  1  2  3   2   120   25
  2  3  4   2   120   25
  3  4  9   2   180   25
  4  3  5   2   180   25
  3  5  6   2   180   25
  5  6  7   2   180   25
  6  7  8   2   180   25
  4  9 10   2   180   25
  9 10 11   2   180   25

[ moleculetype ]
QMPE              3
[ atoms ]
  1  Q1   1     QMPE    NH3  1    0      74
  2  Q5   1     QMPE    PO4  2    0      74
  3  P1   1     QMPE    GL1  3    0      74
  4  P1   1     QMPE    GL2  4    0      74
  5  C1   1     QMPE    C1A  5    0      74
  6  D3   1     QMPE    D2A  6    0      74
  7  C1   1     QMPE    C3A  7    0      74
  8  C1   1     QMPE    C4A  8    0      74
  9  C1   1     QMPE    C1B  9    0      74
 10  C1   1     QMPE    C2B  10   0      74
 11  C1   1     QMPE    C3B  11   0      74
[ bonds ]
  1  2   1   0.47  12500
  2  3   1   0.47  12500
  3  4   1   0.37  25000
  3  5   1   0.47  12500
  5  6   1   0.47   7500
  6  7   1   0.47   7500
  7  8   1   0.47   7500
  4  9   1   0.47  12500
  9 10   1   0.47   7500
 10 11   1   0.47   7500
[ angles ]
  1  2  3   2   120   25
  2  3  4   2   120   25
  3  4  9   2   180   25
  4  3  5   2   180   25
  3  5  6   2   180   25
  5  6  7   2   120   45
  6  7  8   2   180   25
  4  9 10   2   180   25
  9 10 11   2   180   25

[ moleculetype ]
OYPE              3
[ atoms ]
  1  Q1   1     OYPE    NH3  1    0      74
  2  Q5   1     OYPE    PO4  2    0      74
  3  P1   1     OYPE    GL1  3    0      74
  4  P1   1     OYPE    GL2  4    0      74
  5  C1   1     OYPE    C1A  5    0      74
  6  C1   1     OYPE    C2A  6    0      74
  7  D3   1     OYPE    D3A  7    0      74
  8  C1   1     OYPE    C4A  8    0      74
  9  C1   1     OYPE    C5A  9    0      74
 10  C1   1     OYPE    C1B  10   0      74
 11  C1   1     OYPE    C2B  11   0      74
 12  D3   1     OYPE    D3B  12   0      74
 13  C1   1     OYPE    C4B  13   0      74
[ bonds ]
  1  2   1   0.47  12500
  2  3   1   0.47  12500
  3  4   1   0.37  25000
  3  5   1   0.47  12500
  5  6   1   0.47   7500
  6  7   1   0.47   7500
  7  8   1   0.47   7500
  8  9   1   0.47   7500
  4 10   1   0.47  12500
 10 11   1   0.47   7500
 11 12   1   0.47   7500
 12 13   1   0.47   7500
[ angles ]
  1  2  3   2   120   25
  2  3  4   2   120   25
  3  4 10   2   180   25
  4  3  5   2   180   25
  3  5  6   2   180   25
  5  6  7   2   180   25
  6  7  8   2   120   45
  7  8  9   2   180   25
  4 10 11   2   180   25
 10 11 12   2   180   25
 11 12 13   2   120   45

[ moleculetype ]
PMPG              3
[ atoms ]
  1  P4   1     PMPG    GL0  1    0      74
  2  QA   1     PMPG    PO4  2   -1      74
  3  P1   1     PMPG    GL1  3    0      74
  4  P1   1     PMPG    GL2  4    0      74
  5  C1   1     PMPG    C1A  5    0      74
  6  C1   1     PMPG    C2A  6    0      74
  7  C1   1     PMPG    C3A  7    0      74
  8  C1   1     PMPG    C4A  8    0      74
  9  C1   1     PMPG    C1B  9    0      74
 10  C1   1     PMPG    C2B  10   0      74
 11  C1   1     PMPG    C3B  11   0      74
[ bonds ]
  1  2   1   0.47  12500
  2  3   1   0.47  12500
  3  4   1   0.37  25000
  3  5   1   0.47  12500
  5  6   1   0.47   7500
  6  7   1   0.47   7500
  7  8   1   0.47   7500
  4  9   1   0.47  12500
  9 10   1   0.47   7500
 10 11   1   0.47   7500
[ angles ]
  1  2  3   2   120   25
  2  3  4   2   120   25
  3  4  9   2   180   25
  4  3  5   2   180   25
  3  5  6   2   180   25
  5  6  7   2   180   25
  6  7  8   2   180   25
  4  9 10   2   180   25
  9 10 11   2   180   25

[ moleculetype ]
PVCL2             3
[ atoms ]
  1  P4   1     PVCL2   GL0  1    0      74
  2  QA   1     PVCL2   PO41 2   -1      74
  3  P1   1     PVCL2   GL11 3    0      74
  4  P1   1     PVCL2   GL12 4    0      74
  5  QA   1     PVCL2   PO42 5   -1      74
  6  P1   1     PVCL2   GL21 6    0      74
  7  P1   1     PVCL2   GL22 7    0      74
  8  C1   1     PVCL2   C1A  8    0      74
  9  C1   1     PVCL2   C2A  9    0      74
 10  C1   1     PVCL2   C3A  10   0      74
 11  C1   1     PVCL2   C4A  11   0      74
 12  C1   1     PVCL2   C1B  12   0      74
 13  C1   1     PVCL2   C2B  13   0      74
 14  D3   1     PVCL2   D3B  14   0      74
 15  C1   1     PVCL2   C4B  15   0      74
 16  C1   1     PVCL2   C1C  16   0      74
 17  C1   1     PVCL2   C2C  17   0      74
 18  C1   1     PVCL2   C3C  18   0      74
 19  C1   1     PVCL2   C4C  19   0      74
 20  C1   1     PVCL2   C1D  20   0      74
 21  C1   1     PVCL2   C2D  21   0      74
 22  D3   1     PVCL2   D3D  22   0      74
 23  C1   1     PVCL2   C4D  23   0      74
[ bonds ]
  1  2   1   0.47  12500
  2  3   1   0.47  12500
  3  4   1   0.37  25000
  1  5   1   0.47  12500
  5  6   1   0.47  12500
  6  7   1   0.37  25000
  3  8   1   0.47  12500
  8  9   1   0.47   7500
  9 10   1   0.47   7500
 10 11   1   0.47   7500
  4 12   1   0.47  12500
 12 13   1   0.47   7500
 13 14   1   0.47   7500
 14 15   1   0.47   7500
  6 16   1   0.47  12500
 16 17   1   0.47   7500
 17 18   1   0.47   7500
 18 19   1   0.47   7500
  7 20   1   0.47  12500
 20 21   1   0.47   7500
 21 22   1   0.47   7500
 22 23   1   0.47   7500
[ angles ]
  2  1  5   2   120   25
  1  2  3   2   120   25
  2  3  4   2   120   25
  1  5  6   2   120   25
  5  6  7   2   120   25
  3  4 12   2   180   25
  4  3  8   2   180   25
  3  8  9   2   180   25
  8  9 10   2   180   25
  9 10 11   2   180   25
  4 12 13   2   180   25
 12 13 14   2   180   25
 13 14 15   2   120   45
  6  7 20   2   180   25
  7  6 16   2   180   25
  6 16 17   2   180   25
 16 17 18   2   180   25
 17 18 19   2   180   25
  7 20 21   2   180   25
 20 21 22   2   180   25
 21 22 23   2   120   45
```

---

### ⚠️ Critical Evaluation of the Custom Topologies Above

The five ITP blocks provided in Step 3 (PMPE, QMPE, OYPE, PMPG, PVCL2) are **experimental
placeholders** for lipids that have no official Martini 3 parameters. Before using them in
production simulations, the following deficiencies must be understood and, where possible,
corrected.

---

#### 1. Wrong Non-Bonded Bead Types — Most Critical

Every bead in these topologies uses **Martini 2-era bead type names** that are not
equivalent in the Martini 3 non-bonded interaction matrix. Using them with the
`martini_v3.0.0.itp` installed in the Singularity image will silently apply wrong
van der Waals and electrostatic parameters for every bead.

| Structural group | Type used in README | Correct Martini 3 type | Reference |
|---|---|---|---|
| PE headgroup amine (NH3) | `Q1` | `Q4p` (+1) | Pedersen et al. 2025 |
| Phosphate (PO4), PE/PG | `Q5` at **charge 0** | `Q5` at **charge −1** | Pedersen et al. 2025 |
| Glycerol ester linker (GL1/GL2), PE | `P1` | `SN4a` | Pedersen et al. 2025 |
| Glycerol headgroup (GL0), PG | `P4` | `P4r` | Pedersen et al. 2025 |
| Phosphate (PO4), PG | `QA` (−1) | `Q5` (−1) | Pedersen et al. 2025 |
| cis double-bond carbon | `D3` | `C4h` | Pedersen et al. 2025 |

In particular, the PE zwitterion is broken: NH3 and PO4 are both assigned charge 0
(net 0, but no electrostatic dipole), whereas correct Martini 3 PE has NH3 = **+1** and
PO4 = **−1**. The `Q1` bead used for NH3 is a generic polar charged bead in Martini 3
with significantly different σ/ε from the intended `Q4p`.

To verify the correct parameters, compare against POPE in
`martini_v3.0.0_phospholipids_PE_v2.itp` and PYPG in
`martini_v3.0.0_phospholipids_PG_v2.itp`, both installed at
`/usr/local/share/gromacs/top/martini3.ff/` inside the Singularity image
[(Pedersen et al. 2025)](https://doi.org/10.1021/acscentsci.5c00755).

---

#### 2. cy17:0 Acyl Chain is Under-beaded (PMPE, QMPE, PMPG)

All three cyclopropane-containing lipids assign **3 `C1` beads** (C1B–C2B–C3B) to
the cy17:0 sn2 chain. cy17:0 has 17 carbons; the standard Martini 3 4:1 mapping
requires **4 beads** to represent a 16–17-carbon chain.

| Lipid | cy17:0 bead count | Expected | Effective chain carbons represented |
|---|---|---|---|
| PMPE (sn2) | 3 | 4 | ~12 of 17 |
| PMPG (sn2) | 3 | 4 | ~12 of 17 |
| QMPE (sn2) | 3 | 4 | ~12 of 17 |

A 3-bead sn2 chain is ~4 Å shorter than a 16:0 palmitoyl sn2, biasing the bilayer
toward artificial asymmetric thickness.

---

#### 3. Cyclopropane Ring is Entirely Absent

cy17:0 contains a **three-membered carbocycle** (cyclopropane ring at C9–C10 with a
methylene bridge, total 17 carbons). This ring structure:

- **Restricts C–C–C–C dihedral rotation** around C8–C9 and C10–C11, stiffening the chain
  mid-section. There is no dihedral or constraint encoding this in the `[ angles ]` blocks.
- **Adds steric bulk** (methylene bridge projects out of the chain plane), slightly
  increasing effective cross-sectional area relative to a linear 17-carbon chain.
- **Modulates phase behaviour**: in *E. coli*, cyclopropane FA formation from cis-vaccenic
  acid (18:1Δ11) by CFA synthase reduces membrane fluidity at the cost of a shorter
  effective acyl length, aiding acid and stationary-phase stress resistance
  [(Grogan & Cronan 1997)](https://doi.org/10.1128/mr.61.4.429-441.1997).

None of this is captured by three linear C1 beads with 180°/25 kJ mol⁻¹ angles. The
physically closest available Martini 3 alternative for an exploratory study is the
straight 16:0 palmitoyl chain (DPPE/DPPG substitution from the workaround table above),
which at least preserves chain length.

A proper cyclopropane CG model would require treating the ring carbons as an SC-type
(small ring) pseudo-particle following the Martini 3 ring-building rules detailed in
[(Souza et al. 2021)](https://doi.org/10.1038/s41592-021-01098-3), with QM-derived
bonded parameters validated against atomistic reference data — this constitutes original
parameterization work.

---

#### 4. QMPE sn1 Chain: D3 Bead on a Saturated Odd-Chain Acid

QMPE is labeled "15:0/cy17:0" (sn1 = pentadecanoic acid, fully saturated). Yet atom 6
of the sn1 chain carries bead type `D3` (double-bond bead), which is physically incorrect
for a saturated chain. The kink angle 120°/45 kJ mol⁻¹ placed at positions 6–7–8
encodes a cis double bond that does not exist in 15:0. If QMPE's sn1 is truly a
saturated C15:0 chain, all four beads should be `C1` type with 180°/25 kJ mol⁻¹ angles
throughout.

If "Q" in the lipid name instead encodes a *monounsaturated* odd-chain acid (e.g.,
15:1Δ9 or 14:1Δ9-methyl-branched), the D3 placement at bead 2 corresponds to a Δ5–Δ7
double bond, inconsistent with conventional 15:1Δ9 (Δ9 in bead 3 for a 4-bead
15-carbon chain). This ambiguity must be resolved before the topology can be used.

---

#### 5. OYPE sn1 Chain Has 5 Beads (20:1 equivalent), Not 18:1

The OYPE sn1 tail (C1A–C2A–D3A–C4A–C5A, 5 beads with D3 at position 3) corresponds in
4:1 Martini mapping to a **C20:1Δ11** (gondoic acid) chain, not 18:1Δ9 (oleoyl).
For reference, the official POPE in M3-Lipid-Parameters maps the 18:1Δ9 oleoyl sn2
chain to **4 beads** with altail pattern CDCC (D at position 2). If OYPE represents
18:1(sn1)/16:1(sn2)-PE as implied by the name, the sn1 chain should use 4 beads. A
5-bead sn1 inflates the effective chain length, biasing leaflet thickness and area-per-
lipid. If the intent was truly a 20:1Δ11 chain on sn1, the name OYPE should be revised.

The **closest officially available Martini 3 substitute** for OYPE (18:1/16:1-PE) is
`YOPE` (16:1Δ9/18:1Δ9-PE, sn1/sn2 positions swapped), found in
`martini_v3.0.0_phospholipids_PE_v2.itp`. Alternatively, `POPE` (16:0/18:1-PE) is
fully validated.

---

#### 6. PVCL2 Atom-Name Typo

Atom 19 in the PVCL2 topology is named `C4A` but should be **`C4C`** (it belongs to
the C-chain acyl tail). This has been corrected in the topology above; verify that your
working copy reflects the fix before running `gmx grompp`.

---

#### Recommended Path Forward

| Option | Description | Phase risk | Fidelity |
|---|---|---|---|
| **A (wrong for translocation)** | Replace PMPE/QMPE/PMPG → DPPE/DPPG | **GEL at 37°C — false negative** | Low |
| **B (correct minimum)** | Replace cy17:0 lipids → vaccenic-acid equivalents (POPE/DVPE/DVPG/PYPG) with corrected Martini 3 bead types | Fluid ✓ | Medium |
| **C (rigorous)** | Full QM-guided Martini 3 parameterization of cy17:0, validated against CHARMM36 reference MD | Fluid ✓ | High — publishable |

For a first-pass translocation study, **Option B** (vaccenic-acid substitution with correct bead types) is
the minimum scientifically defensible choice. Option A with DPPE/DPPG will produce a
membrane that is too ordered, artificially suppressing translocation.

---

### Step 4: Map Temperature Groups in GROMACS

Create a comprehensive index file combining your peptide and all 6 separate lipid types into a single cohesive temperature-coupling block to ensure uniform heating profiles.

```bash
gmx make_ndx -f ku04amp01_gim.gro -o index.ndx <<EOF
"KU04" | "PMPE" | "POPE" | "QMPE" | "OYPE" | "PMPG" | "PYPG" | "PVCL2"
name 21 KU04_Membrane
q
EOF
```
## Step 5: Adjust MDP Parameters

Audit your .mdp file configuration fields. Ensure the integration variables are altered to match your updated group names:
gmx tc-grps = KU04_Membrane solvent ref-t = 310 310 gen-t = 310 pcoupltype = semiisotropic 

## Step 6: Assemble and Run the Compiler

Compile the target workspace into a runtime binary (.tpr) layout:

```bash 
gmx grompp -f minimization.mdp \ -c ku04amp01_gim.gro \ -p ku04amp01_gim.top \ -n index.ndx \ -o min.tpr 
```

------------------------------

## 🧬 Scientific Parameter Rationale

> **⚠️ Note on bead-type nomenclature**: The rationale below describes the *physical intent*
> behind bead assignments. The actual bead type labels used in the custom ITP blocks
> (Q1, P1, QA, D3) are **Martini 2-era names** and do not map to the correct Martini 3
> interaction levels. For correct Martini 3 simulations with the installed M3-Lipid-Parameters
> [(Pedersen et al. 2025)](https://doi.org/10.1021/acscentsci.5c00755), replace them as shown
> in the Critical Topology Evaluation table above.
> The standard Martini 3 parameters ([Souza et al. 2021](https://doi.org/10.1038/s41592-021-01098-3))
> should be consulted for all bead-type decisions.

The parameters provided above are derived directly from the modular design rules of Martini 3. Rather than guessing values, the force field assigns physical behaviors based on specific structural fragments:
## 1. Choice of Non-Bonded Bead Subtypes (The [ atoms ] block)

* PE Zwitterionic Heads: The primary headgroup uses Q1 (highly polar amine) and Q5 (polar zwitterionic phosphate pair) beads, reproducing a net neutral charge. **Correction**: Martini 3 (Pedersen et al. 2025) uses `Q4p` (+1) for the amine and `Q5` (−1) for the phosphate. Net charge is zero, but the individual bead charges create the correct dipole moment. Using uncharged Q1/Q5 as in the custom ITPs removes this dipole entirely.
* PG Anionic Heads: Phosphatidylglycerols feature an outer neutral glycerol linker (P4) coupled with an explicit charged unshielded phosphate bead (QA), assigning a net negative charge (-1). **Correction**: Martini 3 uses `P4r` for the PG glycerol head and `Q5` (−1) for the phosphate.
* Saturated Tails (Palmitoyl / Myristoyl): Mapped with regular C1 beads representing 4-carbon saturated segments. Palmitoyl (16:0) is assigned exactly 4 beads (C1A-C4A), while Myristoyl (14:0) is assigned 3 beads (C1B-C3B). **Note**: The first tail bead connected to the glycerol linker should be `SC1` (small C1) in Martini 3, as verified in the DPPE and DPPG topologies from M3-Lipid-Parameters.
* Unsaturated Tails (Oleoyl / Vaccenoyl / Palmitoleoyl): Every cis double bond uses a special D3 bead (less attractive and more sterically volume-occupying). The kink location changes by bead order to capture distinct tail interactions (e.g., D2A for Palmitoleoyl vs. D3A for Oleoyl). **Correction**: Martini 3 uses `C4h` instead of `D3` for cis double bond beads. The bead position for ∆9 in 18:1 oleoyl is bead **2** (CDCC pattern in 4 beads), not bead 3; and for 16:1∆9 palmitoleoyl it is bead **3** (CCDC pattern), as confirmed in POPE and DYPE in the M3-Lipid-Parameters ITPs.

## 2. Bonded Constraints (The [ bonds ] and [ angles ] blocks)

* Headgroup Equilibrium Lengths: Backbones use standard Martini phosphate-to-glycerol lengths of 0.47 nm and an aggressive force constant of 12,500 kJ/mol/nm². Linkers bridging the glycerol groups (GL1 to GL2) use a tighter constraint of 0.37 nm and 25,000 kJ/mol/nm² to retain chemical proximity. **Note**: The M3-Lipid-Parameters (Pedersen et al. 2025) uses named bond-type references (e.g., `b_NH3_PO4_def`, `b_GL_GL_glyc`) from `martini_v3.0.0_ffbonded_v2.itp` rather than explicit numerical constants; the values are similar but not identical to those in the custom ITPs.
* Acyl Tail Inter-bead Distance: Saturated and unsaturated tail links use standard Martini spacing weights equal to 0.47 nm with structural tracking constants set to 7,500 kJ/mol/nm².
* Chain Angles: Saturated backbones use a straight equilibrium parameter configuration (180° with a force of 25 kJ/mol) to mimic an extended configuration. Double-bond segments containing the D3 entity use a rigid equilibrium value of 120° with an increased force constraint of 45 kJ/mol to lock in the structural cis kink that drives membrane fluidity. **Note for cyclopropane chains**: The cyclopropane ring in cy17:0 constrains C–C–C–C dihedral rotation at C8–C9 and C10–C11 in a way fundamentally different from a simple linear kink. No standard angle-based representation captures this; the 180°/25 kJ mol⁻¹ angles used for cy17:0 beads in the custom PMPE/PMPG/QMPE topologies treat the ring as if it were a straight saturated chain, which is physically incorrect [(Grogan & Cronan 1997)](https://doi.org/10.1128/mr.61.4.429-441.1997).

---

## References

1. **Martini 3 Force Field**  
   Souza, P.C.T., Alessandri, R., Barnoud, J. et al.  
   *Martini 3: a general purpose force field for coarse-grained molecular dynamics.*  
   Nature Methods **18**, 382–388 (2021). <https://doi.org/10.1038/s41592-021-01098-3>

2. **M3-Lipid-Parameters (Martini 3 expanded lipidome)**  
   Pedersen, K.B., Ingólfsson, H.I., Ramirez-Echemendia, D.P. et al.  
   *The Martini 3 Lipidome: Expanded and Refined Parameters Improve Lipid Phase Behavior.*  
   ACS Central Science (2025). <https://doi.org/10.1021/acscentsci.5c00755>

3. **Cyclopropane fatty acids in *E. coli***  
   Grogan, D.W. & Cronan, J.E.  
   *Cyclopropane ring formation in membrane lipids of bacteria.*  
   Microbiology and Molecular Biology Reviews **61**, 429–441 (1997).  
   <https://doi.org/10.1128/mr.61.4.429-441.1997>

4. **Cyclopropane FA membrane biophysics (CHARMM36 all-atom reference)**  
   Poger, D. & Mark, A.E.  
   *A Ring to Rule Them All: The Effect of Cyclopropane Fatty Acids on the Fluidity of Lipid Bilayers.*  
   Journal of Physical Chemistry B **119**, 5487–5495 (2015).  
   <https://doi.org/10.1021/jp5092717>

5. **ElNeDyn elastic network for Martini**  
   Periole, X., Cavalli, M., Marrink, S.J. & Ceruso, M.A.  
   *Combining an Elastic Network With a Coarse-Grained Molecular Force Field: Structure, Dynamics, and Intermolecular Recognition.*  
   Journal of Chemical Theory and Computation **5**, 2531–2543 (2009).  
   <https://doi.org/10.1021/ct9002114>

6. **martinize2 / vermouth**  
   Kroon, P.C., Grunewald, F., Barnoud, J. et al.  
   *Martinize2 and Vermouth: Unified Framework for Topology Construction.*  
   eLife **12**, RP90627 (2023). <https://doi.org/10.7554/eLife.90627>

7. **COBY membrane builder**  
   *COBY: CG Membrane Builder for GROMACS.*  
   <https://pypi.org/project/COBY/>

8. **cg2at backmapping**  
   Vickery, O.N. & Stansfeld, P.J.  
   *CG2AT2: an Enhanced Fragment-Based Approach for Serial Multi-scale Molecular Dynamics Simulations.*  
   Journal of Chemical Theory and Computation **17**, 6472–6482 (2021).  
   <https://doi.org/10.1021/acs.jctc.1c00652>