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
- **Cyclopropane-chain lipids (cy17:0) have no official MARTINI 3 parameters.** PMPE, QMPE,
  and PMPG (68 of 100 lipids per leaflet) cannot be used directly. See workaround below.
- MARTINI tends to underestimate the free-energy penalty for charged residues crossing the
  membrane core, so translocation barriers may be lower than in reality. The Asp and Glu
  residues on KU04AMP01 are most affected.
- Backbone secondary structure must be constrained via an elastic network (ElNeDyn);
  helix-coil transitions during insertion are artificially restricted.
- Free energies are qualitative — reporting ΔG‡ in kJ/mol requires extensive validation.

---

## The Cyclopropane Problem and Workaround

Three of the seven G-IM lipid species contain **cy17:0** (cyclopropane) acyl chains, for
which MARTINI 3 has no validated bead type. These account for 68 of the 100 lipids per
leaflet:

| G-IM Lipid | cy17:0 species | Replace with | MARTINI 3 name |
|---|---|---|---|
| PMPE (16:0/cy17:0) | cy17:0 | dipalmitoyl (16:0/16:0) | `DPPE` |
| QMPE (15:0/cy17:0) | cy17:0 | dipalmitoyl (approximate) | `DPPE` |
| PMPG (16:0/cy17:0) | cy17:0 | dipalmitoyl (16:0/16:0) | `DPPG` |

Keep POPE, OYPE, PYPG, and PVCL2 as-is — all have validated MARTINI 3 parameters. This
substitution preserves the charge balance (PE-dominated, ~10% PG, ~2% CL) and the
qualitative membrane physics, while sacrificing the specific packing effects of the
cyclopropane ring.

Re-parametrizing cy17:0 from scratch using the MARTINI 3 small-molecule workflow (QM target
data → iterative Boltzmann inversion) is publishable work in itself and is not recommended
unless cyclopropane specificity is a primary scientific question.

---

## Comparison Summary

| Method | Compute cost (V100) | Free-energy accuracy | Reaction coord needed | Handles peptide flexibility | cy17:0 issue |
|---|---|---|---|---|---|
| Steered MD | Very low (1–5 ns) | Poor | Yes (Z) | Poor | None |
| Umbrella sampling | Very high (8–80 GPU-days) | High (if CV is right) | Yes (Z) | Poor for flexible peptides | None |
| Metadynamics | High (0.5–5 µs) | Moderate–high | Yes (2D helps) | Moderate | None |
| **MARTINI CG** | **Low (1–2 GPU-days / 10 µs)** | Qualitative | **No** | **Good** | Substitute DPPE/DPPG |

---

## Recommended Combined Strategy

The methods are complementary. A staged approach minimises total compute cost while
extracting the most mechanistic information:

### Stage 1 — MARTINI CG pathway discovery (3 × 10 µs, ~3 GPU-days)

```bash
# 1. Equilibrate free peptide atomistically to determine dominant secondary structure
gmx mdrun -v -deffnm ku04amp01_free -nt 12   # 5–10 ns

# 2. Convert to CG with ElNeDyn elastic network
martinize2 -f ku04amp01_equil.pdb \
    -o ku04amp01_cg.top -x ku04amp01_cg.pdb \
    --ff martini3001 --elastic --ef 500 --eu 0.9 --el 0.5

# 3. Build simplified G-IM bilayer (CHARMM-GUI MARTINI Membrane Builder or insane.py)
#    Composition: DPPE×46, POPE×13, DPPE×12, OYPE×8, DPPG×10, PYPG×9, PVCL2×2

# 4. Equilibrate CG system (~50 ns NPT)
# 5. Production: 3 independent 10 µs runs
```

Analyse: does KU04AMP01 translocate? What is the dominant orientation at each depth?
Which residues contact PG/CL head groups first?

### Stage 2 — Atomistic backmapping (per key intermediate, ~1 hour each)

```bash
# Extract a representative translocation intermediate from CG trajectory
python3 backward.py -f cg_intermediate.gro -o aa_intermediate.gro \
    -from martini -to charmm36

# Short 50–100 ns CHARMM36 run from backmapped frame
gmx mdrun -v -deffnm aa_refinement -nt 48 -nb gpu -pme gpu
```

### Stage 3 — Metadynamics PMF refinement (optional, ~500 ns–1 µs)

Using the CG-identified translocation pathway as a guide, set up 2D well-tempered
metadynamics (Z distance + tilt angle) with PLUMED, starting from the backmapped
intermediate. This gives an atomistic free-energy profile along the *physically relevant*
pathway rather than along a geometrically imposed one.

```ini
# plumed.dat (excerpt)
z_dist: DISTANCE ATOMS=peptide_com,bilayer_com COMPONENTS
tilt:   ANGLE ATOMS=peptide_N_term,peptide_C_term,bilayer_normal

metad: METAD ARG=z_dist.z,tilt ...
       SIGMA=0.05,0.1 HEIGHT=1.0 PACE=500 BIASFACTOR=15
       TEMP=310.15
```

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

