
# Comparative Analysis of Two Molecular Docking Strategies  

## Report of PandaDock:

This one should have been a better option than the detailed manual workflow given below

Pandadock: https://github.com/pritampanda15/PandaDock

Explore this with test cases. Compare with published protein-ligand crystal forms (1hsg with MK1): https://www.rcsb.org/structure/1HSG

Preliminary results are bad, however. Comparison of 1htvs protein from RCSB docked with ligand MK1 (Indinavir: https://go.drugbank.com/drugs/DB00224) and comparing multiple basic settings of pandadock CPU and GPU algorithms showed poor comparison with the experimentally verified pose in the PDB. AMDock using autoligand, however, showed fairly good agreement with default settings.

One thing left TODO is to test whether changing the default scoring in PandaDock (from physics based) to ML+MMGBSA might be better.

The workflow below approaches AMDock's workflow closest.

Workflow suggestion:
---
**(EasyDock vs. HTVS Pipeline)**
This report compares the scientific merits and limitations of two docking workflows based on the provided documents:

1. **EasyDock** — a distributed high-throughput docking system built on AutoDock Vina / gnina / smina.  
2. **HTVS Pipeline** — a multi-stage GPU-accelerated workflow including pocket detection, docking, ML rescoring, MM-GBSA, and short MD equilibration.


---
## Overview of the Two Strategies

### **EasyDock (CPU-Distributed Docking Framework)**  

A Python-based platform enabling large-scale CPU docking across distributed nodes using Dask, with automated ligand handling, workflow continuation, and support for Vina-family docking programs. Its focus is *throughput and automation*.

### **HTVS Pipeline (GPU-Accelerated Multi-Stage Workflow)**  

A modular pipeline designed for scientific rigor and accuracy: automated pocket detection, docking with GPUs, ML-based rescoring, MM-GBSA, and short MD equilibration. Its focus is *accuracy and physical realism*.


#### 🧬 Manual HTVS Pipeline Details: Pocket Detection → GPU Docking → MM-GBSA Rescoring → GROMACS Equilibration

This repository provides a reproducible, SLURM-ready pipeline for high-throughput virtual screening (HTVS) using open-source tools. It automates:

- Binding pocket detection
- GPU-accelerated docking
- MM-GBSA rescoring
- Short GROMACS equilibration before production MD

> ⚠️ GROMACS production MD is assumed to be handled separately. This pipeline prepares high-confidence, physically plausible starting structures.

---

## 📌 Why This Workflow?

Docking scores are fast but approximate. They often:
- Ignore solvation, entropy, and receptor flexibility
- Misrank poses due to scoring bias
- Produce false positives that collapse in MD

This pipeline improves reliability by:
- Detecting pockets automatically (P2Rank, fpocket, rbcavity, sphgen)
- Running GPU-accelerated docking (PocketVina, Uni-Dock)
- Rescoring top hits with MM-GBSA
- Relaxing complexes in explicit solvent before MD

---

## 🧰 Tools Used

| Task | Tool | Link |
|------|------|------|
| Pocket detection | P2Rank, fpocket, DOCK6 sphgen, rDock rbcavity | [P2Rank](https://github.com/rdk/p2rank), [fpocket](https://github.com/Discngine/fpocket), [DOCK6](http://dock.compbio.ucsf.edu/), [rDock](https://github.com/rickgoss/rDock) |
| GPU docking | PocketVina, Uni-Dock, AutoDock Vina | [PocketVina](https://github.com/DeltaGroupNJUPT/PocketVina), [Uni-Dock](https://github.com/DeepPurpose/Uni-Dock), [AutoDock Vina](http://vina.scripps.edu/) |
| ML rescoring | DiffDock, Gnina | [DiffDock](https://github.com/gkipf/DiffDock), [Gnina](https://github.com/gnina/gnina) |
| MM-GBSA | gmx_MMPBSA, AMBER MMPBSA.py | [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA), [AMBER](https://ambermd.org/) |
| MD engine | GROMACS | [GROMACS](http://www.gromacs.org/) |

---

## 🧪 Pipeline Overview

```text
1. Ligand prep (SMILES → 3D → protonation)
2. Pocket detection (P2Rank/fpocket/sphgen/rbcavity)
3. Grid generation per pocket
4. GPU docking (PocketVina/Uni-Dock)
5. Rescoring (MM-GBSA or ML)
6. Short GROMACS equilibration
7. Production MD (user-defined)
```

---

## 🧮 MM-GBSA Explained

MM-GBSA estimates binding free energy:

\[
\Delta G_{\text{bind}} \approx G_{\text{complex}} - G_{\text{receptor}} - G_{\text{ligand}}
\]

Each term includes:
- Molecular mechanics energy (EMM)
- Solvation energy (GGB + GSA)
- Optional entropy (TS)

Use MM-GBSA for relative ranking, not absolute ΔG.

---

## 🖥️ SLURM Job Templates

### Pocket Detection (P2Rank)

```bash
#!/bin/bash
#SBATCH --job-name=p2rank
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

module load java
java -Xmx8G -jar p2rank.jar predict -f receptor.pdb -o p2rank_out
```

### GPU Docking (PocketVina)

```bash
#!/bin/bash
#SBATCH --job-name=vdock
#SBATCH --gres=gpu:1
#SBATCH --array=1-1000%50

LIB_CHUNK="chunks/chunk_${SLURM_ARRAY_TASK_ID}.smi"
singularity exec --nv pocketvina.sif vina_gpu \
  --receptor receptor.pdbqt --ligand ${LIB_CHUNK} \
  --center_x X --center_y Y --center_z Z \
  --size_x SX --size_y SY --size_z SZ \
  --out out/chunk_${SLURM_ARRAY_TASK_ID}.pdbqt
```

### MM-GBSA Rescoring (gmx_MMPBSA)

```bash
#!/bin/bash
#SBATCH --job-name=mmgbsa
#SBATCH --cpus-per-task=16

parallel -j 16 singularity exec mmgbsa.sif bash -c "
  prepare_complex.sh {1}
  gmx pdb2gmx -f complex.pdb -o complex.gro -p topol.top
  run_mmgbsa.sh complex
" ::: $(cat top_candidates.txt)
```

---

## 📦 Singularity Container Recipes

### Ligand Parameterization (AmberTools + OpenBabel)

```def
Bootstrap: docker
From: ubuntu:20.04

%post
  apt-get update && apt-get install -y build-essential openbabel python3-pip
  pip3 install rdkit-pypi
  wget https://ambermd.org/ambertools.tar.bz2
  tar xjf ambertools.tar.bz2 && cd ambertools && ./configure && make install
```

### MM-GBSA (gmx_MMPBSA)

```def
Bootstrap: docker
From: ubuntu:20.04

%post
  apt-get update && apt-get install -y gromacs python3-pip
  pip3 install gmx_MMPBSA numpy scipy biopython
```

---

## 🧠 Tips for Stability

- Validate ligand protonation and charges
- Treat metals and cofactors explicitly
- Use restrained MD before MM-GBSA
- Average MM-GBSA over multiple snapshots
- Record all versions and seeds for reproducibility

---

## 📚 References

- [EasyDock](https://github.com/ci-lab-cz/easydock)
- Minibaeva, G.; Ivanova, A.; Polishchuk, P., [Journal of Cheminformatics 2023, 15 (1), 102](https://doi.org/10.1186/s13321-023-00772-2).
- [P2Rank](https://github.com/rdk/p2rank)  
- [fpocket](https://github.com/Discngine/fpocket)  
- [DOCK6](http://dock.compbio.ucsf.edu/)  
- [rDock](https://github.com/rickgoss/rDock)  
- [AutoDock Vina](http://vina.scripps.edu/)  
- [Uni-Dock](https://github.com/DeepPurpose/Uni-Dock)  
- [DiffDock](https://github.com/gkipf/DiffDock)  
- [Gnina](https://github.com/gnina/gnina)  
- [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA)  
- [AMBER MMPBSA.py](https://ambermd.org/)  
- [GROMACS](http://www.gromacs.org/)  
- [Singularity](https://apptainer.org/docs/)

---

## 2. Scientific Merits

### **EasyDock Merits**
- **High scalability across CPUs** via Dask, enabling screening of millions of ligands.  
- **Automatic task resumption** with centralized SQLite database.  
- **Multi-engine support** (AutoDock Vina, gnina, smina).  
- **Boron-containing ligand support** with validated substitution protocol.  
- **Predictive runtime model** improves throughput by prioritizing slow ligands.  
- **Straightforward integration** through Python APIs.

### **HTVS Pipeline Merits**
- **End-to-end structure-based drug design workflow**: pocket detection → docking → rescoring → MD relaxation.  
- **GPU docking** (PocketVina, Uni-Dock) offers orders-of-magnitude speed gains.  
- **Higher physical accuracy** via MM-GBSA and MD equilibration.  
- **Reduces false positives** by eliminating unstable complexes that collapse in MD.  
- **Supports ML methods** (DiffDock, gnina) for improved pose prediction and rescoring.  
- **SLURM-ready scripts** for HPC environments and reproducible containerized execution.

---

## 3. Scientific Demerits

### **EasyDock Demerits**
- **Only docking-based scoring** (Vina, gnina) → limited ranking accuracy.  
- **No MM-GBSA or MD refinement**, making false positives more likely.  
- **CPU-only focus** limits throughput relative to modern GPU-based engines.  
- **No built-in pocket detection**, relies on user-specified grids.  
- **No automatic stereoisomer enumeration** unless done manually.  
- **Scaling beyond tested limits (20 nodes)** uncertain due to potential DB bottlenecks.

### **HTVS Pipeline Demerits**
- **Computationally expensive** due to MM-GBSA and MD.  
- **More complex toolchain** with multiple points of failure.  
- **Requires GPU cluster with SLURM**; less flexible than Dask-based scheduling.  
- **No centralized DB or built-in restart system**.  
- **More sensitive to ligand/force-field preparation issues**, which propagate into rescoring/MD.

---

## 4. Appropriate Scientific Use Cases

### **Use EasyDock When:**
- The goal is to screen **10⁶–10⁸ compounds** rapidly.  
- CPU cluster or ad-hoc distributed environment is available.  
- Integration into **generative design loops** or Python workflows is needed.  
- Acceptable to rely on approximate Vina-family scoring.

### **Use HTVS Pipeline When:**
- You need **high accuracy** for ranking ~10³–10⁴ shortlisted ligands.  
- GPU resources are available for docking and MD.  
- You require **physically sound complexes** for MD or ABFE workflows.  
- You want to exploit **ML-based pose prediction** and **MM-GBSA** refinement.

---

## 5. Recommended Combined Workflow
A two-tier strategy is scientifically optimal:

1. **EasyDock** for initial ultra-large docking (CPU distributed).  
2. **HTVS pipeline** for rigorous evaluation of the top ~0.1–1% of ligands.

This balances throughput and accuracy efficiently.

---

## 6. Final Comparison Table

| Criterion | EasyDock (CPU Docking) | HTVS Pipeline (Docking + Physics) |
|----------|-------------------------|-----------------------------------|
| **Primary emphasis** | Throughput, scaling, automation | Physical realism, accuracy |
| **Docking engines** | Vina, smina, gnina | PocketVina, Uni-Dock, Vina |
| **GPU support** | None | Full (GPUs required for speed) |
| **Scalability** | Excellent CPU scaling via Dask | GPU scaling strong but requires SLURM |
| **Pocket detection** | Manual grid specification | Automated (P2Rank, fpocket, sphgen, rbcavity) |
| **Post-processing** | None | MM-GBSA + MD refinement |
| **False positives** | Higher (no physics filters) | Lower (MM-GBSA + MD reduce unstable hits) |
| **Accuracy of ranking** | Moderate | Significantly higher |
| **Ease of deployment** | Simple | Complex, multi-tool |
| **Robustness** | High (single DB, auto-resume) | Moderate (toolchain-dependent) |
| **Best for** | Very large first-pass screens | Final shortlisting and hit triage |

---

## 7. Summary

EasyDock is a robust, scalable, and practical framework for **large-scale CPU-based docking** but is limited by the inherent shortcomings of docking scores. The HTVS pipeline is far more **scientifically rigorous**, combining docking, ML rescoring, MM-GBSA, and MD equilibration, at the cost of substantially higher computational requirements and workflow complexity.  

A **combined two-phase strategy** — high-throughput CPU docking followed by GPU-based physics refinement — provides the best balance of speed and accuracy for modern virtual screening.
