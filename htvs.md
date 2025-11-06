
# 🧬 HTVS Pipeline: Pocket Detection → GPU Docking → MM-GBSA Rescoring → GROMACS Equilibration

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