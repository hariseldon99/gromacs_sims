# Logs for GROMACS Simulation
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/gromacs/badges/version.svg)](https://anaconda.org/conda-forge/gromacs)

[![Static Badge](https://img.shields.io/badge/gromacspy-2.0.3-blue)](https://github.com/samuelmurail/gromacs_py)
## Zebrafish Hemoglobin and Phenyl Hydrazine

New Simulations performed using conda-forge version of GROMACS compiled for nVidia-GPUs (version 2024.2) and [gromacs_py](https://gromacs-py.readthedocs.io/en/latest/index.html)

## Setup

All simulations were carried out on Ubuntu Linux 22.04 or  CentOS Linux 7 (for HPC).

For Windows users, it is recommended that the [WSL (Windows Subsystem for Linux)](https://documentation.ubuntu.com/wsl/en/latest/guides/install-ubuntu-wsl2/) be used.

Setup requires at least the following dependencies:

### Minimum Requirements

1.  [AMDock](https://github.com/Valdes-Tresanco-MS/AMDock) for moleular docking
2.  [acpype](https://acpype.readthedocs.io/en/latest/) to generate ligand topologies using the [AMBER force field](https://ambermd.org/AmberModels.php)
3.  [ambertools](https://ambermd.org/AmberModels.php) for the AMBER foece-field parameters
4.  [cudatoolkit>=11.8](https://anaconda.org/conda-forge/cudatoolkit) for python
5.  [gromacs=2024.2](https://anaconda.org/conda-forge/gromacs) from conda-forge
6.  [gromacs_py](https://gromacs-py.readthedocs.io/en/latest/index.html)
7.  [mdanalysis](https://www.mdanalysis.org/) for post-processing
8.  [nglview](https://github.com/nglviewer/nglview) and [pymol-3](https://github.com/schrodinger/pymol-open-source) for viewing and adjusting molecular data.
   
### Instructions:

* In order to setup a local installation, best to use [anaconda](https://www.anaconda.com/) or [miniforge](https://github.com/conda-forge/miniforge) to setup a conda-python environment using the [mdanalysis_gromacspy_condaenv.yml](config/mdanalysis_gromacspy_condaenv.yml) template. See [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for details.
* Alternatively, a singularity image can be generated from the [gromacs_2024.2-GPU.def](config/gromacs_2024.2-GPU.def) definition file and used in a reasonably portable manner. Better suited for HPC clusters. See [singularity documentation](https://docs.sylabs.io/guides/latest/user-guide/build_a_container.html#building-containers-from-singularityce-definition-files) for details
  
## Logs
The ligand PDB can be built using [rdkit](https://www.rdkit.org/) from the [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) string for the ligand. The relevant python code is in [build_ligand_pdb.py](build_ligand_pdb.py).

Note that the PDB include all hydrogen atoms explicitly in order to work with all-atom force-fields like AMBER or CHARMM. 

This can be done in many standard molecule viewers like pymol or chimera ([click here for instructions](https://kpwulab.com/2021/02/19/pymolchimeraadd-hydrogens-to-your-structure/)).

Docking was performed separately using Autodock Vina via [AMDock](https://github.com/Valdes-Tresanco-MS/AMDock). Docked complex at [zeb_hb_phz_docking_20240912](zeb_hb/zeb_hb_phz_docking_20240912).

The best-docked complex (the one with the lowest bond affinity) was then used as initial configuration for the GROMACS run.

For further information, see the following Jupyter Notebooks:
1. [zeb-hb protein only](only_zeb_hb_Refined.ipynb) 
2. [zeb-hb protein + phz ligand complex](complex_zeb_hb_Refined_and_phz.ipynb)
