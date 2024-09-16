# Logs for GROMACS Simulation

## Zebrafish Hemoglobin and Phenyl Hydrazine

New Simulations performed using conda-forge version of GROMACS compiled for nVidia-GPUs (version 2024.2) and [gromacs_py](https://gromacs-py.readthedocs.io/en/latest/index.html)

## Setup
Setup requires at least the following dependencies

### Minimum Requirements

1.  [acpype](https://acpype.readthedocs.io/en/latest/) to generate ligand topologies using the [AMBER force field](https://ambermd.org/AmberModels.php)
2.  [ambertools](https://ambermd.org/AmberModels.php) for the AMBER foece-field parameters
3.  [cudatoolkit>=11.8](https://anaconda.org/conda-forge/cudatoolkit) for python
4.  [gromacs=2024.2](https://anaconda.org/conda-forge/gromacs) from conda-forge
5.  [gromacs_py](https://gromacs-py.readthedocs.io/en/latest/index.html)
6.  [mdanalysis](https://www.mdanalysis.org/) for post-processing
7.  [nglview](https://github.com/nglviewer/nglview) and [pymol-3](https://github.com/schrodinger/pymol-open-source) for viewing and adjusting molecular data.
   
### Instructions:

* In order to setup a local installation, best to use [anaconda](https://www.anaconda.com/) or [miniforge](https://github.com/conda-forge/miniforge) to setup a conda-python environment using the [mdanalysis_gromacspy_condaenv.yml](config/mdanalysis_gromacspy_condaenv.yml) template. See [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for details.
* Alternatively, a singularity image can be generated from the [gromacs_2024.2-GPU.def](config/gromacs_2024.2-GPU.def) definition file and used in a reasonably portable manner. Better suited for HPC clusters. See [singularity documentation](https://docs.sylabs.io/guides/latest/user-guide/build_a_container.html#building-containers-from-singularityce-definition-files) for details
  
## Logs
The ligand PDB can be built using [rdkit](https://www.rdkit.org/) from the [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) string for the ligand. The relevant python code is in [build_ligand_pdb.py](build_ligand_pdb.py).

Note that the PDB must be fully hydrogenated for all-atom force-fields like AMBER to work. This can be done in many standard molecule viewers like pymol or chimera ([click here for instructions](https://kpwulab.com/2021/02/19/pymolchimeraadd-hydrogens-to-your-structure/)).

For further information, see the following Jupyter Notebooks:
1. [zeb-hb protein only](only_zeb_hb_Refined.ipynb) 
2. [zeb-hb protein + phz ligand complex](complex_zeb_hb_Refined_and_phz.ipynb)
