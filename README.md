# Logs for GROMACS Simulations - Protein-ligand Complexes:
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/gromacs/badges/version.svg)](https://anaconda.org/conda-forge/gromacs)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gromacs_py/badges/version.svg)](https://anaconda.org/bioconda/gromacs_py)
## Test: Lysozyme in Water

1. [Simulation Logs](lysozyme/gromacs_logs.md)
2. [Results](lysozyme/post-processing/EM_pp.ipynb)


## MMP9 Protein in Water

1. [Simulation Logs](MMP9_protein/gromacs_logs.md)

## MMP Protein Simulations: 

There are two control Molecular Dynamics simulations, that of MMP2 in SPC water,
and that of Rolipram in SPC water (separately).
 
1. [MMP2 Simulation Logs](MMP2_Rolipram/MMP2_Protein/gromacs_logs.md)
2. [Rolipram Simulation Logs](MMP2_Rolipram/Rolipram/gromacs_logs.md)
3. [Results](MMP2_Rolipram/post-processing/EM_pp.ipynb)


## Zebrafish Hemoglobin and Phenyl Hydrazine

New Simulations performed using [singularity image](zeb_hb/config/gromacs_2024.2-GPU.def) with the [conda-forge version of GROMACS](https://anaconda.org/conda-forge/gromacs) (version 2024.2) compiled for nVidia-GPUs, as well as [gromacs_py](https://gromacs-py.readthedocs.io/en/latest/index.html)
 and [mdanalysis](https://www.mdanalysis.org/).
 
See [README for zeb_hb](zeb_hb/README.md)

## Docs

1. [MD-Tutorials by Justin A. Lemkul](http://www.mdtutorials.com/)- Living J. Comp. Mol. Sci. [1 (1): 5068 (2018)](https://doi.org/10.33011/livecoms.1.1.5068).
2. [GROMACS on CompChems](https://www.compchems.com/categories/gromacs/)
3. [GROMACS Website](https://www.gromacs.org/about.html)
4. [DNA-ligand simulation in gromacs-py](https://gromacs-py.readthedocs.io/en/latest/notebook/01_dna_ligand_ambertools.html)
5. [ACPYPE-AnteChamber PYthon Parser interfacE
](https://alanwilter.github.io/acpype/)
6. [RMSF Analysis of MD Simulation](https://gromacs.bioexcel.eu/t/rmsf-analysis-of-protein-ligand-md-simulation/4919/4)
7. [How to calculate RMSD of complex - researchgate](https://www.researchgate.net/post/how_to_calculated_RMSD_of_protein-ligand_complex_in_gromacs)
8. [GROMACS-Protein-Ligand](https://github.com/DweipayanG/GROMACS-Protein-Ligand/blob/main/Gromacs%20Codes)
9. [GROMACS: MD Simulation of a Protein-Ligand Complex](https://angeloraymondrossi.github.io/workshop/charmm-gromacs-small-organic-molecules-new.html)
10. [Gromacs Protein Ligand Complex Simulations](https://github.com/leelasd/ligpargen/wiki/Gromacs-Protein-Ligand-Complex-Simulations)
11. [Protein-Ligand Interactions](https://projects.volkamerlab.org/teachopencadd/talktorials/T016_protein_ligand_interactions.html)
12. [How to Compare Proteins](https://physics.aps.org/articles/v17/40) - Physics 17, 40 (2024). 
13. [CHARMM GUI](https://charmm-gui.org/)
14. [SwissParams CHARMM topology generator](http://swissparam.ch/howto_gromacs.php)
15. [Hydrogen Bond Tutorial](https://ouchidekaiseki.com/en/hydrogenbond.php)
16. [Bioinformatics Stackexchange discussion on Hydrogen adding to PDB](https://bioinformatics.stackexchange.com/questions/17916/how-can-i-programmatically-add-hydrogen-to-a-pdb-structure-using-biopython)
17. [Automated Topology Builder](https://atb.uq.edu.au/)
