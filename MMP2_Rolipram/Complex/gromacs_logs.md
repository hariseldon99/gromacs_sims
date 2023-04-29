# Logs for Using GROMACS

# Control: MMP2 + Rolipram Protein-Ligand Complex in Water

* Lysozyme Tutorial: [Lysozyme](http://www.mdtutorials.com/gmx/lysozyme/index.html)
* Protein-Ligand Complex: [Tutorial](http://www.mdtutorials.com/gmx/complex/)
* Another Tutorial: [Molecular Dynamics simulation of a protein in water environment](https://www.compchems.com/gromacs-tutorial-molecular-dynamics-simulation-of-a-protein-in-water-environment/#protein-selection-and-initial-setup)
* Checkpointing: [Restart from Checkpoint or Continue MD-run](https://www.compchems.com/extend-or-continue-a-gromacs-simulation/)
* RMSF Calculation: [Obtain RMSF vs residue number using GROMACS](https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/)


Now, we are running the Full Protein-Ligand Complex after running the Protein and Ligand simulations separately

### Index of files:
 
1. MMP2 protein PDB file: https://files.rcsb.org/view/1CK7.pdb
2. MDP file for generating ions tpr: http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
3. MDP file for energy minimization: 
    * Original: http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
    * Modified: [minim.mdp](minim.mdp)
4. MDP file for NVT equilibriation: http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
5. MDP file for NPT equilibriation: http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
6. MDP file for actual production MD: http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp

## Steps:
1. Copied the original ['processed' GRO file](../MMP2_Protein/1ck7_processed.gro) of protein from the protein-only simulation to a new file called [complex.gro](complex.gro)
2. Prepared a fresh topol file with OPLS/AA force fields.
3. Copy-pasted Ligand topology (basically the co-ordinates from [rol.gro](../Rolipram/rol.gro) from the Ligand simulation) into the complex.gro file below the last line of the protein atoms, and before the box vectors.
4. Edited the [topol.top](topol.top) file manually to include ligand topology after forcefield transclusion & updated "molecules" section with UNL
    
    **Note:** Had to transclude the ligand topology before any definitions, as a ligand that introduces new [atomtypes] must be transcluded
              after the parent force field, and prior to any [moleculetype] definition. Force field-level directives must all appear before
              any molecule-level directives. If not, genion yields error: "Invalid order for directive atomtypes"

5. Added SPCE water model force fields  
4. Added dodecahedral box using 'gmx editconf'
5. Solvated with water using 'gmx solvate'
6. (Mostly) neutralized  spurious charges with Na and Cl ions.
7. Ran Energy minimization with local mods to mdp: Manual addition of define DFLEXIBLE will use flexible water instead of rigid water. Reduces warnings
8. Before running nvt, added posre of ligand to topol.top manually AFTER posre of protein with the command

    ```bash
    ; Ligand position restraints
    #ifdef POSRES_LIG
    #include "../Rolipram/posre_rol.itp"
    #endif
    ```
    Now, to restrain both the protein and the ligand, we would need to specify define = -DPOSRES -DPOSRES_LIG in the .mdp file
