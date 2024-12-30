# Logs for Using GROMACS

# Control: MMP9 + Sb3CT Protein-Ligand Complex in Water

* Lysozyme Tutorial: [Lysozyme](http://www.mdtutorials.com/gmx/lysozyme/index.html)
* Protein-Ligand Complex: [Tutorial](http://www.mdtutorials.com/gmx/complex/)
* Another Tutorial: [Molecular Dynamics simulation of a protein in water environment](https://www.compchems.com/gromacs-tutorial-molecular-dynamics-simulation-of-a-protein-in-water-environment/#protein-selection-and-initial-setup)
* Checkpointing: [Restart from Checkpoint or Continue MD-run](https://www.compchems.com/extend-or-continue-a-gromacs-simulation/)
* RMSF Calculation: [Obtain RMSF vs residue number using GROMACS](https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/)


Now, we are running the Full Protein-Ligand Complex 

### Index of files:
 
1. MMP9 protein PDB file: https://files.rcsb.org/view/1L6J.pdb
2. MDP file for generating ions tpr: http://www.mdtutorials.com/gmx/complex/Files/ions.mdp
3. MDP file for energy minimization: 
    * Original: http://www.mdtutorials.com/gmx/complex/Files/em.mdp
    * Modified: [em.mdp](em.mdp)
4. MDP file for NVT equilibriation: http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp
5. MDP file for NPT equilibriation: http://www.mdtutorials.com/gmx/complex/Files/npt.mdp
6. MDP file for actual production MD: http://www.mdtutorials.com/gmx/complex/Files/md.mdp

## Steps:
1. Copied the original 'processed' GRO file of protein from the protein-only simulation to a new file called [complex.gro](complex.gro)
2. Prepared a fresh topol file with charmm27 force fields.
3. Copy-pasted Ligand topology (basically the co-ordinates from lig.gro) into the complex.gro file below the last line of the protein atoms, and before the box vectors.
4. Edited the [topol.top](topol.top) file manually to include ligand topology after forcefield transclusion & updated "molecules" section with UNL
    
    **Note:** Had to transclude the ligand topology before any definitions, as a ligand that introduces new [atomtypes] must be transcluded
              after the parent force field, and prior to any [moleculetype] definition. Force field-level directives must all appear before
              any molecule-level directives. If not, genion yields error: "Invalid order for directive atomtypes"

5. Added SPCE water model force fields  
6. Added dodecahedral box using 'gmx editconf'
7. Solvated with water using 'gmx solvate'
8. (Mostly) neutralized  spurious charges with Na and Cl ions.
9. Ran Energy minimization with local mods to mdp: Manual addition of define DFLEXIBLE will use flexible water instead of rigid water. Reduces warnings
10. Before running nvt, added posre of ligand to topol.top manually AFTER posre of protein with the command

    ```bash
    ; Ligand position restraints
    #ifdef POSRES_LIG
    #include "posre_LIG.itp"
    #endif
    ```
    Now, to restrain both the protein and the ligand, we would need to specify define = -DPOSRES -DPOSRES_LIG in the .mdp file
11. Since Rolipram and the protein are physically linked very tightly, it is best to consider them as a single entity. That is, UNL is grouped with the protein for the purposes of temperature coupling. In the same way, the few ions we inserted are considered part of the solvent. To do this, we need a special index group that merges the protein and UNL. We accomplish this with:

    ```bash
    $ gmx make_ndx -f em.gro -o index.ndx
    ```
    At the prompt, merge the 'Protein' and 'UNL' groups by entering their numbers with a pipe "|" in between, like "1|13", then 'q' to quit
12. Then, download [this mdp file](http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp) and mod it to add our protein-unl group to temperature control groups, and two posres macros for restraining both protein and ligand.

    ```diff
    2c2
    < define                  = -DPOSRES  ; position restrain the protein and ligand
    ---
    > define                  = -DPOSRES -DPOSRES_LIG  ; position restrain the protein and ligand
    33c33
    < tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
    ---
    > tc-grps                 = Protein_LIG Water_and_ions    ; two coupling groups - more accurate
    ```

