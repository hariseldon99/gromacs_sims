# Logs for Using GROMACS

# Control: Colicinpk + PYC Protein-Ligand Complex in Water

* Protein-Ligand Complex: [Tutorial](http://www.mdtutorials.com/gmx/complex/)
* SwissParam ligand topology creation [http://swissparam.ch/command-line.php]

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
7. Download latest charmm field and unpack in the working directory: https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz

## Steps:
1. Generated the original 'processed' GRO file of protein from the protein-onlly pdb and copied to a new file called [complex.gro](colpk_pyc_complex.gro)
```bash
gmx pdb2gmx -f protein.pdb -o colpk_processed.gro -ignh -vsite none
```
Choose the downloaded charmm36 ff field and the corresponding charmm-compatible water model.

2. Prepared a fresh topol file with charmm27 force fields. First protonated with avogadro2, then saved as mol2, then used the Lemkul perl script to fix bond orders, then fixed residue name and atom name(s) as per the mdtutorials website, and ran

```bash
curl -F "myMol2=@pyc_fix.mol2" "swissparam.ch:5678/startparam?approach=both" 
curl "swissparam.ch:5678/retrievesession?sessionNumber=48182684" -o results.tar.gz
```
3. Overwrite old ligand pdb with pdb in results archive, since it is protonated and 'fixed'. Do not forget to get the .itp files from the topology results and change the name LIG back to PYC throughout. Convert pdb to gro using obgui

4. Copy-pasted Ligand topology (basically the co-ordinates from pyc.gro) into the [complex.gro](colpk_pyc_complex.gro) file below the last line of the protein atoms, and before the box vectors. Also, update the number of atoms in the first line by adding the number of newly added atoms.

5. Edited the [topol.top](topol.top) file manually to include ligand topology after forcefield transclusion & updated "molecules" section with PYC
    
    **Note:** Had to transclude the ligand topology before any definitions but AFTER charmm forcefield transclusion, 
	      as a ligand that introduces new [atomtypes] must be transcluded after the parent force field, and prior to any [moleculetype] definition. Force field-level 
              directives must all appear before any molecule-level directives. If not, genion yields error: "Invalid order for directive atomtypes"

6. Added SPCE water model force fields  
7. Added dodecahedral box using 'gmx editconf'
8. Solvated with water using 'gmx solvate'
9. (Mostly) neutralized  spurious charges with Na and Cl ions.
10. Ran Energy minimization:
```bash
gmx grompp -f em.mdp -c colpk_pyc_complex_box_solv_ions.gro -p topol.top -o em.tprs -maxwarn 10
gmx mdrun -nt 12 -pin on -gpu_id "0" -nb gpu -v -deffnm em
```
11. Before running nvt, [generate restraints config for lig](http://www.mdtutorials.com/gmx/complex/06_equil.html) added posre of ligand to topol.top manually AFTER  ligand itp transclusion

    ```bash
    ; Ligand position restraints
    #ifdef POSRES_PYC
    #include "posre_PYC.itp"
    #endif
    ```
    Now, to restrain both the protein and the ligand, we would need to specify define = -DPOSRES -DPOSRES_PYC in the .mdp file
12. Since pyc and the protein are physically linked very tightly, it is best to consider them as a single entity. That is, PYC is grouped with the protein for the purposes of temperature coupling. In the same way, the few ions we inserted are considered part of the solvent. To do this, we need a special index group that merges the protein and UNL. We accomplish this with:

    ```bash
    $ gmx make_ndx -f em.gro -o index.ndx
    ```
    At the prompt, merge the 'Protein' and 'PYC' groups by entering their numbers with a pipe "|" in between, like "1|13", then 'q' to quit
13. Then, download [this mdp file](http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp) and mod it to add our protein-unl group to temperature control groups, and two posres macros for restraining both protein and ligand.

    ```diff
    2c2
    < define                  = -DPOSRES  ; position restrain the protein and ligand
    ---
    > define                  = -DPOSRES -DPOSRES_PYC  ; position restrain the protein and ligand
    33c33
    < tc-grps                 = Protein_JZ4 Water_and_ions    ; two coupling groups - more accurate
    ---
    > tc-grps                 = Protein_PYC Water_and_ions    ; two coupling groups - more accurate
    ```
13. Rest of the instructions are similar to the Lemkul ones in mdtutorials.com, except set the `maxwarn` level to 10. The equilibration and production run scripts are [here](gromacs_equi.pbs) and  [here](gromacs_prod.pbs), respectively. GPU options are added as command flags in the scripts.

