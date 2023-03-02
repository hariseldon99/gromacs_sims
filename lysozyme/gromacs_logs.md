# Logs for Using GROMACS

# Test: Lysozyme in Water

Tutorial: [Lysozyme](http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html)


## Steps:


### Initial Steps @ 20230302:1247

1. Downloaded pdb file: [1AKI](https://www.rcsb.org/structure/1AKI)
2. Stripped out all the crystal water: 
```bash
grep -v HOH 1aki.pdb > 1AKI_clean.pdb
```
3. Executing pdb2gmx to generate topology
```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```

   - Chose OPLS-AA/L all-atom force field (option 15)
   - Chose SPC water as solvent
   - Output: 
       + 1AKI_processed.gro: Post processed structure file with pot info. Can also be outputted to another PDB
       + topol.top: Topology file for the molecule + pot 
       + posre.itp: File that has info to constrain heavy atoms
4. Executing pdb2gmx again to build a PDB (just as an exercise)

```bash
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.pdb -water spce
```


   - Chose OPLS-AA/L all-atom force field (option 15)
   - Output:
       + 1AKI_processed.pdb: Post processed structure file with pot info. Diff with gro file from prev step
       + topol.top: Topology file for mol + pot. No diffs with previous top file except in creation date/time metadata
       + posre.itp: Identical to the one in the previous step. The previoius one was autobacked up as \#posre.itp\#  
       

**Notes:** 

1. GROMACS can output mol + pot structure into gro files or pdb files
2. In the [SPC (Simple Point Charge)  model of water](https://en.wikipedia.org/wiki/Water_model), the water molecules are rigid and have Lennard-Jones pots for OO interaction, but not OH bonds, presumably because they are held rigidly by *holonomic constraints* (remember those from Goldstein?). The length of OH and angle of HOH are set in the model. HOH is assumed to be at an ideal tetrahedral angle, so the water is essentially crystalline (crystal drops). Charges are specified for calculating coulombic interactions

### Solvation step