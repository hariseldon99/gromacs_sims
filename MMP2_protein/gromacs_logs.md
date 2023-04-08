# Logs for Using GROMACS

# Test: Lysozyme in Water

* Tutorial: [Lysozyme](http://www.mdtutorials.com/gmx/lysozyme/index.html)

* Another Tutorial: [GROMACS Tutorial: Molecular Dynamics simulation of a protein in water environment](https://www.compchems.com/gromacs-tutorial-molecular-dynamics-simulation-of-a-protein-in-water-environment/#protein-selection-and-initial-setup)


We have mostly followed the first tutorial, and the second is just for reference purposes.

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

### Solvation step @ 20230313:2313

1. First, define the container in which the protein molecule will be placed. This is a shape that will be replicated using PBC. This will be done with "gmx editconf"
```bash
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```

- The protein is centered in the box with '-c'
- It is placed a minimum distance of 1 nm from the container edge (the edge could have a complex shape) with '-d 1.0'
- The shape is set to cubic with '-bt cubic'. This means that the cube should be about 2 nm big
- The main issue is that the box should be big enough that imposing PBC does not appreciably affect the dynamics of the complex. The pots should die out or be cut-off before they hit the complex in adjacent boxes. 
- Most pot configs have builtin cutoffs that are less than a nanometer.
- Output:
    + 1AKI_newbox.gro : gromacs file with proper box for PBC

2. Solvate with water using "gmx solvate" and pass the previous config file using flag '-cp'
```bash
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```
- Option '-cs' chooses the solvent model, in this case, [a standard equilibriated 3-point solvent model](https://github.com/gromacs/gromacs/blob/main/share/top/spc216.gro). Use this with any 3-point water model. The coordinates of the solvent molecules are programmed by the choice of model.
- The previous "topol.top" file is provided with '-p' flag for modification (adding the solvent)
```bash
$ diff topol.top \#topol.top.2\# 
18404c18404
< LYSOZYME in water
---
> LYSOZYME
18409d18408
< SOL             10644
```
- Output is the gromacs file that has protein and solvate.

### Neutralize spurious charges by adding ions @ 20230314:0030

Last line of [atoms] directive in 'topol.top' reads
```json
; residue 129 LEU rtp LEU  q -1.0
  1941   opls_238    129    LEU      N    677       -0.5    14.0067
  1942   opls_241    129    LEU      H    677        0.3      1.008
  1943   opls_283    129    LEU     CA    677       0.04     12.011
  1944   opls_140    129    LEU     HA    677       0.06      1.008
  1945   opls_136    129    LEU     CB    678      -0.12     12.011
  1946   opls_140    129    LEU    HB1    678       0.06      1.008
  1947   opls_140    129    LEU    HB2    678       0.06      1.008
  1948   opls_137    129    LEU     CG    679      -0.06     12.011
  1949   opls_140    129    LEU     HG    679       0.06      1.008
  1950   opls_135    129    LEU    CD1    680      -0.18     12.011
  1951   opls_140    129    LEU   HD11    680       0.06      1.008
  1952   opls_140    129    LEU   HD12    680       0.06      1.008
  1953   opls_140    129    LEU   HD13    680       0.06      1.008
  1954   opls_135    129    LEU    CD2    681      -0.18     12.011
  1955   opls_140    129    LEU   HD21    681       0.06      1.008
  1956   opls_140    129    LEU   HD22    681       0.06      1.008
  1957   opls_140    129    LEU   HD23    681       0.06      1.008
  1958   opls_271    129    LEU      C    682        0.7     12.011
  1959   opls_272    129    LEU     O1    682       -0.8    15.9994
  1960   opls_272    129    LEU     O2    682       -0.8    15.9994   ; qtot 8

```
This is a molecule with total charge 8e? need to add ions to neutralize th(e)s(e) charge(s)

To create ion data, we need a [simple mdp file](http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp) (ordinarily used for the the steepest descent or dynamics part). The 'grompp' command will assemble the parameters specified in the .mdp file with the coordinates and topology information to generate a binary .tpr file with ion data, as well as atomic-level description of the protein, water etc.

```bash
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

...

NOTE 1 [file ions.mdp]:
  With Verlet lists the optimal nstlist is >= 10, with GPUs >= 20. Note
  that with the Verlet scheme, nstlist has no effect on the accuracy of
  your simulation.

...

NOTE 2 [file topol.top, line 18409]:
  System has non-zero total charge: 8.000000
  Total charge should normally be an integer. See
  http://www.gromacs.org/Documentation/Floating_Point_Arithmetic
  for discussion on how close it should be to an integer.
  
...

NOTE 3 [file ions.mdp]:
  You are using a plain Coulomb cut-off, which might produce artifacts.
  You might want to consider using PME electrostatics.



```

The grompp is the GROMACS pre-processor, reads a molecular topology file, checks the validity of the file, expands the topology from a molecular description to an atomic description. This goes into the tpr file.


Now, use the 'genion' command

```bash
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

...

Will try to add 0 NA ions and 8 CL ions.
Select a continuous group of solvent molecules
Group     0 (         System) has 33892 elements
Group     1 (        Protein) has  1960 elements
Group     2 (      Protein-H) has  1001 elements
Group     3 (        C-alpha) has   129 elements
Group     4 (       Backbone) has   387 elements
Group     5 (      MainChain) has   517 elements
Group     6 (   MainChain+Cb) has   634 elements
Group     7 (    MainChain+H) has   646 elements
Group     8 (      SideChain) has  1314 elements
Group     9 (    SideChain-H) has   484 elements
Group    10 (    Prot-Masses) has  1960 elements
Group    11 (    non-Protein) has 31932 elements
Group    12 (          Water) has 31932 elements
Group    13 (            SOL) has 31932 elements
Group    14 (      non-Water) has  1960 elements
Select a group: 13
Selected 13: 'SOL'
Number of (3-atomic) solvent molecules: 10644

Processing topology
Replacing 8 solute molecules in topology file (topol.top)  by 0 NA and 8 CL ions.

Back Off! I just backed up topol.top to ./#topol.top.3#
Using random seed -30511106.
Replacing solvent molecule 9742 (atom 31186) with CL
Replacing solvent molecule 5926 (atom 19738) with CL
Replacing solvent molecule 7124 (atom 23332) with CL
Replacing solvent molecule 5805 (atom 19375) with CL
Replacing solvent molecule 389 (atom 3127) with CL
Replacing solvent molecule 7492 (atom 24436) with CL
Replacing solvent molecule 1740 (atom 7180) with CL
Replacing solvent molecule 6568 (atom 21664) with CL

```

- Choose "13" (SOL)

```bash
$diff topol.top \#topol.top.3\# 
18409,18410c18409
< SOL         10636
< CL               8
---
> SOL             10644

```

- Ion positions are random

### Steepest Descent for Energy Minimization @ 20230314:1130
Generic input [MDP param file](http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp) for EM

1. Download MDP file and build input file for EM 

```bash
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```
- Input:
   - minim.mdp is generic mdp parameter file
   - 1AKI_solv_ions.gro 1aki + ions gro file
   - topol.top latest topology file
   
- Output: em.tpr is binary topology file with all the Initial Condition data for running EM
2. Did actual EM (INCL PREV STEP) using SLURM batch script

```bash
#!/bin/bash

#SBATCH --job-name=lysozyme-EM
#This sets the name of the job

#SBATCH --partition=CPU

#SBATCH --ntasks=32
#This sets the number of processes to 10.

#SBATCH --cpus-per-task=1
#This allocates the number of cpus per tasks. 

#SBATCH --time=01:00:00 
#This allocates the walltime to 1 day. The program will not run for longer.

#SBATCH --qos=normal 
#This sets the quality of service to 'normal'


#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703472

export GROFILE=1AKI_solv_ions.gro
export MDPFILE=minim.mdp
export TOPOL_FILE=topol.top

export EMFILE=em

#Preprocessing 
gmx grompp -f $MDPFILE -c $GROFILE -p $TOPOL_FILE -o ${EMFILE}.tpr

#Actual EM
gmx mdrun -ntmpi ${SLURM_NTASKS} -pin on -ntomp 1 -v -deffnm $EMFILE
# -deffnm means "Default filename", GROMACS looks for files with this name and
# chooses the right extension for the right data
```

Ran in a few seconds. Output files are

1. em.edr - Portable energy file. The energies are stored here.
2. em.log - Text log files of the whole NVT simulation.
3. em.trr - Full trajectory of the run, all Coordinates and velocities, forces and energies in binary form.
4. em.gro - Final structure at end of sim

### Actual Molecular Dynamics (MD)

####  Equilibriation

Now, the initial conditions are set from the energy minimizing position by randomly choosing velocities drawn from a Maxwell-Boltzmann distribution at T=300K. Now, all the topology data needs to be put into the .tpr binary.

The first stage is to reach thermal equilibrium in the standard canonical ensemble (NVT - isothermal-isochoric). The corresponding [mdp parameter file is here](http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp).

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```

MDP options:

1. gen_vel = yes: Initiates velocity generation.
2. tcoupl = V-rescale: The velocity rescaling thermostat.
3. pcoupl = no: Pressure coupling is not applied. Pressure evolves dynamically. Volume and temp are pre-set.

Once the tpr file is pre-processed, the MD can be run from script (the batch file above can be adapted):

```bash
gmx mdrun -deffnm nvt
```

This is for a short sim time - 100 ps (see mdp file). Output files are:

1.  nvt.cpt - Checkpoint file of coordinates and velocities and topology at the time when sim ended.
2.  nvt.edr - Energy data for the full simulation
3.  nvt.gro - The final structure
4.  nvt.trr - The full phase space trajectory of the simulation, frame by frame, stored in binary to machine precision
5.  nvt.log - Simulation log in text

Check for equilibriation by extracting temperature-time data (temp is computed from avg. KE) from edr file

```bash
gmx energy -f nvt.edr -o temperature.xvg
```

Type "16 0" at the prompt to select the temperature. Result is a 'xvg' file (text file with column data).


Technically, equilibriation should be done if temp stabilizes dynamically. However, since experimental conditions are at constant pressure rather than volume, pressure and density are what should also be stabilized. Thus, make another run from the current sim time in a NPT ensemble. Start another mdrun from the previous NVT checkpoint. As before, generate the tpr initcond file with pre-processing and launch another mdrun, this time with a [new mdp file for NPT ensemble](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp). Do this in the batch script:

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -deffnm npt
```

This is also for a short sim time - 100 ps (see mdp file). Output files are:

1.  npt.cpt - Checkpoint file of coordinates and velocities and topology at the time when sim ended.
2.  npt.edr - Energy data for the full simulation
3.  npt.gro - The final structure
4.  npt.trr - The full phase space trajectory of the simulation, frame by frame, stored in binary to machine precision
5.  npt.log - Simulation log in text

Check for equilibriation by extracting pressure-time data and density-time from edr file:

```bash
printf "18 0" | gmx_mpi energy -f npt.edr -o npt_pressure.xvg
printf "24 0" |gmx energy -f npt.edr -o density.xvg
```
Select 18 0 for pressure and 24 0 for density.

####  Production MD Run

Now that the equilibriation coordinates have been checkpointed, together with potentials, constraints etc., running an actual MD is straightforward as long as the mdp file is given. [Get it from here](http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp). First, prepare the tpr file for the MD run from the npt checkpoint and final config. This is for a 1 nanosecond simulation.


```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
```

Once the .tpr file is ready, just run the mdrun command with the tpr file name as the default filename. Convention is to set is as 'md_\<init-time\>_\<final-time\>'.

```bash
gmx mdrun -deffnm md_0_1
```

For parallel systems with GPU, we have benchmarked that OpenMP threads are faster than thread-MPI for a single node computation. No benchmarks have been done for multinode runs. The mdrun script contains:  

```bash
Load basic OHPC tools
module load ohpc
#Load cuda
module load cuda
#Load singularity
module load singularity
export SIFPATH=$HOME/images
export SIFIMG=gromacs_2022.3.sif

echo "Starting"
echo '---------------------------------------------'

export MDNAME=md_0_1
export MPI_NUM_PROCS=1
export OMP_NUM_THREADS=<whatever>

LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx mdrun -ntmpi $MPI_NUM_PROCS -nb gpu -pin on -v -ntomp $OMP_NUM_THREADS -deffnm $MDNAME
```

GROMACS for GPU is installed in a singularity container [prepared from the docker image @ the nvidia ngc catalog](https://catalog.ngc.nvidia.com/orgs/hpc/containers/gromacs).

**Note for long sims** The .trr binary trajectory files will be very large. Instead, you can put parameters in the mdp file that outputs the trajectories (co-ordinates only) in lower precision floats as a .xtc file. The params in mdp are:

```json
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
```

Similar the sims above, the output files are:

1.  md_0_1.cpt - Checkpoint file of coordinates and velocities and topology at the time when sim ended.
2.  md_0_1.edr - Energy data for the full simulation
3.  md_0_1.gro - The final structure
4.  md_0_1.xtc - The frame-by-frame coordinates for trajectory of the simulation stored in lower precision floats
5.  md_0_1.log - Simulation log in text


### Post-Processing

The first is trjconv, which is used as a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc). For this exercise, we will use trjconv to account for any periodicity in the system. The protein will diffuse through the unit cell, and may appear "broken" or may "jump" across to the other side of the box. To account for such behaviors, issue the following:

```bash
printf "1 0" |gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```

The choice of 1 0 selects the protein to be centered and the full system as output. The modified trajectories are saved to 'md_0_1_noPBC.xtc'. 

Let's look at structural stability first. GROMACS has a built-in utility for RMSD calculations called rms. To use rms, issue this command:

```bash
printf "4 4" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```
Choose 4 ("Backbone") for both the least-squares fit and the group for RMSD calculation. The -tu flag will output the results in terms of ns, even though the trajectory was written in ps. This is done for clarity of the output (especially if you have a long simulation - 1e+05 ps does not look as nice as 100 ns). The output plot will show the RMSD relative to the structure present in the minimized, equilibrated system:

If we wish to calculate RMSD relative to the crystal structure, we could issue the following:

```bash
gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
```

To obtain radius of gyration of the molecule as a function of time, do:

```bash
printf "1" | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
```
 The choice of 1 puts out data for the protein only.
 
 ### Index of files:
 
1. Lysozyme protein PDB file: https://files.rcsb.org/view/1AKI.pdb
2. MDP file for generating ions tpr: http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
3. MDP file for energy minimization: http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
4. MDP file for NVT equilibriation: http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
5. MDP file for NPT equilibriation: http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
6. MDP file for actual production MD: http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp