# Logs for Using GROMACS

# Control: MMP2 in Water

* Lysozyme Tutorial: [Lysozyme](http://www.mdtutorials.com/gmx/lysozyme/index.html)

* Another Tutorial: [Molecular Dynamics simulation of a protein in water environment](https://www.compchems.com/gromacs-tutorial-molecular-dynamics-simulation-of-a-protein-in-water-environment/#protein-selection-and-initial-setup)
* Checkpointing: [Restart from Checkpoint or Continue MD-run](https://www.compchems.com/extend-or-continue-a-gromacs-simulation/)
* RMSF Calculation: [Obtain RMSF vs residue number using GROMACS](https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/)


After our [successful run for Lyzosyme](../lysozyme/gromacs_logs.md), we ran a simulation of the MMP9 protein using the same prescriptions and steps as we did for lyzozyme.

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
1. Downloaded pdb file
2. Stripped out all the crystal water
3. Added SPCE water model force fields (OPLS-AA/L all-atom). 
4. Added cubic box using 'gmx editconf'
5. Solvated with water using 'gmx solvate'
6. Neutralized spurious charges with Na and Cl ions.
7. Ran steepest descent minimization of energy. Our first deviation from the lysozyme sim occurred here.
8. Ran NVT and NPT equilibriation
9. Production MD from NPT equilibrium. Had to extend the simulation way beyond the 1 ns time that was set for Lyzosyme.
10. Post processing: Obtained RMSD and Radius of Gyration. New quantity (RMSF) as a function of residue number was also obtained.

## Deviations from Lysozyme:

#### Energy Minimization

During energy minimization, GROMACS spat out a warning with the output for the following command:

```bash
$ gmx grompp -f minim.mdp -c 1ck7_solv_ions.gro -p topol.top -o em.tpr

...

WARNING 1 [file minim.mdp]:
  The largest distance between excluded atoms is 1.694 nm, which is larger
  than the cut-off distance. This will lead to missing long-range
  corrections in the forces and energies. If you expect that minimization
  will bring such distances within the cut-off, you can ignore this warning.
```

This is what ChatGPT said after asking it about this issue:

This error message is indicating that the distance between two atoms that are excluded from the nonbonded interactions (such as van der Waals and electrostatic interactions) is greater than the cut-off distance specified in your simulation parameters.

In GROMACS, the nonbonded interactions are typically truncated beyond a certain distance, which is specified by the user. This is done to reduce the computational cost of the simulation. However, if the distance between excluded atoms is greater than the specified cut-off distance, the simulation will miss the long-range corrections in the forces and energies, which can affect the accuracy of the simulation.

In the context of energy minimization, the warning is suggesting that the minimization may not converge to the true energy minimum of the system because the excluded atoms are not interacting correctly. If you expect that the minimization will bring the excluded atoms within the cut-off distance, you can ignore the warning. However, if you are concerned about the accuracy of your simulation, you may want to consider increasing the cut-off distance or modifying the exclusion settings to ensure that all relevant atoms are included in the nonbonded interactions.

To increase the cut-off distance in GROMACS, you need to modify the value of the rcoulomb and rvdw parameters in your simulation input file. These parameters specify the cut-off distance for the Coulomb (electrostatic) and van der Waals interactions, respectively.

You can modify the cut-off distance by adding the following lines to your .mdp file:

```makefile
rlist = <desired value> 
rcoulomb = <desired value> 
rvdw = <desired value>
```

Here, rlist specifies the distance at which the neighbor list is constructed, and rcoulomb and rvdw specify the cut-off distance for Coulomb and van der Waals interactions, respectively. You should set the value of rlist to be slightly larger than the value of rcoulomb and rvdw.

Note that increasing the cut-off distance will increase the computational cost of the simulation, as more atoms will need to be included in the nonbonded interactions. Therefore, you should choose a reasonable value for the cut-off distance based on the accuracy you need and the available computational resources.

I tried to do this to the mdp file, but noticed that there is an 'nstlist' parameter (Frequency to update the neighbor list and long range forces) but no 'rlist'.

The nstlist parameter in the .mdp file specifies the frequency with which the neighbor list is updated, but it does not directly control the cut-off distance for nonbonded interactions.

To specify the cut-off distance for nonbonded interactions in GROMACS, you need to set the rvdw and rcoulomb parameters. These parameters define the cut-off distance for van der Waals and electrostatic interactions, respectively.

However, the neighbor list distance, which is defined by the rlist parameter, should be larger than the cut-off distance to ensure that all interacting pairs of atoms are included in the neighbor list. If you do not specify a value for rlist, GROMACS will use a default value of rvdw + 0.1 nm.

Therefore, if you want to increase the cut-off distance for nonbonded interactions in your simulation, you should set the rcoulomb and rvdw parameters to your desired values and ensure that rlist is larger than both of these values.

So I modified the 'minim.mdp' file as follows:

```bash
$ diff minim.mdp minim.mdp.works
9c9
< nstlist         = 1         ; Frequency to update the neighbor list and long range forces
---
> nstlist         = 20        ; Frequency to update the neighbor list and long range forces: Updated manually after gromacs spat out a note
13,14c13,15
< rcoulomb        = 1.0       ; Short-range electrostatic cut-off
< rvdw            = 1.0       ; Short-range Van der Waals cut-off
---
> rcoulomb        = 2.0       ; Short-range electrostatic cut-off
> rvdw            = 2.0       ; Short-range Van der Waals cut-off
> rlist		= 2.5       ; ChatGPT says that I shuld add this for this particular protein
```


Then, the em ran with no serious warnings.

####  Production MD Run

The production MD went fine. We ran it from t=0 to t=30 ns. However, the RMSD's were rather large, so we decided to extend the simulation to 30 ns, then to 100 ns. In order to extend the simulation past the last checkpoint, we followed [this tutorial](https://www.compchems.com/extend-or-continue-a-gromacs-simulation/).

By default, GROMACS writes a checkpoint file of your system in the cpt format (md.cpt) every 15 minutes. It also automatically backs up the previous checkpoint as md_prev.cpt.

The checkpoint files contain a complete description of the system storing the coordinates and velocities of the system at full precision. Therefore, you can always continue the simulation from the last checkpoint exactly as if the interruption never occurred. To do that, you only need to pass it to the gmx mdrun with the -cpi flag.

In some cases, you may run a simulation for a certain period of time but once it finishes you decide that you want to proceed further. To accomplish this you can use the built-in gmx convert-tpr module.

This module allows you to edit the tpr files you assembled via the gmx grompp module in various ways. Among others, we can decide to extend the duration of the simulation for a certain tpr file.

This procedure involved two steps:

1. First, we used the 'convert-tpr' module to extend the simulation time without creating a new tpr file. In other words, the name of the tpr file that we created via the -o flag needed to be the same as the one we already used.

    ```bash
    $ gmx convert-tpr -s md_0_30.tpr -extend 100000 -o md_0_30.tpr
    ```
    To append to the previous files we will need to have all the files generated by the previous simulation with the same prefix in the current directory (md_0_30.xtc, md_0_30.edr, md_0_30.log, md_0_30.trr, md_0_30.cpt).

    In the next step you will only need to run the new tpr file starting from the last checkpoint (as explained above).

2. Finally, we ran the simulation starting from the last cpt file:
    ```bash
    $ gmx mdrun -v -deffnm md_0_30 -cpi md_0_30.cpt
    ```
    Thus, even though the filenames say 'md_0_30', the data that they contain are from a t=0-100 ns simulation. This took several hours (more than a day).

#### Post-Processing
 During NVT equilibriation, the temperature-time data was obtained the same way as lysozyme. In addition, the 'gmx energy' command was used to obtain Total energy E(t). Although the actual equilibriation ***should*** be seen with the Helmholtz free energy  H(t) = E(t) - T S(t) where 'T' is the equilibriated (time-average) temperature and S(t) is the entropy, getting the entropy involces using [Schlitter’s Formula on the covariance matrix](https://doi.org/10.1021/jp046022f), and large cov matrices have to be diagonalized at each time frame. This is a slow memory intensive computation (see [This researchgate discussion](https://www.researchgate.net/post/How_can_I_calculate_configurational_entropy_along_the_complete_trajectory_of_gromacs)) that does not yield anything terribly interesting, as entropy equilibriates pretty fast. So equilibriation was seen with merely the total energy.

During NPT equilibriation, the free energy is the enthalpy, easily computed since H = E + PV and pressure can already be obtained from 'gmx energy' and, in any case, 'gmx energy' has an option for enthalpy.

For the full MD post equilibriation, In addition to plotting the RMSD and Radius of Gyration, we also obtained RMSF plots [based on this tutorial](https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/). This is different from RMSD in two ways.
 
 1. The RMSD is either the instantaneous RMS deviation of all the chosen atoms (usually just the backbone) from a reference configuration, either the initial configuration (after equilibriating) or from the crystal structure (after just EM). 
 2. The RMSF is the RMS of the **temporal** RMS deviation of a group of atoms from the **temporal** average.
 

We elected to plot the RMSF over all simulation time of each C-alpha residue as a function of the residue number.

```bash
$ gmx rmsf -f md_0_30.xtc -s md_0_30.tpr -o rmsf.xvg -res
```
Description of flags:
   * '-s' chooses the reference configuration for the RMSF
   * '-f' chooses the trajectory file. In this case, xtc was prepared during MD
   * '-res' is used to compute the RMSF for each residue. 
