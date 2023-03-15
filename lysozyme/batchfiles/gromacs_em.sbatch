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
#DO NOT USE 'srun' as it launches multiple independent jobs
