#!/bin/bash

#SBATCH --job-name=lysozyme-NPT
#This sets the name of the job

#SBATCH --partition=CPU

#SBATCH --ntasks=16
#This sets the number of processes to 16

#SBATCH --cpus-per-task=1
#This allocates the number of cpus per tasks. 

#SBATCH --time=12:00:00 
#This allocates the walltime to 1/2 day. The program will not run for longer.

#SBATCH --qos=normal 
#This sets the quality of service to 'normal'


#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703472

export NPT=npt
export GROFILE=nvt.gro
export MDPFILE=npt.mdp
export CPTFILE=nvt.cpt
export TOPOL_FILE=topol.top

#Preprocessing 
gmx grompp -f $MDPFILE -c $GROFILE -r $GROFILE -t $CPTFILE -p $TOPOL_FILE -o ${NPT}.tpr

#Actual EM
gmx mdrun -ntmpi ${SLURM_NTASKS} -pin on -ntomp 1 -v -deffnm $NPT
#DO NOT USE 'srun' as it launches multiple independent jobs
