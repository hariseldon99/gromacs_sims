#!/bin/bash
#----------------------------------------------------------
# Job name
#PBS -N singularity_job
# queue select
#PBS -q workq
#PBS -l select=1:mpiprocs=24:host=kuhpcgn2
# stdout output file
#PBS -o singularity_job.out
#PBS -j oe
#PBS -l walltime=48:00:00 
#----------------------------------------------------------

export PBS_NCPUS=$(wc -l $PBS_NODEFILE | awk '{print $1}')

#Load basic OHPC tools
module load ohpc
#Load cuda
module load cuda
#Load singularity
module load singularity
export SIFPATH=$HOME/images
export SIFIMG=gromacs_2024.2-GPU.sif

unset OMP_NUM_THREADS #PBS sets this and it creates problems with GROMACS
#Executing command
LD_LIBRARY_PATH="" singularity shell --nv -B$HOME/gitrepos/gromacs_sims:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG