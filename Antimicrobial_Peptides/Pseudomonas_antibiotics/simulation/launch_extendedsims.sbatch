#!/bin/bash
#SBATCH --job-name=mdsim_ext
#This sets the name of the job

#SBATCH --partition=GPU
#This sets the partition to the GPU partition. Important for GPU jobs

#SBATCH --gres=gpu:1
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs

#SBATCH --ntasks=1 
#This sets the number of processes to 1.

#SBATCH --cpus-per-task=24
#This allocates the number of cpus per tasks. 

#SBATCH --time=168:00:00 
#This allocates the walltime to 5 minutes. The program will not run for longer.

#SBATCH --qos=elevated 
#This sets the quality of service to 'normal'

#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472

source /usr/local/condaenv/bin/activate
source /usr/local/gromacs/bin/GMXRC
CONDA_ENV_NAME="gromacs_compiled"

export PYFILE="./resume_mdsims.py"

conda run -n $CONDA_ENV_NAME python $PYFILE 
