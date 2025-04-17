#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=colpm_pyc_prod
#SBATCH --partition=GPU
#SBATCH --gres=gpu
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
# This allocates the walltime to 4  days. The program will not run for longer.
#SBATCH --time=96:00:00
# This sets the quality of service to 'elevated'
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703
#----------------------------------------------------------
IPYNB="colpmpyocynanin.ipynb"

#Start time
start=`date +%s.%N`

echo "Starting"
echo '---------------------------------------------'
num_proc=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))
echo 'num_proc='$num_proc
echo '---------------------------------------------'

export CUDA_VISIBLE_DEVICES=0
export GMX_ENABLE_DIRECT_GPU_COMM=1

export GMX_IMGDIR=$HOME/SIFDIR/gromacs
export GMX_IMG=gromacs_2024.4-GPU.sif
SINGULARITY="singularity exec --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${GMX_IMGDIR}/${GMX_IMG}"

#Actual MD Dynamics: 
LD_LIBRARY_PATH="" ${SINGULARITY} jupyter nbconvert --to notebook --execute $IPYNB

#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
