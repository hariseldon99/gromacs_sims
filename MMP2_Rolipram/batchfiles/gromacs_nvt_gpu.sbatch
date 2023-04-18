#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=roli-nvt
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --ntasks=12
# This allocates the number of OpenMP threads per MPI-Thread
#SBATCH --cpus-per-task=2
# This allocated the number of MPI Threads
# This allocates the walltime to 1/2 day. The program will not run for longer.
#SBATCH --time=12:00:00 
# This sets the quality of service to 'normal'
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:5545394160
#----------------------------------------------------------
#Use openmp. Set to false to switch to MPI
export USE_OPENMP=true

#Start time
start=`date +%s.%N`

echo "Starting"
echo '---------------------------------------------'
num_proc=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))
echo 'num_proc='$num_proc
echo '---------------------------------------------'

export NVT=nvt

export OMP_NUM_THREADS=$SLURM_NTASKS
export MPI_NUM_PROCS=$SLURM_CPUS_PER_TASK

export GMX_IMGDIR=${SIFDIR}/gromacs/
export GMX_IMG=gromacs-2022.3_20230206.sif
export GMX_ENABLE_DIRECT_GPU_COMM=1
SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${GMX_IMGDIR}/${GMX_IMG}"

#Actual NVT Dynamics
${SINGULARITY} gmx mdrun -ntmpi $MPI_NUM_PROCS -ntomp $OMP_NUM_THREADS -nb gpu -bonded gpu -pin on -v -deffnm $NVT
#DO NOT USE 'srun' as it launches multiple independent jobs
#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
