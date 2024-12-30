#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=1l6j-sb3ct
#SBATCH --partition=GPU
#SBATCH --gres=mps:75
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --ntasks=32
# This allocates the number of OpenMP threads per MPI-Thread
#SBATCH --cpus-per-task=1
# This allocates the walltime to 1/2 day. The program will not run for longer.
#SBATCH --time=96:00:00
# This sets the quality of service to 'elevated'
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703
#----------------------------------------------------------

#Start time
start=`date +%s.%N`

echo "Starting"
echo '---------------------------------------------'
NUM_PROCS=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))
echo 'num_procs='$NUM_PROCS
echo '---------------------------------------------'

export CUDA_VISIBLE_DEVICES=0
export GMX_ENABLE_DIRECT_GPU_COMM=1


# Change these to something relevant
export MDNAME=md500
export MDP_FILE=md.mdp
export CPT_FILE=$MDNAME.cpt
unset OMP_NUM_THREADS

source /usr/local/gromacs/bin/GMXRC
#Actual MD Dynamics: Continue from checkpoint 
gmx mdrun -ntmpi 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -v -ntomp $NUM_PROCS -deffnm $MDNAME -cpi $CPT_FILE -nstlist 400

#End time
end=`date +%s.%N`

export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
