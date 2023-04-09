#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=MMP2-resume
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
# This sets the number of processes to 24
#SBATCH --ntasks=24
# This allocates the number of cpus per tasks.
#SBATCH --cpus-per-task=1
# This allocates the walltime to 2 days. The program will not run for longer.
#SBATCH --time=48:00:00 
# This sets the quality of service to 'normal'
#SBATCH --qos=elevated 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703472
#----------------------------------------------------------

#Start time
start=`date +%s.%N`

echo starting
echo '---------------------------------------------'
num_proc=${SLURM_NTASKS}
echo 'num_proc='$num_proc
echo '---------------------------------------------'

export MDNAME=md_0_30
export CPTFILE=md_0_30.cpt

#Actual MD
export OMP_NUM_THREADS=${SLURM_NTASKS}
export MPI_NUM_PROCS=${SLURM_CPUS_PER_TASK} 


export SIFPATH=$SIFDIR/gromacs
export SIFIMG=gromacs-2022.3_20230206.sif

LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx mdrun -ntmpi ${MPI_NUM_PROCS} -pin on -ntomp ${OMP_NUM_THREADS} -v -deffnm $MDNAME -cpi ${CPTFILE}
#DO NOT USE 'srun' as it launches multiple independent jobs

#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
