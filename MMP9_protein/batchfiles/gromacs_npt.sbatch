#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=lysozyme-NPT
#SBATCH --partition=CPU
# This sets the number of processes to 24
#SBATCH --ntasks=24
# This allocates the number of cpus per tasks.
#SBATCH --cpus-per-task=1
# This allocates the walltime to 1/2 day. The program will not run for longer.
#SBATCH --time=12:00:00 
# This sets the quality of service to 'normal'
#SBATCH --qos=normal 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:-1001215703472
#----------------------------------------------------------

#Use openmp. Set to false to switch to MPI
export USE_OPENMP=false

#Start time
start=`date +%s.%N`

echo starting
echo '---------------------------------------------'
num_proc=${SLURM_NTASKS}
echo 'num_proc='$num_proc
echo '---------------------------------------------'

export NPT=npt
# Output from NVT
export GROFILE=nvt.gro
export MDPFILE=npt.mdp
# Checkpoint from NVT 
export CPTFILE=nvt.cpt
export TOPOL_FILE=topol.top

#Preprocessing 
#gmx grompp -f $MDPFILE -c $GROFILE -r $GROFILE -t $CPTFILE -p $TOPOL_FILE -o ${NPT}.tpr

if [ "$USE_OPENMP" = true ]
then
    export OMP_NUM_THREADS=$num_proc
    export MPI_NUM_PROCS=1
else
    export OMP_NUM_THREADS=1
    export MPI_NUM_PROCS=$num_proc
fi

#Actual NPT
gmx mdrun -ntmpi ${MPI_NUM_PROCS} -pin on -ntomp ${OMP_NUM_THREADS} -v -deffnm $NPT
#DO NOT USE 'srun' as it launches multiple independent jobs

#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
