#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=lysozyme-extend
#SBATCH --partition=GPU
#SBATCH --gres=mps:20
# This sets the number of processes to 24
#SBATCH --ntasks=24
# This allocates the number of cpus per tasks.
#SBATCH --cpus-per-task=1
# This allocates the walltime to 1/2 day. The program will not run for longer.
#SBATCH --time=182:00:00 
#SBATCH --qos=elevated
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:5545394160
#----------------------------------------------------------


#Start time
start=`date +%s.%N`

echo starting
echo '---------------------------------------------'
num_proc=${SLURM_NTASKS}
echo 'num_proc='$num_proc
echo '---------------------------------------------'

export SIFPATH=$SIFDIR/gromacs
export SIFIMG=gromacs-2022.3_20230206.sif
export MDNAME=md_0_1
export CPTFILE=${MDNAME}.cpt

export OMP_NUM_THREADS=$num_proc
export MPI_NUM_PROCS=1

#Actual MD
LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx mdrun -ntmpi $MPI_NUM_PROCS -nb gpu -pin on -v -ntomp $OMP_NUM_THREADS -deffnm $MDNAME -cpi ${CPTFILE}

#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
