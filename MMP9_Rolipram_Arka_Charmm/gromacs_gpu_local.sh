#!/bin/bash

export num_proc=64

#Use openmp. Set to false to switch to MPI
export USE_OPENMP=true

#Start time
start=`date +%s.%N`

#Load basic OHPC tools
module load ohpc
#Load cuda
module load cuda
#Load singularity
module load singularity
export SIFPATH=$HOME/.config/sifdir
export SIFIMG=gromacs_2022.3.sif

echo "Starting"
echo '---------------------------------------------'


#Preprocessing 
LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx hbond -s md100_rescale.tpr -f md100_center.xtc -num hb.xvg -tu ns

if [ "$USE_OPENMP" = true ]
then
    export OMP_NUM_THREADS=$num_proc
    export MPI_NUM_PROCS=1
else
    export OMP_NUM_THREADS=1
    export MPI_NUM_PROCS=$num_proc
fi

#Actual MD Dynamics: 
#LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd $SIFPATH/$SIFIMG gmx mdrun -ntmpi $MPI_NUM_PROCS -nb gpu -pin on -v -ntomp $OMP_NUM_THREADS  -deffnm <ADD OPTIONS>


#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
