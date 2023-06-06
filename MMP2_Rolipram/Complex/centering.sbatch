#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=gmx_centering
#SBATCH --partition=CPU
# This sets the number of processes to 24
#SBATCH --ntasks=1
# This allocates the number of cpus per tasks.
#SBATCH --cpus-per-task=1
# This allocates the walltime to 1 hour. The program will not run for longer.
#SBATCH --time=02:00:00 
# This sets the quality of service to 'normal'
#SBATCH --qos=normal 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=telegram:5545394160
#----------------------------------------------------------


#Start time
start=`date +%s.%N`

echo "Starting"
echo '---------------------------------------------'
num_proc=${SLURM_NTASKS}
echo 'num_proc='$num_proc
echo '---------------------------------------------'

LD_LIBRARY_PATH="" printf "4 0"|singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ~/SIFDIR/gromacs/gromacs-2022.3_20230206.sif gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
#End time
end=`date +%s.%N`

echo "OMP_NUM_THREADS= "$OMP_NUM_THREADS", MPI_NUM_PROCS= "$MPI_NUM_PROCS
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
