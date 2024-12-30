#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=1l6j-sb3ct
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs
#SBATCH --partition=GPU
#SBATCH --gres=mps:50
#SBATCH --ntasks=16
# This allocates the number of OpenMP threads per MPI-Thread
#SBATCH --cpus-per-task=1
# This allocates the walltime to 1/2 day. The program will not run for longer.
#SBATCH --time=4:00:00
# This sets the quality of service to 'elevated'
#SBATCH --qos=elevated
#----------------------------------------------------------

#Start time
start=`date +%s.%N`

export CUDA_VISIBLE_DEVICES=0
export GMX_ENABLE_DIRECT_GPU_COMM=1

export NVT=nvt
export GROFILE=em.gro
export MDPFILE=nvt.mdp
export TOPOL_FILE=topol.top
export INDEX_FILE=index.ndx
unset OMP_NUM_THREADS

source /usr/local/gromacs/bin/GMXRC

#Preprocessing 
gmx grompp -f $MDPFILE -c $GROFILE -r $GROFILE -p $TOPOL_FILE -n $INDEX_FILE -o ${NVT}.tpr -maxwarn 2

#Actual NVT Dynamics: Note that, in ku cluster, gmx is compiled such that mpirun needs to be called explicitly
gmx mdrun -pin on -ntmpi 0 -nb gpu -pme gpu -update gpu -gpu_id 0 -ntomp ${SLURM_NTASKS} -v -deffnm $NVT

#End time
end=`date +%s.%N`

echo "NTHREADS= ${SLURM_NTASKS}", 
export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
