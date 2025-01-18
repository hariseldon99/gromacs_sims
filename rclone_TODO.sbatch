#!/bin/bash
#SBATCH --job-name=rclone-sync
#This sets the name of the job

#SBATCH --ntasks=1
#This sets the number of processes to 8. Change if needed

#SBATCH --cpus-per-task=16
#This allocates the number of cpus per tasks. If the number of tasks is 4, and cpus per task is 1, then slurm will assign 4 X 1 = 4 cpus to the job

#SBATCH --time=3:00:00 
#This allocates the walltime to 2 day. The program will not run for longer.


#Start time
start=`date +%s.%N`

let nprocs=${SLURM_NTASKS}*${SLURM_CPUS_PER_TASK}


#Start time
start=`date +%s.%N`
export SRC=$HOME/gitrepos/gromacs_sims/Colpm_pyocyanin
export DEST=gdrive://data/GROMACS/gromacs_sims/Colpm_pyocyanin


echo "Starting"
echo '---------------------------------------------'

rclone --checkers ${SLURM_NTASKS} --transfers ${SLURM_CPUS_PER_TASK} sync -P $SRC $DEST --exclude ".git/" 

#End time
end=`date +%s.%N`

export RUNTIME=$( echo "$end - $start" | bc -l )
echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
