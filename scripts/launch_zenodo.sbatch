#!/bin/bash
#SBATCH --job-name=zenodo_upload_ndxgyra2007C
#This sets the name of the job

#SBATCH --partition=CPU
#This sets the partition to the GPU partition. Important for GPU jobs

#SBATCH --ntasks=1 
#This sets the number of processes to 1.

#SBATCH --cpus-per-task=1
#This allocates the number of cpus per tasks. 

#SBATCH --time=168:00:00 
#This allocates the walltime to 5 minutes. The program will not run for longer.

#SBATCH --qos=elevated 
#This sets the quality of service to 'normal'

#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472

export PYFILE="scripts/zenodo_upload.py"
export CONDA_ENV_NAME="zenodo"
export PYTHONWARNINGS="ignore::SyntaxWarning"
export ZENODO_API_TOKEN=<YOUR_ZENODO_API_TOKEN_HERE>
# Run the biobb command with SLURM_CPUS_PER_TASK
conda run -n $CONDA_ENV_NAME $PYFILE

#End time
end=`date +%s.%N`

RUNTIME=$( echo "$end - $start" | bc -l )

echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
