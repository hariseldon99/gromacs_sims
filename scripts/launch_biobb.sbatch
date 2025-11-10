#!/bin/bash
#SBATCH --job-name=biobb_pseudo
#This sets the name of the job

#SBATCH --partition=GPU
#This sets the partition to the GPU partition. Important for GPU jobs

#SBATCH --gres=mps:50
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs

#SBATCH --ntasks=1 
#This sets the number of processes to 1.

#SBATCH --cpus-per-task=32
#This allocates the number of cpus per tasks. 

#SBATCH --time=168:00:00 
#This allocates the walltime to 5 minutes. The program will not run for longer.

#SBATCH --qos=elevated 
#This sets the quality of service to 'normal'

#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472

source /usr/local/condaenv/bin/activate
source /usr/local/gromacs/bin/GMXRC
CONDA_ENV_NAME="gromacs_compiled"

export PYFILE=$HOME/gitrepos/gromacs_sims/scripts/biobb_protein_ligand_simulation.py
export PDB_FILE=$HOME/gitrepos/gromacs_sims/Antimicrobial_Peptides/Pseudomonas_antibiotics/simulation/new_simulation_complexes/pbppen2015C.pdb
export LIGAND_CODE='PEN'
export BASENAME='pbppen2015C'
echo "Starting"
echo '---------------------------------------------'


echo "NUM_CPUS= ${SLURM_CPUS_PER_TASK}"
    
export PYTHONWARNINGS="ignore::SyntaxWarning"

# Run the biobb command with SLURM_CPUS_PER_TASK
conda run -n $CONDA_ENV_NAME $PYFILE --input_structure "$PDB_FILE" \
        --ligand_code "$LIGAND_CODE" \
        --outdir "$BASENAME" \
        --nprocs $SLURM_CPUS_PER_TASK \
        --usegpu --protonated

#End time
end=`date +%s.%N`

RUNTIME=$( echo "$end - $start" | bc -l )

echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
