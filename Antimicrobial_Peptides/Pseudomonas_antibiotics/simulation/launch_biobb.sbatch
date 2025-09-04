#!/bin/bash
#SBATCH --job-name=biobb_pseudo
#This sets the name of the job

#SBATCH --partition=GPU
#This sets the partition to the GPU partition. Important for GPU jobs

#SBATCH --gres=gpu:1
#This allocates 1 GPU as a Global Resource (gres). Important for GPU jobs

#SBATCH --ntasks=1 
#This sets the number of processes to 1.

#SBATCH --cpus-per-task=24
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

echo "Starting"
echo '---------------------------------------------'

# Define the directory with PDB files
PDB_DIR="$HOME/gitrepos/gromacs_sims/Antimicrobial_Peptides/Pseudomonas_antibiotics/simulation/new_simulation_complexes"

# Iterate over all PDB files
for pdb_file in $PDB_DIR/*.pdb; do
    # Extract just the filename for output directory
    filename=$(basename "$pdb_file")
    basename="${filename%.pdb}"
    
    # Try to extract ligand code from PDB file
    # This gets the first 3-letter code after HET records
    ligand_code=$(grep -m 1 "HETATM\|HET" "$pdb_file" | awk '{print $4}')
    
    # If ligand extraction fails, use default LIG
    if [[ -z "$ligand_code" || ${#ligand_code} != 3 ]]; then
        echo "Warning: Could not extract valid ligand code from $pdb_file. Using default LIG."
        ligand_code="LIG"
    fi
    
    echo "Processing $pdb_file with ligand code $ligand_code"
    
    # Run the biobb command with SLURM_CPUS_PER_TASK
    conda run -n $CONDA_ENV_NAME $PYFILE --input_structure "$pdb_file" \
        --ligand_code "$ligand_code" \
        --outdir "$basename" \
        --nprocs $SLURM_CPUS_PER_TASK \
        --usegpu --protonated
done

#End time
end=`date +%s.%N`

echo "NUM_CPUS= ${SLURM_CPUS_PER_TASK}"
RUNTIME=$( echo "$end - $start" | bc -l )

echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'