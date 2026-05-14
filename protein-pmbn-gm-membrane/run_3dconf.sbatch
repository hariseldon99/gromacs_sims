#!/bin/bash
#----------------------------------------------------------
# This sets the name of the job
#SBATCH --job-name=run3d-conf
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
# This allocates the walltime to 4 days. The program will not run for longer.
#SBATCH --time=96:00:00
#SBATCH --output=run3d-conf.out
#SBATCH --error=run3d-conf.out
#SBATCH --qos=elevated 
#This sets the quality of service to 'normal'

#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472
#----------------------------------------------------------

RUNSCRIPT="../scripts/generate_3d_conformer.py"
CONDA_ENV="htvs"
INPUT_SDF="pmbn_pubchem.sdf"
OUTPUT_PREFIX="polymyxin_b_3d"
TEMPLATE_PDB="polymyxin_b_pubchem_residues.pdb"
NUM_CONFORMERS=200
NUM_THREADS=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))

#Start time
start=`date +%s.%N`

echo "Starting"
echo '---------------------------------------------'
source /usr/local/condaenv/bin/activate
conda run -n ${CONDA_ENV} python ${RUNSCRIPT} \
		--input ${INPUT_SDF} \
		--output-prefix ${OUTPUT_PREFIX} \
		--template-pdb ${TEMPLATE_PDB} \
		--keep-existing-hydrogens \
		--num-conformers ${NUM_CONFORMERS} \
		--num-threads ${NUM_THREADS}

#End time
end=`date +%s.%N`

echo "NUM_THREADS= ${NUM_THREADS}"
RUNTIME=$( echo "$end - $start" | bc -l )

echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'

