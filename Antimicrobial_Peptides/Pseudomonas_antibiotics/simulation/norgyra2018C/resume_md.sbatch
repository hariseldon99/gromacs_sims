#!/bin/bash
#SBATCH --job-name=biobb_pseudo
#This sets the name of the job

#SBATCH --partition=GPU
#This sets the partition to the GPU partition. Important for GPU jobs

#SBATCH --gres=gpu:1
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

echo "Starting"
echo '---------------------------------------------'

# Define the directory with PDB files

echo "NUM_CPUS= ${SLURM_CPUS_PER_TASK}"
export MD_DIR=/home/admin/gitrepos/gromacs_sims/Antimicrobial_Peptides/Pseudomonas_antibiotics/simulation/norgyra2018C/sandbox_e29b3078-14e8-4ffc-879e-c5fb85811538
export MD_NAME=prot_norgyra2018C_NOR
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

gmx -nobackup -nocopyright mdrun -o ${MD_DIR}/${MD_NAME}_md.trr -s ${MD_DIR}/${MD_NAME}_gppmd.tpr -c ${MD_DIR}/${MD_NAME}_md.gro -e ${MD_DIR}/${MD_NAME}_md.edr -g ${MD_DIR}/${MD_NAME}_md.log -x ${MD_DIR}/${MD_NAME}_md.xtc -cpi ${MD_DIR}/${MD_NAME}_md.cpt -cpo ${MD_DIR}/${MD_NAME}_md.cpt

echo "20 20"|gmx trjconv -f ${MD_DIR}/${MD_NAME}_md.xtc -s ${MD_DIR}/${MD_NAME}_gppmd.tpr -fit none -o ${MD_DIR}/${MD_NAME}_cluster_traj.xtc -n ${MD_DIR}/../${MD_NAME}_index.ndx -nocenter -pbc cluster -ur compact
echo "20 20"|gmx trjconv -f ${MD_DIR}/${MD_NAME}_cluster_traj.xtc -s ${MD_DIR}/${MD_NAME}_gppmd.tpr -fit none -o ${MD_DIR}/${MD_NAME}_cluster_center_traj.xtc -n ${MD_DIR}/../${MD_NAME}_index.ndx -center -pbc none -ur compact
echo "20 20"|gmx trjconv -f ${MD_DIR}/${MDNAME}_cluster_center_traj.xtc -s ${MD_DIR}/${MD_NAME}_gppmd.tpr -dump 0 -o ${MD_DIR}/${MD_NAME}_md_start_dry.pdb -n ${MD_DIR}/../${MD_NAME}_index.ndx
#End time
end=`date +%s.%N`

RUNTIME=$( echo "$end - $start" | bc -l )

echo '---------------------------------------------'
echo "Runtime: "$RUNTIME" sec"
echo '---------------------------------------------'
