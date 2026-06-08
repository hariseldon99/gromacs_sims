#!/bin/bash
#SBATCH --job-name=resume_zenodo
#This sets the name of the job

#SBATCH --ntasks=1 
#This sets the number of processes to 1.

#SBATCH --cpus-per-task=1
#This allocates the number of cpus per tasks. 

#SBATCH --time=11:59:00 
#This allocates the walltime to 5 minutes. The program will not run for longer.

<<<<<<< HEAD
export ZENODO_API_TOKEN="tWEcQkYyyVoDAyyBEZemxpJ4RGHRqL1avVP3kaZBxHfVA7jYmyknCcsxYDaR"  # replace with your actual token or set in environment
=======
export ZENODO_API_TOKEN=""  # replace with your actual token or set in environment
>>>>>>> 53faccd6d99d75a164144f9b22cc251006877f8c

export DEPOSITION="20525267"
for CURRENT_ZIPPATH in /tmp/pbpmez1999C.zip.*; do
	python ${HOME}/gitrepos/gromacs_sims/scripts/zenodo_upload.py -d ${DEPOSITION} -f $CURRENT_ZIPPATH
<<<<<<< HEAD
done
=======
done
>>>>>>> 53faccd6d99d75a164144f9b22cc251006877f8c
