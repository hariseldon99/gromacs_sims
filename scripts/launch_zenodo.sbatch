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

export ZENODO_TOKEN=<ZENODO_API_TOKEN>

export PATH="$HOME/gitrepos/zenodo_upload:$PATH"
export ZIPPATH=$HOME/gitrepos/gromacs_sims/Antimicrobial_Peptides/Pseudomonas_antibiotics/simulation/ndxgyra2007C

# create destination directory and target zip path (use directory name of ZIPPATH)
DEST_DIR="/tmp"
mkdir -p "$DEST_DIR"
DEST_ZIP="$DEST_DIR/$(basename "$ZIPPATH").zip"

# remove any existing archive, ensure ZIPPATH exists, then zip recursively excluding *.trr files
rm -f "$DEST_ZIP"
cd "$ZIPPATH" || { echo "ZIPPATH not found: $ZIPPATH" >&2; exit 1; }

zip -r "$DEST_ZIP" . -x '*.trr' '*.trr/*' || { echo "zip failed" >&2; exit 1; }

# export the created archive path so later code can use it
export UPLOAD_FILE="$DEST_ZIP"


# Start timer (integer seconds only)
start=$(date +%s)

# create a new deposition and capture output
create_out=$(zenodo_create.sh 2>&1) || { echo "zenodo_create.sh failed: $create_out" >&2; exit 1; }

# extract numeric deposition id (first number in output) and export it
DEPOSITION=$(printf '%s' "$create_out" | grep -oE '[0-9]+' | head -n1)
if [ -z "$DEPOSITION" ]; then
    echo "Failed to parse deposition id from zenodo_create.sh output:" >&2
    printf '%s\n' "$create_out" >&2
    exit 1
fi
export DEPOSITION

zenodo_upload.sh $DEPOSITION $UPLOAD_FILE -v
# End timer
end=$(date +%s)
runtime=$((end - start))

echo "---------------------------------------------"
echo "Runtime: ${runtime} sec"
echo "---------------------------------------------"