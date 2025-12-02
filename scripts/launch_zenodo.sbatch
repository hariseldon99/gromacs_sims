#!/bin/bash
#SBATCH --job-name=zenodo_upload_ndxgyra2007C
#This sets the name of the job

#SBATCH --partition=CPU
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00 
#SBATCH --qos=elevated 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472

# Require a real token in the environment (do not hardcode)
if [ -z "${ZENODO_TOKEN:-}" ]; then
    echo "ERROR: ZENODO_TOKEN is not set. Export it before submitting the job." >&2
    exit 1
fi

# Use the actual repo directory (hyphen) and ensure scripts are found
export PATH="$HOME/gitrepos/zenodo-upload:$PATH"
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

## Create a new empty deposition inline (no external script)
echo "Creating new deposition..."
create_out=$(curl --silent --show-error --fail \
    --retry 5 --retry-delay 5 \
    -H "Content-Type: application/json" \
    -X POST --data '{}' \
    "https://zenodo.org/api/deposit/depositions?access_token=${ZENODO_TOKEN}" 2>&1) || {
    echo "ERROR: create deposition failed: $create_out" >&2
    exit 1
}
DEPOSITION=$(printf '%s' "$create_out" | jq -r .id 2>/dev/null || echo "")
if ! printf '%s' "$DEPOSITION" | grep -qE '^[0-9]+$'; then
    echo "ERROR: failed to obtain deposition id; API output:" >&2
    printf '%s\n' "$create_out" >&2
    exit 1
fi
export DEPOSITION
echo "Deposition ID: $DEPOSITION"

# Query bucket URL and verify
ZENODO_ENDPOINT="https://zenodo.org"
BUCKET=$(curl -sS "${ZENODO_ENDPOINT}/api/deposit/depositions/${DEPOSITION}?access_token=${ZENODO_TOKEN}" | jq -r .links.bucket // empty)
if [ -z "$BUCKET" ] || [ "$BUCKET" = "null" ]; then
    echo "ERROR: bucket URL is empty for deposition ${DEPOSITION}." >&2
    echo "API response:" >&2
    curl -sS "${ZENODO_ENDPOINT}/api/deposit/depositions/${DEPOSITION}?access_token=${ZENODO_TOKEN}" >&2
    exit 1
fi
echo "Bucket URL: $BUCKET"

# Upload file
echo "Uploading file: $UPLOAD_FILE"
curl --progress-bar --retry 5 --retry-delay 5 -o /dev/null --upload-file "$UPLOAD_FILE" "${BUCKET}/$(basename "$UPLOAD_FILE")?access_token=${ZENODO_TOKEN}" || {
    echo "ERROR: upload failed" >&2
    exit 1
}

# End timer
end=$(date +%s)
runtime=$((end - start))

echo "---------------------------------------------"
echo "Runtime: ${runtime} sec"
echo "---------------------------------------------"