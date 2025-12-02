#!/bin/bash
#SBATCH --job-name=zenodo_upload_levparc2007C
#This sets the name of the job
#SBATCH --partition=CPU
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00 
#SBATCH --qos=elevated 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=telegram:-1001215703472

export ZENODO_TOKEN="XXXXXXXXXXX"  # replace with your actual token or set in environment
# user-supplied metadata (set these before submitting or keep defaults)
TITLE="${TITLE:-MD Simulation of DNA topoisomerase IV subunit C (parC) protein of P. aeruginosa, with mutations recorded in 2007, each in complex with Levofloxacin (LEV) ligand molecule}"
DESCRIPTION="${DESCRIPTION:-Mutated sequences of DNA topoisomerase IV (subunit C) protein (mutations recorded in 1999 and 2018) associated with P. aeruginosa are modeled and individually docked with a molecule of Nalidixic Acid. The data from subsequent Molecular Dynamics simulation is presented.\n The data from the Molecular Modelling and Docking Simulations can be found here:\n Roy, A. (2025). Docking of several proteins associated with strains of Pseudomonas aeruginosa with various antibiotic molecules [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15383904}"

# Path to the directory to be zipped and uploaded
export ZIPPATH=$HOME/workspace/levparc2007C_20251106_AMDock

# Ensure endpoint variable is defined (allows overrides)
ZENODO_ENDPOINT="${ZENODO_ENDPOINT:-https://zenodo.org}"

# Require a real token in the environment (do not hardcode)
if [ -z "${ZENODO_TOKEN:-}" ]; then
    echo "ERROR: ZENODO_TOKEN is not set. Export it before submitting the job." >&2
    exit 1
fi

# create destination directory and target zip path (use directory name of ZIPPATH)
DEST_DIR="/tmp"
mkdir -p "$DEST_DIR"
# intercept zip to automatically split the created archive into 500MB parts
# and intercept the later upload curl to upload parts individually.
PART_SIZE=500M

# Will hold the original archive path and the list of part files
ARCHIVE_TO_SPLIT=""
PART_FILES=()

zip() {
    # call the real zip (use type -P to get the filesystem binary, not the shell function)
    local real_zip
    real_zip=$(type -P zip) || { echo "zip not found" >&2; return 1; }
    "$real_zip" "$@" || return $?

    # determine archive name: first non-option argument
    local archive=""
    for a in "$@"; do
        case "$a" in
            -*) continue ;;
            *) archive="$a"; break ;;
        esac
    done
    [ -z "$archive" ] && return 0

    ARCHIVE_TO_SPLIT="$archive"

    # remove any old parts for this archive
    rm -f "${ARCHIVE_TO_SPLIT}."* 2>/dev/null || true

    # split into 500MB parts using numeric suffixes (001, 002, ...)
    # fallback if split doesn't support --numeric-suffixes -> use default suffixes
    if split --version >/dev/null 2>&1; then
        split --numeric-suffixes=1 --suffix-length=3 -b "$PART_SIZE" "$ARCHIVE_TO_SPLIT" "${ARCHIVE_TO_SPLIT}."
    else
        split -b "$PART_SIZE" "$ARCHIVE_TO_SPLIT" "${ARCHIVE_TO_SPLIT}."
    fi

    # collect parts
    PART_FILES=()
    for p in "${ARCHIVE_TO_SPLIT}."*; do
        [ -f "$p" ] || continue
        PART_FILES+=("$p")
    done

    # if splitting produced parts and there is more than one part, keep parts and (optionally) remove original
    if [ "${#PART_FILES[@]}" -gt 1 ]; then
        # remove the original full archive to avoid accidental upload of the single large file
        rm -f "$ARCHIVE_TO_SPLIT" 2>/dev/null || true
    else
        # no splitting occurred (single part), normalize PART_FILES to the single file
        if [ -f "${PART_FILES[0]:-}" ]; then
            :
        else
            # if split didn't produce a named part (older split behavior), fallback to the archive itself
            PART_FILES=("$ARCHIVE_TO_SPLIT")
        fi
    fi

    return 0
}

## Do not override curl; perform uploads explicitly below using the real curl binary.
REAL_CURL=$(command -v curl) || { echo "curl not found" >&2; exit 1; }

DEST_ZIP="$DEST_DIR/$(basename "$ZIPPATH").zip"

# remove any existing archive, ensure ZIPPATH exists, then zip recursively excluding *.trr files
rm -f "$DEST_ZIP"
cd "$ZIPPATH" || { echo "ZIPPATH not found: $ZIPPATH" >&2; exit 1; }


zip -r "$DEST_ZIP" . -x '*.trr' '*.trr/*' || { echo "zip failed" >&2; exit 1; }

# export the created archive path so later code can use it
export UPLOAD_FILE="$DEST_ZIP"

# Start timer (integer seconds only)
start=$(date +%s)

CREATOR_NAME="${CREATOR_NAME:-Roy, Analabha}"
CREATOR_AFF="${CREATOR_AFF:-Department of Physics, The University of Burdwan, India}"
CREATOR_ORCID="${CREATOR_ORCID:-0000-0002-4797-0624}"
CREATOR_ORCID_URL="https://orcid.org/${CREATOR_ORCID##*/}"
 
# Additional metadata: version and contributors
VERSION="${VERSION:-1.0}"
CONTRIB_SUPERVISOR_NAME="${CONTRIB_SUPERVISOR_NAME:-Bandopadhyay, Rajib}"
CONTRIB_SUPERVISOR_AFF="${CONTRIB_SUPERVISOR_AFF:-Department of Physics, The University of Burdwan, India}"
CONTRIB_SUPERVISOR_ORCID="${CONTRIB_SUPERVISOR_ORCID:-0000000247970624}"
# normalize to full ORCID url
CONTRIB_SUPERVISOR_ORCID_URL="https://orcid.org/${CONTRIB_SUPERVISOR_ORCID##*/}"

echo "Creating empty deposition..."
# create empty deposition (first step)
create_res="$(mktemp)"
http_code=$("$REAL_CURL" --silent --show-error --write-out "%{http_code}" \
    -H "Content-Type: application/json" -X POST --data '{}' \
    "${ZENODO_ENDPOINT}/api/deposit/depositions?access_token=${ZENODO_TOKEN}" \
    -o "$create_res" 2>&1) || true
create_body="$(cat "$create_res" 2>/dev/null || true)"
rm -f "$create_res"
if [ "${http_code:-0}" -ge 400 ] || [ -z "$create_body" ]; then
    echo "ERROR: failed to create deposition (HTTP ${http_code})" >&2
    printf '%s\n' "$create_body" >&2
    exit 1
fi
DEPOSITION=$(printf '%s' "$create_body" | jq -r .id 2>/dev/null || echo "")
if ! printf '%s' "$DEPOSITION" | grep -qE '^[0-9]+$'; then
    echo "ERROR: could not parse deposition id from create response:" >&2
    printf '%s\n' "$create_body" >&2
    exit 1
fi
echo "Deposition ID: $DEPOSITION"

# now update the deposition with metadata (include version + contributors)
update_json=$(jq -n \
--arg title "$TITLE" \
--arg description "$DESCRIPTION" \
--arg upload_type "dataset" \
--arg version "$VERSION" \
--arg name "$CREATOR_NAME" \
--arg affiliation "$CREATOR_AFF" \
--arg orcid "$CREATOR_ORCID_URL" \
--arg contrib_name "$CONTRIB_SUPERVISOR_NAME" \
--arg contrib_aff "$CONTRIB_SUPERVISOR_AFF" \
--arg contrib_orcid "$CONTRIB_SUPERVISOR_ORCID_URL" \
'{
    metadata: {
    title: $title,
    upload_type: $upload_type,
    version: $version,
    description: $description,
    creators: [
        { name: $name, affiliation: $affiliation, orcid: $orcid }
    ],
    contributors: [
        { name: $name, affiliation: $affiliation, orcid: $orcid, type: "Researcher" },
        { name: $contrib_name, affiliation: $contrib_aff, orcid: $contrib_orcid, type: "Supervisor" }
    ]
    }
}')

echo "update_json (preview):"
printf '%s\n' "$update_json" | jq . || echo "update_json not valid JSON"

echo "Updating deposition ${DEPOSITION} with metadata ..."
upd_res="$(mktemp)"
upd_code=$("$REAL_CURL" --silent --show-error --write-out "%{http_code}" \
    -H "Content-Type: application/json" -X PUT --data "$update_json" \
    "${ZENODO_ENDPOINT}/api/deposit/depositions/${DEPOSITION}?access_token=${ZENODO_TOKEN}" \
    -o "$upd_res" 2>&1) || true
upd_body="$(cat "$upd_res" 2>/dev/null || true)"
rm -f "$upd_res"
if [ "${upd_code:-0}" -ge 400 ]; then
    echo "ERROR: update deposition returned HTTP ${upd_code}" >&2
    printf '%s\n' "$upd_body" >&2
    exit 1
fi
echo "Update succeeded."

# Parse bucket URL directly from the creation output to avoid an extra API call
BUCKET=$(printf '%s' "$create_body" | jq -r '.links.bucket // empty' 2>/dev/null || echo "")
if [ -z "$BUCKET" ] || [ "$BUCKET" = "null" ]; then
    # fallback: query the API explicitly and show response for debugging
    echo "Warning: bucket not found in creation output; querying deposition endpoint..." >&2
    api_out=$(curl -sS "${ZENODO_ENDPOINT}/api/deposit/depositions/${DEPOSITION}?access_token=${ZENODO_TOKEN}" 2>&1) || {
        echo "ERROR: failed to query deposition endpoint: $api_out" >&2
        exit 1
    }
    BUCKET=$(printf '%s' "$api_out" | jq -r '.links.bucket // empty' 2>/dev/null || echo "")
    if [ -z "$BUCKET" ] || [ "$BUCKET" = "null" ]; then
        echo "ERROR: bucket URL is empty for deposition ${DEPOSITION}." >&2
        echo "API response:" >&2
        printf '%s\n' "$api_out" >&2
        exit 1
    fi
fi
echo "Bucket URL: $BUCKET"

# Upload file (handle split parts if present)
echo "Uploading file(s)..."
if [ "${#PART_FILES[@]}" -gt 0 ]; then
    for part in "${PART_FILES[@]}"; do
        echo "Uploading part: $part"
        target="${BUCKET}/$(basename "$part")?access_token=${ZENODO_TOKEN}"
        "$REAL_CURL" --progress-bar --retry 5 --retry-delay 5 -o /dev/null --upload-file "$part" "$target" || {
            echo "ERROR: upload of $part failed" >&2
            exit 1
        }
    done
else
    target="${BUCKET}/$(basename "$UPLOAD_FILE")?access_token=${ZENODO_TOKEN}"
    "$REAL_CURL" --progress-bar --retry 5 --retry-delay 5 -o /dev/null --upload-file "$UPLOAD_FILE" "$target" || {
        echo "ERROR: upload failed" >&2
        exit 1
    }
fi

# End timer
end=$(date +%s)
runtime=$((end - start))

echo "---------------------------------------------"
echo "Runtime: ${runtime} sec"
echo "---------------------------------------------"