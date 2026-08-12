#!/bin/bash
# Step 3 — In-vacuo ligand simulation CLI wrapper (delegates to Python API)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${SCRIPT_DIR}/03_ligand_invacuo.py" "$@"
