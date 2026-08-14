#!/bin/bash
# Step 2 — Run acpype CLI wrapper (delegates to Python API)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${SCRIPT_DIR}/02_run_acpype.py" "$@"
