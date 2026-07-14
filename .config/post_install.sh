#!/bin/bash
# post_install.sh
# Run AFTER: conda env create -f gromacs2024-gpu-environment.yml && conda activate gromacs2024-gpu
#
# Installs everything that isn't a plain conda/pip package:
#   - gmx_MMPBSA                 (pip, --no-deps -- see PIP_CONSTRAINT note below)
#   - CHARMM36 (jul2022) force field, dropped into $CONDA_PREFIX/share/gromacs/top
#   - Martini 3 force field files, dropped into $CONDA_PREFIX/share/gromacs/top/martini3.ff
#   - cgenff_charmm2gmx.py       -> $CONDA_PREFIX/bin
#   - insane_compchems.py        -> $CONDA_PREFIX/bin (py2->py3 converted)
#
# Mirrors the equivalent %post steps from gromacs_2024.5-GPU.def.

set -euo pipefail

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: no active conda environment detected. Run 'conda activate gromacs2024-gpu' first." >&2
    exit 1
fi

echo ">>> Using CONDA_PREFIX=${CONDA_PREFIX}"
GMX_TOP="${CONDA_PREFIX}/share/gromacs/top"
mkdir -p "${GMX_TOP}" "${CONDA_PREFIX}/bin" "${CONDA_PREFIX}/etc"

# --- 0. Build-isolation constraint file (referenced by PIP_CONSTRAINT, see environment.yml) ---
# pandas==1.5.3 (hard-pinned by gmx_MMPBSA's setup.py) has no cp312 wheel, so pip
# builds it from sdist; that legacy setup.py imports pkg_resources, which was
# removed from setuptools as of 83.0. Constraining the isolated build env to
# setuptools<81 keeps pkg_resources available.
echo "setuptools<81" > "${CONDA_PREFIX}/etc/pip-build-constraints.txt"
export PIP_CONSTRAINT="${CONDA_PREFIX}/etc/pip-build-constraints.txt"
echo ">>> Wrote ${PIP_CONSTRAINT}"

# --- 1. gmx_MMPBSA, no-deps (runtime deps already satisfied by the conda env) ---
echo ">>> Installing gmx_MMPBSA (--no-deps)"
pip install --no-cache-dir --no-deps gmx_MMPBSA

# --- 2. CHARMM36 (jul2022) ---
echo ">>> Fetching CHARMM36 jul2022 force field"
TMPDIR_FF=$(mktemp -d)
wget -q "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz" \
    -O "${TMPDIR_FF}/charmm36.tgz"
tar -xzf "${TMPDIR_FF}/charmm36.tgz" -C "${TMPDIR_FF}"
cp -r "${TMPDIR_FF}/charmm36-jul2022.ff" "${GMX_TOP}/"
# MacKerell lab tarballs often use CRLF -- fix in-place so GROMACS can parse forcefield.doc
find "${GMX_TOP}/charmm36-jul2022.ff" -type f | xargs dos2unix -q
# Guarantee forcefield.doc exists with a one-line description (required by pdb2gmx).
# Must end with a newline -- GROMACS' line-reader skips unterminated lines.
printf 'CHARMM36m all-atom force field (July 2022)\n' \
    > "${GMX_TOP}/charmm36-jul2022.ff/forcefield.doc"

# cgenff script for small-molecule parameters
wget -q "https://raw.githubusercontent.com/Lemkul-Lab/cgenff_charmm2gmx/refs/heads/main/cgenff_charmm2gmx.py" \
    -O "${CONDA_PREFIX}/bin/cgenff_charmm2gmx.py"
chmod +x "${CONDA_PREFIX}/bin/cgenff_charmm2gmx.py"

# --- 3. Martini 3 ---
echo ">>> Fetching Martini 3 force field"
git clone --depth 1 https://github.com/marrink-lab/martini-forcefields "${TMPDIR_FF}/martini3" || true
mkdir -p "${GMX_TOP}/martini3.ff"
if [ -d "${TMPDIR_FF}/martini3/martini_forcefields" ]; then
    cp -r "${TMPDIR_FF}/martini3/martini_forcefields/regular/v3.0.0/gmx_files/"* "${GMX_TOP}/martini3.ff/"
    cp -r "${TMPDIR_FF}/martini3/martini_forcefields/regular/v3.0.0/gmx_files_contributed/"* "${GMX_TOP}/martini3.ff/"
else
    cp -r "${TMPDIR_FF}/martini3/martini_forcefields/regular/v3.0.0/gmx_files/"* "${GMX_TOP}/martini3.ff/" || true
    cp -r "${TMPDIR_FF}/martini3/martini_forcefields/regular/v3.0.0/gmx_files_contributed/"* "${GMX_TOP}/martini3.ff/" || true
fi

# --- 4. insane_compchems.py (py2 source, converted to py3 in-place) ---
echo ">>> Fetching and converting insane_compchems.py"
wget -q "https://www.compchems.com/gromacs_cg/insane.py" -O "${CONDA_PREFIX}/bin/insane_compchems.py"
2to3 -w "${CONDA_PREFIX}/bin/insane_compchems.py"
rm -f "${CONDA_PREFIX}/bin/insane_compchems.py.bak"
chmod +x "${CONDA_PREFIX}/bin/insane_compchems.py"

# --- cleanup ---
rm -rf "${TMPDIR_FF}"

echo ">>> DONE. Force fields installed under: ${GMX_TOP}"
echo ">>> Helper scripts installed under:     ${CONDA_PREFIX}/bin"
echo ">>> PIP_CONSTRAINT will auto-export on every 'conda activate gromacs2024-gpu' from now on."
