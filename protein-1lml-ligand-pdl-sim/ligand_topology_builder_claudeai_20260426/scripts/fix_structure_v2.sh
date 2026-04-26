#!/bin/bash
# fix_structure_v2.sh — corrected version of fix_structure_and_rerun.sh
# Builds PDL_final_v2.pdb directly (no 00_prep_ligand.py), then runs
# antechamber + MCPB.py step 1 on the corrected 65-atom structure.
set -euo pipefail

echo "================================================================"
echo " PDL structure fix v2 — corrected 65-atom rebuild"
echo "================================================================"

# ── 0. Require the corrected PDB ────────────────────────────────────
if [[ ! -f PDL_prep_H_corrected.pdb ]]; then
    echo "[ERROR] PDL_prep_H_corrected.pdb not found."
    echo "        Download it from the Claude session outputs and place it here."
    exit 1
fi

NATOMS=$(grep -c "^HETATM" PDL_prep_H_corrected.pdb)
echo "PDL_prep_H_corrected.pdb: ${NATOMS} atoms"
[[ $NATOMS -eq 65 ]] || { echo "[ERROR] Expected 65 atoms, got ${NATOMS}"; exit 1; }

# ── 1. Backup originals ──────────────────────────────────────────────
mkdir -p backup_wrong_structure
for f in PDL_prep_H.pdb PDL_final_v2.pdb PDL.mol2 mcpb.in ligand.pdb; do
    [[ -f "$f" ]] && cp "$f" "backup_wrong_structure/$f" && echo "  backed up $f"
done

# ── 2. Install corrected PDB as PDL_prep_H.pdb ───────────────────────
cp PDL_prep_H_corrected.pdb PDL_prep_H.pdb
echo "[OK] PDL_prep_H.pdb ← PDL_prep_H_corrected.pdb (65 atoms)"

# ── 3. Build PDL_final_v2.pdb directly ──────────────────────────────
# Split into metal residue (ZN, serial 1) and ligand residue (PDL, serials 2-65).
# Do NOT call 00_prep_ligand.py — that tool is for the raw docked PDB only.
echo ""
echo "[1/3] Building PDL_final_v2.pdb..."
python3 - << 'PYEOF'
from pathlib import Path

lines_in = Path("PDL_prep_H_corrected.pdb").read_text().splitlines()
out = []
new_serial = 0

for line in lines_in:
    if not line.startswith("HETATM"):
        continue
    new_serial += 1
    elem = line[76:78].strip()

    if elem in ("Pd", "PD", "pd"):
        res_name  = "ZN "
        chain     = "A"
        res_seq   = 1
        atom_name = " ZN "
        elem_out  = "Zn"
    else:
        res_name  = "PDL"
        chain     = "A"
        res_seq   = 2
        raw       = f"{elem}{new_serial}"
        atom_name = f"{raw:<4s}"[:4]
        elem_out  = elem

    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])

    rec = (f"HETATM{new_serial:5d} {atom_name:<4s} {res_name:<3s} {chain}{res_seq:4d}"
           f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem_out:>2}\n")
    out.append(rec)

out.append("END\n")
Path("PDL_final_v2.pdb").write_text("".join(out))

heavy = sum(1 for l in out if l.startswith("HETATM") and l[76:78].strip() != "H")
h_n   = sum(1 for l in out if l.startswith("HETATM") and l[76:78].strip() == "H")
total = heavy + h_n
print(f"  Written PDL_final_v2.pdb: {total} atoms ({heavy} heavy + {h_n} H)")
assert total == 65, f"Expected 65, got {total}"
# Verify column alignment
s = [l for l in out if l.startswith("HETATM")]
assert s[0][17:20].strip() == "ZN",  f"Pd residue name wrong: {s[0][17:20]!r}"
assert s[1][17:20].strip() == "PDL", f"PDL residue name wrong: {s[1][17:20]!r}"
print("  [OK] Atom count and column alignment verified")
PYEOF

# ── 4. Extract ligand-only PDB for antechamber ───────────────────────
grep "^HETATM" PDL_final_v2.pdb | grep "PDL" > ligand.pdb
NLIG=$(wc -l < ligand.pdb)
echo "  ligand.pdb: ${NLIG} atoms (expected 64)"
[[ $NLIG -eq 64 ]] || { echo "[ERROR] Expected 64 ligand atoms, got ${NLIG}"; exit 1; }

# ── 5. Run antechamber ───────────────────────────────────────────────
echo ""
echo "[2/3] Running antechamber (GAFF2 atom typing)..."
rm -f ANTECHAMBER_AC.AC ANTECHAMBER_AC.AC0 ANTECHAMBER_BOND_TYPE.AC \
       ANTECHAMBER_BOND_TYPE.AC0 ATOMTYPE.INF sqm.in sqm.out sqm.pdb

antechamber \
    -i  ligand.pdb \
    -fi pdb \
    -o  PDL.mol2 \
    -fo mol2 \
    -at gaff2 \
    -rn PDL \
    -c  wc \
    -cf 0 \
    -nc -2 \
    -j  5 \
    -dr no \
    2>&1 | tee logs/antechamber_corrected.log

sed -i 's/MOL/PDL/g' PDL.mol2
echo "  PDL.mol2 generated"

# Check atom count in mol2
NMOL2=$(awk '/@<TRIPOS>ATOM/{found=1;next} /@<TRIPOS>/{found=0} found && NF>=6' PDL.mol2 | wc -l)
echo "  mol2 atom count: ${NMOL2} (expected 64)"

# ── 6. Patch carbonyl O types: oh → o ───────────────────────────────
echo ""
echo "[3/3] Patching carbonyl O types in PDL.mol2 (oh → o)..."
python3 - << 'PYEOF'
import re, numpy as np
from pathlib import Path

content = Path("PDL.mol2").read_text()

# Parse atoms
atoms = {}
in_atom = False
for line in content.splitlines():
    if "@<TRIPOS>ATOM" in line:  in_atom = True;  continue
    if "@<TRIPOS>"     in line:  in_atom = False; continue
    if in_atom and line.strip():
        p = line.split()
        if len(p) >= 9:
            atoms[int(p[0])] = {
                "name": p[1], "type": p[5],
                "pos": np.array([float(p[2]), float(p[3]), float(p[4])]),
                "line": line,
            }

# Parse bonds
bonds = {}
in_bond = False
for line in content.splitlines():
    if "@<TRIPOS>BOND" in line: in_bond = True;  continue
    if "@<TRIPOS>"     in line and in_bond: in_bond = False; continue
    if in_bond and line.strip():
        p = line.split()
        if len(p) >= 3:
            a, b = int(p[1]), int(p[2])
            bonds.setdefault(a, []).append(b)
            bonds.setdefault(b, []).append(a)

patched = 0
for idx, info in atoms.items():
    if info["type"] not in ("oh", "o-"):
        continue
    nbrs = bonds.get(idx, [])
    # Any H neighbours? If yes → genuine OH, leave alone
    h_nbrs = [n for n in nbrs
               if atoms[n]["type"].startswith("h") or atoms[n]["type"].startswith("H")]
    if h_nbrs:
        continue
    # C neighbour with d < 1.26 Å → carbonyl
    for c in nbrs:
        if atoms[c]["name"][0].upper() != "C":
            continue
        d = np.linalg.norm(info["pos"] - atoms[c]["pos"])
        if d < 1.26:
            old_type = info["type"]
            info["type"] = "o"
            old_line = info["line"]
            # Rebuild the line with patched type
            p = old_line.split()
            p[5] = "o"
            # Reconstruct with fixed-width fields matching mol2 format
            new_line = (f"    {p[0]:>4s} {p[1]:<8s} {p[2]:>10s} {p[3]:>10s} "
                        f"{p[4]:>10s} {p[5]:<8s} {p[6]:>3s} {p[7]:<8s} {p[8]}")
            content = content.replace(old_line, new_line, 1)
            info["line"] = new_line
            patched += 1
            print(f"  Patched {info['name']} ({old_type} → o, C-O = {d:.3f} Å)")
            break

Path("PDL.mol2").write_text(content)
print(f"  {patched} atom(s) patched. Remaining 'oh': "
      + str(content.count(" oh ")) + " (should be 0 for this ligand)")
PYEOF

# ── 7. Verify mol2 O types ───────────────────────────────────────────
echo "  O atom types in final PDL.mol2:"
awk '/@<TRIPOS>ATOM/{in_a=1;next} /@<TRIPOS>/{in_a=0} in_a && /^[ ]*[0-9]+ O/' PDL.mol2 \
  | awk '{printf "    %-8s type=%-8s\n", $2, $6}'

# ── 8. Update mcpb.in ────────────────────────────────────────────────
sed -i 's/^original_pdb.*/original_pdb         PDL_final_v2.pdb/' mcpb.in
echo ""
echo "mcpb.in — current key settings:"
grep "original_pdb\|naa_mol2_charge\|ion_mol2_charge\|total_charge\|spin_mult" mcpb.in

# ── 9. Clean stale MCPB step-1 outputs ──────────────────────────────
echo ""
echo "Removing stale MCPB.py step-1 outputs..."
rm -f PDL_small.pdb PDL_large.pdb PDL_standard.pdb PDL_small.res \
       PDL_standard.fingerprint PDL_large.fingerprint \
       PDL_small_opt.com PDL_small_fc.com PDL_large_opt.com PDL_large_mk.com

# ── 10. Run MCPB.py step 1 ──────────────────────────────────────────
echo ""
echo "Running MCPB.py step 1..."
MCPB.py -i mcpb.in -s 1 2>&1 | tee logs/mcpb_step1_corrected.log

echo ""
echo "================================================================"
echo " Generated Gaussian input files:"
ls -lh PDL_small_opt.com PDL_small_fc.com PDL_large_mk.com 2>/dev/null \
  || { echo "  [WARNING] .com files not found — check logs/mcpb_step1_corrected.log"; exit 1; }

echo ""
echo " Next steps:"
echo "  1. Restore Pd in the generated .com files:"
echo "       sed -i 's/Zn/Pd/g; s/ZN/Pd/g'   PDL_small_opt.com PDL_small_fc.com PDL_large_mk.com"
echo "       sed -i 's/^2  1/0  1/'            PDL_small_opt.com PDL_small_fc.com"
echo "  2. Verify route lines match expected (see README Stage 1)"
echo "  3. Submit: sbatch 02a_gaussian_sm_opt_restart.sh"
echo "================================================================"
