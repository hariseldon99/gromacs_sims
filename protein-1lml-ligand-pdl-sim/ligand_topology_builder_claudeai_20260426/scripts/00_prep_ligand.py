#!/usr/bin/env python3
"""
00_prep_ligand.py
=================
Prepare best_docked_ligandonly.pdb for AMBER topology build.

Actions:
  1. Strip AutoDock virtual/dummy BH atoms (d < 1.0 Å to a real atom).
  2. Fix the element column (cols 77-78) — BH → B, enforce title-case.
  3. Assign unique 4-char atom names compatible with AMBER / antechamber.
  4. Rebuild CONECT records based on covalent radii.
  5. Report coordination environment of Pd and write PDL_prep.pdb.

Usage:
    python3 00_prep_ligand.py [input.pdb]
    Default input: best_docked_ligandonly.pdb
"""

import sys, os, math, itertools
from collections import defaultdict

# --------------------------------------------------------------------------
# Covalent radii (Å) for bond detection
# --------------------------------------------------------------------------
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Pd": 1.39, "Zn": 1.22, "Cu": 1.32, "Fe": 1.52,
}
BOND_TOLERANCE = 0.40   # Å extra on sum of radii

def cov_radius(elem):
    return COVALENT_RADII.get(elem, 0.77)

def is_bonded(a, b):
    d = distance(a["xyz"], b["xyz"])
    r = cov_radius(a["elem"]) + cov_radius(b["elem"]) + BOND_TOLERANCE
    return d < r

def distance(p, q):
    return math.sqrt(sum((pi-qi)**2 for pi,qi in zip(p,q)))

# --------------------------------------------------------------------------
# PDB parser (HETATM / ATOM only)
# --------------------------------------------------------------------------
def parse_pdb(path):
    """
    Robust PDB parser that handles AMDock/AutoDock column quirks:
      - 2-char element names (Pd) can start at position 11 or 12
      - AutoDock replaces occupancy/B-factor with partial charges
      - Element field may sit at [75:77] or [76:78]
    """
    atoms = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            rec = line[:6].strip()
            if rec not in ("HETATM", "ATOM"):
                continue
            serial = int(line[6:11])

            # --- atom name: try standard [12:16] first, then [11:15] ---
            name_std   = line[12:16].strip() if len(line) > 15 else ""
            name_shift = line[11:15].strip() if len(line) > 14 else ""
            # Prefer whichever starts with an upper-case letter (real element)
            if name_std and name_std[0].isupper():
                name = name_std
            elif name_shift and name_shift[0].isupper():
                name = name_shift
            else:
                name = name_std or name_shift

            resname = line[17:20].strip() if len(line) > 19 else "UNK"
            chain   = line[21].strip()    if len(line) > 21 else "A"
            resseq  = int(line[22:26])    if len(line) > 25 else 1
            x       = float(line[30:38])
            y       = float(line[38:46])
            z       = float(line[46:54])

            # --- element field ---
            # AMDock writes element at [75:77] (line length 77); standard PDB is [76:78].
            # Strategy: try both windows and pick the one that yields a known symbol.
            KNOWN_2CHAR = {"Pd","Fe","Cu","Zn","Mg","Ca","Na","Cl","Br","Al",
                           "Si","Mn","Co","Ni","Mo","Ru","Rh","Os","Ir","Pt",
                           "Au","Ag","Cd","Hg","Pb","Cr","Ti","V","W","Nb"}
            elem_raw = ""
            for start in (75, 76, 74, 77):
                candidate = line[start:start+2].strip() if len(line) >= start+1 else ""
                if not candidate or not candidate.replace("-","").isalpha():
                    continue
                # Prefer 2-char symbols explicitly listed, then any alpha
                cand_title = candidate[0].upper() + (candidate[1].lower() if len(candidate)>1 else "")
                if len(cand_title) == 2 and cand_title in KNOWN_2CHAR:
                    elem_raw = cand_title
                    break
                if not elem_raw:  # keep as fallback
                    elem_raw = candidate
            if not elem_raw:
                elem_raw = ""

            atoms.append({
                "serial":   serial,
                "name":     name,
                "resname":  resname,
                "chain":    chain,
                "resseq":   resseq,
                "xyz":      (x, y, z),
                "elem_raw": elem_raw,
            })
    return atoms

# --------------------------------------------------------------------------
# Element cleaning
# --------------------------------------------------------------------------
def clean_element(elem_raw, name):
    """Return a valid 1-2 char element symbol."""
    # Known problematic cases from AutoDock/AMDock
    mapping = {"BH": "B", "Pd": "Pd", "ZN": "Zn", "MG": "Mg",
               "FE": "Fe", "CU": "Cu", "CA": "Ca", "NA": "Na"}
    if elem_raw in mapping:
        return mapping[elem_raw]
    if elem_raw and len(elem_raw) <= 2 and elem_raw.isalpha():
        return elem_raw.capitalize() if len(elem_raw)==1 else elem_raw[0].upper()+elem_raw[1].lower()
    # Guess from atom name
    for e in ["Pd","Cl","Br","Fe","Cu","Zn","Mg","Ca","Na"]:
        if name.upper().startswith(e.upper()):
            return e
    return name[0].upper()

# --------------------------------------------------------------------------
# Detect and strip virtual/dummy atoms
# --------------------------------------------------------------------------
def strip_virtual_atoms(atoms, threshold=1.05):
    """
    Strip AutoDock virtual/dummy atoms.

    Criteria:
      - Distance < threshold to ANY atom with a DIFFERENT element, AND
      - The candidate to strip is a non-standard element (BH, DU, ...) OR
        is an H that is implausibly close (< 0.95 Å) to a non-H heavy atom.

    Real N-H bonds are ~1.01 Å; AutoDock places all dummy H/BH atoms at
    exactly ~0.86 Å from their parent heavy atom, which is below the N-H
    minimum of ~0.95 Å. We use 0.96 Å as the cutoff for stripping H.
    """
def strip_virtual_atoms(atoms, threshold=1.05):
    """
    Strip AutoDock/AMDock virtual atoms by examining raw element strings and
    inter-element distances.

    An atom is virtual if it satisfies ANY of:
      (a) elem_raw is in NON_STANDARD  (BH, DU, EP, LP, PS)
      (b) elem_raw == "H" and d < 0.95 Å to a different-element neighbour
          (AutoDock NH2 pseudo-H placed at 0.86 Å; real X-H bonds ≥ 0.95 Å)
    """
    NON_STANDARD = {"BH", "DU", "EP", "LP", "PS", "BHDOCK"}
    H_CUTOFF     = 0.95   # Å — real X-H bonds are ≥ ~0.95 Å

    stripped   = []
    keep_flags = [True] * len(atoms)

    for i, a in enumerate(atoms):
        if not keep_flags[i]:
            continue
        a_raw = (a.get("elem_raw") or a.get("elem") or "").strip().upper()
        for j, b in enumerate(atoms):
            if i == j or not keep_flags[j]:
                continue
            b_raw = (b.get("elem_raw") or b.get("elem") or "").strip().upper()
            if a_raw == b_raw:
                continue
            d = distance(a["xyz"], b["xyz"])
            if d >= threshold:
                continue

            # Determine which is the virtual atom
            a_is_virtual = (a_raw in NON_STANDARD or
                            (a_raw == "H" and d < H_CUTOFF))
            b_is_virtual = (b_raw in NON_STANDARD or
                            (b_raw == "H" and d < H_CUTOFF))

            if a_is_virtual and keep_flags[i]:
                keep_flags[i] = False
                stripped.append((a, b, d))
            elif b_is_virtual and keep_flags[j]:
                keep_flags[j] = False
                stripped.append((b, a, d))
            # If neither is unambiguously virtual (shouldn't happen with real
            # PDB data), leave both and let the user decide.

    kept = [a for a, f in zip(atoms, keep_flags) if f]
    return kept, stripped

    kept = [a for a, f in zip(atoms, keep_flags) if f]
    return kept, stripped

# --------------------------------------------------------------------------
# Unique atom name assignment
# --------------------------------------------------------------------------
def assign_unique_names(atoms):
    """
    Assign unique 4-char PDB atom names: element + index within element.
    E.g., C1, C2, ..., N1, N2, ..., Pd1
    """
    counters = defaultdict(int)
    for a in atoms:
        e = a["elem"]
        counters[e] += 1
        n = counters[e]
        # Build name: elem (1-2 chars) + number — total ≤ 4 chars
        raw = f"{e}{n}"
        if len(raw) > 4:
            raw = raw[:4]
        a["name_new"] = f"{raw:<4s}"
    return atoms

# --------------------------------------------------------------------------
# CONECT record builder
# --------------------------------------------------------------------------
def build_conect(atoms):
    bonds = defaultdict(list)
    for i, a in enumerate(atoms):
        for j, b in enumerate(atoms):
            if j <= i:
                continue
            if is_bonded(a, b):
                bonds[a["serial_new"]].append(b["serial_new"])
                bonds[b["serial_new"]].append(a["serial_new"])
    return bonds

# --------------------------------------------------------------------------
# PDB writer
# --------------------------------------------------------------------------
def write_pdb(atoms, bonds, outpath):
    with open(outpath, "w") as fh:
        fh.write("REMARK Prepared by 00_prep_ligand.py for MCPB.py / AMBER\n")
        fh.write("REMARK Residue: PDL  |  Metal: Pd(II)  |  Coordination: square-planar N4\n")
        for a in atoms:
            x, y, z = a["xyz"]
            e = a["elem"]
            elem_col = f"{e:>2s}"
            line = (
                f"HETATM{a['serial_new']:5d} {a['name_new']:4s} "
                f"{a['resname']:<3s} {a['chain']:1s}"
                f"{a['resseq']:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          {elem_col}\n"
            )
            fh.write(line)
        fh.write("END\n")
        # CONECT records
        for serial, neighbors in sorted(bonds.items()):
            # Split into lines of 4 neighbors max
            for i in range(0, len(neighbors), 4):
                chunk = neighbors[i:i+4]
                conn_str = "".join(f"{n:5d}" for n in chunk)
                fh.write(f"CONECT{serial:5d}{conn_str}\n")
    print(f"  → Written: {outpath}")

# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def report_coordination(atoms, metal_elem="Pd", cutoff=2.8):
    metals = [a for a in atoms if a["elem"] == metal_elem]
    if not metals:
        print(f"  [!] No {metal_elem} atom found.")
        return
    for m in metals:
        print(f"\n  Metal: {m['name_new'].strip()} (serial {m['serial_new']})  "
              f"xyz = ({m['xyz'][0]:.3f}, {m['xyz'][1]:.3f}, {m['xyz'][2]:.3f})")
        coord = []
        for a in atoms:
            if a is m: continue
            d = distance(a["xyz"], m["xyz"])
            if d < cutoff:
                coord.append((d, a))
        coord.sort()
        print(f"  Coordination environment (d < {cutoff} Å):")
        for d, a in coord:
            print(f"    {a['name_new'].strip():6s} ({a['elem']:2s})  serial={a['serial_new']:3d}  d={d:.3f} Å")
        print(f"\n  Ion serial for mcpb.in:  ion_ids {metals[0]['serial_new']}")

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else "best_docked_ligandonly.pdb"
    outfile = "PDL_prep.pdb"

    print(f"\n{'='*60}")
    print(f"  Ligand PDB Prep: {infile}")
    print(f"{'='*60}\n")

    atoms = parse_pdb(infile)
    print(f"  Parsed {len(atoms)} atoms from {infile}")

    # ── Step 1: strip virtual atoms using RAW element strings ──────────────
    # Must happen before clean_element because BH→B conversion would make
    # the virtual atom indistinguishable from a real boron.
    # We use elem_raw for stripping decisions; real atoms are not BH/DU/EP/LP
    # and do not have inter-element distances below 0.96 Å in reality.
    for a in atoms:
        a["elem"] = a["elem_raw"]   # temporary: use raw for stripping

    atoms, stripped = strip_virtual_atoms(atoms, threshold=1.05)
    if stripped:
        print(f"\n  Stripped {len(stripped)} virtual/dummy atoms:")
        for v, neighbor, d in stripped:
            print(f"    Removed: serial={v['serial']} name={v['name']} "
                  f"elem_raw='{v['elem_raw']}' — d={d:.3f} Å to "
                  f"{neighbor['name']} ({neighbor['elem_raw']})")
    else:
        print("  No virtual atoms detected.")
    print(f"  Remaining atoms: {len(atoms)}")

    # ── Step 2: assign cleaned elements now that dummies are gone ──────────
    for a in atoms:
        a["elem"] = clean_element(a["elem_raw"], a["name"])
    changed = [(a['serial'], a['name'], a['elem_raw'], a['elem'])
               for a in atoms if a['elem_raw'] != a['elem']]
    if changed:
        print(f"\n  Element cleaning:")
        for s, n, raw, clean in changed:
            print(f"    Atom {s:3d} {n:6s}: '{raw}' → '{clean}'")

    # Renumber serials
    for idx, a in enumerate(atoms, start=1):
        a["serial_new"] = idx
    for a in atoms:
        a["resname"] = "PDL"
        a["chain"]   = "A"
        a["resseq"]  = 1

    # Unique names
    atoms = assign_unique_names(atoms)

    # Coordination report
    report_coordination(atoms, metal_elem="Pd", cutoff=2.8)

    # CONECT
    bonds = build_conect(atoms)
    n_bonds = sum(len(v) for v in bonds.values()) // 2
    print(f"\n  Detected {n_bonds} covalent bonds "
          f"(cov-radii + {BOND_TOLERANCE:.2f} Å tolerance)")

    write_pdb(atoms, bonds, outfile)

    print(f"""
{'='*60}
  CHECKLIST before proceeding:
  1. Verify PDL_prep.pdb visually in VMD/PyMOL.
  2. Confirm total charge & spin for QM:
       Current assumption: Pd(II) square-planar, total_charge=0,
       spin_multiplicity=1 (diamagnetic d8).
       Adjust mcpb.in if ligand is cationic/anionic.
  3. Add missing H atoms:
       reduce PDL_prep.pdb > PDL_prep_H.pdb  (AmberTools)
       Then set original_pdb PDL_prep_H.pdb in mcpb.in.
  4. ion_ids in mcpb.in should match Pd serial in PDL_prep.pdb.
{'='*60}
""")

if __name__ == "__main__":
    main()
