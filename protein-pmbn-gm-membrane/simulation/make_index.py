#!/usr/bin/env python3
"""
Generate index.ndx with MEMB and SOLV groups for the PMB + G-OM membrane system.

MEMB = all lipid/LPS/PMB atoms (ECLIPA sugar residues + inner leaflet lipids + PMB)
SOLV = water + ions

Run from the simulation directory:
    python3 make_index.py
Then pass -n index.ndx to gmx grompp.
"""

# All residue names that belong to the membrane / solute group.
# ECLIPA (LPSA in topology) is split into multiple glycan residue types in the GRO.
MEMB_RES = {
    "ECLIP", "ARHM", "AGAL", "BMANN", "BGLCN",   # LPS outer core / Lipid A
    "AHEP", "AKDO", "AGLC", "BGLC",               # LPS inner core
    "LPSA",                                         # fallback if written as one residue
    "PPPE", "PVPG", "PVCL2", "PVCL",              # inner leaflet lipids
    "PMB",                                          # Polymyxin B
}

SOLV_RES = {"CAL", "SOD", "CLA", "TIP3"}


def write_group(f, name, atoms):
    f.write(f"[ {name} ]\n")
    for i, a in enumerate(atoms):
        f.write(f"{a:7d}")   # 7-wide: guarantees a space separator even for 6-digit atom numbers
        if (i + 1) % 10 == 0:
            f.write("\n")
    if atoms and len(atoms) % 10 != 0:
        f.write("\n")
    f.write("\n")


def main(gro="system.gro", ndx="index.ndx"):
    memb_atoms = []
    solv_atoms = []
    unknown = set()

    with open(gro) as f:
        f.readline()                        # title
        natoms = int(f.readline().strip())
        for i in range(natoms):
            line = f.readline()
            resname = line[5:10].strip()
            atom_nr = i + 1
            if resname in MEMB_RES:
                memb_atoms.append(atom_nr)
            elif resname in SOLV_RES:
                solv_atoms.append(atom_nr)
            else:
                unknown.add(resname)

    with open(ndx, "w") as f:
        write_group(f, "MEMB", memb_atoms)
        write_group(f, "SOLV", solv_atoms)

    print(f"MEMB: {len(memb_atoms)} atoms")
    print(f"SOLV: {len(solv_atoms)} atoms")
    print(f"Total: {len(memb_atoms) + len(solv_atoms)} / {natoms}")
    if unknown:
        print(f"WARNING - unassigned residues (not in MEMB or SOLV): {unknown}")
        print("Add them to MEMB_RES or SOLV_RES in this script.")


if __name__ == "__main__":
    main()
