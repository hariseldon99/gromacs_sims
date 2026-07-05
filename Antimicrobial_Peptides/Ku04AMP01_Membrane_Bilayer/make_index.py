#!/usr/bin/env python3
"""
Generate index.ndx with SOLU, MEMB, and SOLV groups for the KU04 + G-IM system.

SOLU = KU04 Peptide (Standard Amino Acid Residues)
MEMB = Lipids
SOLV = Water + Ions
SOLU_MEMB = Peptide + Lipids (Integrated group for COM removal)
"""

# Standard amino acid residue names for the peptide
SOLU_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

# G-IM Membrane lipids
MEMB_RES = {
    "PMPE", "POPE", "QMPE", "OYPE",   # PE lipids
    "PMPG", "PYPG",                    # PG lipids
    "PVCL2",                           # cardiolipin
}

# Solvent and ions
SOLV_RES = {"SOD", "CLA", "TIP3"}

def write_group(f, name, atoms):
    f.write(f"[ {name} ]\n")
    for i, a in enumerate(atoms):
        f.write(f"{a:7d}")
        if (i + 1) % 10 == 0:
            f.write("\n")
    if atoms and len(atoms) % 10 != 0:
        f.write("\n")
    f.write("\n")

def main(gro="system.gro", ndx="index.ndx"):
    solu_atoms = []
    memb_atoms = []
    solv_atoms = []
    unknown = set()

    with open(gro) as f:
        f.readline() # Title
        line2 = f.readline().strip()
        if not line2: return
        natoms = int(line2)
        
        for i in range(natoms):
            line = f.readline()
            resname = line[5:10].strip()
            atom_nr = i + 1
            
            if resname in SOLU_RES:
                solu_atoms.append(atom_nr)
            elif resname in MEMB_RES:
                memb_atoms.append(atom_nr)
            elif resname in SOLV_RES:
                solv_atoms.append(atom_nr)
            else:
                unknown.add(resname)

    with open(ndx, "w") as f:
        write_group(f, "SOLU", solu_atoms)
        write_group(f, "MEMB", memb_atoms)
        write_group(f, "SOLV", solv_atoms)
        # Required for comm_grps in mdp
        write_group(f, "SOLU_MEMB", solu_atoms + memb_atoms)

    print(f"Generated {ndx}:")
    print(f"  SOLU: {len(solu_atoms)} atoms")
    print(f"  MEMB: {len(memb_atoms)} atoms")
    print(f"  SOLV: {len(solv_atoms)} atoms")
    print(f"  SOLU_MEMB: {len(solu_atoms) + len(memb_atoms)} atoms")
    
    if unknown:
        print(f"WARNING: Unassigned residues: {unknown}")

if __name__ == "__main__":
    main()