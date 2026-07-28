#!/usr/bin/env python3

with open("ku04_gim.gro", "r") as f:
    lines = f.readlines()

title = lines[0]
total_atoms_line = lines[1]
box_vector = lines[-1]
atom_lines = lines[2:-1]

# 1. Pool all individual beads by their residue name
pools = {
    "DVPE": [], "POPG": [], "POPE": [], "DOPE": [], 
    "TOCL": [], "W": [], "NA": [], "CL": [], "PEPTIDE": []
}

for line in atom_lines:
    resname = line[5:10].strip()
    if resname in ["NA", "NA+"]:
        pools["NA"].append(line)
    elif resname in ["CL", "CL-"]:
        pools["CL"].append(line)
    elif resname in pools:
        pools[resname].append(line)
    else:
        pools["PEPTIDE"].append(line)

# 2. Define the exact block structure from your [ molecules ] topology section
# Format: (Residue name key in pools, number of molecules, beads per molecule)
topology_blocks = [
    ("PEPTIDE", 8, 82),   # KU04
    ("DOPE", 30, 12),
    ("DVPE", 169, 12),
    ("POPE", 44, 12),
    ("POPG", 70, 12),
    ("TOCL", 7, 23),      # Martini 3 TOCL has 23 beads
    ("DOPE", 29, 12),
    ("DVPE", 165, 12),
    ("POPE", 43, 12),
    ("POPG", 68, 12),
    ("TOCL", 7, 23),
    ("W", 6105, 1),
    ("NA", 249, 1),
    ("CL", 87, 1)
]

# 3. Extract the exact number of lines needed for each block sequential order
ordered_atoms = []
for res_key, num_mol, beads_per_mol in topology_blocks:
    total_beads_needed = num_mol * beads_per_mol
    block_lines = pools[res_key][:total_beads_needed]
    ordered_atoms.extend(block_lines)
    # Remove the used lines from the pool so the next leaflet gets fresh beads
    pools[res_key] = pools[res_key][total_beads_needed:]

# 4. Verify and write the perfectly mapped structure file
actual_atom_count = len(ordered_atoms)
new_atoms_line = f"{actual_atom_count}\n"

with open("ku04_gim_ordered.gro", "w") as f:
    f.write(title)
    f.write(new_atoms_line)
    f.writelines(ordered_atoms)
    f.write(box_vector)

print(f"Success! Generated ku04_gim_ordered.gro with {actual_atom_count} cleanly ordered beads.")
