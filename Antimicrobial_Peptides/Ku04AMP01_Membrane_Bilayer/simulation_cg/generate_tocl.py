# Fixed-width GROMACS format string builder
lines = [
    "Single TOCL molecule for COBY",
    "23"
]

# (Residue ID, Residue Name, Atom Name, Atom ID, X, Y, Z)
atoms_data = [
    (1, "TOCL", "GLC", 1, 2.5, 2.5, 2.5),
    (1, "TOCL", "PO41", 2, 2.5, 2.5, 2.4),
    (1, "TOCL", "GL11", 3, 2.4, 2.4, 2.3),
    (1, "TOCL", "GL21", 4, 2.4, 2.6, 2.3),
    (1, "TOCL", "C1A1", 5, 2.3, 2.3, 2.2),
    (1, "TOCL", "D2A1", 6, 2.2, 2.2, 2.1),
    (1, "TOCL", "C3A1", 7, 2.1, 2.1, 2.0),
    (1, "TOCL", "C4A1", 8, 2.0, 2.0, 1.9),
    (1, "TOCL", "C1B1", 9, 2.3, 2.7, 2.2),
    (1, "TOCL", "D2B1", 10, 2.2, 2.8, 2.1),
    (1, "TOCL", "C3B1", 11, 2.1, 2.9, 2.0),
    (1, "TOCL", "C4B1", 12, 2.0, 3.0, 1.9),
    (1, "TOCL", "PO42", 13, 2.5, 2.5, 2.6),
    (1, "TOCL", "GL12", 14, 2.6, 2.4, 2.7),
    (1, "TOCL", "GL22", 15, 2.6, 2.6, 2.7),
    (1, "TOCL", "C1A2", 16, 2.7, 2.3, 2.8),
    (1, "TOCL", "D2A2", 17, 2.8, 2.2, 2.9),
    (1, "TOCL", "C3A2", 18, 2.9, 2.1, 3.0),
    (1, "TOCL", "C4A2", 19, 3.0, 2.0, 3.1),
    (1, "TOCL", "C1B2", 20, 2.7, 2.7, 2.8),
    (1, "TOCL", "D2B2", 21, 2.8, 2.8, 2.9),
    (1, "TOCL", "C3B2", 22, 2.9, 2.9, 3.0),
    (1, "TOCL", "C4B2", 23, 3.0, 3.0, 3.1)
]

for atom in atoms_data:
    # Standard GROMACS format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
    line = f"{atom[0]:5d}{atom[1]:<5s}{atom[2]:>5s}{atom[3]:5d}{atom[4]:8.3f}{atom[5]:8.3f}{atom[6]:8.3f}"
    lines.append(line)

lines.append("   5.00000   5.00000   5.00000")

with open("tocl_single.gro", "w") as f:
    f.write("\n".join(lines) + "\n")

print("Created strict formatted tocl_single.gro successfully!")