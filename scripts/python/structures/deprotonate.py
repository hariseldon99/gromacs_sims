#!/usr/bin/env python3
import sys

def deprotonate_pdb(input_file, output_file):
    """
    Remove hydrogen atoms from a PDB file.

    Parameters:
        input_file (str): Path to the input PDB file.
        output_file (str): Path where the deprotonated PDB file will be saved.

    The function reads the input PDB file line by line. For lines that start with "ATOM" 
    or "HETATM", it checks the element column (columns 77-78) to see if it is a hydrogen (H).
    In cases where the element field is missing or ambiguous, the atom name (columns 13-16)
    is examined. If a hydrogen is detected, the line is skipped; otherwise, the line is written 
    to the output file unchanged. Note that the PDB format is fixed-width, so this script assumes 
    that the file follows standard formatting.
    """
    with open(input_file, "r") as f:
        lines = f.readlines()

    with open(output_file, "w") as outf:
        for line in lines:
            # Only check lines that define atoms.
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the element symbol (columns 77-78) if possible.
                element = line[76:78].strip() if len(line) >= 78 else ""
                # If the element symbol is missing, fallback to the atom name (columns 13-16).
                atom_name = line[12:16].strip() if len(line) >= 16 else ""
                # Determine if this atom is hydrogen:
                if element.upper() == "H" or (not element and atom_name and atom_name[0].upper() == "H"):
                    continue  # Skip writing hydrogen atoms.
            # Write all non-hydrogen lines (and non-ATOM/HETATM lines) to the output.
            outf.write(line)

def main():
    if len(sys.argv) != 3:
        print("Usage: python deprotonate.py input.pdb output.pdb")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    deprotonate_pdb(input_file, output_file)
    print(f"Deprotonated PDB file written to {output_file}")

if __name__ == "__main__":
    main()
