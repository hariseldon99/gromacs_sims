#!/usr/bin/env python3
import sys
import argparse

def split_protein_ligand(input_file, protein_output, ligand_output, ligand_resname):
    """
    Split a PDBQT file into separate files for the protein and ligand based on the ligand residue name.

    Parameters:
        input_file (str): Path to the input PDBQT file.
        protein_output (str): Path where the protein PDBQT file will be saved.
        ligand_output (str): Path where the ligand PDBQT file will be saved.
        ligand_resname (str): Residue name of the ligand to identify it in the file.

    The function writes all lines corresponding to the ligand (matching the residue name)
    to the ligand output file, and the rest to the protein output file.
    """
    with open(input_file, "r") as f:
        lines = f.readlines()

    with open(protein_output, "w") as protein_out, open(ligand_output, "w") as ligand_out:
        for line in lines:
            # Check if the line corresponds to the ligand based on the residue name (columns 18-20).
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == ligand_resname:
                ligand_out.write(line)
            else:
                protein_out.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Split a PDBQT file into separate files for the protein and ligand based on the ligand residue name."
    )
    parser.add_argument(
        "input_file", 
        type=str, 
        help="Path to the input PDBQT file."
    )
    parser.add_argument(
        "protein_output", 
        type=str, 
        help="Path where the protein PDBQT file will be saved."
    )
    parser.add_argument(
        "ligand_output", 
        type=str, 
        help="Path where the ligand PDBQT file will be saved."
    )
    parser.add_argument(
        "ligand_resname", 
        type=str, 
        help="Residue name of the ligand to identify it in the file."
    )
    
    args = parser.parse_args()
    
    split_protein_ligand(
        args.input_file, 
        args.protein_output, 
        args.ligand_output, 
        args.ligand_resname
    )
    print(f"Protein PDBQT file written to {args.protein_output}")
    print(f"Ligand PDBQT file written to {args.ligand_output}")

if __name__ == "__main__":
    main()