#!/usr/bin/env python3
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AddHs
from rdkit.Chem.rdmolfiles import MolToPDBBlock
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray


def fix_protonate_pdb(input_file, output_file, ligand_code):
    # Read the PDB file using Biotite
    pdb_file = PDBFile.read(input_file)
    atom_array = pdb_file.get_structure()
    ligand_atoms  = []
    for structure in atom_array:
        for atom in structure:
            print(atom, atom.res_name, atom.element)
            if atom.res_name == ligand_code:
                ligand_atoms.append(atom)

    # Extract bond information for the ligand
    ligand_bonds = pdb_file.get_bonds()
    ligand_bonds = ligand_bonds[ligand_atoms.array_index]
    # Remove hydrogen atoms from the ligand
    non_hydrogen_atoms = ligand_atoms[ligand_atoms.element != "H"]

    # Convert the deprotonated ligand to an RDKit molecule
    ligand_mol = Chem.MolFromPDBBlock(MolToPDBBlock(non_hydrogen_atoms))
    if ligand_mol is None:
        raise ValueError("Failed to create RDKit molecule from ligand atoms.")

    # Re-protonate the ligand using RDKit
    protonated_ligand_mol = AddHs(ligand_mol)

    # Convert the protonated ligand back to a PDB block
    protonated_ligand_pdb = Chem.MolToPDBBlock(protonated_ligand_mol)

    # Parse the protonated ligand PDB block into a Biotite AtomArray
    protonated_ligand_atoms = PDBFile.read(protonated_ligand_pdb).get_structure()

    # Combine the protein atoms with the re-protonated ligand
    atom_array = atom_array[atom_array.res_name != ligand_code]
    atom_array = atom_array + protonated_ligand_atoms

    # Write the updated structure back to the output file
    pdb_file.set_structure(atom_array)
    pdb_file.write(output_file)
    

def main():
    parser = argparse.ArgumentParser(description="Deprotonate a PDB file by removing hydrogen atoms.")
    parser.add_argument("input_pdb", help="Input PDB file name")
    parser.add_argument("output_pdb", help="Output PDB file name")
    parser.add_argument("ligand_code", help="Three-letter code for the ligand")
    args = parser.parse_args()

    input_file = args.input_pdb
    output_file = args.output_pdb
    ligand_code = args.ligand_code
   
    fix_protonate_pdb(input_file, output_file, ligand_code)
    print(f"PDB file with ligand protonation fixed is written to {output_file}")

if __name__ == "__main__":
    main()