#!/usr/bin/env python3
import subprocess
import os
import argparse
from Bio.PDB import PDBParser
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel

def protonate_protein(protein_pdb):
    """
    Protonate the protein using OpenBabel.
    """

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")

    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, protein_pdb):
        raise ValueError(f"Could not read protein PDB file: {protein_pdb}")

    # Add hydrogens
    mol.AddHydrogens()

    protonated_pdb = protein_pdb.replace('.pdb', '_h.pdb')
    if not obConversion.WriteFile(mol, protonated_pdb):
        raise ValueError(f"Could not write protonated protein PDB file: {protonated_pdb}")

    print(f"Protein protonated successfully: {protonated_pdb}")
    return protonated_pdb

def protonate_ligand(ligand_pdb):
        """
        Protonate the ligand using RDKit.
        """
        print(ligand_pdb)
        ligand_mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
        if ligand_mol is None:
            raise ValueError(f"Could not parse ligand PDB file: {ligand_pdb}")
        
        # Add hydrogens
        ligand_mol = Chem.AddHs(ligand_mol)
        
        # Perform geometry optimization
        AllChem.EmbedMolecule(ligand_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(ligand_mol)
        # Get the residue name from the ligand PDB
        with open(ligand_pdb, 'r') as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    residue_name = line[17:20].strip()
                    print(f"Residue name found: {residue_name}")
                    break
                else:
                    raise ValueError("Could not determine residue name from ligand PDB file.")

        
        # Write the protonated ligand back to PDB format
        protonated_pdb = ligand_pdb.replace('.pdb', '_h.pdb')
        with open(protonated_pdb, 'w') as f:
            f.write(Chem.MolToPDBBlock(ligand_mol))
        
        print(f"Ligand protonated successfully: {protonated_pdb}")
        # Load the protonated PDB and replace residues named "UNL" with the residue name
        lines = []
        with open(protonated_pdb, 'r') as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    if line[17:20].strip() == "UNL":
                        line = line[:17] + f"{residue_name:<3}" + line[20:]
                lines.append(line)
        
        with open(protonated_pdb, 'w') as f:
            f.writelines(lines)
        return protonated_pdb

def prepare_ligand(ligand_pdb, protonate=True):
    """
    Prepare the ligand using ADFR.
    """
    # Protonate the ligand before preparing it
    if protonate:
        ligand_pdb = protonate_ligand(ligand_pdb)
    ligand_pdbqt = ligand_pdb.replace('.pdb', '.pdbqt')
    try:
        subprocess.run(['prepare_ligand4', 
                        '-l', ligand_pdb, 
                        '-o', ligand_pdbqt], check=True)
        print(f"Ligand prepared successfully: {ligand_pdbqt}")
        return ligand_pdbqt
    except subprocess.CalledProcessError as e:
        print(f"Error preparing ligand: {e}")
        raise

def prepare_receptor(receptor_pdb, protonate=True):
    """
    Prepare the receptor using ADFR.
    """
    if protonate:
        receptor_pdb = protonate_protein(receptor_pdb)
    receptor_pdbqt = receptor_pdb.replace('.pdb', '.pdbqt')
    try:
        subprocess.run(['prepare_receptor4', 
                        '-r', receptor_pdb, 
                        '-o', receptor_pdbqt], check=True)
        print(f"Receptor prepared successfully: {receptor_pdbqt}")
        return receptor_pdbqt
    except subprocess.CalledProcessError as e:
        print(f"Error preparing receptor: {e}")
        raise

def calculate_center_of_mass(receptor_pdb, geometry=False):
    """
    Calculate the center of mass of the protein using Biopython.
    The geometry parameter determines whether to use the geometry or the mass of the atoms.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', receptor_pdb)
    atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
    
    total_mass = 0.0
    weighted_sum = [0.0, 0.0, 0.0]
    for atom in atoms:
        mass = 1.0 if geometry else atom.mass
        total_mass += mass
        weighted_sum[0] += atom.coord[0] * mass
        weighted_sum[1] += atom.coord[1] * mass
        weighted_sum[2] += atom.coord[2] * mass
    com =  np.array([coord / total_mass for coord in weighted_sum])
    return com

def estimate_box_size(receptor_pdb):
    """
    Estimate the docking box size based on the receptor structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', receptor_pdb)
    atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
    
    min_coords = np.min([atom.coord for atom in atoms], axis=0)
    max_coords = np.max([atom.coord for atom in atoms], axis=0)
    
    box_size = max_coords - min_coords
    return box_size

def run_docking(receptor_pdbqt, 
                ligand_pdbqt, 
                output_pdbqt, 
                center_x, center_y, center_z, 
                size_x, size_y, size_z, 
                nprocs=12):
    """
    Run the docking using AutoDock Vina.
    """
    try:
        subprocess.run([
            'vina',
            '--receptor', receptor_pdbqt, 
            '--ligand', ligand_pdbqt,
            '--exhaustiveness', '32', 
            '--num_modes', '1', 
            '--cpu', str(nprocs),
            '--center_x', str(center_x), '--center_y', str(center_y), '--center_z', str(center_z),
            '--size_x', str(size_x), '--size_y', str(size_y), '--size_z', str(size_z),
            '--out', output_pdbqt
        ], check=True)
        print(f"Docking completed successfully: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        print(f"Error during docking: {e}")
        raise

def docking(protein_pdb, ligand_pdb, output_pdbqt, nprocs, verbose=False):

    if verbose:
        print("Verbose mode enabled.")
        print("Preparing receptor ...")

    receptor_pdbqt = prepare_receptor(protein_pdb)
    
    if verbose:
        print("Verbose mode enabled.")
        print("Preparing ligand ...")
    ligand_pdbqt = prepare_ligand(ligand_pdb)

    if verbose:
        print("Calculating center of mass and box size ...")

    center_x, center_y, center_z = calculate_center_of_mass(protein_pdb)
    box_size = estimate_box_size(protein_pdb)
    size_x, size_y, size_z = box_size

    if verbose:
        print(f"Center of mass: {center_x}, {center_y}, {center_z}")
        print(f"Estimated box size: {box_size}")

    if verbose:
        print("Running docking ...")
        run_docking(receptor_pdbqt, ligand_pdbqt, output_pdbqt, 
                center_x, center_y, center_z, 
                size_x, size_y, size_z, 
                nprocs=nprocs)

def batch_docking(receptor_pdbqt, ligand_pdbqtdir, 
                  center_x, center_y, center_z, 
                  size_x, size_y, size_z,
                  nprocs, verbose=False):
    """
    Perform batch docking for multiple ligands.
    """
    receptor_name = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
    output_dir = os.path.join(os.getcwd(), f'{receptor_name}_docking_results')
    os.makedirs(output_dir)
    if verbose:
        print("Verbose mode enabled.")
    ligand_pdbqts = [os.path.join(ligand_pdbqtdir, f) for f in os.listdir(ligand_pdbqtdir) if f.endswith('.pdbqt')]
    command_list = ['vina', '--receptor', receptor_pdbqt,
                    '--exhaustiveness', '32', 
                    '--num_modes', '1', 
                    '--cpu', str(nprocs),
                    '--center_x', str(center_x), '--center_y', str(center_y), '--center_z', str(center_z),
                    '--size_x', str(size_x), '--size_y', str(size_y), '--size_z', str(size_z),'--dir', output_dir]
    for ligand_pdbqt in ligand_pdbqts:
        command_list.extend(['--batch', ligand_pdbqt])
    try:
        subprocess.run(command_list, check=True)
        print(f"Docking completed successfully @: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error during docking: {e}")
        raise


def main(nprocs = 12):
    parser = argparse.ArgumentParser(description='Blind Molecular Docking using AutoDock Vina and AutoLigand.')
    parser.add_argument('receptor_pdb', type=str, help='Path to the receptor PDB file.')
    parser.add_argument('ligand_pdb', type=str, help='Path to the ligand PDB file.')
    parser.add_argument('--output', type=str, default='output.pdbqt', help='Output PDBQT file for the docking results.')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output.')
    args = parser.parse_args()

    
    args = parser.parse_args()
    if args.verbose:
        print("Verbose mode enabled.")
        print("Preparing receptor ...")

    receptor_pdbqt = prepare_receptor(args.receptor_pdb)
    
    if args.verbose:
        print("Verbose mode enabled.")
        print("Preparing ligand ...")
    ligand_pdbqt = prepare_ligand(args.ligand_pdb)

    if args.verbose:
        print("Calculating center of mass and box size ...")

    center_x, center_y, center_z = calculate_center_of_mass(args.receptor_pdb, geometry=True)
    box_size = estimate_box_size(args.receptor_pdb)
    size_x, size_y, size_z = box_size

    if args.verbose:
        print(f"Center of mass: {center_x}, {center_y}, {center_z}")
        print(f"Estimated box size: {box_size}")

    if args.verbose:
        print("Running docking ...")
    run_docking(receptor_pdbqt, ligand_pdbqt, args.output, 
                center_x, center_y, center_z, 
                size_x, size_y, size_z, 
                nprocs=nprocs)

if __name__ == "__main__":
    main()