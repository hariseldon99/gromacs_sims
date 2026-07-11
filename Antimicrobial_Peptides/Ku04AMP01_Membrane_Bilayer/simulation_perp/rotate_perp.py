import argparse
import numpy as np
import MDAnalysis as mda

def main():
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(
        description="Rotate a peptide around its Center of Mass (COM) using MDAnalysis."
    )
    
    # Required positional arguments for file paths
    parser.add_argument("input_file", help="Path to the input file (e.g., input.gro, input.pdb)")
    parser.add_argument("output_file", help="Path to save the rotated output file (e.g., output.gro)")
    
    # Optional arguments to easily tweak the rotation
    parser.add_argument("-a", "--angle", type=float, default=90.0, help="Rotation angle in degrees (default: 90.0)")
    parser.add_argument("-x", "--axis", choices=['x', 'y', 'z'], default='y', help="Axis to rotate around (default: y)")
    
    args = parser.parse_args()

    # 1. Load the universe
    u = mda.Universe(args.input_file)
    
    # 2. Select the peptide atoms
    peptide = u.select_atoms("protein")
    if len(peptide) == 0:
        raise ValueError("No protein atoms found in the selection! Check your input file formatting.")

    # 3. Calculate Center of Mass
    peptide_com = peptide.center_of_mass()
    
    # 4. Map the chosen axis letter to a directional vector
    axis_mapping = {
        'x': [1, 0, 0],
        'y': [0, 1, 0],
        'z': [0, 0, 1]
    }
    direction_vector = axis_mapping[args.axis]

    # 5. Perform the in-place rotation
    peptide.rotateby(angle=args.angle, axis=direction_vector, point=peptide_com)
    
    # 6. Save the new coordinates
    peptide.write(args.output_file)
    print(f"Successfully rotated {args.input_file} by {args.angle}° around the {args.axis}-axis via its COM.")
    print(f"Output saved to: {args.output_file}")

if __name__ == "__main__":
    main()
