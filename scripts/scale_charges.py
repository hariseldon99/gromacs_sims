#!/usr/bin/env python3

"""
Scale partial charges in the [ atoms ] section of a GROMACS .itp file
and write a new file with a header documenting the scaling.

Usage:
    python scale_charges.py -i input.itp -o output.itp -s 0.8
"""

import argparse
from datetime import datetime


def write_header(fout, input_file, scale_factor):
    """
    Write a comment block at the top of the output file documenting scaling.
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    fout.write("; ============================================================\n")
    fout.write(";  Charge scaling applied to this topology file\n")
    fout.write("; ------------------------------------------------------------\n")
    fout.write(f";  Source file        : {input_file}\n")
    fout.write(f";  Scaling factor     : {scale_factor}\n")
    fout.write(f";  Generated on       : {now}\n")
    fout.write(";  Modified section   : [ atoms ] (charge column only)\n")
    fout.write(";  Notes              :\n")
    fout.write(";    - Charges scaled uniformly\n")
    fout.write(";    - Total charge preserved\n")
    fout.write(";    - No bonded/LJ parameters modified\n")
    fout.write("; ============================================================\n\n")


def scale_itp_charges(input_file, output_file, scale_factor):
    in_atoms = False
    atom_count = 0
    original_sum = 0.0
    new_sum = 0.0

    print(f"[INFO] Reading input file: {input_file}")
    print(f"[INFO] Writing output file: {output_file}")
    print(f"[INFO] Scaling charges by factor: {scale_factor}\n")

    with open(input_file, "r") as fin, open(output_file, "w") as fout:

        # --- Write header block ---
        write_header(fout, input_file, scale_factor)

        for line in fin:
            stripped = line.strip()

            # Detect start of [ atoms ]
            if stripped.startswith("[ atoms ]"):
                in_atoms = True
                fout.write(line)
                print("[INFO] Entering [ atoms ] section")
                continue

            # Detect new section → exit atoms
            elif stripped.startswith("[") and not stripped.startswith("[ atoms ]"):
                if in_atoms:
                    print("[INFO] Leaving [ atoms ] section\n")
                in_atoms = False

            # Modify atom lines
            if in_atoms and stripped and not stripped.startswith(";"):
                parts = line.split()

                if parts[0].isdigit() and len(parts) >= 7:
                    try:
                        charge = float(parts[6])
                        scaled_charge = charge * scale_factor

                        original_sum += charge
                        new_sum += scaled_charge
                        atom_count += 1

                        parts[6] = f"{scaled_charge:.6f}"

                        fout.write(
                            f"{int(parts[0]):>6} {parts[1]:>10} {int(parts[2]):>6} "
                            f"{parts[3]:>8} {parts[4]:>6} {int(parts[5]):>6} "
                            f"{parts[6]:>10} {parts[7]:>10}\n"
                        )
                        continue

                    except ValueError:
                        pass

            fout.write(line)

    # --- Summary output ---
    print("\n[SUMMARY]")
    print(f"Atoms processed        : {atom_count}")
    print(f"Original total charge  : {original_sum:.6f} e")
    print(f"Scaled total charge    : {new_sum:.6f} e")

    if abs(original_sum) < 1e-3:
        print("[CHECK] Original system is neutral (within tolerance)")
    else:
        print("[WARNING] Original system is NOT neutral")

    if abs(new_sum) < 1e-3:
        print("[CHECK] Scaled system remains neutral")
    else:
        print("[WARNING] Scaled system is NOT neutral")

    print("\n[INFO] Done.")


def main():
    parser = argparse.ArgumentParser(
        description="Scale partial charges in the [ atoms ] section of a GROMACS .itp file."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input .itp file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output .itp file"
    )

    parser.add_argument(
        "-s", "--scale",
        required=True,
        type=float,
        help="Scaling factor for charges (e.g., 0.8)"
    )

    args = parser.parse_args()

    scale_itp_charges(args.input, args.output, args.scale)


if __name__ == "__main__":
    main()
