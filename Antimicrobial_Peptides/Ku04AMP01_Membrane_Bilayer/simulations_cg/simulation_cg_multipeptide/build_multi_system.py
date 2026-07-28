#!/usr/bin/env python
import COBY

input_pdb = "KU04AMP01_cg.pdb"

# 1. Read coordinates from single peptide
atoms = []
with open(input_pdb, "r") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            atoms.append({
                "line": line,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54])
            })

# Calculate geometric center
cx = sum(a["x"] for a in atoms) / len(atoms)
cy = sum(a["y"] for a in atoms) / len(atoms)
cz = sum(a["z"] for a in atoms) / len(atoms)

# COBY internal frame is centered at (0, 0, 0): [-L/2, +L/2]
# For a 140x140x120 Å box:
# X and Y boundaries: -70.0 Å to +70.0 Å
# Z boundaries: -60.0 Å to +60.0 Å (membrane center at Z = 0)
x_grid = [-40.0, 0.0, 40.0]
y_grid = [-45.0, -15.0, 15.0, 45.0]
z_target = 32.0  # +3.2 nm in COBY frame (solvent phase above upper leaflet at ~ +2.0 nm)

grid_positions = [(x, y) for x in x_grid for y in y_grid]

# 2. Generate 12 individual single-peptide PDB files in COBY's centered frame
protein_arg_list = []
for idx, (gx, gy) in enumerate(grid_positions):
    filename = f"ku04_copy_{idx}.pdb"
    out_lines = []
    atom_id = 1
    
    for a in atoms:
        # Rotate 90 degrees around Y-axis (long helix axis aligned along X)
        rx = a["z"] - cz
        ry = a["y"] - cy
        rz = -(a["x"] - cx)
        
        # Translate to centered grid position
        nx = rx + gx
        ny = ry + gy
        nz = rz + z_target
        
        line = a["line"]
        new_line = (
            f"{line[:6]}{atom_id:5d}{line[11:30]}"
            f"{nx:8.3f}{ny:8.3f}{nz:8.3f}"
            f"{line[54:]}"
        )
        out_lines.append(new_line)
        atom_id += 1
        
    with open(filename, "w") as f:
        f.writelines(out_lines)
        f.write("END\n")
        
    protein_arg_list.append(f"file:{filename} moleculetypes:KU04 center_protein:False")

print("Generated 12 centered PDB files for COBY.")

# 3. Build system in COBY
sysname = "ku04_gim_multi"

COBY.COBY(
    box=[14, 14, 12],
    membrane="lipid:DVPE:46 lipid:POPG:19 lipid:POPE:25 lipid:DOPE:8 lipid:TOCL:2 apl:0.61",
    protein=protein_arg_list,
    solvation="default",
    itp_input="file:top_for_COBY.itp",
    molecule_import="file:tocl_single.gro moleculetypes:TOCL",
    out_sys=sysname,
    out_top=sysname + ".top",
    out_log=sysname + ".log",
    sn=sysname,
)