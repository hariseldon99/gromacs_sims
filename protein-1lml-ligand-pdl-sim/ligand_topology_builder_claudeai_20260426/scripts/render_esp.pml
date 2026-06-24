# ==============================================================================
# MULTIWFN INPUT PREPARATION REFERENCE (For Manuscript Methodology)
# ==============================================================================
# 1. Load formatted checkpoint: `multiwfn PDL_small_opt.fchk`
# 2. Enter Menu 12 (Quantitative molecular surface analysis), choose Option 0.
# 3. From Post-processing menu:
#    - Type `13` to export mapped ESP grid -> saves as 'mapfunc.cub'
#    - Type `-2` to export 0.001 a.u. spatial density surface -> saves as 'surf.cub'
#    - Type `5` to export matching atom geometry -> saves as 'PDL_final_multiwfn.pdb'
# ==============================================================================

# 1. Load the molecular structure and Multiwfn grid outputs
load PDL_final_multiwfn.pdb, ligand
load surf.cub, density_grid

# Pre-scale ESP cube from a.u. to kcal/mol so the ramp bar shows kcal/mol units
python
def _scale_cube(infile, outfile, factor=627.5095):
    with open(infile, 'r') as f:
        lines = f.readlines()
    natom = abs(int(lines[2].split()[0]))
    header_end = 6 + natom
    with open(outfile, 'w') as f:
        f.writelines(lines[:header_end])
        for line in lines[header_end:]:
            vals = line.split()
            if vals:
                f.write('  '.join(f'{float(v)*factor:14.6E}' for v in vals) + '\n')
            else:
                f.write(line)
_scale_cube('mapfunc.cub', 'mapfunc_kcal.cub')
python end
load mapfunc_kcal.cub, esp_grid

# 2. Automatically form the missing Pd-N coordination bonds
bond element Pd, element N within 2.5 of element Pd

# 3. Automatic Ball-and-Stick Styling
hide everything
show sticks, ligand
show spheres, ligand
set stick_radius, 0.12, ligand
set sphere_scale, 0.25, ligand

# 4. Color the elements conventionally
color green, element C and ligand
color blue, element N and ligand
color red, element O and ligand
color orange, element Pd and ligand
color white, element H and ligand

# 5. Generate a smooth, solid 3D surface at the vdW edge (0.001 a.u.)
isosurface esp_surface, density_grid, 0.001

# 6. Color mapping ramp using data grid thresholds (31.4 kcal/mol calculated from 0.05 a.u. limit in Multiwfn)
ramp_new esp_colors, esp_grid, [-31.4, 0.0, 31.4], [red, white, blue]
color esp_colors, esp_surface

# 7. Set surface transparency
set transparency, 0.4, esp_surface

# 8. Enable the screen-anchored color scale bar (values in kcal/mol)
enable esp_colors

# 9. Snap the camera target cleanly onto the complex
zoom ligand
