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
load mapfunc.cub, esp_grid

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

# 6. Color mapping ramp using data grid thresholds
ramp_new esp_colors, esp_grid, [-0.05, 0.00, 0.05], [red, white, blue]
color esp_colors, esp_surface

# 7. Set surface transparency
set transparency, 0.4, esp_surface

# 8. Hide the automatic raw data scale bar
disable esp_colors

# 9. TOP CENTER TEXT OVERLAY: Tuned Angstrom offset from central Pd core
set label_color, black
set label_font_id, 7
set label_size, 24
set depth_cue, 0

# Offset syntax in Angstroms: [X, Y, Z]
# -4 shifts the start slightly left so the long string balances perfectly in the middle
# 10 pushes it higher up, completely clearing the top edge of the surface
set label_position, [0, 9.5, 0]

# Print the title using the Palladium atom as the core anchor
label element Pd and ligand, "ESP Surface Range: -31.4 to +31.4 kcal/mol (Red to Blue)"

# 10. Snap the camera target cleanly onto the complex
zoom ligand
