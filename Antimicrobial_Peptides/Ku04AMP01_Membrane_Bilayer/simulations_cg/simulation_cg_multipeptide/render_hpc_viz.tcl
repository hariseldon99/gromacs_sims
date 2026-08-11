# ==============================================================================
# Remote Headless HPC VMD Visualization & Movie Script (v6 - M       Performance Soft Viz)
# Optimized for remote HPC nodes with NVIDIA RTX graphics (TachyonLOptiXInternal)
# Gracefully falls back to TachyonInternal on machines without RTX GPUs
# ==============================================================================

puts "======================================================="
puts " Starting Headless VMD Visualization (Soft Palette & Fast) "
puts "======================================================="

set script_dir [file dirname [info script]]
if {$script_dir == "" || $script_dir == "."} {
    set script_dir [pwd]
}
set output_dir [file join $script_dir "viz_output"]
set frames_dir [file join $output_dir "frames"]

# Clean old frames
file delete -force $output_dir
file mkdir $output_dir
file mkdir $frames_dir

# 1. Load Base Visualization State
set viz_state_file [file join $script_dir "prod_dry_viz.vmd"]
if {[file exists $viz_state_file]} {
    puts "Loading visualization state from: $viz_state_file"
    source $viz_state_file
} else {
    puts "ERROR: Could not find $viz_state_file"
    exit 1
}

set top_mol [molinfo top]
if {$top_mol < 0} {
    puts "ERROR: No molecule loaded!"
    exit 1
}

set num_frames [molinfo $top_mol get numframes]
puts "Successfully loaded molecule ID $top_mol with $num_frames frames."
axes location off
# 2. Dynamic Camera Centering on Peptides & Wide Field of View
proc apply_centered_wide_view {} {
    set pep [atomselect top "resname ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"]
    set c [measure center $pep]
    $pep delete
    
    set cx [lindex $c 0]
    set cy [lindex $c 1]
    set cz [lindex $c 2]
    
    # Center camera directly on the 12 peptides
    set c_mat [list [list 1.0 0.0 0.0 [expr {-$cx}]] [list 0.0 1.0 0.0 [expr {-$cy}]] [list 0.0 0.0 1.0 [expr {-$cz}]] [list 0.0 0.0 0.0 1.0]]
    molinfo top set center_matrix [list $c_mat]
    
    # Set wide perspective scale to include most of the membrane plane (s = 0.0038)
    set s 0.0038
    set s_mat [list [list $s 0.0 0.0 0.0] [list 0.0 $s 0.0 0.0] [list 0.0 0.0 $s 0.0] [list 0.0 0.0 0.0 1.0]]
    molinfo top set scale_matrix [list $s_mat]
}

# 3. Soft, Highly Aesthetic Color Scheme Setup (No Garish Red)
puts "Applying Soft Aesthetic Color Palette..."

# Peptides: Glowing High-Luminance Custom Colors
color change rgb 4 1.00 0.85 0.00   ; # Hydrophobic: Glowing Gold
color change rgb 0 0.00 0.95 1.00   ; # Basic: Electric Cyan
color change rgb 7 1.00 0.20 0.70   ; # Acidic/Polar: Vivid Magenta
color change rgb 3 1.00 0.45 0.00   ; # Aromatic: Bright Amber

mol modmaterial 0 top AOShiny
mol modcolor 0 top ResType

# Membrane: Softer Slate-Teal Palette (Replaces garish red)
color scale colors BlkW {0.22 0.28 0.35} {0.45 0.55 0.65} {0.72 0.82 0.90}
color scale method BlkW

# Soft AOChalky material with mild translucency for depth
material change opacity AOChalky 0.88
mol modmaterial 1 top AOChalky

# Set camera centering on peptides at middle frame for optimal trajectory centering
animate goto [expr {int(($num_frames - 1) / 2)}]
apply_centered_wide_view

# 4. HPC RTX Renderer Auto-Detection
set render_list [render list]
if {[lsearch -exact $render_list "TachyonLOptiXInternal"] != -1} {
    set selected_renderer "TachyonLOptiXInternal"
    puts "Selected Renderer: TachyonLOptiXInternal (NVIDIA RTX Hardware Ray Tracing)"
} elseif {[lsearch -exact $render_list "TachyonInternal"] != -1} {
    set selected_renderer "TachyonInternal"
    puts "Selected Renderer: TachyonInternal (CPU Multi-threaded Ray Tracing)"
} else {
    set selected_renderer "snapshot"
    puts "Selected Renderer: snapshot"
}

# 5. Save Soft-Palette Centered Snapshots (Beginning, Middle, End)
set f_beg 0
set f_mid [expr {int(($num_frames - 1) / 2)}]
set f_end [expr {$num_frames - 1}]

puts "\n--- Generating Soft Aesthetic Snapshots ---"

# Beginning Snapshot
animate goto $f_beg
display update ui
set snap_beg_ppm [file join $output_dir "snapshot_beginning.ppm"]
puts "Rendering Beginning Snapshot (frame $f_beg) -> $snap_beg_ppm"
render $selected_renderer $snap_beg_ppm

# Middle Snapshot
animate goto $f_mid
display update ui
set snap_mid_ppm [file join $output_dir "snapshot_middle.ppm"]
puts "Rendering Middle Snapshot (frame $f_mid) -> $snap_mid_ppm"
render $selected_renderer $snap_mid_ppm

# End Snapshot
animate goto $f_end
display update ui
set snap_end_ppm [file join $output_dir "snapshot_end.ppm"]
puts "Rendering End Snapshot (frame $f_end) -> $snap_end_ppm"
render $selected_renderer $snap_end_ppm

# 6. Slower Trajectory Movie Generation (600 Rendered Frames)
set target_movie_frames 600
set step_size [expr {int(ceil(double($num_frames) / $target_movie_frames))}]
if {$step_size < 1} { set step_size 1 }

puts "\n--- Rendering Trajectory Movie Frames ---"
puts "  Total Simulation Frames: $num_frames"
puts "  Frame Step Size: $step_size"
puts "  Total Frames to Render: [expr {int(ceil(double($num_frames) / $step_size))}]"

set frame_idx 0

for {set i 0} {$i < $num_frames} {incr i $step_size} {
    animate goto $i
    display update ui
    set frame_ppm [file join $frames_dir [format "frame_%05d.ppm" $frame_idx]]
    if {$frame_idx % 100 == 0} {
        puts "Rendering frame [format "%5d" $frame_idx] / [format "%5d" $i] ($selected_renderer)..."
    }
    render $selected_renderer $frame_ppm
    incr frame_idx
}

puts "\nSuccessfully rendered $frame_idx frames to $frames_dir"
puts "Completed VMD rendering script."
exit 0
