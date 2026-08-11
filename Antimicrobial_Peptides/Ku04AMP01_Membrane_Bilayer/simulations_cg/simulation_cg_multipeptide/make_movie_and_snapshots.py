#!/usr/bin/env python3
"""
Automation script for Headless HPC VMD Rendering & Video Encoding (v3).
Soft aesthetic membrane color palette (no garish red), dynamic peptide centering,
wide membrane view, OptiX GPU ray tracing for HPC nodes, and slow 60 fps video encoding.
"""

import os
import sys
import subprocess

def run_cmd(cmd, cwd=None):
    print(f"[RUNNING]: {cmd}")
    res = subprocess.run(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(res.stdout)
    if res.returncode != 0:
        print(f"[ERROR]: Command failed with exit code {res.returncode}")
        sys.exit(res.returncode)
    return res.stdout

def main():
    work_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(work_dir, "viz_output")
    frames_dir = os.path.join(output_dir, "frames")
    
    # 1. Run VMD rendering script
    vmd_cmd = "vmd -dispdev text -e render_hpc_viz.tcl"
    print("=== Step 1: Running VMD Headless OptiX GPU Ray Tracing ===")
    run_cmd(vmd_cmd, cwd=work_dir)
    
    # 2. Convert Snapshots PPM to PNG
    print("=== Step 2: Converting Soft-Palette Snapshots to PNG ===")
    for snap_name in ["snapshot_beginning", "snapshot_middle", "snapshot_end"]:
        ppm_file = os.path.join(output_dir, f"{snap_name}.ppm")
        png_file = os.path.join(output_dir, f"{snap_name}.png")
        if os.path.exists(ppm_file):
            convert_cmd = f"ffmpeg -y -i {ppm_file} {png_file}"
            run_cmd(convert_cmd, cwd=work_dir)
            print(f"-> Saved PNG snapshot: {png_file}")
            
    # 3. Encode Movie Frames at 60 fps (Standard & Slow Playback)
    print("=== Step 3: Encoding 60 FPS Trajectory Movies ===")
    frame_pattern = os.path.join(frames_dir, "frame_%05d.ppm")
    
    # Video A: 60 FPS Standard Video
    mp4_60fps = os.path.join(output_dir, "trajectory_movie_60fps.mp4")
    ffmpeg_60_cmd = (
        f"ffmpeg -y -r 60 -i {frame_pattern} "
        f"-c:v libx264 -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -pix_fmt yuv420p -crf 18 {mp4_60fps}"
    )
    run_cmd(ffmpeg_60_cmd, cwd=work_dir)
    print(f"-> Generated 60 FPS Video (10s duration): {mp4_60fps}")
    
    # Video B: 60 FPS Slowed-Down Video (20s duration, zero choppiness)
    mp4_60fps_slow = os.path.join(output_dir, "trajectory_movie_60fps_slow.mp4")
    ffmpeg_slow_cmd = (
        f"ffmpeg -y -i {mp4_60fps} "
        f"-filter:v 'setpts=2.0*PTS,fps=60' "
        f"-c:v libx264 -pix_fmt yuv420p -crf 18 {mp4_60fps_slow}"
    )
    run_cmd(ffmpeg_slow_cmd, cwd=work_dir)
    print(f"-> Generated Slow 60 FPS Video (20s duration, zero choppiness): {mp4_60fps_slow}")

    # Summary
    print("\n=== Summary of Generated Artifacts ===")
    for name in ["snapshot_beginning.png", "snapshot_middle.png", "snapshot_end.png"]:
        path = os.path.join(output_dir, name)
        if os.path.exists(path):
            size = os.path.getsize(path) / 1024.0
            print(f"  - Snapshot: {name} ({size:.1f} KB)")
            
    for vname in ["trajectory_movie_60fps.mp4", "trajectory_movie_60fps_slow.mp4"]:
        vpath = os.path.join(output_dir, vname)
        if os.path.exists(vpath):
            vsize = os.path.getsize(vpath) / (1024.0 * 1024.0)
            print(f"  - Video: {vname} ({vsize:.2f} MB)")

if __name__ == "__main__":
    main()
