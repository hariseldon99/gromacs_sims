#!/usr/bin/env python
# nglview.contrib.movie.MovieMaker requires a live Jupyter kernel to take widget
# screenshots and cannot be used in a standalone script.  This version renders
# each frame with matplotlib (Agg backend, no display required) and writes the
# MP4 directly with imageio-ffmpeg.
import argparse
import os

import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import MDAnalysis as mda
from MDAnalysis import transformations as trans
import numpy as np

# ── CLI arguments ─────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Render a GROMACS trajectory as an MP4 movie (no display needed)."
)
parser.add_argument("-t", "--topology", required=True,
                    help="Topology file (e.g. .tpr or .gro)")
parser.add_argument("-x", "--trajectory", required=True,
                    help="Trajectory file (e.g. .xtc or .trr)")
parser.add_argument("-o", "--output", required=True,
                    help="Output MP4 filename")
parser.add_argument("--lps-sel",
                    default="resname AGAL AGLC AHEP AKDO ARHM BGLC BGLCNA BMANNA",
                    help="MDAnalysis selection for the LPS / membrane sugar group "
                        "(default: GM LPS resnames)")
parser.add_argument("--lipid-sel",
                    default="resname ECLIPA PPPE PVPG PVCL2",
                    help="MDAnalysis selection for the lipid group "
                        "(default: GM inner-membrane lipid resnames)")
parser.add_argument("--protein-sel",
                    default="resname PMB",
                    help="MDAnalysis selection for the protein / peptide / ligand "
                        "of interest (default: 'resname PMB')")
parser.add_argument("--lps-label",   default="LPS sugars",
                    help="Legend label for the LPS group (default: 'LPS sugars')")
parser.add_argument("--lipid-label", default="Lipids",
                    help="Legend label for the lipid group (default: 'Lipids')")
parser.add_argument("--protein-label", default="PMB",
                    help="Legend label for the protein/peptide/ligand (default: 'PMB')")
parser.add_argument("--title", default="",
                    help="Static title prefix shown above every frame, "
                        "followed by the current simulation time (default: none)")
args = parser.parse_args()

tpr    = args.topology
xtc    = args.trajectory
output = args.output

# ── CPU count for ffmpeg threading ───────────────────────────────────────────
# Set NPROCS (or SLURM_CPUS_PER_TASK / OMP_NUM_THREADS) in your job script.
_n_cpus = (
    os.environ.get("PBS_NCPUS")
    or os.environ.get("SLURM_CPUS_PER_TASK")
    or os.environ.get("OMP_NUM_THREADS")
)
N_CPUS = int(_n_cpus) if _n_cpus else os.cpu_count() or 1
print(f"Using {N_CPUS} CPU thread(s) for ffmpeg encoding.", flush=True)

u = mda.Universe(tpr, xtc)

# ── atom selections ───────────────────────────────────────────────────────────
lps      = u.select_atoms(args.lps_sel)
lipids   = u.select_atoms(args.lipid_sel)
pmb      = u.select_atoms(args.protein_sel)
membrane = u.select_atoms(f"({args.lps_sel}) or ({args.lipid_sel})")
all_atoms = u.select_atoms("all")

print(f"LPS group    : {len(lps)} atoms  [{args.lps_sel}]", flush=True)
print(f"Lipid group  : {len(lipids)} atoms  [{args.lipid_sel}]", flush=True)
print(f"Protein/mol  : {len(pmb)} atoms  [{args.protein_sel}]", flush=True)

# ── PBC transformations ───────────────────────────────────────────────────────
workflow = [
    trans.unwrap(all_atoms),
    trans.center_in_box(membrane, wrap=True),
    trans.wrap(all_atoms, compound="residues"),
]
u.trajectory.add_transformations(*workflow)

# ── rendering parameters ──────────────────────────────────────────────────────
FPS     = 30
STRIDE  = 1          # increase to skip frames and speed up rendering
DPI     = 100
FIGSIZE = (16, 9)   # → 1600×900 px, both divisible by macro_block_size=16

# colour map
C_LPS   = "#00aa00"
C_LIPID = "#ff8800"
C_PMB   = "#ff2222"
BG      = "black"

# ── set up figure (reused every frame) ───────────────────────────────────────
fig, (ax_top, ax_side) = plt.subplots(1, 2, figsize=FIGSIZE, facecolor=BG)
for ax in (ax_top, ax_side):
    ax.set_facecolor(BG)
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")

legend_patches = [
    mpatches.Patch(color=C_LPS,   label=args.lps_label),
    mpatches.Patch(color=C_LIPID, label=args.lipid_label),
    mpatches.Patch(color=C_PMB,   label=args.protein_label),
]

n_frames = len(u.trajectory[::STRIDE])
print(f"Rendering {n_frames} frames → {output}", flush=True)

# ── one-time layout: reserve space for suptitle and legend ───────────────────
fig.subplots_adjust(left=0.07, right=0.97, bottom=0.10, top=0.93, wspace=0.35)
leg = fig.legend(handles=legend_patches, loc="lower center", ncol=3,
                facecolor=BG, labelcolor="white", framealpha=0.4,
                bbox_to_anchor=(0.5, 0.0))
title_txt = fig.suptitle("", color="white")
_title_prefix = (args.title + "  ") if args.title else ""

# ── write movie frame-by-frame (memory-efficient) ────────────────────────────
with imageio.get_writer(output, fps=FPS, quality=8,
                        output_params=["-threads", str(N_CPUS)]) as writer:
    for i, ts in enumerate(u.trajectory[::STRIDE]):
        if i % max(1, n_frames // 20) == 0:
            print(f"  frame {i}/{n_frames}  t={ts.time:.0f} ps", flush=True)

        box = ts.dimensions[:3]   # (Lx, Ly, Lz) in Å

        lps_pos   = lps.positions
        lip_pos   = lipids.positions
        pmb_pos   = pmb.positions

        # ── top-down view (XY) ────────────────────────────────────────────────
        ax_top.cla()
        ax_top.set_facecolor(BG)
        ax_top.scatter(lps_pos[:, 0],  lps_pos[:, 1],  c=C_LPS,   s=2,  alpha=0.6, linewidths=0, rasterized=True)
        ax_top.scatter(lip_pos[:, 0],  lip_pos[:, 1],  c=C_LIPID, s=2,  alpha=0.6, linewidths=0, rasterized=True)
        ax_top.scatter(pmb_pos[:, 0],  pmb_pos[:, 1],  c=C_PMB,   s=12, alpha=1.0, linewidths=0, zorder=5)
        ax_top.set_xlim(0, box[0])
        ax_top.set_ylim(0, box[1])
        ax_top.set_aspect("equal")
        ax_top.set_title("Top view (XY)", color="white")
        ax_top.set_xlabel("x (Å)", color="white")
        ax_top.set_ylabel("y (Å)", color="white")
        ax_top.tick_params(colors="white")
        for spine in ax_top.spines.values():
            spine.set_color("white")

        # ── side view (XZ) ────────────────────────────────────────────────────
        ax_side.cla()
        ax_side.set_facecolor(BG)
        ax_side.scatter(lps_pos[:, 0],  lps_pos[:, 2],  c=C_LPS,   s=2,  alpha=0.6, linewidths=0, rasterized=True)
        ax_side.scatter(lip_pos[:, 0],  lip_pos[:, 2],  c=C_LIPID, s=2,  alpha=0.6, linewidths=0, rasterized=True)
        ax_side.scatter(pmb_pos[:, 0],  pmb_pos[:, 2],  c=C_PMB,   s=12, alpha=1.0, linewidths=0, zorder=5)
        ax_side.set_xlim(0, box[0])
        ax_side.set_ylim(0, box[2])
        ax_side.set_aspect("equal")
        ax_side.set_title("Side view (XZ)", color="white")
        ax_side.set_xlabel("x (Å)", color="white")
        ax_side.set_ylabel("z (Å)", color="white")
        ax_side.tick_params(colors="white")
        for spine in ax_side.spines.values():
            spine.set_color("white")

        title_txt.set_text(f"{_title_prefix}t = {ts.time:.0f} ps")

        # Draw to the fixed-size canvas — guarantees identical dimensions every frame.
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        # buffer_rgba() → (h, w, 4) RGBA; drop alpha for imageio (RGB only).
        img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)[:, :, :3]
        writer.append_data(img)

plt.close(fig)
print(f"Movie saved to {output}")