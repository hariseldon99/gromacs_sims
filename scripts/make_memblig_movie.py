#!/usr/bin/env python
"""make_memblig_movie – render GROMACS membrane/ligand trajectories.

Can be used as a **CLI tool** or imported as a **Python module**::

    from make_memblig_movie import build_universe, make_figure, draw_frame
    from make_memblig_movie import render_movie, render_single_frame

Public API
----------
build_universe(topology, trajectory, lps_sel, lipid_sel, protein_sel)
    Build an MDAnalysis Universe with PBC corrections applied.
    Returns (universe, lps, lipids, protein) atom groups.

make_figure(lps_label, lipid_label, protein_label, figsize, bg, c_lps, c_lipid, c_protein)
    Create and configure the matplotlib figure.
    Returns (fig, ax_top, ax_side).

draw_frame(fig, ax_top, ax_side, ts, lps, lipids, protein, title_prefix, bg, c_lps, c_lipid, c_protein)
    Draw one trajectory timestep onto the figure.
    Returns an (H, W, 3) uint8 RGB numpy array.

render_single_frame(topology, trajectory, output, frame=0, ...)
    Render one frame from a trajectory and save it as an image file.

render_movie(topology, trajectory, output, ...)
    Render the full trajectory as an MP4 movie.

CLI usage
---------
Run ``python make_memblig_movie.py --help`` for CLI options.
Pass ``--frame N`` to render a single frame instead of the whole movie.
"""
# nglview.contrib.movie.MovieMaker requires a live Jupyter kernel to take widget
# screenshots and cannot be used in a standalone script.  This version renders
# each frame with matplotlib (Agg backend, no display required) and writes the
# MP4 directly with imageio-ffmpeg.

import argparse
import os

import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")          # safe headless default; call set_backend() to override
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import MDAnalysis as mda
from MDAnalysis import transformations as trans
import numpy as np


def set_backend(backend: str) -> None:
    """Switch the matplotlib rendering backend.

    Must be called **before** :func:`make_figure` (i.e. before the first
    pyplot window is created).  Has no effect if the same backend is already
    active.

    Parameters
    ----------
    backend:
        Any backend name accepted by matplotlib, e.g.

        * ``"Agg"``        – off-screen / headless (module default)
        * ``"TkAgg"``      – interactive Tk window
        * ``"Qt5Agg"``     – interactive Qt5 window
        * ``"WxAgg"``      – interactive wx window
        * ``"module://matplotlib_inline.backend_inline"`` – Jupyter inline

    Example
    -------
    >>> import make_memblig_movie as mmm
    >>> mmm.set_backend("TkAgg")          # show figures in a Tk window
    >>> fig, ax_top, ax_side = mmm.make_figure()
    """
    plt.switch_backend(backend)

# ── Default selection strings ─────────────────────────────────────────────────
DEFAULT_LPS_SEL     = "resname AGAL AGLC AHEP AKDO ARHM BGLC BGLCNA BMANNA"
DEFAULT_LIPID_SEL   = "resname ECLIPA PPPE PVPG PVCL2"
DEFAULT_PROTEIN_SEL = "resname PMB"

# ── Default rendering parameters ─────────────────────────────────────────────
DEFAULT_FPS     = 30
DEFAULT_STRIDE  = 1
DEFAULT_DPI     = 100
DEFAULT_FIGSIZE = (16, 9)   # → 1600×900 px, both divisible by macro_block_size=16

# ── Default colours ───────────────────────────────────────────────────────────
DEFAULT_C_LPS     = "#00aa00"
DEFAULT_C_LIPID   = "#ff8800"
DEFAULT_C_PROTEIN = "#ff2222"
DEFAULT_BG        = "black"


# ─────────────────────────────────────────────────────────────────────────────
# Public helpers
# ─────────────────────────────────────────────────────────────────────────────

def _n_cpus_from_env() -> int:
    """Return the CPU count advertised by the job scheduler, or os.cpu_count()."""
    _val = (
        os.environ.get("PBS_NCPUS")
        or os.environ.get("SLURM_CPUS_PER_TASK")
        or os.environ.get("OMP_NUM_THREADS")
    )
    return int(_val) if _val else (os.cpu_count() or 1)


def build_universe(
    topology: str,
    trajectory: str,
    lps_sel: str = DEFAULT_LPS_SEL,
    lipid_sel: str = DEFAULT_LIPID_SEL,
    protein_sel: str = DEFAULT_PROTEIN_SEL,
):
    """Build an MDAnalysis Universe with PBC unwrap/center/wrap applied.

    Parameters
    ----------
    topology:
        Topology file path (.tpr, .gro, …).
    trajectory:
        Trajectory file path (.xtc, .trr, …).
    lps_sel:
        MDAnalysis selection string for the LPS / membrane sugar group.
    lipid_sel:
        MDAnalysis selection string for the lipid group.
    protein_sel:
        MDAnalysis selection string for the molecule of interest.

    Returns
    -------
    u : mda.Universe
    lps : AtomGroup
    lipids : AtomGroup
    protein : AtomGroup
    """
    u = mda.Universe(topology, trajectory)

    lps      = u.select_atoms(lps_sel)
    lipids   = u.select_atoms(lipid_sel)
    protein  = u.select_atoms(protein_sel)
    membrane = u.select_atoms(f"({lps_sel}) or ({lipid_sel})")
    all_atoms = u.select_atoms("all")

    print(f"LPS group    : {len(lps)} atoms  [{lps_sel}]", flush=True)
    print(f"Lipid group  : {len(lipids)} atoms  [{lipid_sel}]", flush=True)
    print(f"Protein/mol  : {len(protein)} atoms  [{protein_sel}]", flush=True)

    workflow = [
        trans.unwrap(all_atoms),
        trans.center_in_box(membrane, wrap=True),
        trans.wrap(all_atoms, compound="residues"),
    ]
    u.trajectory.add_transformations(*workflow)

    return u, lps, lipids, protein


def make_figure(
    lps_label: str = "LPS sugars",
    lipid_label: str = "Lipids",
    protein_label: str = "PMB",
    figsize: tuple = DEFAULT_FIGSIZE,
    bg: str = DEFAULT_BG,
    c_lps: str = DEFAULT_C_LPS,
    c_lipid: str = DEFAULT_C_LIPID,
    c_protein: str = DEFAULT_C_PROTEIN,
):
    """Create and configure the dual-panel matplotlib figure.

    Returns
    -------
    fig : Figure
    ax_top : Axes   (XY top-down view)
    ax_side : Axes  (XZ side view)
    """
    fig, (ax_top, ax_side) = plt.subplots(1, 2, figsize=figsize, facecolor=bg)
    for ax in (ax_top, ax_side):
        ax.set_facecolor(bg)
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_color("white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")

    legend_patches = [
        mpatches.Patch(color=c_lps,     label=lps_label),
        mpatches.Patch(color=c_lipid,   label=lipid_label),
        mpatches.Patch(color=c_protein, label=protein_label),
    ]
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.10, top=0.93, wspace=0.35)
    fig.legend(handles=legend_patches, loc="lower center", ncol=3,
            facecolor=bg, labelcolor="white", framealpha=0.4,
            bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("", color="white")

    return fig, ax_top, ax_side


def draw_frame(
    fig,
    ax_top,
    ax_side,
    ts,
    lps,
    lipids,
    protein,
    title_prefix: str = "",
    bg: str = DEFAULT_BG,
    c_lps: str = DEFAULT_C_LPS,
    c_lipid: str = DEFAULT_C_LIPID,
    c_protein: str = DEFAULT_C_PROTEIN,
) -> np.ndarray:
    """Draw one trajectory timestep onto *fig* and return an RGB image array.

    The figure is reused across calls; both axes are cleared and redrawn.

    Parameters
    ----------
    fig, ax_top, ax_side:
        Objects returned by :func:`make_figure`.
    ts:
        Current MDAnalysis ``Timestep`` (the trajectory must already be
        positioned at this frame).
    lps, lipids, protein:
        Atom groups returned by :func:`build_universe`.
    title_prefix:
        Optional text prepended to the ``t = … ps`` suptitle.

    Returns
    -------
    img : np.ndarray, shape (H, W, 3), dtype uint8
        RGB pixel data suitable for imageio or matplotlib ``imshow``.
    """
    box = ts.dimensions[:3]  # (Lx, Ly, Lz) in Å

    lps_pos  = lps.positions
    lip_pos  = lipids.positions
    pmb_pos  = protein.positions

    # ── top-down view (XY) ────────────────────────────────────────────────────
    ax_top.cla()
    ax_top.set_facecolor(bg)
    ax_top.scatter(lps_pos[:, 0], lps_pos[:, 1], c=c_lps,     s=2,  alpha=0.6, linewidths=0, rasterized=True)
    ax_top.scatter(lip_pos[:, 0], lip_pos[:, 1], c=c_lipid,   s=2,  alpha=0.6, linewidths=0, rasterized=True)
    ax_top.scatter(pmb_pos[:, 0], pmb_pos[:, 1], c=c_protein, s=12, alpha=1.0, linewidths=0, zorder=5)
    ax_top.set_xlim(0, box[0])
    ax_top.set_ylim(0, box[1])
    ax_top.set_aspect("equal")
    ax_top.set_title("Top view (XY)", color="white")
    ax_top.set_xlabel("x (Å)", color="white")
    ax_top.set_ylabel("y (Å)", color="white")
    ax_top.tick_params(colors="white")
    for spine in ax_top.spines.values():
        spine.set_color("white")

    # ── side view (XZ) ────────────────────────────────────────────────────────
    ax_side.cla()
    ax_side.set_facecolor(bg)
    ax_side.scatter(lps_pos[:, 0], lps_pos[:, 2], c=c_lps,     s=2,  alpha=0.6, linewidths=0, rasterized=True)
    ax_side.scatter(lip_pos[:, 0], lip_pos[:, 2], c=c_lipid,   s=2,  alpha=0.6, linewidths=0, rasterized=True)
    ax_side.scatter(pmb_pos[:, 0], pmb_pos[:, 2], c=c_protein, s=12, alpha=1.0, linewidths=0, zorder=5)
    ax_side.set_xlim(0, box[0])
    ax_side.set_ylim(0, box[2])
    ax_side.set_aspect("equal")
    ax_side.set_title("Side view (XZ)", color="white")
    ax_side.set_xlabel("x (Å)", color="white")
    ax_side.set_ylabel("z (Å)", color="white")
    ax_side.tick_params(colors="white")
    for spine in ax_side.spines.values():
        spine.set_color("white")

    _prefix = (title_prefix.rstrip() + "  ") if title_prefix else ""
    fig.texts[0].set_text(f"{_prefix}t = {ts.time:.0f} ps")

    # Draw to the fixed-size canvas — guarantees identical dimensions every frame.
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    # buffer_rgba() → (h, w, 4) RGBA; drop alpha for imageio (RGB only).
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)[:, :, :3]
    return img


def render_single_frame(
    topology: str,
    trajectory: str,
    output: str,
    frame: int = 0,
    lps_sel: str = DEFAULT_LPS_SEL,
    lipid_sel: str = DEFAULT_LIPID_SEL,
    protein_sel: str = DEFAULT_PROTEIN_SEL,
    lps_label: str = "LPS sugars",
    lipid_label: str = "Lipids",
    protein_label: str = "PMB",
    title: str = "",
    figsize: tuple = DEFAULT_FIGSIZE,
    bg: str = DEFAULT_BG,
    c_lps: str = DEFAULT_C_LPS,
    c_lipid: str = DEFAULT_C_LIPID,
    c_protein: str = DEFAULT_C_PROTEIN,
    dpi: int = DEFAULT_DPI,
):
    """Render one trajectory frame and save it as an image file.

    Parameters
    ----------
    topology, trajectory:
        Input files passed to :func:`build_universe`.
    output:
        Output image path.  The format is inferred from the extension
        (e.g. ``.png``, ``.jpg``, ``.pdf``).
    frame:
        Zero-based frame index to render.
    *:
        All remaining keyword arguments are forwarded to
        :func:`build_universe` / :func:`make_figure` / :func:`draw_frame`.
    """
    u, lps, lipids, protein = build_universe(
        topology, trajectory, lps_sel, lipid_sel, protein_sel
    )
    u.trajectory[frame]
    ts = u.trajectory.ts

    fig, ax_top, ax_side = make_figure(
        lps_label, lipid_label, protein_label, figsize, bg, c_lps, c_lipid, c_protein
    )
    img = draw_frame(fig, ax_top, ax_side, ts, lps, lipids, protein, title,
                    bg, c_lps, c_lipid, c_protein)
    imageio.imwrite(output, img)
    plt.close(fig)
    print(f"Frame {frame} (t={ts.time:.0f} ps) saved to {output}", flush=True)


def render_movie(
    topology: str,
    trajectory: str,
    output: str,
    lps_sel: str = DEFAULT_LPS_SEL,
    lipid_sel: str = DEFAULT_LIPID_SEL,
    protein_sel: str = DEFAULT_PROTEIN_SEL,
    lps_label: str = "LPS sugars",
    lipid_label: str = "Lipids",
    protein_label: str = "PMB",
    title: str = "",
    figsize: tuple = DEFAULT_FIGSIZE,
    fps: int = DEFAULT_FPS,
    stride: int = DEFAULT_STRIDE,
    dpi: int = DEFAULT_DPI,
    bg: str = DEFAULT_BG,
    c_lps: str = DEFAULT_C_LPS,
    c_lipid: str = DEFAULT_C_LIPID,
    c_protein: str = DEFAULT_C_PROTEIN,
):
    """Render the full trajectory as an MP4 movie.

    Parameters
    ----------
    topology, trajectory, output:
        Input topology/trajectory and output MP4 path.
    stride:
        Step between frames (1 = every frame).
    fps:
        Frames per second in the output video.
    *:
        All remaining keyword arguments are forwarded to
        :func:`build_universe` / :func:`make_figure` / :func:`draw_frame`.
    """
    n_cpus = _n_cpus_from_env()
    print(f"Using {n_cpus} CPU thread(s) for ffmpeg encoding.", flush=True)

    u, lps, lipids, protein = build_universe(
        topology, trajectory, lps_sel, lipid_sel, protein_sel
    )
    fig, ax_top, ax_side = make_figure(
        lps_label, lipid_label, protein_label, figsize, bg, c_lps, c_lipid, c_protein
    )

    n_frames = len(u.trajectory[::stride])
    print(f"Rendering {n_frames} frames → {output}", flush=True)

    with imageio.get_writer(output, fps=fps, quality=8,
                            output_params=["-threads", str(n_cpus)]) as writer:
        for i, ts in enumerate(u.trajectory[::stride]):
            if i % max(1, n_frames // 20) == 0:
                print(f"  frame {i}/{n_frames}  t={ts.time:.0f} ps", flush=True)
            img = draw_frame(fig, ax_top, ax_side, ts, lps, lipids, protein, title,
                            bg, c_lps, c_lipid, c_protein)
            writer.append_data(img)

    plt.close(fig)
    print(f"Movie saved to {output}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Render a GROMACS trajectory as an MP4 movie or a single image "
            "(no display needed)."
        )
    )
    parser.add_argument("-t", "--topology", required=True,
                        help="Topology file (e.g. .tpr or .gro)")
    parser.add_argument("-x", "--trajectory", required=True,
                        help="Trajectory file (e.g. .xtc or .trr)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file (.mp4 for movie, .png/.jpg for a single frame)")
    parser.add_argument("--frame", type=int, default=None, metavar="N",
                        help="Render only frame N (0-based index) as an image instead of "
                            "the full movie. Ignored when --output ends in .mp4 unless "
                            "this flag is given explicitly.")
    parser.add_argument("--lps-sel",
                        default=DEFAULT_LPS_SEL,
                        help="MDAnalysis selection for the LPS / membrane sugar group "
                            "(default: GM LPS resnames)")
    parser.add_argument("--lipid-sel",
                        default=DEFAULT_LIPID_SEL,
                        help="MDAnalysis selection for the lipid group "
                            "(default: GM inner-membrane lipid resnames)")
    parser.add_argument("--protein-sel",
                        default=DEFAULT_PROTEIN_SEL,
                        help="MDAnalysis selection for the protein / peptide / ligand "
                            "of interest (default: 'resname PMB')")
    parser.add_argument("--lps-label",     default="LPS sugars",
                        help="Legend label for the LPS group (default: 'LPS sugars')")
    parser.add_argument("--lipid-label",   default="Lipids",
                        help="Legend label for the lipid group (default: 'Lipids')")
    parser.add_argument("--protein-label", default="PMB",
                        help="Legend label for the protein/peptide/ligand (default: 'PMB')")
    parser.add_argument("--title", default="",
                        help="Static title prefix shown above every frame, "
                            "followed by the current simulation time (default: none)")
    parser.add_argument("--fps", type=int, default=DEFAULT_FPS,
                        help=f"Frames per second for the output movie (default: {DEFAULT_FPS})")
    parser.add_argument("--stride", type=int, default=DEFAULT_STRIDE,
                        help=f"Step between trajectory frames (default: {DEFAULT_STRIDE})")
    parser.add_argument("--figsize", type=float, nargs=2, default=list(DEFAULT_FIGSIZE),
                        metavar=("W", "H"),
                        help=f"Figure width and height in inches "
                            f"(default: {DEFAULT_FIGSIZE[0]} {DEFAULT_FIGSIZE[1]})")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                        help=f"Resolution in dots per inch (default: {DEFAULT_DPI})")
    args = parser.parse_args()

    common_kw = dict(
        lps_sel=args.lps_sel,
        lipid_sel=args.lipid_sel,
        protein_sel=args.protein_sel,
        lps_label=args.lps_label,
        lipid_label=args.lipid_label,
        protein_label=args.protein_label,
        title=args.title,
        figsize=tuple(args.figsize),
        dpi=args.dpi,
    )

    if args.frame is not None:
        render_single_frame(
            args.topology, args.trajectory, args.output,
            frame=args.frame,
            **common_kw,
        )
    else:
        render_movie(
            args.topology, args.trajectory, args.output,
            fps=args.fps,
            stride=args.stride,
            **common_kw,
        )


if __name__ == "__main__":
    main()