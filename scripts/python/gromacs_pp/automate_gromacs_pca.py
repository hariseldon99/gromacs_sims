#!/usr/bin/env python3
"""
Automate the GROMACS PCA/FEL workflow used in the 1LML-PDL post-processing notebook.

The script:
  1. Builds a non-interactive GROMACS index with MDAnalysis.
  2. Optionally concatenates trajectories with gmx trjcat.
  3. Runs gmx covar to compute eigenvalues/eigenvectors.
  4. Runs gmx anaeig for requested principal components.
  5. Combines requested PC pairs, runs gmx sham, and converts FES XPM to DAT.
  6. Plots eigenvalues, PC projections, and free energy landscapes.
"""
from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
from pathlib import Path

import numpy as np


LOG = logging.getLogger("gromacs-pca")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run GROMACS PCA, PC projections, and free energy landscape generation."
    )
    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=Path,
        help="Input structure/topology for GROMACS and MDAnalysis, usually .gro, .pdb, or .tpr.",
    )
    parser.add_argument(
        "-f",
        "--trajectory",
        required=True,
        nargs="+",
        type=Path,
        help="Input trajectory file(s), for example one or more .xtc files.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=Path("pca"),
        type=Path,
        help="Output directory. Default: pca",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix. Default: structure stem.",
    )
    parser.add_argument(
        "--gmx",
        default="gmx",
        help="GROMACS executable. Default: gmx",
    )
    parser.add_argument(
        "--selection",
        default="backbone or resname LIG",
        help=(
            "MDAnalysis atom selection for PCA. "
            "Default: 'backbone or resname LIG'."
        ),
    )
    parser.add_argument(
        "--fit-selection",
        default=None,
        help=(
            "Optional MDAnalysis selection for least-squares fitting in gmx covar. "
            "Default: use --selection."
        ),
    )
    parser.add_argument(
        "--group-name",
        default="PCA_SELECTION",
        help="Name of the PCA index group written to the .ndx file.",
    )
    parser.add_argument(
        "--fit-group-name",
        default="PCA_FIT",
        help="Name of the optional fit index group when --fit-selection is used.",
    )
    parser.add_argument(
        "--pcs",
        default="1,2",
        help="Comma/range list of PCs to project, e.g. '1,2,3' or '1-3'. Default: 1,2",
    )
    parser.add_argument(
        "--pairs",
        default="1,2",
        help=(
            "Comma-separated PC pairs for FEL, e.g. '1,2' or '1,2;2,3;3,1'. "
            "Default: 1,2"
        ),
    )
    parser.add_argument(
        "--sham-extra",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra arguments passed to gmx sham, e.g. --sham-extra -ngrid 50 50",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep .ndx, eigenvectors.trr, combined.xtc, FES.xpm, and GROMACS logs.",
    )
    parser.add_argument(
        "--keep-xpm",
        action="store_true",
        help="Keep gmx sham FES .xpm files even when other intermediates are cleaned.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write data files only; skip PNG plots.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def parse_pc_list(text: str) -> list[int]:
    pcs: set[int] = set()
    for token in text.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            start, end = int(left), int(right)
            if start > end:
                raise ValueError(f"Invalid descending PC range: {token}")
            pcs.update(range(start, end + 1))
        else:
            pcs.add(int(token))
    if not pcs or min(pcs) < 1:
        raise ValueError("PC numbers must be positive integers.")
    return sorted(pcs)


def parse_pairs(text: str) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for token in text.split(";"):
        token = token.strip()
        if not token:
            continue
        parts = [part.strip() for part in token.split(",")]
        if len(parts) != 2:
            raise ValueError(f"Invalid PC pair '{token}'. Use forms like '1,2;2,3'.")
        pair = (int(parts[0]), int(parts[1]))
        if pair[0] < 1 or pair[1] < 1:
            raise ValueError("PC pair numbers must be positive integers.")
        pairs.append(pair)
    if not pairs:
        raise ValueError("At least one PC pair is required.")
    return pairs


def run_command(
    cmd: list[str],
    cwd: Path,
    *,
    stdin: str | None = None,
    log_file: Path | None = None,
) -> None:
    LOG.info("Running: %s", " ".join(cmd))
    completed = subprocess.run(
        cmd,
        cwd=cwd,
        input=stdin,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if log_file:
        log_file.write_text(completed.stdout or "")
    if completed.returncode != 0:
        if completed.stdout:
            LOG.error(completed.stdout.strip())
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {' '.join(cmd)}")
    if completed.stdout:
        LOG.debug(completed.stdout.strip())


def gmx_command(gmx: str, tool: str) -> list[str]:
    return [*shlex.split(gmx), "-nobackup", tool]


def write_ndx_group(handle, name: str, indices_zero_based: np.ndarray) -> None:
    handle.write(f"[ {name} ]\n")
    indices = [str(int(index) + 1) for index in indices_zero_based]
    for start in range(0, len(indices), 15):
        handle.write(" ".join(indices[start : start + 15]) + "\n")
    handle.write("\n")


def create_index(
    structure: Path,
    output_ndx: Path,
    selection: str,
    group_name: str,
    fit_selection: str | None,
    fit_group_name: str,
) -> tuple[str, str]:
    try:
        import MDAnalysis as mda
    except ImportError as exc:  # pragma: no cover - exercised by user environment
        raise SystemExit(
            "MDAnalysis is required for non-interactive index creation. "
            "Install it in this environment, then rerun."
        ) from exc

    LOG.info("Creating PCA index with MDAnalysis selection: %s", selection)
    universe = mda.Universe(str(structure))
    pca_atoms = universe.select_atoms(selection)
    if len(pca_atoms) == 0:
        raise ValueError(f"PCA selection matched no atoms: {selection}")

    fit_group = group_name
    with output_ndx.open("w") as handle:
        if fit_selection:
            LOG.info("Creating fit index with MDAnalysis selection: %s", fit_selection)
            fit_atoms = universe.select_atoms(fit_selection)
            if len(fit_atoms) == 0:
                raise ValueError(f"Fit selection matched no atoms: {fit_selection}")
            write_ndx_group(handle, fit_group_name, fit_atoms.indices)
            fit_group = fit_group_name
        write_ndx_group(handle, group_name, pca_atoms.indices)

    LOG.info("Wrote %s with %d PCA atoms", output_ndx, len(pca_atoms))
    return fit_group, group_name


def read_xvg_xy(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@")) or stripped == "&":
                continue
            rows.append([float(part) for part in stripped.split()])
    if not rows:
        raise ValueError(f"No numeric data found in {path}")
    return np.array(rows, dtype=float)


def combine_pc_pair(pc_a: Path, pc_b: Path, output: Path) -> None:
    data_a = read_xvg_xy(pc_a)
    data_b = read_xvg_xy(pc_b)
    if len(data_a) != len(data_b):
        raise ValueError(f"PC files have different lengths: {pc_a} and {pc_b}")
    if not np.allclose(data_a[:, 0], data_b[:, 0]):
        raise ValueError(f"PC files have different time axes: {pc_a} and {pc_b}")
    combined = np.column_stack([data_a[:, 0], data_a[:, 1], data_b[:, 1]])
    np.savetxt(output, combined, fmt="%12.5f")
    LOG.info("Wrote PC pair file: %s", output)


def xpm_to_dat(xpm_file: Path, dat_file: Path) -> None:
    x_axis: list[float] = []
    y_axis: list[float] = []
    xpm_rows: list[str] = []
    symbol_to_value: dict[str, float] = {}

    with xpm_file.open() as handle:
        for line in handle:
            if line.startswith("/* x-axis"):
                x_axis = [float(value) for value in line.split()[2:-2]]
                continue
            if line.startswith("/* y-axis"):
                y_axis = [float(value) for value in line.split()[2:-2]]
                continue
            if line.startswith('"') and len(line.split()) > 4:
                parts = line.split()
                symbol = parts[0][1:]
                try:
                    symbol_to_value[symbol] = float(parts[-2][1:-1])
                except ValueError:
                    continue
            if line.startswith('"') and x_axis and y_axis:
                xpm_rows.insert(0, line.strip().strip(",")[1:-1])

    if not x_axis or not y_axis or not xpm_rows or not symbol_to_value:
        raise ValueError(f"Could not parse GROMACS XPM data from {xpm_file}")

    values: list[tuple[float, float, float]] = []
    for y_index, row in enumerate(xpm_rows):
        for x_index, symbol in enumerate(row):
            values.append((x_axis[x_index], y_axis[y_index], symbol_to_value[symbol]))

    with dat_file.open("w") as handle:
        for x_value, y_value, z_value in values:
            handle.write(f"{x_value:12.5f}\t{y_value:12.5f}\t{z_value:12.5f}\n")
    LOG.info("Converted %s to %s", xpm_file.name, dat_file.name)


def plot_eigenvalues(eigenvalues_xvg: Path, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    data = read_xvg_xy(eigenvalues_xvg)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data[:, 0], data[:, 1], marker="o", linewidth=1.5, markersize=3.5)
    ax.set_xlabel("Eigenvector index")
    ax.set_ylabel("Eigenvalue (nm$^2$)")
    ax.set_title("PCA Eigenvalues")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)
    LOG.info("Wrote plot: %s", output_png)


def plot_pc_projection(pc_xvg: Path, pc_number: int, output_png: Path) -> None:
    import matplotlib.pyplot as plt

    data = read_xvg_xy(pc_xvg)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data[:, 0], data[:, 1], linewidth=1.4)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(f"PC{pc_number} projection (nm)")
    ax.set_title(f"PC{pc_number} Projected Trajectory")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)
    LOG.info("Wrote plot: %s", output_png)


def plot_fel(dat_file: Path, pair: tuple[int, int], output_png: Path) -> None:
    import matplotlib.pyplot as plt

    data = np.loadtxt(dat_file)
    x_values = np.unique(data[:, 0])
    y_values = np.unique(data[:, 1])
    x_lookup = {value: index for index, value in enumerate(x_values)}
    y_lookup = {value: index for index, value in enumerate(y_values)}
    z = np.full((len(y_values), len(x_values)), np.nan)
    for x_value, y_value, energy in data:
        z[y_lookup[y_value], x_lookup[x_value]] = energy / 4.184
    z = np.ma.masked_invalid(z)

    fig, ax = plt.subplots(figsize=(7, 6))
    contour = ax.contourf(x_values, y_values, z, levels=30, cmap="viridis")
    colorbar = fig.colorbar(contour, ax=ax)
    colorbar.set_label("Free energy (kcal/mol)")
    ax.set_xlabel(f"PC{pair[0]}")
    ax.set_ylabel(f"PC{pair[1]}")
    ax.set_title(f"Free Energy Landscape PC{pair[0]}-PC{pair[1]}")
    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)
    LOG.info("Wrote plot: %s", output_png)


def remove_if_exists(path: Path) -> None:
    try:
        path.unlink()
        LOG.debug("Removed intermediate: %s", path)
    except FileNotFoundError:
        pass


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    structure = args.structure.resolve()
    trajectories = [path.resolve() for path in args.trajectory]
    for path in [structure, *trajectories]:
        if not path.exists():
            raise FileNotFoundError(path)

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or structure.stem

    pcs = parse_pc_list(args.pcs)
    pairs = parse_pairs(args.pairs)
    required_pcs = sorted(set(pcs).union(pc for pair in pairs for pc in pair))
    if required_pcs != pcs:
        LOG.info("Adding PCs needed by --pairs: %s", ",".join(map(str, required_pcs)))
        pcs = required_pcs

    index_file = outdir / f"{prefix}_pca.ndx"
    eigenvalues = outdir / f"{prefix}_eigenvalues.xvg"
    eigenvectors = outdir / f"{prefix}_eigenvectors.trr"
    average_pdb = outdir / f"{prefix}_average.pdb"
    covar_log = outdir / f"{prefix}_covar.log"
    combined_traj = outdir / f"{prefix}_combined.xtc"
    transient: list[Path] = [index_file, eigenvectors, average_pdb, covar_log]

    fit_group, pca_group = create_index(
        structure,
        index_file,
        args.selection,
        args.group_name,
        args.fit_selection,
        args.fit_group_name,
    )

    if len(trajectories) == 1:
        trajectory_for_gmx = trajectories[0]
    else:
        trajectory_for_gmx = combined_traj
        transient.append(combined_traj)
        run_command(
            [*gmx_command(args.gmx, "trjcat"), "-f", *map(str, trajectories), "-o", str(combined_traj)],
            cwd=outdir,
            log_file=outdir / f"{prefix}_trjcat.log",
        )
        transient.append(outdir / f"{prefix}_trjcat.log")

    run_command(
        [
            *gmx_command(args.gmx, "covar"),
            "-s",
            str(structure),
            "-f",
            str(trajectory_for_gmx),
            "-o",
            str(eigenvalues),
            "-v",
            str(eigenvectors),
            "-av",
            str(average_pdb),
            "-n",
            str(index_file),
        ],
        cwd=outdir,
        stdin=f"{fit_group}\n{pca_group}\n",
        log_file=covar_log,
    )

    pc_files: dict[int, Path] = {}
    for pc in pcs:
        pc_file = outdir / f"{prefix}_pc{pc}.xvg"
        pc_log = outdir / f"{prefix}_anaeig_pc{pc}.log"
        run_command(
            [
                *gmx_command(args.gmx, "anaeig"),
                "-f",
                str(trajectory_for_gmx),
                "-s",
                str(structure),
                "-v",
                str(eigenvectors),
                "-first",
                str(pc),
                "-last",
                str(pc),
                "-proj",
                str(pc_file),
                "-n",
                str(index_file),
            ],
            cwd=outdir,
            stdin=f"{pca_group}\n",
            log_file=pc_log,
        )
        pc_files[pc] = pc_file
        transient.append(pc_log)

    for pair in pairs:
        pair_label = f"{pair[0]}{pair[1]}"
        pair_file = outdir / f"{prefix}_pc{pair[0]}pc{pair[1]}.xvg"
        fes_xpm = outdir / f"{prefix}_FES{pair_label}.xpm"
        fel_dat = outdir / f"{prefix}_FEL{pair_label}.dat"
        sham_log = outdir / f"{prefix}_sham{pair_label}.log"
        sham_prob = outdir / f"{prefix}_prob{pair_label}.xpm"
        sham_enthalpy = outdir / f"{prefix}_enthalpy{pair_label}.xpm"
        sham_entropy = outdir / f"{prefix}_entropy{pair_label}.xpm"
        sham_dist = outdir / f"{prefix}_ener{pair_label}.xvg"
        combine_pc_pair(pc_files[pair[0]], pc_files[pair[1]], pair_file)
        run_command(
            [
                *gmx_command(args.gmx, "sham"),
                "-f",
                str(pair_file),
                "-ls",
                str(fes_xpm),
                "-lp",
                str(sham_prob),
                "-lsh",
                str(sham_enthalpy),
                "-lss",
                str(sham_entropy),
                "-dist",
                str(sham_dist),
                *args.sham_extra,
            ],
            cwd=outdir,
            log_file=sham_log,
        )
        xpm_to_dat(fes_xpm, fel_dat)
        transient.extend([sham_log, sham_prob, sham_enthalpy, sham_entropy, sham_dist])
        if not args.keep_xpm:
            transient.append(fes_xpm)

    if not args.no_plots:
        plot_eigenvalues(eigenvalues, outdir / f"{prefix}_eigenvalues.png")
        for pc, pc_file in pc_files.items():
            plot_pc_projection(pc_file, pc, outdir / f"{prefix}_pc{pc}.png")
        for pair in pairs:
            pair_label = f"{pair[0]}{pair[1]}"
            plot_fel(outdir / f"{prefix}_FEL{pair_label}.dat", pair, outdir / f"{prefix}_FEL{pair_label}.png")

    if not args.keep_intermediates:
        for path in transient:
            remove_if_exists(path)

    LOG.info("Done. Results written to %s", outdir)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        LOG.exception("PCA automation failed: %s", exc)
        raise SystemExit(1) from exc
