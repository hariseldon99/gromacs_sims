#!/usr/bin/env python3
"""Generate 3D conformers from a 2D molecular file with RDKit.

This script is intended for small molecules and nonstandard peptide-like
ligands such as polymyxin B. It reads an SDF/MOL/MOL2/PDB input structure,
embeds multiple 3D conformers with RDKit ETKDG, optimizes them with MMFF94s
when possible and UFF otherwise, then writes the lowest-energy conformer as
SDF and PDB.

Parallelization is handled through RDKit's native worker threads. Use
``--num-threads`` to control the number of threads used during conformer
embedding and force-field optimization; ``--num-threads 0`` tells RDKit to use
all CPU cores visible to the process.

For molecules that already have a residue-labelled PDB mapping, use
``--template-pdb``. The template PDB must contain the same atom ordering as the
input molecule, or ``REMARK 200`` lines of the form produced by
``polymyxin_b_pubchem_residues.pdb``:

    REMARK 200 PDB_ATOM     1 SDF_ATOM    73 C    MOA    1

When a template is supplied, output PDB atom names, residue names, residue
numbers, and CONECT records are written by this script instead of RDKit's
generic PDB writer.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

Chem = None
AllChem = None


@dataclass(frozen=True)
class TemplateAtom:
    """PDB identity fields to reuse for one atom."""

    atom_name: str
    res_name: str
    chain_id: str
    res_seq: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a 3D conformer with RDKit ETKDG and write optimized SDF "
            "and PDB outputs."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input molecule file. Supported extensions: .sdf, .mol, .mol2, .pdb",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        help=(
            "Output prefix. Defaults to the input path stem beside the input "
            "file, with '_3d' appended."
        ),
    )
    parser.add_argument(
        "--template-pdb",
        type=Path,
        help=(
            "Optional residue-labelled PDB template used to preserve atom names, "
            "residue names, residue numbers, chain ID, and CONECT records."
        ),
    )
    parser.add_argument(
        "--num-conformers",
        type=int,
        default=200,
        help="Number of conformers to embed before optimization. Default: 200",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=0,
        help=(
            "Number of RDKit worker threads for embedding/optimization. "
            "Use 0 for all available cores. Default: 0"
        ),
    )
    parser.add_argument(
        "--max-optimization-iters",
        type=int,
        default=2000,
        help="Maximum MMFF/UFF optimization iterations per conformer. Default: 2000",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260514,
        help="Random seed for reproducible embedding. Default: 20260514",
    )
    parser.add_argument(
        "--prune-rms-thresh",
        type=float,
        default=0.5,
        help="ETKDG RMS pruning threshold in angstrom. Default: 0.5",
    )
    parser.add_argument(
        "--keep-existing-hydrogens",
        action="store_true",
        help=(
            "Keep explicit hydrogens present in the input. By default RDKit "
            "removes and regenerates hydrogens before embedding."
        ),
    )
    parser.add_argument(
        "--no-sdf",
        action="store_true",
        help="Only write PDB output; skip the optimized SDF file.",
    )
    parser.add_argument(
        "--pdb-resname",
        default="LIG",
        help="Residue name to use when no template PDB is supplied. Default: LIG",
    )
    parser.add_argument(
        "--chain-id",
        default="A",
        help="Chain ID to use when no template PDB is supplied. Default: A",
    )
    return parser.parse_args()


def load_rdkit() -> None:
    """Import RDKit after argparse so --help works without RDKit installed."""

    global Chem, AllChem
    try:
        from rdkit import Chem as rdkit_chem
        from rdkit.Chem import AllChem as rdkit_all_chem
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "RDKit is required to generate conformers, but it is not importable "
            "from this Python environment. Activate the conda/venv environment "
            "that contains RDKit, or install it with conda-forge:\n"
            "  conda install -c conda-forge rdkit"
        ) from exc
    Chem = rdkit_chem
    AllChem = rdkit_all_chem


def read_molecule(input_path: Path, keep_existing_hydrogens: bool) -> Chem.Mol:
    suffix = input_path.suffix.lower()
    remove_hs = not keep_existing_hydrogens

    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(input_path), removeHs=remove_hs, sanitize=True)
        mol = next((m for m in supplier if m is not None), None)
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(input_path), removeHs=remove_hs, sanitize=True)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(input_path), removeHs=remove_hs, sanitize=True)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(input_path), removeHs=remove_hs, sanitize=True)
    else:
        raise ValueError(f"Unsupported input extension: {input_path.suffix}")

    if mol is None:
        raise ValueError(f"RDKit could not parse input molecule: {input_path}")

    return Chem.AddHs(mol, addCoords=True)


def optimize_conformers(
    mol: Chem.Mol,
    num_conformers: int,
    num_threads: int,
    max_optimization_iters: int,
    seed: int,
    prune_rms_thresh: float,
) -> tuple[Chem.Mol, int, str, float]:
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = prune_rms_thresh
    params.useRandomCoords = True
    # RDKit does the parallel work internally. This avoids pickling molecule
    # objects through Python multiprocessing and works cleanly on HPC nodes.
    if hasattr(params, "numThreads"):
        params.numThreads = num_threads

    conformer_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params))
    if not conformer_ids:
        raise RuntimeError("RDKit failed to embed any 3D conformers")

    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if mmff_props is not None and AllChem.MMFFHasAllMoleculeParams(mol):
        force_field = "MMFF94s"
        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            numThreads=num_threads,
            maxIters=max_optimization_iters,
            mmffVariant="MMFF94s",
        )
    else:
        force_field = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(
            mol,
            numThreads=num_threads,
            maxIters=max_optimization_iters,
        )

    best_conf_id = None
    best_energy = None
    for conformer_id, (_not_converged, energy) in zip(conformer_ids, results):
        if best_energy is None or energy < best_energy:
            best_conf_id = conformer_id
            best_energy = energy

    if best_conf_id is None or best_energy is None:
        raise RuntimeError("Optimization finished without an energy result")

    return mol, best_conf_id, force_field, best_energy


def template_atoms_from_pdb(template_pdb: Path, input_atom_count: int) -> list[TemplateAtom]:
    pdb_lines = template_pdb.read_text().splitlines()
    atom_records = [line for line in pdb_lines if line.startswith(("ATOM  ", "HETATM"))]
    if len(atom_records) != input_atom_count:
        raise ValueError(
            f"Template atom count ({len(atom_records)}) does not match embedded "
            f"molecule atom count ({input_atom_count})."
        )

    ordered_records = atom_records
    remark_map: dict[int, int] = {}
    for line in pdb_lines:
        if not line.startswith("REMARK 200"):
            continue
        parts = line.split()
        if len(parts) >= 7 and parts[2] == "PDB_ATOM" and parts[4] == "SDF_ATOM":
            pdb_serial = int(parts[3])
            sdf_atom = int(parts[5])
            remark_map[sdf_atom] = pdb_serial

    if remark_map:
        ordered_records = []
        for sdf_atom in range(1, input_atom_count + 1):
            pdb_serial = remark_map.get(sdf_atom)
            if pdb_serial is None:
                raise ValueError(f"Template is missing REMARK 200 mapping for SDF atom {sdf_atom}")
            try:
                ordered_records.append(atom_records[pdb_serial - 1])
            except IndexError as exc:
                raise ValueError(f"Template PDB serial {pdb_serial} is out of range") from exc

    atoms = []
    for line in ordered_records:
        atoms.append(
            TemplateAtom(
                atom_name=line[12:16].strip(),
                res_name=line[17:20].strip() or "LIG",
                chain_id=(line[21].strip() or "A"),
                res_seq=int(line[22:26]),
            )
        )
    return atoms


def pdb_atom_name(name: str, element: str) -> str:
    if len(name) >= 4:
        return name[:4]
    if len(element) == 1:
        return f" {name:<3}"[:4]
    return f"{name:<4}"[:4]


def write_template_pdb(
    mol: Chem.Mol,
    conf_id: int,
    output_pdb: Path,
    template_atoms: list[TemplateAtom] | None,
    force_field: str,
    energy: float,
    fallback_resname: str,
    fallback_chain_id: str,
) -> None:
    conf = mol.GetConformer(conf_id)
    records = [
        "REMARK   1 Generated by scripts/generate_3d_conformer.py",
        f"REMARK   2 RDKit ETKDG conformer optimized with {force_field}; energy {energy:.6f}",
    ]

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        serial = idx + 1
        element = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        if template_atoms is None:
            atom_name = f"{element}{serial}"
            res_name = fallback_resname[:3].upper()
            chain_id = fallback_chain_id[:1] or "A"
            res_seq = 1
        else:
            tmpl = template_atoms[idx]
            atom_name = tmpl.atom_name
            res_name = tmpl.res_name[:3].upper()
            chain_id = tmpl.chain_id[:1] or "A"
            res_seq = tmpl.res_seq

        records.append(
            f"HETATM{serial:5d} {pdb_atom_name(atom_name, element):4s} {res_name:>3s} "
            f"{chain_id}{res_seq:4d}    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
        )

    records.append("TER")
    records.extend(conect_records(mol))
    records.append("END")
    output_pdb.write_text("\n".join(records) + "\n")


def conect_records(mol: Chem.Mol) -> Iterable[str]:
    connections: dict[int, list[int]] = {}
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx() + 1
        end = bond.GetEndAtomIdx() + 1
        order = max(1, int(round(bond.GetBondTypeAsDouble())))
        connections.setdefault(begin, []).extend([end] * order)
        connections.setdefault(end, []).extend([begin] * order)

    for serial in sorted(connections):
        targets = connections[serial]
        for offset in range(0, len(targets), 4):
            chunk = targets[offset : offset + 4]
            yield "CONECT" + f"{serial:5d}" + "".join(f"{target:5d}" for target in chunk)


def write_sdf(mol: Chem.Mol, conf_id: int, output_sdf: Path, force_field: str, energy: float) -> None:
    mol.SetProp("RDKit_3D_force_field", force_field)
    mol.SetProp("RDKit_3D_energy", f"{energy:.6f}")
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(mol, confId=conf_id)
    finally:
        writer.close()


def main() -> None:
    args = parse_args()
    load_rdkit()
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = args.input.with_name(f"{args.input.stem}_3d")

    mol = read_molecule(args.input, args.keep_existing_hydrogens)
    mol, best_conf_id, force_field, energy = optimize_conformers(
        mol,
        args.num_conformers,
        args.num_threads,
        args.max_optimization_iters,
        args.seed,
        args.prune_rms_thresh,
    )

    output_pdb = output_prefix.with_suffix(".pdb")
    output_sdf = output_prefix.with_suffix(".sdf")
    template_atoms = None
    if args.template_pdb is not None:
        template_atoms = template_atoms_from_pdb(args.template_pdb, mol.GetNumAtoms())

    write_template_pdb(
        mol,
        best_conf_id,
        output_pdb,
        template_atoms,
        force_field,
        energy,
        args.pdb_resname,
        args.chain_id,
    )
    if not args.no_sdf:
        write_sdf(mol, best_conf_id, output_sdf, force_field, energy)

    print(f"Input: {args.input}")
    print(f"Embedded conformers: {args.num_conformers}")
    print(f"RDKit worker threads: {args.num_threads} (0 means all visible cores)")
    print(f"Selected conformer ID: {best_conf_id}")
    print(f"Optimization force field: {force_field}")
    print(f"Best energy: {energy:.6f}")
    print(f"Wrote PDB: {output_pdb}")
    if not args.no_sdf:
        print(f"Wrote SDF: {output_sdf}")


if __name__ == "__main__":
    main()
