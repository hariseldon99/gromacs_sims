#!/usr/bin/env python3
"""Protonate polymyxin B Dab side-chain amines for CHARMM-GUI.

The PubChem PMB SDF is neutral, with the five free Dab side-chain amines
represented as neutral -NH2 groups. For the usual membrane-binding PMB model,
those five side-chain amines should be protonated to -NH3+ while the lactam
ring nitrogen remains neutral.

This script preserves the original atom order, appends one hydrogen to each
target side-chain nitrogen, and writes an SDF with total formal charge +5.
The target atom IDs are the 1-based SDF atom IDs used in
protein-pmbn-gm-membrane/polymyxin_b_pubchem_residues.pdb REMARK 200 records.
"""

from __future__ import annotations

import argparse
from pathlib import Path

Chem = None


TARGET_SIDE_CHAIN_NITROGENS = {
    29: "DAB 2 ND",
    27: "DAB 4 ND",
    26: "DAB/DDB 6 ND",
    25: "DAB 9 ND",
    28: "DAB 10 ND",
}
LACTAM_NITROGEN = 22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Protonate the five free Dab side-chain amines in polymyxin B "
            "from neutral -NH2 to charged -NH3+."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input PMB SDF with explicit hydrogens, e.g. polymyxin_b_3d.sdf.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output protonated SDF, e.g. polymyxin_b_3d_plus5.sdf.",
    )
    parser.add_argument(
        "--target-sdf-atoms",
        default=",".join(str(i) for i in TARGET_SIDE_CHAIN_NITROGENS),
        help=(
            "Comma-separated 1-based SDF atom IDs to protonate. Default: "
            "29,27,26,25,28."
        ),
    )
    return parser.parse_args()


def load_rdkit() -> None:
    global Chem
    try:
        from rdkit import Chem as rdkit_chem
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "RDKit is required, but it is not importable from this Python "
            "environment. Activate the conda/venv environment that contains "
            "RDKit, or install it with conda-forge:\n"
            "  conda install -c conda-forge rdkit"
        ) from exc
    Chem = rdkit_chem


def parse_target_atoms(raw_targets: str) -> list[int]:
    targets = []
    for value in raw_targets.split(","):
        value = value.strip()
        if not value:
            continue
        targets.append(int(value))
    if not targets:
        raise ValueError("At least one target atom ID is required")
    return targets


def read_sdf(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError(f"RDKit could not parse input SDF: {path}")
    if not mol.GetNumConformers():
        raise ValueError(f"Input SDF has no conformer coordinates: {path}")
    return mol


def explicit_hydrogen_neighbors(atom) -> int:
    return sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)


def protonate_targets(mol: Chem.Mol, target_atom_ids: list[int]) -> Chem.Mol:
    editable = Chem.RWMol(mol)
    zero_based_targets = [atom_id - 1 for atom_id in target_atom_ids]

    for atom_id, atom_idx in zip(target_atom_ids, zero_based_targets):
        if atom_idx < 0 or atom_idx >= editable.GetNumAtoms():
            raise ValueError(f"Target atom ID {atom_id} is outside the molecule")
        atom = editable.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "N":
            raise ValueError(
                f"Target atom ID {atom_id} should be nitrogen, found {atom.GetSymbol()}"
            )
        atom.SetFormalCharge(1)
        atom.SetNoImplicit(False)

    charged = editable.GetMol()
    charged.UpdatePropertyCache(strict=False)
    protonated = Chem.AddHs(charged, addCoords=True, onlyOnAtoms=zero_based_targets)
    Chem.SanitizeMol(protonated)

    for atom_id in target_atom_ids:
        atom = protonated.GetAtomWithIdx(atom_id - 1)
        h_count = explicit_hydrogen_neighbors(atom)
        if atom.GetFormalCharge() != 1 or h_count != 3:
            raise ValueError(
                f"Expected atom {atom_id} to be -NH3+, but formal charge is "
                f"{atom.GetFormalCharge()} and explicit H count is {h_count}"
            )

    lactam = protonated.GetAtomWithIdx(LACTAM_NITROGEN - 1)
    if lactam.GetFormalCharge() != 0:
        raise ValueError(f"Lactam nitrogen atom {LACTAM_NITROGEN} was charged unexpectedly")

    total_charge = Chem.GetFormalCharge(protonated)
    if total_charge != len(target_atom_ids):
        raise ValueError(
            f"Expected total formal charge +{len(target_atom_ids)}, got {total_charge}"
        )

    return protonated


def write_sdf(mol: Chem.Mol, output_path: Path, target_atom_ids: list[int]) -> None:
    mol.SetProp("PMB_PROTONATED_DAB_SIDECHAIN_SDF_ATOMS", ",".join(map(str, target_atom_ids)))
    mol.SetProp("PMB_LACTAM_NITROGEN_LEFT_NEUTRAL_SDF_ATOM", str(LACTAM_NITROGEN))
    mol.SetProp("PMB_EXPECTED_TOTAL_FORMAL_CHARGE", str(Chem.GetFormalCharge(mol)))
    if mol.HasProp("PUBCHEM_TOTAL_CHARGE"):
        mol.SetProp("PUBCHEM_TOTAL_CHARGE", str(Chem.GetFormalCharge(mol)))

    writer = Chem.SDWriter(str(output_path))
    try:
        writer.write(mol)
    finally:
        writer.close()


def main() -> None:
    args = parse_args()
    load_rdkit()
    target_atom_ids = parse_target_atoms(args.target_sdf_atoms)
    mol = read_sdf(args.input)
    protonated = protonate_targets(mol, target_atom_ids)
    write_sdf(protonated, args.output, target_atom_ids)

    print(f"Wrote: {args.output}")
    print(f"Input atom count: {mol.GetNumAtoms()}")
    print(f"Output atom count: {protonated.GetNumAtoms()}")
    print(f"Total formal charge: {Chem.GetFormalCharge(protonated):+d}")
    for atom_id in target_atom_ids:
        label = TARGET_SIDE_CHAIN_NITROGENS.get(atom_id, "custom target")
        atom = protonated.GetAtomWithIdx(atom_id - 1)
        print(
            f"  SDF atom {atom_id:3d} {label:12s}: "
            f"formal charge {atom.GetFormalCharge():+d}, "
            f"explicit H {explicit_hydrogen_neighbors(atom)}"
        )


if __name__ == "__main__":
    main()
