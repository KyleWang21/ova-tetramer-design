#!/usr/bin/env python3
"""Complete an AF3 homodimer into geometrically exact D2 tetramers.

The validated A/B interface is kept rigid.  Candidate C/D chains are generated
by 180-degree rotations around axes perpendicular to the fitted A/B C2 axis.
This is a geometry generator only; sequence design and AF3 remain independent.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from scipy.spatial import cKDTree


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sx = source - source.mean(0)
    tx = target - target.mean(0)
    u, _, vt = np.linalg.svd(sx.T @ tx)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = vt.T @ u.T
    translation = target.mean(0) - rotation @ source.mean(0)
    return rotation, translation


def chain_arrays(chain):  # noqa: ANN001
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    ca = np.stack([np.asarray(r["CA"].coord, dtype=float) for r in residues])
    atoms = []
    atom_residue = []
    for residue in residues:
        for atom in residue:
            if atom.element != "H" and not atom.name.startswith("H"):
                atoms.append(np.asarray(atom.coord, dtype=float))
                atom_residue.append(int(residue.id[1]))
    return residues, ca, np.stack(atoms), np.asarray(atom_residue, dtype=int)


def rotation_about(axis: np.ndarray) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    return 2.0 * np.outer(axis, axis) - np.eye(3)


def transformed(coords: np.ndarray, center: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    return (rotation @ (coords - center).T).T + center


def pdb_atom_line(serial: int, atom, residue, chain_id: str, xyz: np.ndarray) -> str:  # noqa: ANN001
    name = atom.name[:4].center(4)
    altloc = atom.altloc if atom.altloc not in {"", "?"} else " "
    resname = residue.resname[:3]
    icode = residue.id[2] if residue.id[2] not in {"", "?"} else " "
    occupancy = atom.occupancy if atom.occupancy is not None else 1.0
    bfactor = atom.bfactor if atom.bfactor is not None else 0.0
    element = (atom.element or atom.name[0]).strip()[:2].rjust(2)
    return (
        f"ATOM  {serial:5d} {name}{altloc}{resname:>3s} {chain_id}{int(residue.id[1]):4d}{icode}   "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element}\n"
    )


def write_pose(path: pathlib.Path, chains, center: np.ndarray, rotation: np.ndarray) -> None:  # noqa: ANN001
    serial = 1
    with path.open("w") as handle:
        for source_chain, output_id, apply_rotation in (
            (chains[0], "A", False), (chains[1], "B", False),
            (chains[0], "C", True), (chains[1], "D", True),
        ):
            for residue in source_chain:
                if residue.id[0] != " ":
                    continue
                for atom in residue:
                    xyz = np.asarray(atom.coord, dtype=float)
                    if apply_rotation:
                        xyz = transformed(xyz[None], center, rotation)[0]
                    handle.write(pdb_atom_line(serial, atom, residue, output_id, xyz))
                    serial += 1
            handle.write("TER\n")
        handle.write("END\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dimer", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--samples", type=int, default=360)
    ap.add_argument("--offset-min", type=float, default=-80.0)
    ap.add_argument("--offset-max", type=float, default=80.0)
    ap.add_argument("--offset-steps", type=int, default=41)
    ap.add_argument("--heavy-evaluate", type=int, default=300)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--clash-cutoff", type=float, default=2.0)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True) if args.dimer.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    model = next(parser.get_structure("dimer", str(args.dimer)).get_models())
    chains = [chain for chain in model if sum(1 for r in chain if r.id[0] == " ") >= 380]
    if len(chains) != 2:
        raise SystemExit(f"expected exactly two full protein chains, found {len(chains)}")
    r_a, ca_a, atoms_a, atom_res_a = chain_arrays(chains[0])
    r_b, ca_b, atoms_b, atom_res_b = chain_arrays(chains[1])
    if len(ca_a) != len(ca_b):
        raise SystemExit("chain lengths differ")

    fitted, _ = kabsch(ca_a, ca_b)
    angle = math.degrees(math.acos(float(np.clip((np.trace(fitted) - 1.0) / 2.0, -1.0, 1.0))))
    values, vectors = np.linalg.eig(fitted)
    z_axis = np.real(vectors[:, np.argmin(np.abs(values - 1.0))])
    z_axis /= np.linalg.norm(z_axis)
    center = 0.5 * (ca_a.mean(0) + ca_b.mean(0))
    probe = np.array([1.0, 0.0, 0.0])
    if abs(float(probe @ z_axis)) > 0.8:
        probe = np.array([0.0, 1.0, 0.0])
    u_axis = probe - (probe @ z_axis) * z_axis
    u_axis /= np.linalg.norm(u_axis)
    v_axis = np.cross(z_axis, u_axis)

    ab_atoms = np.concatenate([atoms_a, atoms_b])
    ab_atom_labels = np.concatenate([
        np.stack([np.zeros(len(atoms_a), dtype=int), atom_res_a], axis=1),
        np.stack([np.ones(len(atoms_b), dtype=int), atom_res_b], axis=1),
    ])
    coarse_rows = []
    ab_ca = np.concatenate([ca_a, ca_b])
    tree_ab_ca = cKDTree(ab_ca)
    offsets = np.linspace(args.offset_min, args.offset_max, args.offset_steps)
    for index in range(args.samples):
        phi = 180.0 * index / args.samples
        axis = math.cos(math.radians(phi)) * u_axis + math.sin(math.radians(phi)) * v_axis
        rotation = rotation_about(axis)
        for axial_offset in offsets:
            pose_center = center + axial_offset * z_axis
            cd_ca = transformed(ab_ca, pose_center, rotation)
            tree_cd_ca = cKDTree(cd_ca)
            ca_clashes = sum(map(len, tree_ab_ca.query_ball_tree(tree_cd_ca, 3.5)))
            ca_contacts = sum(map(len, tree_ab_ca.query_ball_tree(tree_cd_ca, 8.0)))
            coarse_rows.append({
                "index": index,
                "phi_deg": phi,
                "axial_offset_a": float(axial_offset),
                "fitted_dimer_rotation_deg": angle,
                "ca_clash_pairs_lt3_5": ca_clashes,
                "ca_pairs_lt8": ca_contacts,
                "rotation": rotation,
                "center": pose_center,
            })
    coarse_ranked = sorted(
        coarse_rows,
        key=lambda row: (
            row["ca_clash_pairs_lt3_5"] > 0,
            row["ca_clash_pairs_lt3_5"],
            -row["ca_pairs_lt8"],
        ),
    )[: args.heavy_evaluate]
    rows = []
    for coarse in coarse_ranked:
        rotation = coarse["rotation"]
        pose_center = coarse["center"]
        cd_atoms = transformed(ab_atoms, pose_center, rotation)
        tree_ab = cKDTree(ab_atoms)
        tree_cd = cKDTree(cd_atoms)
        clash_lists = tree_ab.query_ball_tree(tree_cd, args.clash_cutoff)
        clash_atoms = sum(map(len, clash_lists))
        contact_lists = tree_ab.query_ball_tree(tree_cd, args.contact_cutoff)
        residue_pairs = set()
        chain_pairs = set()
        for left_atom, partners in enumerate(contact_lists):
            left_chain, left_res = ab_atom_labels[left_atom]
            for right_atom in partners:
                right_chain, right_res = ab_atom_labels[right_atom]
                residue_pairs.add((int(left_chain), int(left_res), int(right_chain) + 2, int(right_res)))
                chain_pairs.add((int(left_chain), int(right_chain) + 2))
        rows.append({
            **coarse,
            "min_interdimer_heavy_a": float(tree_ab.query(tree_cd.data, k=1)[0].min()),
            "clash_atom_pairs": clash_atoms,
            "residue_contact_pairs": len(residue_pairs),
            "chain_pair_contacts": len(chain_pairs),
        })
    ranked = sorted(
        rows,
        key=lambda row: (
            row["clash_atom_pairs"] > 0,
            row["clash_atom_pairs"],
            -row["chain_pair_contacts"],
            -row["residue_contact_pairs"],
        ),
    )
    selected = ranked[: args.top]
    for rank, row in enumerate(selected, 1):
        path = args.out / f"d2_complete_{rank:02d}_phi{row['phi_deg']:.1f}.pdb"
        write_pose(path, chains, row["center"], row["rotation"])
        row["pdb"] = str(path)
        row["rank"] = rank
    fields = [
        "rank", "index", "phi_deg", "axial_offset_a", "fitted_dimer_rotation_deg",
        "min_interdimer_heavy_a", "clash_atom_pairs", "ca_pairs_lt8",
        "ca_clash_pairs_lt3_5", "residue_contact_pairs", "chain_pair_contacts", "pdb",
    ]
    with (args.out / "d2_completion_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(selected)
    best = selected[0]
    print(
        f"fitted_AB_rotation={angle:.2f}deg best_phi={best['phi_deg']:.1f} "
        f"clashes={best['clash_atom_pairs']} contacts={best['residue_contact_pairs']} "
        f"chain_pairs={best['chain_pair_contacts']}"
    )


if __name__ == "__main__":
    main()
