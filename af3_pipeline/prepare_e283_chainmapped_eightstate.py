#!/usr/bin/env python3
"""Prepare four AF3 backbones after canonical contact-map chain remapping.

Homotetramer chain identifiers are arbitrary across independent AF3 models.
This script maps every physical interface to one of the two reference contact
maps before building the eight tied-AF2 dimer states.  It prevents an unrelated
chain pair from being optimized merely because it happened to be labelled AD.
"""

from __future__ import annotations

import csv
import itertools
import pathlib
import shutil

import numpy as np


PROTECTED = {12, 31, 74, 121, 368, 383} | set(range(258, 266))
CONTACT_MAP_CUTOFF = 5.3
RESTRAINT_CUTOFF = 8.0
MINIMUM_JACCARD = 0.70


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_positions(path: pathlib.Path) -> set[int]:
    return {
        int(value)
        for value in path.read_text().replace("\n", ",").split(",")
        if value.strip()
    }


def write_positions(path: pathlib.Path, values: set[int]) -> None:
    path.write_text(",".join(map(str, sorted(values))) + "\n")


def pdb_atoms(path: pathlib.Path) -> tuple[list[str], dict[str, list[tuple[int, np.ndarray]]]]:
    lines = path.read_text().splitlines()
    chains: dict[str, list[tuple[int, np.ndarray]]] = {}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element == "H":
            continue
        chain = line[21]
        residue = int(line[22:26])
        xyz = np.asarray(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=np.float32,
        )
        chains.setdefault(chain, []).append((residue, xyz))
    if set(chains) != set("ABCD"):
        raise RuntimeError(f"expected chains ABCD in {path}, found {sorted(chains)}")
    return lines, chains


def residue_contact_map(
    chains: dict[str, list[tuple[int, np.ndarray]]], pair: str, cutoff: float
) -> set[tuple[int, int]]:
    left = chains[pair[0]]
    right = chains[pair[1]]
    left_residue = np.asarray([item[0] for item in left], dtype=np.int16)
    right_residue = np.asarray([item[0] for item in right], dtype=np.int16)
    left_xyz = np.stack([item[1] for item in left])
    right_xyz = np.stack([item[1] for item in right])
    cutoff2 = cutoff * cutoff
    contacts: set[tuple[int, int]] = set()
    for start in range(0, len(left_xyz), 256):
        distance2 = np.square(
            left_xyz[start : start + 256, None, :] - right_xyz[None, :, :]
        ).sum(axis=2)
        left_index, right_index = np.where(distance2 <= cutoff2)
        contacts.update(
            (int(left_residue[start + i]), int(right_residue[j]))
            for i, j in zip(left_index, right_index, strict=True)
        )
    return contacts


def contact_positions(
    chains: dict[str, list[tuple[int, np.ndarray]]], pair: str, cutoff: float
) -> set[int]:
    contacts = residue_contact_map(chains, pair, cutoff)
    return {position for contact in contacts for position in contact}


def jaccard(left: set[tuple[int, int]], right: set[tuple[int, int]]) -> float:
    return len(left & right) / len(left | right) if left | right else 0.0


def best_orientation_jaccard(
    candidate: set[tuple[int, int]], reference: set[tuple[int, int]]
) -> float:
    reversed_candidate = {(right, left) for left, right in candidate}
    return max(jaccard(candidate, reference), jaccard(reversed_candidate, reference))


def extract_pair(lines: list[str], pair: str, path: pathlib.Path) -> None:
    output: list[str] = []
    for chain in pair:
        output.extend(
            line for line in lines
            if line.startswith("ATOM") and line[21] == chain
        )
        output.append("TER")
    output.append("END")
    path.write_text("\n".join(output) + "\n")


def main() -> None:
    project = pathlib.Path(__file__).resolve().parents[1]
    source = project / "experiments/e268_c4_e251_softbasin_stabilization_tied_design"
    old = project / "experiments/e280_c4_eightstate_multibackbone_tied_design"
    out = project / "experiments/e283_c4_chainmapped_eightstate_tied_design"
    if (out / "INPUT_FROZEN").exists():
        print(f"{out} already frozen")
        return
    if out.exists():
        raise RuntimeError(f"refusing to overwrite partially prepared {out}")
    out.mkdir(parents=True)
    target_dir = out / "targets"
    target_dir.mkdir()
    shutil.copy2(old / "anchors.fasta", out / "anchors.fasta")
    shutil.copy2(old / "anchor_manifest.tsv", out / "anchor_manifest.tsv")
    shutil.copy2(old / "surface_positions.txt", out / "surface_positions.txt")

    contexts = read_tsv(source / "target_manifest.tsv")
    if len(contexts) != 4:
        raise RuntimeError(f"expected four backbone contexts, found {len(contexts)}")
    parsed: list[tuple[list[str], dict[str, list[tuple[int, np.ndarray]]]]] = []
    pair_maps: list[dict[str, set[tuple[int, int]]]] = []
    for context in contexts:
        lines, chains = pdb_atoms(project / context["target_c4"])
        parsed.append((lines, chains))
        pair_maps.append({
            "".join(pair): residue_contact_map(chains, "".join(pair), CONTACT_MAP_CUTOFF)
            for pair in itertools.combinations("ABCD", 2)
        })

    reference_pairs = (contexts[0]["pair1"], contexts[0]["pair2"])
    reference_maps = tuple(pair_maps[0][pair] for pair in reference_pairs)
    surface = read_positions(out / "surface_positions.txt")
    design_positions: set[int] = set()
    for row in read_tsv(out / "anchor_manifest.tsv"):
        for mutation in row["mutations"].split(","):
            if mutation:
                design_positions.add(int(mutation[1:-1]))

    state_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    for context_index, (context, (lines, chains), maps) in enumerate(
        zip(contexts, parsed, pair_maps, strict=True), 1
    ):
        legacy_pairs = (context["pair1"], context["pair2"])
        selected_pairs: list[str] = []
        for class_index, (reference_pair, reference_map, legacy_pair) in enumerate(
            zip(reference_pairs, reference_maps, legacy_pairs, strict=True), 1
        ):
            scored = sorted(
                (
                    best_orientation_jaccard(contact_map, reference_map),
                    pair.startswith("A"),
                    pair,
                )
                for pair, contact_map in maps.items()
            )
            score, _, selected_pair = scored[-1]
            if score < MINIMUM_JACCARD:
                raise RuntimeError(
                    f"context {context_index} class {class_index} has no canonical mapping: {score:.3f}"
                )
            if selected_pair in selected_pairs:
                raise RuntimeError(f"context {context_index} mapped both classes to {selected_pair}")
            selected_pairs.append(selected_pair)
            positions = contact_positions(chains, selected_pair, RESTRAINT_CUTOFF)
            restraint = positions - PROTECTED
            design_positions.update((positions & surface) - PROTECTED)
            stem = f"context{context_index}_class{class_index}_{selected_pair.lower()}"
            target = target_dir / f"{stem}.pdb"
            restraint_path = target_dir / f"{stem}_restraint.txt"
            extract_pair(lines, selected_pair, target)
            write_positions(restraint_path, restraint)
            state_rows.append({
                "state_index": len(state_rows) + 1,
                "context_index": context_index,
                "source_cif": context["source_cif"],
                "source_iptm": context["source_iptm"],
                "interface": selected_pair,
                "target_pdb": str(target.relative_to(project)),
                "restraint_positions": str(restraint_path.relative_to(project)),
                "default_weight": 3.0,
            })
            audit_rows.append({
                "context_index": context_index,
                "canonical_class": class_index,
                "reference_pair": reference_pair,
                "legacy_pair": legacy_pair,
                "legacy_jaccard": best_orientation_jaccard(maps[legacy_pair], reference_map),
                "selected_pair": selected_pair,
                "selected_jaccard": score,
                "contact_pairs_5p3": len(maps[selected_pair]),
                "restraint_positions_8p0": len(restraint),
                "mapping_pass": int(score >= MINIMUM_JACCARD),
            })

    design_positions -= PROTECTED
    if not set(range(1, 9)) == {int(row["state_index"]) for row in state_rows}:
        raise RuntimeError("expected exactly eight mapped states")
    write_positions(out / "joint_design.txt", design_positions)
    write_tsv(out / "state_manifest.tsv", state_rows)
    write_tsv(out / "chain_mapping_audit.tsv", audit_rows)
    (out / "INPUT_FROZEN").touch()
    print(
        f"prepared E283 states={len(state_rows)} design_positions={len(design_positions)} "
        f"minimum_mapping_jaccard={min(float(row['selected_jaccard']) for row in audit_rows):.3f}"
    )


if __name__ == "__main__":
    main()
