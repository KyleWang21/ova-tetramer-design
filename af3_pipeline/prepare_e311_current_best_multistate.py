#!/usr/bin/env python3
"""Build E311 around four independent AF3 backbones of the current best hit."""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
PROTECTED = NATIVE_CYS | set(range(258, 266))
ORIGINAL_HOTSPOTS = {99, 155, 157, 182, 184, 334, 336, 338}


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


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    active = False
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if active:
                return "".join(chunks)
            active = line[1:].split()[0] == name
        elif active:
            chunks.append(line.strip())
    raise ValueError(f"missing FASTA record {name!r} in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-ensemble", type=pathlib.Path, required=True)
    parser.add_argument("--seed1-summary", type=pathlib.Path, required=True)
    parser.add_argument("--multiseed-summary", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-system", action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    project = pathlib.Path(__file__).resolve().parents[1]
    out = args.out.resolve()
    if (out / "INPUT_FROZEN").exists():
        print(f"{out} already frozen")
        return
    if out.exists():
        raise RuntimeError(f"refusing to overwrite partial output: {out}")

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    surface = read_positions(args.surface_positions)
    strict = [
        row for row in read_tsv(args.seed1_summary)
        if int(row["seed1_strict_pass"]) == 1
    ]
    for path in args.multiseed_summary:
        for row in read_tsv(path):
            if int(row.get("accepted_af3", 0)) != 1:
                continue
            strict.append({
                "system": row.get("source_system") or row["system"],
                "mean_iptm": row["mean_iptm"],
                "min_iptm": row["min_iptm"],
                "sequence": row["sequence"],
                "seed1_strict_pass": "1",
                "accepted_af3": "1",
            })
    strict.sort(key=lambda row: (
        -int(row.get("accepted_af3", 0)),
        -float(row["mean_iptm"]), -float(row["min_iptm"]), row["system"]
    ))
    anchors: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in strict:
        if row["system"] in set(args.exclude_system):
            continue
        sequence = row["sequence"]
        if sequence in seen:
            continue
        mutations = [
            index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
            if old != new
        ]
        surface_fraction = sum(index in surface for index in mutations) / len(mutations)
        if (
            len(mutations) > 20 or surface_fraction < 0.80
            or {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS
            or sequence[257:265] != "SIINFEKL"
        ):
            continue
        seen.add(sequence)
        anchors.append({
            "anchor": f"E311-A{len(anchors) + 1:02d}",
            "task_index": len(anchors),
            "source_system": row["system"],
            "source_seed1_mean_iptm": row["mean_iptm"],
            "source_seed1_min_iptm": row["min_iptm"],
            "n_mutations": len(mutations),
            "surface_mutation_fraction": surface_fraction,
            "mutations": ",".join(
                f"{reference[index - 1]}{index}{sequence[index - 1]}" for index in mutations
            ),
            "sequence": sequence,
        })
        if len(anchors) == 8:
            break
    if len(anchors) != 8:
        raise RuntimeError(f"expected eight strict AF3 anchors, found {len(anchors)}")

    target_rows = read_tsv(args.target_ensemble / "target_manifest.tsv")
    if len(target_rows) != 4:
        raise RuntimeError(f"expected four AF3 target backbones, found {len(target_rows)}")
    out.mkdir(parents=True)
    shutil.copytree(args.target_ensemble / "targets", out / "targets")
    shutil.copy2(args.surface_positions, out / "surface_positions.txt")
    state_rows: list[dict[str, object]] = []
    joint_design: set[int] = set()
    first_interface_hotspots: set[int] | None = None
    for context_index, row in enumerate(target_rows, 1):
        joint_design |= read_positions(pathlib.Path(row["joint_design"]))
        for interface_index in (1, 2):
            pair = row[f"pair{interface_index}"]
            target = out / "targets" / pathlib.Path(row[f"target{interface_index}"]).name
            restraint = out / "targets" / pathlib.Path(row[f"restraint{interface_index}"]).name
            state_rows.append({
                "state_index": len(state_rows) + 1,
                "context_index": context_index,
                "source_cif": row["source_cif"],
                "source_iptm": row["source_iptm"],
                "interface": pair,
                "target_pdb": str(target.relative_to(project)),
                "restraint_positions": str(restraint.relative_to(project)),
                "default_weight": 3.0,
            })
            if interface_index == 1:
                observed = read_positions(restraint) & ORIGINAL_HOTSPOTS
                first_interface_hotspots = (
                    observed if first_interface_hotspots is None
                    else first_interface_hotspots & observed
                )
    if len(state_rows) != 8 or first_interface_hotspots is None:
        raise RuntimeError("failed to create exactly eight AF3 interface states")
    joint_design -= PROTECTED
    write_positions(out / "joint_design.txt", joint_design)
    write_positions(out / "joint_design_surface.txt", joint_design)
    write_positions(out / "joint_design_hotspot.txt", joint_design | first_interface_hotspots)
    write_tsv(out / "state_manifest.tsv", state_rows)
    write_tsv(out / "anchor_manifest.tsv", anchors)
    with (out / "anchors.fasta").open("w") as handle:
        for row in anchors:
            handle.write(
                f">{row['anchor']} source={row['source_system']}\n{row['sequence']}\n"
            )
    (out / "hotspot_expansion_audit.tsv").write_text(
        "original_hotspots\tall_four_first_interface_hotspots\tnewly_enabled_hotspots\n"
        f"{','.join(map(str, sorted(ORIGINAL_HOTSPOTS)))}\t"
        f"{','.join(map(str, sorted(first_interface_hotspots)))}\t"
        f"{','.join(map(str, sorted(first_interface_hotspots - joint_design)))}\n"
    )
    (out / "INPUT_FROZEN").touch()
    print(
        f"prepared E311 anchors=8 AF3_states=8 design_positions={len(joint_design)} "
        f"best_seed1={float(anchors[0]['source_seed1_mean_iptm']):.3f}"
    )


if __name__ == "__main__":
    main()
