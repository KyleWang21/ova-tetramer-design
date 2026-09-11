#!/usr/bin/env python3
"""Prepare the next eight-state Zhongwei tied-AF2 generation.

The source batch is an input generator only.  Candidate acceptance remains
template-free AF3 5-seed x 5-sample validation followed by v1.1 core screening.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
PROTECTED = NATIVE_CYS | set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty TSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


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


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def candidate_row(
    reference: str,
    sequence: str,
    surface: set[int],
    source: str,
    score: float,
    priority: int,
) -> dict[str, object] | None:
    if len(sequence) != len(reference):
        return None
    mutations = {
        i for i, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
        if old != new
    }
    if not mutations or len(mutations) > 24:
        return None
    if len(mutations & surface) / len(mutations) < 0.80:
        return None
    if {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
        return None
    if sequence[257:265] != "SIINFEKL":
        return None
    if any(i in PROTECTED and reference[i - 1] != sequence[i - 1] for i in mutations):
        return None
    return {
        "source": source,
        "sequence": sequence,
        "score": score,
        "source_priority": priority,
        "n_mutations": len(mutations),
        "surface_mutation_fraction": len(mutations & surface) / len(mutations),
        "mutations": ",".join(
            f"{reference[i - 1]}{i}{sequence[i - 1]}" for i in sorted(mutations)
        ),
    }


def add_diverse(
    pool: list[dict[str, object]], selected: list[dict[str, object]], count: int
) -> None:
    ordered = sorted(
        pool,
        key=lambda row: (-int(row["source_priority"]), -float(row["score"]), str(row["source"])),
    )
    for row in ordered:
        if len(selected) >= count:
            return
        sequence = str(row["sequence"])
        if any(sequence == str(old["sequence"]) for old in selected):
            continue
        if selected and min(hamming(sequence, str(old["sequence"])) for old in selected) < 3:
            continue
        selected.append(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument(
        "--pool-dir",
        type=pathlib.Path,
        default=None,
        help="Candidate-pool directory; defaults to SOURCE/candidate_pool_novel.",
    )
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--anchor-prefix", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--registry", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if (args.out / "INPUT_FROZEN").exists():
        print(f"{args.out} already frozen")
        return
    if args.out.exists():
        raise RuntimeError(f"refusing to overwrite partial output: {args.out}")

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    surface = {
        int(token)
        for token in (args.source / "surface_positions.txt").read_text().replace("\n", ",").split(",")
        if token.strip()
    }
    pool_dir = args.pool_dir if args.pool_dir is not None else args.source / "candidate_pool_novel"
    novel: list[dict[str, object]] = []
    for path, score_key in (
        (pool_dir / "hamming_ai_ranking_current_weighted.tsv", "robust_prediction"),
        (pool_dir / "generation_shortlist.tsv", "af2_iptm"),
    ):
        for row in read_tsv(path):
            try:
                score = float(row.get(score_key, "nan"))
            except ValueError:
                continue
            if score != score:
                score = -1.0
            prepared = candidate_row(
                reference, row.get("sequence", ""), surface,
                f"{args.source.name}:{row.get('name', '')}", score, 3,
            )
            if prepared is not None:
                novel.append(prepared)

    accepted: list[dict[str, object]] = []
    for row in read_tsv(args.registry):
        try:
            score = float(row.get("mean_iptm", "nan"))
        except ValueError:
            score = -1.0
        if score != score:
            score = -1.0
        prepared = candidate_row(
            reference, row.get("sequence", ""), surface,
            f"Accepted_AF3:{row.get('candidate', '')}", score, 2,
        )
        if prepared is not None:
            accepted.append(prepared)

    unique: dict[str, dict[str, object]] = {}
    for row in novel + accepted:
        sequence = str(row["sequence"])
        old = unique.get(sequence)
        if old is None or (int(row["source_priority"]), float(row["score"])) > (
            int(old["source_priority"]), float(old["score"])
        ):
            unique[sequence] = row
    novel = [row for row in unique.values() if int(row["source_priority"]) == 3]
    accepted = [row for row in unique.values() if int(row["source_priority"]) == 2]

    selected: list[dict[str, object]] = []
    add_diverse(novel, selected, 4)
    add_diverse(accepted, selected, 8)
    add_diverse(list(unique.values()), selected, 8)
    if len(selected) != 8:
        raise RuntimeError(f"could only select {len(selected)} diverse anchors")

    args.out.mkdir(parents=True)
    shutil.copytree(args.source / "targets", args.out / "targets")
    for name in ("surface_positions.txt", "joint_design.txt", "joint_design_surface.txt", "joint_design_hotspot.txt"):
        shutil.copy2(args.source / name, args.out / name)
    states: list[dict[str, object]] = []
    try:
        out_rel = args.out.resolve().relative_to(pathlib.Path.cwd().resolve())
    except ValueError:
        out_rel = pathlib.Path("experiments") / args.out.name
    for row in read_tsv(args.source / "state_manifest.tsv"):
        updated = dict(row)
        updated["target_pdb"] = str(out_rel / "targets" / pathlib.Path(row["target_pdb"]).name)
        updated["restraint_positions"] = str(out_rel / "targets" / pathlib.Path(row["restraint_positions"]).name)
        states.append(updated)
    write_tsv(args.out / "state_manifest.tsv", states)

    anchors: list[dict[str, object]] = []
    with (args.out / "anchors.fasta").open("w") as fasta:
        for index, row in enumerate(selected, 1):
            anchor = f"{args.anchor_prefix}-A{index:02d}"
            anchors.append({
                "anchor": anchor, "task_index": index - 1,
                "source": row["source"], "score": row["score"],
                "n_mutations": row["n_mutations"],
                "surface_mutation_fraction": row["surface_mutation_fraction"],
                "mutations": row["mutations"], "sequence": row["sequence"],
            })
            fasta.write(f">{anchor} source={row['source']}\n{row['sequence']}\n")
    write_tsv(args.out / "anchor_manifest.tsv", anchors)
    (args.out / "selection_notes.txt").write_text(
        "Next Zhongwei batch: four novel source-batch anchors plus four Accepted_AF3 anchors where possible.\n"
        "All anchors satisfy native Cys, SIINFEKL, <=24 mutations and surface mutation fraction >=0.80;\n"
        "greedy sequence Hamming distance is at least 3 where possible.\n"
    )
    (args.out / "INPUT_FROZEN").touch()
    print(f"prepared {args.out} anchors=8 states={len(states)} novel={sum(int(row['source_priority']) == 3 for row in selected)}")


if __name__ == "__main__":
    main()
