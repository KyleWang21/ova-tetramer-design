#!/usr/bin/env python3
"""Prepare a diverse second tied-AF2 generation from E311 and AF3 hits."""

from __future__ import annotations

import csv
import pathlib
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
PROTECTED = NATIVE_CYS | set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty TSV: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


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
    source_priority: int,
) -> dict[str, object] | None:
    if len(sequence) != len(reference):
        return None
    mutations = {
        index
        for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
        if old != new
    }
    if not mutations or len(mutations) > 24:
        return None
    surface_fraction = len(mutations & surface) / len(mutations)
    if surface_fraction < 0.80:
        return None
    if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
        return None
    if sequence[257:265] != "SIINFEKL":
        return None
    return {
        "source": source,
        "sequence": sequence,
        "score": score,
        "source_priority": source_priority,
        "n_mutations": len(mutations),
        "surface_mutation_fraction": surface_fraction,
        "mutations": ",".join(
            f"{reference[index - 1]}{index}{sequence[index - 1]}"
            for index in sorted(mutations)
        ),
    }


def select_diverse(
    pool: list[dict[str, object]], selected: list[dict[str, object]], count: int
) -> None:
    remaining = sorted(
        pool,
        key=lambda row: (-int(row["source_priority"]), -float(row["score"]), str(row["source"])),
    )
    for row in remaining:
        if len(selected) >= count:
            return
        if any(row["sequence"] == chosen["sequence"] for chosen in selected):
            continue
        if selected and min(hamming(str(row["sequence"]), str(chosen["sequence"])) for chosen in selected) < 3:
            continue
        selected.append(row)


def main() -> None:
    project = pathlib.Path(__file__).resolve().parents[1]
    source = project / "experiments/e311_c4_e294_detbest_continuation_tied_design"
    pool_dir = source / "candidate_pool_novel"
    reference_path = project / "OVA_P3-13R_四聚体候选_AA.fasta"
    reference = fasta_sequence(reference_path, "D0-P3-13R")
    surface = {
        int(token)
        for token in (source / "surface_positions.txt").read_text().replace("\n", ",").split(",")
        if token.strip()
    }
    out = project / "experiments/e318_c4_e311_diverse_tied_design"
    if (out / "INPUT_FROZEN").exists():
        print(f"{out} already frozen")
        return
    if out.exists():
        raise RuntimeError(f"refusing to overwrite partial output: {out}")

    novel: list[dict[str, object]] = []
    for path, score_key in [
        (pool_dir / "hamming_ai_ranking_current_weighted.tsv", "robust_prediction"),
        (pool_dir / "generation_shortlist.tsv", "af2_iptm"),
    ]:
        for row in read_tsv(path):
            sequence = row["sequence"]
            try:
                score = float(row.get(score_key, "nan"))
            except ValueError:
                continue
            if score != score:
                score = -1.0
            prepared = candidate_row(
                reference, sequence, surface, f"E311:{row['name']}", score, 3
            )
            if prepared is not None:
                novel.append(prepared)

    accepted: list[dict[str, object]] = []
    registry = project / "experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
    for row in read_tsv(registry):
        prepared = candidate_row(
            reference,
            row["sequence"],
            surface,
            f"Accepted_AF3:{row['candidate']}",
            float(row["mean_iptm"]),
            2,
        )
        if prepared is not None:
            accepted.append(prepared)

    unique: dict[str, dict[str, object]] = {}
    for row in novel + accepted:
        old = unique.get(str(row["sequence"]))
        if old is None or (int(row["source_priority"]), float(row["score"])) > (
            int(old["source_priority"]), float(old["score"])
        ):
            unique[str(row["sequence"])] = row
    novel = [row for row in unique.values() if str(row["source"]).startswith("E311:")]
    accepted = [row for row in unique.values() if str(row["source"]).startswith("Accepted_AF3:")]

    selected: list[dict[str, object]] = []
    select_diverse(novel, selected, 4)
    select_diverse(accepted, selected, 8)
    select_diverse(list(unique.values()), selected, 8)
    if len(selected) != 8:
        raise RuntimeError(f"could only select {len(selected)} diverse anchors")

    out.mkdir(parents=True)
    shutil.copytree(source / "targets", out / "targets")
    for name in (
        "surface_positions.txt",
        "joint_design.txt",
        "joint_design_surface.txt",
        "joint_design_hotspot.txt",
    ):
        shutil.copy2(source / name, out / name)

    source_states = read_tsv(source / "state_manifest.tsv")
    state_rows: list[dict[str, object]] = []
    for row in source_states:
        updated = dict(row)
        updated["target_pdb"] = str(pathlib.Path("experiments/e318_c4_e311_diverse_tied_design/targets") / pathlib.Path(row["target_pdb"]).name)
        updated["restraint_positions"] = str(pathlib.Path("experiments/e318_c4_e311_diverse_tied_design/targets") / pathlib.Path(row["restraint_positions"]).name)
        state_rows.append(updated)
    write_tsv(out / "state_manifest.tsv", state_rows)

    anchor_rows: list[dict[str, object]] = []
    with (out / "anchors.fasta").open("w") as fasta:
        for index, row in enumerate(selected, 1):
            anchor = f"E318-A{index:02d}"
            anchor_rows.append({
                "anchor": anchor,
                "task_index": index - 1,
                "source": row["source"],
                "score": row["score"],
                "n_mutations": row["n_mutations"],
                "surface_mutation_fraction": row["surface_mutation_fraction"],
                "mutations": row["mutations"],
                "sequence": row["sequence"],
            })
            fasta.write(f">{anchor} source={row['source']}\n{row['sequence']}\n")
    write_tsv(out / "anchor_manifest.tsv", anchor_rows)
    (out / "selection_notes.txt").write_text(
        "E318 uses four E311 trajectory anchors followed by four high-mean Accepted_AF3 anchors;\n"
        "greedy sequence Hamming distance is at least 3 where possible. All anchors satisfy\n"
        "native Cys, SIINFEKL, <=24 mutations and surface mutation fraction >=0.80.\n"
    )
    (out / "INPUT_FROZEN").touch()
    print(
        f"prepared {out} anchors=8 states={len(state_rows)} "
        f"novel={sum(str(row['source']).startswith('E311:') for row in selected)}"
    )


if __name__ == "__main__":
    main()
