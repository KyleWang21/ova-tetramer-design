#!/usr/bin/env python3
"""Build eight continuation tasks from the best deterministic E283/E286 states."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-experiment", action="append", type=pathlib.Path)
    parser.add_argument(
        "--template-experiment", type=pathlib.Path,
        default=pathlib.Path("experiments/e286_c4_chainmapped_bottomk_tied_design"),
    )
    parser.add_argument(
        "--out", type=pathlib.Path,
        default=pathlib.Path("experiments/e294_c4_detbest_ensemble_bottomk_tied_design"),
    )
    parser.add_argument("--require-marker", type=pathlib.Path)
    parser.add_argument("--anchor-prefix", default="E294")
    args = parser.parse_args()
    project = pathlib.Path(__file__).resolve().parents[1]
    def rooted(path: pathlib.Path) -> pathlib.Path:
        return path if path.is_absolute() else project / path

    source_experiments = [
        rooted(path) for path in (
            args.source_experiment or [
                pathlib.Path("experiments/e283_c4_chainmapped_eightstate_tied_design"),
                pathlib.Path("experiments/e286_c4_chainmapped_bottomk_tied_design"),
            ]
        )
    ]
    template = rooted(args.template_experiment)
    out = rooted(args.out)
    if (out / "INPUT_FROZEN").exists():
        print(f"{out} already frozen")
        return
    require_marker = rooted(args.require_marker) if args.require_marker else (
        source_experiments[-1] / "REMOTE_ARCHIVE_RECOVERED"
    )
    if not require_marker.exists():
        raise RuntimeError(f"required source marker is absent: {require_marker}")
    if out.exists():
        raise RuntimeError(f"refusing to overwrite partial {out}")

    state_pattern = re.compile(r"state(\d+)_iptm$")
    unique: dict[str, dict[str, object]] = {}
    for experiment in source_experiments:
        for path in sorted(experiment.glob("run_*/**/af2diff_candidates.tsv")):
            for row in read_tsv(path):
                sequence = row["sequence"]
                if (
                    len(sequence) != 386
                    or int(row["n_mutations"]) > 20
                    or float(row["surface_mutation_fraction"]) < 0.80
                    or {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS
                ):
                    continue
                state_columns = sorted(
                    (int(match.group(1)), field)
                    for field in row
                    if (match := state_pattern.fullmatch(field))
                )
                if len(state_columns) != 8:
                    raise RuntimeError(f"expected eight state columns in {path}")
                iptms = [float(row[field]) for _, field in state_columns]
                plddts = [float(row[f"state{index}_plddt"]) for index, _ in state_columns]
                record: dict[str, object] = {
                    "origin": str(path.relative_to(project)),
                    "source_candidate": row["candidate"],
                    "selected_stage": row["selected_stage"],
                    "selected_iteration": int(row["selected_iteration"]),
                    "n_mutations": int(row["n_mutations"]),
                    "surface_mutation_fraction": float(row["surface_mutation_fraction"]),
                    "mutations": row["mutations"],
                    "deterministic_min_iptm": min(iptms),
                    "deterministic_mean_iptm": sum(iptms) / len(iptms),
                    "deterministic_min_plddt": min(plddts),
                    "sequence": sequence,
                }
                old = unique.get(sequence)
                score = (
                    float(record["deterministic_min_iptm"]),
                    float(record["deterministic_mean_iptm"]),
                    float(record["deterministic_min_plddt"]),
                    -int(record["n_mutations"]),
                )
                if old is None or score > (
                    float(old["deterministic_min_iptm"]),
                    float(old["deterministic_mean_iptm"]),
                    float(old["deterministic_min_plddt"]),
                    -int(old["n_mutations"]),
                ):
                    unique[sequence] = record

    ranked = sorted(unique.values(), key=lambda row: (
        -float(row["deterministic_min_iptm"]),
        -float(row["deterministic_mean_iptm"]),
        -float(row["deterministic_min_plddt"]),
        int(row["n_mutations"]),
        str(row["sequence"]),
    ))
    selected: list[dict[str, object]] = []
    for minimum_distance in (3, 2, 1):
        for row in ranked:
            if row in selected:
                continue
            if all(hamming(str(row["sequence"]), str(old["sequence"])) >= minimum_distance for old in selected):
                selected.append(row)
            if len(selected) == 4:
                break
        if len(selected) == 4:
            break
    if len(selected) != 4:
        raise RuntimeError(f"could not select four diverse deterministic anchors: {len(selected)}")

    out.mkdir(parents=True)
    for name in (
        "surface_positions.txt", "joint_design.txt", "joint_design_surface.txt",
        "joint_design_hotspot.txt", "hotspot_expansion_audit.tsv",
    ):
        shutil.copy2(template / name, out / name)
    shutil.copytree(template / "targets", out / "targets")
    state_rows = read_tsv(template / "state_manifest.tsv")
    old_prefix = f"{template.relative_to(project)}/targets/"
    new_prefix = f"{out.relative_to(project)}/targets/"
    for row in state_rows:
        row["target_pdb"] = row["target_pdb"].replace(old_prefix, new_prefix)
        row["restraint_positions"] = row["restraint_positions"].replace(old_prefix, new_prefix)
    write_tsv(out / "state_manifest.tsv", state_rows)

    output_rows: list[dict[str, object]] = []
    with (out / "anchors.fasta").open("w") as handle:
        for task_index in range(8):
            source_rank = task_index % 4
            source = selected[source_rank]
            mode = "surface_bottom2" if task_index < 4 else "hotspot_bottom3"
            output = {
                "anchor": f"{args.anchor_prefix}-A{task_index + 1:02d}",
                "task_index": task_index,
                "source_rank": source_rank + 1,
                "continuation_mode": mode,
                **source,
            }
            output_rows.append(output)
            handle.write(
                f">{args.anchor_prefix}-A{task_index + 1:02d} rank={source_rank + 1} mode={mode}\n"
                f"{source['sequence']}\n"
            )
    write_tsv(out / "anchor_manifest.tsv", output_rows)
    write_tsv(out / "source_selection.tsv", selected)
    (out / "INPUT_FROZEN").touch()
    print(
        f"prepared {args.anchor_prefix} unique_final={len(unique)} anchors=8 sources=4 "
        f"best_min8={float(selected[0]['deterministic_min_iptm']):.3f}"
    )


if __name__ == "__main__":
    main()
