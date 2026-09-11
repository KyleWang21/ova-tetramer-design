#!/usr/bin/env python3
"""Collect diverse <25-mutation discrete states from the e115 AF2 trajectories."""

from __future__ import annotations

import argparse
import csv
import pathlib


NATIVE_CYS = {12, 31, 74, 121, 368, 383}


def read_fasta(path: pathlib.Path) -> str:
    chunks: list[str] = []
    seen_header = False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
        elif seen_header:
            chunks.append(line.strip())
    if not chunks:
        raise ValueError(f"no FASTA record in {path}")
    return "".join(chunks)


def number(value: str, default: float = -999.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--reference", type=pathlib.Path, required=True)
    ap.add_argument("--baseline", type=pathlib.Path, action="append", default=[])
    ap.add_argument(
        "--exclude-fasta",
        type=pathlib.Path,
        action="append",
        default=[],
        help="drop sequences already evaluated in an earlier AF3 campaign",
    )
    ap.add_argument(
        "--out-dir",
        type=pathlib.Path,
        help="write the shortlist here instead of modifying the trajectory experiment",
    )
    ap.add_argument("--top", type=int, default=64)
    ap.add_argument("--per-run", type=int, default=6)
    args = ap.parse_args()
    reference = read_fasta(args.reference)
    out_dir = args.out_dir or args.experiment
    out_dir.mkdir(parents=True, exist_ok=True)
    excluded: set[str] = set()
    for path in args.exclude_fasta:
        header_seen = False
        chunks: list[str] = []
        for line in path.read_text().splitlines() + [">"]:
            if line.startswith(">"):
                if header_seen and chunks:
                    excluded.add("".join(chunks))
                header_seen = True
                chunks = []
            elif header_seen:
                chunks.append(line.strip())

    rows: list[dict[str, object]] = []
    for path in sorted(args.experiment.glob("run_*/stage*/trajectory.tsv")):
        run = path.parents[1].name
        stage = path.parent.name
        for row in csv.DictReader(path.open(), delimiter="\t"):
            sequence = row.get("sequence", "")
            if len(sequence) != 386:
                continue
            nmut = sum(a != b for a, b in zip(reference, sequence))
            if nmut >= 25 or sequence[257:265] != reference[257:265]:
                continue
            if {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
                continue
            rows.append({
                "source": f"{run}/{stage}/{row.get('stage','')}/{row.get('iteration','')}",
                "run": run,
                "stage": stage,
                "af2_stage": row.get("stage", ""),
                "n_mutations": nmut,
                "mutations": ",".join(
                    f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b
                ),
                "af2_iptm": number(row.get("iptm", "")),
                "af2_plddt": number(row.get("plddt", "")),
                "target_contact_loss": number(row.get("first_target_con", ""), 999.0),
                "loss": number(row.get("loss", ""), 999.0),
                "sequence": sequence,
            })

    for index, path in enumerate(args.baseline, 1):
        sequence = read_fasta(path)
        nmut = sum(a != b for a, b in zip(reference, sequence))
        rows.append({
            "source": f"baseline/{path.stem}", "run": f"baseline_{index}", "stage": "baseline",
            "af2_stage": "baseline", "n_mutations": nmut,
            "mutations": ",".join(f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b),
            "af2_iptm": -1.0, "af2_plddt": -1.0, "target_contact_loss": 999.0,
            "loss": 999.0, "sequence": sequence,
        })

    unique: dict[str, dict[str, object]] = {}
    for row in rows:
        seq = str(row["sequence"])
        if seq in excluded:
            continue
        rank = (
            -float(row["af2_iptm"]),
            float(row["target_contact_loss"]),
            0 if row["af2_stage"] in {"soft", "hard"} else 1,
            -int(row["n_mutations"]),
        )
        if seq not in unique or rank < unique[seq]["_rank"]:
            row["_rank"] = rank
            unique[seq] = row

    by_run: dict[str, list[dict[str, object]]] = {}
    for row in unique.values():
        by_run.setdefault(f"{row['run']}/{row['stage']}", []).append(row)
    selected: list[dict[str, object]] = []
    for run_rows in by_run.values():
        run_rows.sort(key=lambda row: row["_rank"])
        selected.extend(run_rows[: args.per_run])
    selected.sort(key=lambda row: row["_rank"])
    selected = selected[: args.top]
    if not selected:
        raise SystemExit("no eligible trajectory states")

    fields = [key for key in selected[0] if key != "_rank"]
    with (out_dir / "trajectory_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows({key: row[key] for key in fields} for row in selected)
    with (out_dir / "trajectory_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            handle.write(
                f">OVA-IPTM80-{index:03d} nmut={row['n_mutations']} source={row['source']} "
                f"af2iptm={float(row['af2_iptm']):.3f}\n{row['sequence']}\n"
            )
    print(
        f"excluded={len(excluded)} eligible_unique_novel={len(unique)} "
        f"selected={len(selected)} out={out_dir}"
    )


if __name__ == "__main__":
    main()
