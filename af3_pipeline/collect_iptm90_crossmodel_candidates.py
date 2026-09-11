#!/usr/bin/env python3
"""Merge joint-C4 AF2 trajectories and cross-model surrogate proposals.

This script only chooses sequences for a new template-free AF3 screen.  It does
not treat AF2 or the surrogate as validation evidence.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re

import numpy as np

from run_final20_failclosed_filter import glyco_sites, reference_rsasa


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def fasta_records(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def fasta_sequence(path: pathlib.Path, prefix: str | None = None) -> str:
    for header, sequence in fasta_records(path):
        if prefix is None or header.startswith(prefix):
            return sequence
    raise ValueError(f"FASTA record {prefix!r} not found in {path}")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: object, default: float = math.nan) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def finite_extreme(values: list[object], *, minimum: bool) -> float:
    parsed = [number(value) for value in values]
    parsed = [value for value in parsed if math.isfinite(value)]
    if not parsed:
        return math.nan
    return min(parsed) if minimum else max(parsed)


def average_ranks(values: np.ndarray, higher: bool) -> np.ndarray:
    transformed = -values if higher else values
    order = np.argsort(transformed, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--reference-cif", type=pathlib.Path, required=True)
    parser.add_argument("--surrogate-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument(
        "--surrogate-rank-column", default="acquisition_average_rank",
        help="lower-is-better ranking column in surrogate TSV inputs",
    )
    parser.add_argument("--mpnn-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--exclude-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--exclude-fasta", type=pathlib.Path, action="append", default=[])
    parser.add_argument(
        "--exclude-glob",
        action="append",
        default=[],
        help="Workspace-relative glob of TSV/FASTA files whose sequences are excluded.",
    )
    parser.add_argument(
        "--exclude-ignore-path",
        type=pathlib.Path,
        action="append",
        default=[],
        help="Exact output path to ignore while expanding --exclude-glob (for idempotent reruns).",
    )
    parser.add_argument(
        "--run-glob",
        default="run_*",
        help="Experiment-relative design-run glob; use e.g. run_073v7_* for one wave.",
    )
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top-af2", type=int, default=64)
    parser.add_argument("--top-surrogate", type=int, default=64)
    parser.add_argument("--top-mpnn", type=int, default=0)
    parser.add_argument("--minimum-distance", type=int, default=2)
    parser.add_argument("--min-mutations", type=int, default=1)
    parser.add_argument("--max-mutations", type=int, default=24)
    args = parser.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa P3-13R reference, found {len(reference)}")
    rsasa = reference_rsasa(args.reference_cif, "A", reference)
    excluded: set[str] = set()
    ignored_exclusion_paths = {path.resolve() for path in args.exclude_ignore_path}
    for path in args.exclude_tsv:
        excluded.update(row.get("sequence", "") for row in read_tsv(path))
    for path in args.exclude_fasta:
        excluded.update(sequence for _, sequence in fasta_records(path))
    for pattern in args.exclude_glob:
        for path in sorted(pathlib.Path.cwd().glob(pattern)):
            if not path.is_file():
                continue
            if path.resolve() in ignored_exclusion_paths:
                continue
            if path.suffix.lower() in {".fa", ".faa", ".fasta"}:
                excluded.update(sequence for _, sequence in fasta_records(path))
            elif path.suffix.lower() == ".tsv":
                excluded.update(row.get("sequence", "") for row in read_tsv(path))

    def validate(sequence: str) -> tuple[bool, int, float, str]:
        if len(sequence) != len(reference):
            return False, 0, 0.0, "length"
        changed = [i for i, (old, new) in enumerate(zip(reference, sequence), 1) if old != new]
        nmut = len(changed)
        if not args.min_mutations <= nmut <= args.max_mutations:
            return False, nmut, 0.0, "mutation_budget"
        if {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            return False, nmut, 0.0, "cysteine"
        if any(sequence[pos - 1] != reference[pos - 1] for pos in EPITOPE):
            return False, nmut, 0.0, "SIINFEKL"
        if glyco_sites(sequence) != glyco_sites(reference):
            return False, nmut, 0.0, "glycosylation"
        surface_fraction = sum(rsasa.get(pos, 0.0) >= 0.20 for pos in changed) / nmut
        if surface_fraction < 0.80:
            return False, nmut, surface_fraction, "surface_fraction"
        if any(
            rsasa.get(pos, 0.0) < 0.20
            and (reference[pos - 1] in "PG" or sequence[pos - 1] in "PG")
            for pos in changed
        ):
            return False, nmut, surface_fraction, "buried_pro_gly"
        if sequence in excluded:
            return False, nmut, surface_fraction, "previously_evaluated"
        return True, nmut, surface_fraction, ""

    rejected: dict[str, int] = {}
    unique: dict[str, dict[str, object]] = {}
    af2_paths = sorted(args.experiment.glob(f"{args.run_glob}/joint_*/trajectory.tsv"))
    af2_paths += sorted(args.experiment.glob(f"{args.run_glob}/joint_*/af2diff_candidates.tsv"))
    for path in af2_paths:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            valid, nmut, surface_fraction, reason = validate(sequence)
            if not valid:
                rejected[reason] = rejected.get(reason, 0) + 1
                continue
            stage = row.get("stage", "final")
            # Late logits states are already nearly discrete and sometimes give
            # a better two-interface compromise than the annealed endpoint.
            if path.name == "trajectory.tsv":
                if stage not in {"logits", "soft", "hard"}:
                    continue
                if stage == "logits" and number(row.get("soft"), 0.0) < 0.70:
                    continue
            def state_values(suffix: str, *, parity: int | None = None) -> list[object]:
                values: list[tuple[int, object]] = []
                for key, value in row.items():
                    match = re.fullmatch(r"state(\d+)_" + suffix, key)
                    if match:
                        index = int(match.group(1))
                        if parity is None or index % 2 == parity:
                            values.append((index, value))
                return [value for _, value in sorted(values)]

            candidate = {
                "method": "joint_tied_af2",
                "source": str(path.relative_to(args.experiment)) + f"/{stage}/{row.get('iteration', '')}",
                "sequence": sequence,
                "n_mutations": nmut,
                "surface_mutation_fraction": surface_fraction,
                "mutations": mutation_string(reference, sequence),
                "af2_iptm": finite_extreme(
                    [row.get("iptm"), *state_values("iptm")],
                    minimum=True,
                ),
                "af2_plddt": finite_extreme(
                    [row.get("plddt"), *state_values("plddt")],
                    minimum=True,
                ),
                "af2_interface_pae": finite_extreme(
                    [row.get("interface_pae"), *state_values("interface_pae")],
                    minimum=False,
                ),
                "af2_fold_target_con": finite_extreme(
                    [row.get("fold_target_con"), *state_values("fold_target_con")],
                    minimum=False,
                ),
                "af2_first_target_con": finite_extreme(
                    [row.get("first_target_con"), *state_values("interface_target_con", parity=1)],
                    minimum=False,
                ),
                "af2_second_target_con": finite_extreme(
                    [row.get("second_target_con"), *state_values("interface_target_con", parity=0)],
                    minimum=False,
                ),
                "surrogate_acquisition_rank": math.nan,
                "predicted_af3": math.nan,
                "predicted_protenix": math.nan,
                "mpnn_score": math.nan,
            }
            old = unique.get(sequence)
            old_key = -math.inf if old is None else number(old.get("af2_iptm"), -math.inf)
            if old is None or number(candidate["af2_iptm"], -math.inf) > old_key:
                unique[sequence] = candidate

    for path in args.surrogate_tsv:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            valid, nmut, surface_fraction, reason = validate(sequence)
            if not valid:
                rejected[reason] = rejected.get(reason, 0) + 1
                continue
            candidate = {
                "method": "crossmodel_surrogate",
                "source": str(path) + f"/{row.get('name', '')}",
                "sequence": sequence,
                "n_mutations": nmut,
                "surface_mutation_fraction": surface_fraction,
                "mutations": mutation_string(reference, sequence),
                "af2_iptm": math.nan,
                "af2_plddt": math.nan,
                "af2_interface_pae": math.nan,
                "af2_fold_target_con": math.nan,
                "af2_first_target_con": math.nan,
                "af2_second_target_con": math.nan,
                "surrogate_acquisition_rank": number(
                    row.get(args.surrogate_rank_column), math.inf
                ),
                "predicted_af3": number(row.get("predicted_af3")),
                "predicted_protenix": number(row.get("predicted_protenix")),
                "mpnn_score": math.nan,
            }
            old = unique.get(sequence)
            if old is None or old["method"] != "joint_tied_af2":
                unique[sequence] = candidate

    for path in args.mpnn_tsv:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            valid, nmut, surface_fraction, reason = validate(sequence)
            if not valid:
                rejected[reason] = rejected.get(reason, 0) + 1
                continue
            candidate = {
                "method": "tied_c4_proteinmpnn",
                "source": str(path) + f"/rank{row.get('rank', '')}",
                "sequence": sequence,
                "n_mutations": nmut,
                "surface_mutation_fraction": surface_fraction,
                "mutations": mutation_string(reference, sequence),
                "af2_iptm": math.nan,
                "af2_plddt": math.nan,
                "af2_interface_pae": math.nan,
                "af2_fold_target_con": math.nan,
                "af2_first_target_con": math.nan,
                "af2_second_target_con": math.nan,
                "surrogate_acquisition_rank": math.nan,
                "predicted_af3": math.nan,
                "predicted_protenix": math.nan,
                "mpnn_score": number(row.get("mpnn_score")),
                "generation_average_rank": number(row.get("generation_average_rank"), math.inf),
            }
            old = unique.get(sequence)
            if (
                old is None
                or old["method"] not in {"joint_tied_af2", "tied_c4_proteinmpnn"}
                or (
                    old["method"] == "tied_c4_proteinmpnn"
                    and number(candidate["generation_average_rank"], math.inf)
                    < number(old.get("generation_average_rank"), math.inf)
                )
            ):
                unique[sequence] = candidate

    by_method = {
        method: [row for row in unique.values() if row["method"] == method]
        for method in ("joint_tied_af2", "crossmodel_surrogate", "tied_c4_proteinmpnn")
    }
    af2 = by_method["joint_tied_af2"]
    if af2:
        metrics = [
            ("af2_iptm", True, 0.30), ("af2_plddt", True, 0.10),
            ("af2_interface_pae", False, 0.15), ("af2_fold_target_con", False, 0.10),
            ("af2_first_target_con", False, 0.175), ("af2_second_target_con", False, 0.175),
        ]
        total = np.zeros(len(af2))
        used_weight = 0.0
        for field, higher, weight in metrics:
            values = np.asarray([number(row[field]) for row in af2])
            finite = np.isfinite(values)
            if not finite.all():
                continue
            total += weight * average_ranks(values, higher)
            used_weight += weight
        if used_weight == 0:
            raise RuntimeError("joint AF2 states have no complete ranking metrics")
        for row, rank in zip(af2, total / used_weight):
            row["generation_average_rank"] = rank
    surrogate = by_method["crossmodel_surrogate"]
    for row in surrogate:
        row["generation_average_rank"] = row["surrogate_acquisition_rank"]
    mpnn = by_method["tied_c4_proteinmpnn"]

    selected: list[dict[str, object]] = []
    selected_sequences: list[str] = []
    for method, quota in (
        ("joint_tied_af2", args.top_af2),
        ("crossmodel_surrogate", args.top_surrogate),
        ("tied_c4_proteinmpnn", args.top_mpnn),
    ):
        if quota <= 0:
            continue
        ordered = sorted(by_method[method], key=lambda row: (
            number(row["generation_average_rank"], math.inf), int(row["n_mutations"]), str(row["sequence"])
        ))
        for minimum in (args.minimum_distance, 1):
            for row in ordered:
                if row in selected:
                    continue
                sequence = str(row["sequence"])
                if all(hamming(sequence, old) >= minimum for old in selected_sequences):
                    selected.append(row)
                    selected_sequences.append(sequence)
                if sum(item["method"] == method for item in selected) >= quota:
                    break
            if sum(item["method"] == method for item in selected) >= quota:
                break

    if not selected:
        raise SystemExit("no novel sequence passed the frozen generation constraints")
    selected.sort(key=lambda row: (
        {"joint_tied_af2": 0, "crossmodel_surrogate": 1, "tied_c4_proteinmpnn": 2}[str(row["method"])],
        number(row["generation_average_rank"], math.inf),
    ))
    args.out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for index, row in enumerate(selected, 1):
        output = dict(row)
        output["name"] = f"OVA-C4-IPTM90-{index:03d}"
        rows.append(output)
    fields = ["name"] + [field for field in rows[0] if field != "name"]
    with (args.out / "generation_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "generation_shortlist.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} method={row['method']} nmut={row['n_mutations']} "
                f"generation_rank={float(row['generation_average_rank']):.3f}\n{row['sequence']}\n"
            )
    with (args.out / "rejection_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reason", "count"])
        writer.writerows(sorted(rejected.items()))
    print(
        f"joint_AF2_paths={len(af2_paths)} valid_unique={len(af2)}; "
        f"surrogate_valid_unique={len(surrogate)}; mpnn_valid_unique={len(mpnn)}; "
        f"selected={len(rows)}; out={args.out}"
    )


if __name__ == "__main__":
    main()
