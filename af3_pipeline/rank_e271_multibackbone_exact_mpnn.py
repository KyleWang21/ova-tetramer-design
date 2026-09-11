#!/usr/bin/env python3
"""Aggregate eight exact ProteinMPNN scores and freeze a diverse AF3 shortlist."""

from __future__ import annotations

import argparse
import csv
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


AA20 = set("ARNDCQEGHILKMFPSTWYV")


def average_ranks(values: list[float], *, higher: bool = False) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index], reverse=higher)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def hamming(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming comparison length mismatch")
    return sum(a != b for a, b in zip(left, right))


def fasta_sequences(path: Path) -> set[str]:
    sequences: set[str] = set()
    chunks: list[str] = []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if chunks:
                sequence = "".join(chunks)
                if len(sequence) == 386 and set(sequence) <= AA20:
                    sequences.add(sequence)
            chunks = []
        elif line:
            chunks.append(line)
    if chunks:
        sequence = "".join(chunks)
        if len(sequence) == 386 and set(sequence) <= AA20:
            sequences.add(sequence)
    return sequences


def tabular_sequences(path: Path) -> set[str]:
    sequences: set[str] = set()
    try:
        with path.open(errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                for field in ("sequence", "monomer_sequence", "candidate_sequence"):
                    sequence = (row.get(field) or "").strip()
                    if len(sequence) == 386 and set(sequence) <= AA20:
                        sequences.add(sequence)
    except (csv.Error, UnicodeDecodeError):
        return set()
    return sequences


def load_exclusions(project: Path, patterns: list[str]) -> tuple[set[str], list[str]]:
    sequences: set[str] = set()
    paths: list[str] = []
    for pattern in patterns:
        for raw_path in glob.glob(str(project / pattern), recursive=True):
            path = Path(raw_path)
            if not path.is_file():
                continue
            paths.append(str(path.relative_to(project)))
            if path.suffix.lower() in {".fa", ".faa", ".fasta"}:
                sequences |= fasta_sequences(path)
            elif path.suffix.lower() in {".tsv", ".txt"}:
                sequences |= tabular_sequences(path)
    return sequences, sorted(set(paths))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--source-ranking", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--top", type=int, default=32)
    parser.add_argument("--minimum-distance", type=int, default=3)
    parser.add_argument("--exclude-glob", action="append", default=[])
    args = parser.parse_args()

    score_paths = sorted(args.experiment.glob("task_*/candidate_scores.tsv"))
    if len(score_paths) != 8:
        raise ValueError(f"expected eight E271 score tables, found {len(score_paths)}")
    with args.source_ranking.open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != 512 or len({row["sequence"] for row in source_rows}) != 512:
        raise ValueError("source ranking is not the frozen 512-sequence library")
    source_by_sequence = {row["sequence"]: row for row in source_rows}

    scores: dict[str, list[dict[str, str]]] = {}
    contexts: set[str] = set()
    for path in score_paths:
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) != 512 or len({row["sequence"] for row in rows}) != 512:
            raise ValueError(f"{path} is not a complete 512-sequence score table")
        context_values = {row["context"] for row in rows}
        if len(context_values) != 1:
            raise ValueError(f"{path} contains multiple contexts")
        context = next(iter(context_values))
        if context in contexts:
            raise ValueError(f"duplicate context {context}")
        contexts.add(context)
        for row in rows:
            if row["sequence"] not in source_by_sequence:
                raise ValueError(f"score sequence absent from source ranking: {row['source_name']}")
            scores.setdefault(row["sequence"], []).append(row)
    if set(scores) != set(source_by_sequence) or any(len(rows) != 8 for rows in scores.values()):
        raise ValueError("incomplete cross-backbone score matrix")

    excluded, exclusion_paths = load_exclusions(args.project, args.exclude_glob)
    rows: list[dict[str, object]] = []
    for sequence, source in source_by_sequence.items():
        if sequence in excluded:
            continue
        context_rows = scores[sequence]
        values = np.asarray([float(row["mean_design_nll"]) for row in context_rows])
        order_stds = np.asarray([float(row["std_order_nll"]) for row in context_rows])
        chain_spreads = np.asarray([float(row["mean_chain_spread"]) for row in context_rows])
        ai11 = np.asarray([
            float(row["mean_design_nll"]) for row in context_rows if row["family"] == "AI11_AF3"
        ])
        ova073 = np.asarray([
            float(row["mean_design_nll"]) for row in context_rows if row["family"] == "OVA073_ROSETTA"
        ])
        if ai11.size != 4 or ova073.size != 4:
            raise ValueError("expected four contexts in each backbone family")
        rows.append({
            "source_name": source["name"],
            "sequence": sequence,
            "n_mutations": int(source["n_mutations"]),
            "surface_mutation_fraction": float(source["surface_mutation_fraction"]),
            "mutations": source["mutations"],
            "predicted_af3_mean_iptm": float(source["predicted_af3_mean_iptm"]),
            "predicted_zero_clash_fraction": float(source["predicted_zero_clash_fraction"]),
            "af3_strict_surrogate_rank_score": float(source["af3_strict_surrogate_rank_score"]),
            "generation_average_rank": float(source["generation_average_rank"]),
            "mean_context_nll": float(values.mean()),
            "worst_context_nll": float(values.max()),
            "std_context_nll": float(values.std(ddof=0)),
            "mean_order_std_nll": float(order_stds.mean()),
            "mean_chain_spread_nll": float(chain_spreads.mean()),
            "ai11_af3_mean_nll": float(ai11.mean()),
            "ova073_rosetta_mean_nll": float(ova073.mean()),
            "backbone_family_gap_nll": float(abs(ai11.mean() - ova073.mean())),
            "context_nlls": ",".join(
                f"{row['context']}:{float(row['mean_design_nll']):.8f}"
                for row in sorted(context_rows, key=lambda item: item["context"])
            ),
        })
    if len(rows) < args.top:
        raise ValueError(f"only {len(rows)} novel candidates after historical exclusion")

    mean_rank = average_ranks([float(row["mean_context_nll"]) for row in rows])
    worst_rank = average_ranks([float(row["worst_context_nll"]) for row in rows])
    context_sd_rank = average_ranks([float(row["std_context_nll"]) for row in rows])
    family_gap_rank = average_ranks([float(row["backbone_family_gap_nll"]) for row in rows])
    ai_rank = average_ranks([float(row["af3_strict_surrogate_rank_score"]) for row in rows])
    mutation_rank = average_ranks([float(row["n_mutations"]) for row in rows])
    for row, mean_r, worst_r, sd_r, gap_r, ai_r, mut_r in zip(
        rows, mean_rank, worst_rank, context_sd_rank, family_gap_rank, ai_rank, mutation_rank
    ):
        robust_mpnn = 0.35 * mean_r + 0.45 * worst_r + 0.10 * sd_r + 0.10 * gap_r
        row["mpnn_mean_rank"] = mean_r
        row["mpnn_worst_rank"] = worst_r
        row["mpnn_context_sd_rank"] = sd_r
        row["mpnn_family_gap_rank"] = gap_r
        row["robust_mpnn_borda"] = robust_mpnn
        row["hamming_ai_rank"] = ai_r
        row["mutation_count_rank"] = mut_r
        row["combined_borda"] = 0.55 * robust_mpnn + 0.35 * ai_r + 0.10 * mut_r
    rows.sort(key=lambda row: (
        float(row["combined_borda"]), float(row["worst_context_nll"]),
        -float(row["predicted_af3_mean_iptm"]), int(row["n_mutations"]), str(row["sequence"]),
    ))
    for rank, row in enumerate(rows, 1):
        row["overall_rank"] = rank

    selected: list[dict[str, object]] = []
    used_threshold = args.minimum_distance
    for threshold in range(args.minimum_distance, 0, -1):
        for row in rows:
            if row in selected:
                continue
            if all(hamming(str(row["sequence"]), str(old["sequence"])) >= threshold for old in selected):
                selected.append(row)
            if len(selected) == args.top:
                used_threshold = threshold
                break
        if len(selected) == args.top:
            break
    if len(selected) != args.top:
        raise ValueError(f"selected only {len(selected)}/{args.top}")
    for rank, row in enumerate(selected, 1):
        row["e271_selection_rank"] = rank

    args.out.mkdir(parents=True, exist_ok=True)
    ranking_fields = ["overall_rank"] + [field for field in rows[0] if field != "overall_rank"]
    with (args.out / "multibackbone_exact_mpnn_ranking.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, ranking_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    shortlist_fields = ["e271_selection_rank"] + [
        field for field in selected[0] if field != "e271_selection_rank"
    ]
    with (args.out / "multibackbone_exact_mpnn_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, shortlist_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
    with (args.out / "multibackbone_exact_mpnn_shortlist.fasta").open("w") as handle:
        for row in selected:
            handle.write(
                f">OVA-C4-MBMPNN-{int(row['e271_selection_rank']):03d} "
                f"nmut={int(row['n_mutations'])} surface={float(row['surface_mutation_fraction']):.3f} "
                f"mean_nll={float(row['mean_context_nll']):.4f} "
                f"worst_nll={float(row['worst_context_nll']):.4f} "
                f"pred_af3={float(row['predicted_af3_mean_iptm']):.3f}\n"
                f"{row['sequence']}\n"
            )
    with (args.out / "exclusion_audit.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["excluded_unique_sequences", len(excluded)])
        writer.writerow(["source_candidates", len(source_rows)])
        writer.writerow(["novel_ranked_candidates", len(rows)])
        writer.writerow(["selected_candidates", len(selected)])
        writer.writerow(["minimum_hamming_distance_used", used_threshold])
        writer.writerow(["exclusion_files", len(exclusion_paths)])

    selected_sequences = {str(row["sequence"]) for row in selected}
    figure, axis = plt.subplots(figsize=(8.2, 6.2))
    scatter = axis.scatter(
        [float(row["mean_context_nll"]) for row in rows],
        [float(row["predicted_af3_mean_iptm"]) for row in rows],
        c=[float(row["backbone_family_gap_nll"]) for row in rows],
        s=18, alpha=0.55, cmap="viridis_r", linewidths=0,
    )
    chosen = [row for row in rows if str(row["sequence"]) in selected_sequences]
    axis.scatter(
        [float(row["mean_context_nll"]) for row in chosen],
        [float(row["predicted_af3_mean_iptm"]) for row in chosen],
        facecolors="none", edgecolors="#d62728", s=58, linewidths=1.2,
        label="E272 shortlist",
    )
    for row in selected[:8]:
        axis.annotate(
            str(row["e271_selection_rank"]),
            (float(row["mean_context_nll"]), float(row["predicted_af3_mean_iptm"])),
            fontsize=7, xytext=(3, 2), textcoords="offset points",
        )
    axis.set_xlabel("Eight-backbone mean exact ProteinMPNN NLL (lower is better)")
    axis.set_ylabel("AF3-labelled Hamming predicted mean ipTM")
    axis.set_title("E271 cross-backbone exact-sequence consensus")
    axis.legend(loc="best", frameon=False)
    colorbar = figure.colorbar(scatter, ax=axis)
    colorbar.set_label("|AI11-AF3 mean NLL − OVA073-Rosetta mean NLL|")
    figure.tight_layout()
    figure.savefig(args.out / "multibackbone_exact_mpnn_selection.png", dpi=180)
    figure.savefig(args.out / "multibackbone_exact_mpnn_selection.pdf")
    print(
        f"ranked={len(rows)} excluded={len(source_rows) - len(rows)} selected={len(selected)} "
        f"minimum_distance={used_threshold}"
    )


if __name__ == "__main__":
    main()
