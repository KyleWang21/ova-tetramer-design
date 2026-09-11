#!/usr/bin/env python3
"""Build one snapshot table for all currently traceable OVA C4 samples.

The table deliberately separates generation/AF2 evidence from template-free
AF3 acceptance.  It includes the current Accepted_AF3 registry, the frozen
E312 seed-1 queue, every valid unique E318 partial-pool sequence, and E321
anchor roles.  Historical negative controls are included when their sequence
and summary are available.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
from collections import defaultdict


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path: pathlib.Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    result: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                result.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        result.append((header, "".join(chunks)))
    return result


def sequence_from_fasta(path: pathlib.Path, prefix: str) -> str:
    for header, sequence in fasta_records(path):
        if header.startswith(prefix):
            return sequence
    return ""


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def safe_float(value: object, default: str = "") -> str:
    try:
        return f"{float(value):.4f}"
    except (TypeError, ValueError):
        return default


def relpath(value: str, root: pathlib.Path) -> str:
    if not value:
        return ""
    try:
        return str(pathlib.Path(value).resolve().relative_to(root.resolve()))
    except ValueError:
        return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    reference_path = root / "OVA_P3-13R_四聚体候选_AA.fasta"
    reference = sequence_from_fasta(reference_path, "D0-P3-13R")
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa reference, found {len(reference)}")

    fields = [
        "sample_id", "af3_rank", "generation_rank", "e321_task", "aliases",
        "source_batch", "role", "stage",
        "n_mutations", "surface_mutation_fraction", "mutations_vs_p3",
        "af3_status", "af3_models", "af3_seed_count", "af3_samples_per_seed",
        "af3_mean_iptm", "af3_min_iptm", "af3_interface_match",
        "v1_status", "v1_failed_required_stages", "af2_generation_iptm",
        "predicted_af3_mean_iptm", "predicted_zero_clash_fraction",
        "source_detail", "cif", "sequence", "notes",
    ]
    by_sequence: dict[str, dict[str, str]] = {}

    def add(row: dict[str, str]) -> None:
        sequence = row.get("sequence", "")
        if not sequence:
            return
        old = by_sequence.get(sequence)
        if old is None:
            by_sequence[sequence] = {field: row.get(field, "") for field in fields}
            return
        for field in fields:
            value = row.get(field, "")
            if not value:
                continue
            if field in {"aliases", "role", "source_batch", "stage", "notes", "source_detail"}:
                old_value = old.get(field, "")
                values = [item for item in (old_value + ";" + value).split(";") if item]
                old[field] = ";".join(dict.fromkeys(values))
            elif not old.get(field):
                old[field] = value

    # Current registry: authoritative complete 5-seed x 5-sample AF3 rows.
    registry = root / "experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
    for item in read_tsv(registry):
        complete = item.get("v1_screen_complete") == "1"
        accepted_v1 = item.get("accepted_v1") == "1"
        if accepted_v1:
            v1_status = "Accepted_v1.1"
        elif complete:
            v1_status = "v1.1 complete—failed required gates"
        else:
            v1_status = "v1.1 pending"
        add({
            "sample_id": item.get("candidate", ""),
            "af3_rank": item.get("rank", ""),
            "aliases": item.get("candidate", ""),
            "source_batch": pathlib.Path(item.get("source_experiment", "")).name,
            "role": "Accepted_AF3",
            "stage": "AF3 complete",
            "n_mutations": item.get("n_mutations", ""),
            "surface_mutation_fraction": safe_float(item.get("surface_mutation_fraction")),
            "mutations_vs_p3": "",
            "af3_status": "Accepted_AF3 (25/25)",
            "af3_models": item.get("n_models", "25"),
            "af3_seed_count": "5",
            "af3_samples_per_seed": "5",
            "af3_mean_iptm": safe_float(item.get("mean_iptm")),
            "af3_min_iptm": safe_float(item.get("min_iptm")),
            "af3_interface_match": "yes (25/25)" if item.get("all_model_design_interface_position_count") else "yes",
            "v1_status": v1_status,
            "v1_failed_required_stages": item.get("failed_required_stages", ""),
            "source_detail": relpath(item.get("source_experiment", ""), root),
            "cif": relpath(item.get("representative_af3_cif", ""), root),
            "sequence": item.get("sequence", ""),
            "notes": "zero cross-chain clash/C4/two designed interfaces recorded in registry",
        })

    # E312: eleven candidates frozen for AF3 seed1, no AF3 models yet.
    e312_manifest = root / "experiments/e312_c4_e311_detbest_af3_seed1/manifest.tsv"
    for item in read_tsv(e312_manifest):
        description = item.get("description", "")
        pred = re.search(r"predicted_af3=([0-9.]+)", description)
        zero = re.search(r"predicted_zero_clash=([0-9.]+)", description)
        source_rank = re.search(r"source_ai_rank=([0-9]+)", description)
        sequence = item.get("sequence", "")
        nmut = sum(a != b for a, b in zip(reference, sequence))
        add({
            "sample_id": item.get("name", ""),
            "generation_rank": source_rank.group(1) if source_rank else "",
            "aliases": description.split(" source_")[0],
            "source_batch": "E312",
            "role": "frozen AF3 seed1 candidate",
            "stage": "AF3 seed1 queued",
            "n_mutations": str(nmut),
            "surface_mutation_fraction": "",
            "mutations_vs_p3": mutation_string(reference, sequence),
            "af3_status": "not started (seed1 pending Volc free>10)",
            "af3_models": "0/5",
            "af3_seed_count": "1 planned",
            "af3_samples_per_seed": "5",
            "v1_status": "not started",
            "af2_generation_iptm": "",
            "predicted_af3_mean_iptm": pred.group(1) if pred else "",
            "predicted_zero_clash_fraction": zero.group(1) if zero else "",
            "source_detail": relpath(item.get("input_json", ""), root),
            "sequence": sequence,
            "notes": "frozen 11-candidate pool; AF3 prediction fields are surrogate only",
        })

    # E318 pool: prefer the complete eight-shard archive when available; keep
    # the seven-shard partial pool as a live fallback while the last shard runs.
    e318_root = root / "experiments/e318_c4_e311_diverse_tied_design"
    complete_pool = e318_root / "candidate_pool_novel"
    partial_pool = e318_root / "candidate_pool_partial7"
    e318_pool = complete_pool if (complete_pool / "generation_shortlist.tsv").exists() else partial_pool
    e318_label = "E318 (8/8 shards)" if e318_pool == complete_pool else "E318 (7/8 shards)"
    e318_ranked = {
        item.get("name", ""): item
        for item in read_tsv(e318_pool / "hamming_ai_ranking_current_weighted.tsv")
    }
    for item in read_tsv(e318_pool / "generation_shortlist.tsv"):
        name = item.get("name", "")
        rank = e318_ranked.get(name, {})
        sequence = item.get("sequence", "")
        add({
            "sample_id": name,
            "generation_rank": rank.get("ai_rank", ""),
            "aliases": name,
            "source_batch": e318_label,
            "role": "valid unique generation candidate",
            "stage": "tied-AF2 generation complete",
            "n_mutations": item.get("n_mutations", ""),
            "surface_mutation_fraction": safe_float(item.get("surface_mutation_fraction")),
            "mutations_vs_p3": item.get("mutations", ""),
            "af3_status": "not started (awaiting full E318 archive)",
            "af3_models": "0/25",
            "af3_seed_count": "5 planned",
            "af3_samples_per_seed": "5",
            "v1_status": "not started",
            "af2_generation_iptm": safe_float(item.get("af2_iptm")),
            "predicted_af3_mean_iptm": safe_float(rank.get("predicted_af3_mean_iptm")),
            "predicted_zero_clash_fraction": safe_float(rank.get("predicted_zero_clash_fraction")),
            "source_detail": item.get("source", ""),
            "sequence": sequence,
            "notes": f"{len(read_tsv(e318_pool / 'generation_shortlist.tsv'))} valid unique sequences from E318 archive; not AF3 evidence",
        })

    # E321 anchor roles; merge by sequence so the total table remains one row/sample.
    e321_manifest = root / "experiments/e321_c4_e318_next_tied_design/anchor_manifest.tsv"
    for item in read_tsv(e321_manifest):
        sequence = item.get("sequence", "")
        task = int(item.get("task_index", "0"))
        add({
            "sample_id": item.get("anchor", ""),
            "e321_task": str(task),
            "aliases": item.get("anchor", ""),
            "source_batch": "E321",
            "role": f"E321 anchor / task {task}",
            "stage": "tied-AF2 running" if task < 7 else "tied-AF2 pending GPU",
            "n_mutations": item.get("n_mutations", ""),
            "surface_mutation_fraction": safe_float(item.get("surface_mutation_fraction")),
            "mutations_vs_p3": item.get("mutations", ""),
            "af3_status": "not started",
            "af3_models": "0/25",
            "af3_seed_count": "5 planned",
            "af3_samples_per_seed": "5",
            "v1_status": "not started",
            "af2_generation_iptm": safe_float(item.get("score")),
            "source_detail": item.get("source", ""),
            "sequence": sequence,
            "notes": "anchor role merged with the same sequence if already present in E318/Accepted_AF3",
        })

    # Historical key controls, included only when sequence records are present.
    controls = [
        ("P3-13R_C4_baseline", "FINAL_CANDIDATE_COMPARISON.tsv", "negative control"),
        ("OVA-BODY-C4-C3HOT-AI01", "OVA-BODY-C4-C3HOT-AI01_AA.fasta", "legacy C4-forming control"),
        ("OVA-BODY-C4-SS3-AI03", "OVA-BODY-C4-SS3-AI03_AA.fasta", "legacy C4-forming control"),
        ("OVA-BODY-C4-P208F-AI01", "experiments/e42_ai_p208f_reseed/manifest.tsv", "legacy C4-forming control"),
    ]
    comparison = {
        item.get("candidate", ""): item
        for item in read_tsv(root / "FINAL_CANDIDATE_COMPARISON.tsv")
    }
    for name, source, role in controls:
        item = comparison.get(name, {})
        if name == "P3-13R_C4_baseline":
            sequence = reference
        elif source.endswith("manifest.tsv"):
            records = read_tsv(root / source)
            sequence = records[0].get("sequence", "") if records else ""
        else:
            sequence = sequence_from_fasta(root / source, name)
        if not sequence:
            continue
        add({
            "sample_id": name,
            "aliases": name,
            "source_batch": "historical baseline/control",
            "role": role,
            "stage": "historical AF3 result",
            "n_mutations": item.get("mutations_vs_p3", "0"),
            "af3_status": "negative control / legacy result",
            "af3_models": item.get("c4_models", ""),
            "af3_mean_iptm": item.get("median_c4_iptm", ""),
            "af3_min_iptm": "",
            "v1_status": "not in current v1 registry",
            "source_detail": source,
            "sequence": sequence,
            "notes": item.get("decision", "historical control") + "; current surface-mask fraction not recalculated",
        })

    output = list(by_sequence.values())
    surface_file = root / "experiments/e318_c4_e311_diverse_tied_design/surface_positions.txt"
    surface_positions = {
        int(token)
        for token in surface_file.read_text().replace("\n", ",").split(",")
        if token.strip()
    } if surface_file.exists() else set()
    for row in output:
        sequence = row.get("sequence", "")
        if sequence and not row.get("mutations_vs_p3"):
            row["mutations_vs_p3"] = mutation_string(reference, sequence)
        if sequence and not row.get("n_mutations"):
            row["n_mutations"] = str(sum(a != b for a, b in zip(reference, sequence)))
        if (
            sequence
            and row.get("source_batch") != "historical baseline/control"
            and not row.get("surface_mutation_fraction")
            and surface_positions
        ):
            changed = [
                index for index, (old, new) in enumerate(zip(reference, sequence), 1)
                if old != new
            ]
            if changed:
                row["surface_mutation_fraction"] = f"{sum(index in surface_positions for index in changed) / len(changed):.4f}"
    def category(row: dict[str, str]) -> int:
        role = row.get("role", "")
        if "Accepted_AF3" in role:
            return 0
        if "frozen AF3 seed1 candidate" in role:
            return 1
        if "valid unique generation candidate" in role:
            return 2
        return 3

    def numeric(value: str, default: int = 10**9) -> int:
        try:
            return int(value)
        except (TypeError, ValueError):
            return default

    output.sort(key=lambda row: (
        category(row), numeric(row.get("af3_rank", "")),
        numeric(row.get("generation_rank", "")), row["sample_id"],
    ))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    print(f"wrote {len(output)} unique sequences and {len(fields)} columns to {args.out}")


if __name__ == "__main__":
    main()
