#!/usr/bin/env python3
"""Freeze exact dropout-free E283 endpoints that have never reached AF3.

The ordinary generation collector deliberately excludes every sequence seen in
an older candidate pool.  That is useful for library novelty, but it can hide a
previously generated sequence when new multibackbone evidence makes it newly
informative.  E292 is a small mechanistic quota for those exact endpoints; it
still excludes anything already scored by AF3 or already present in an AF3
manifest.
"""

from __future__ import annotations

import csv
import math
import pathlib

from run_final20_failclosed_filter import glyco_sites, reference_rsasa


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_sequence(path: pathlib.Path, prefix: str) -> str:
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None and name.startswith(prefix):
                return "".join(chunks)
            name, chunks = line[1:], []
        elif name is not None:
            chunks.append(line.strip())
    raise ValueError(f"missing FASTA record {prefix} in {path}")


def sequences_from(paths: list[pathlib.Path]) -> set[str]:
    sequences: set[str] = set()
    for path in paths:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            if len(sequence) == 386:
                sequences.add(sequence)
    return sequences


def main() -> None:
    project = pathlib.Path(__file__).resolve().parents[1]
    source = project / "experiments/e283_c4_chainmapped_eightstate_tied_design"
    out = project / "experiments/e292_c4_e283_deterministic_final_af3_seed1"
    pool = source / "candidate_pool_deterministic_final"
    if (out / "CANDIDATES_FROZEN").exists():
        print(f"{out} already frozen")
        return
    if not (source / "REMOTE_ARCHIVE_RECOVERED").exists():
        raise RuntimeError("E283 archive is not recovered")
    if out.exists() or pool.exists():
        raise RuntimeError("refusing to overwrite partial E292 preparation")

    reference = fasta_sequence(project / "OVA_P3-13R_四聚体候选_AA.fasta", "D0-P3-13R")
    rsasa = reference_rsasa(project / "references/1OVA.cif", "A", reference)

    scored_patterns = (
        "experiments/e*/body_candidate_ranking.tsv",
        "experiments/e*/**/seed1_strict_summary.tsv",
        "experiments/e*/**/full25_summary.tsv",
        "experiments/e*/**/combined_25model_summary.tsv",
        "experiments/e*/**/af3_candidate_recomputed_summary_with_sequence.tsv",
        "experiments/current_accepted_af3_registry/accepted_af3_registry.tsv",
    )
    scored_paths = sorted({path for pattern in scored_patterns for path in project.glob(pattern)})
    scored = sequences_from(scored_paths)
    manifest_paths = sorted(project.glob("experiments/e*/**/manifest.tsv"))
    scheduled = sequences_from(manifest_paths)

    rows: list[dict[str, object]] = []
    for path in sorted(source.glob("run_chainmapped_*/**/af2diff_candidates.tsv")):
        for raw in read_tsv(path):
            sequence = raw["sequence"]
            changed = [
                index
                for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
                if old != new
            ]
            surface_fraction = sum(rsasa.get(index, 0.0) >= 0.20 for index in changed) / len(changed)
            state_iptm = [float(raw[f"state{index}_iptm"]) for index in range(1, 9)]
            state_plddt = [float(raw[f"state{index}_plddt"]) for index in range(1, 9)]
            reasons: list[str] = []
            if not 1 <= len(changed) <= 20:
                reasons.append("mutation_budget")
            if surface_fraction < 0.80:
                reasons.append("surface_fraction")
            if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
                reasons.append("cysteine")
            if any(sequence[index - 1] != reference[index - 1] for index in EPITOPE):
                reasons.append("SIINFEKL")
            if glyco_sites(sequence) != glyco_sites(reference):
                reasons.append("glycosylation")
            if sequence in scored:
                reasons.append("already_af3_scored")
            if sequence in scheduled:
                reasons.append("already_af3_scheduled")
            rows.append({
                "source": str(path.relative_to(project)),
                "candidate": raw["candidate"],
                "selected_stage": raw["selected_stage"],
                "selected_iteration": raw["selected_iteration"],
                "n_mutations": len(changed),
                "surface_mutation_fraction": surface_fraction,
                "mutations": raw["mutations"],
                "deterministic_min_iptm": min(state_iptm),
                "deterministic_mean_iptm": sum(state_iptm) / len(state_iptm),
                "deterministic_min_plddt": min(state_plddt),
                "prior_af3_scored": int(sequence in scored),
                "prior_af3_scheduled": int(sequence in scheduled),
                "eligible": int(not reasons),
                "exclusion_reasons": ";".join(reasons),
                "sequence": sequence,
            })

    if len(rows) != 8:
        raise RuntimeError(f"expected eight deterministic E283 endpoints, found {len(rows)}")
    rows.sort(key=lambda row: (
        -int(row["eligible"]),
        -float(row["deterministic_min_iptm"]),
        -float(row["deterministic_mean_iptm"]),
        -float(row["deterministic_min_plddt"]),
        int(row["n_mutations"]),
        str(row["sequence"]),
    ))
    selected = [row for row in rows if int(row["eligible"])]
    if not selected:
        raise RuntimeError("no untested deterministic E283 endpoint remains for E292")

    pool.mkdir(parents=True)
    write_tsv(pool / "selection_audit.tsv", rows)
    write_tsv(pool / "deterministic_final_shortlist.tsv", selected)
    with (pool / "deterministic_final_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            handle.write(
                f">E292-DET-{index:02d} nmut={row['n_mutations']} "
                f"min8={float(row['deterministic_min_iptm']):.3f} "
                f"mean8={float(row['deterministic_mean_iptm']):.3f}\n"
                f"{row['sequence']}\n"
            )
    print(
        f"E292 deterministic endpoints total={len(rows)} eligible={len(selected)} "
        f"best_min8={max(float(row['deterministic_min_iptm']) for row in selected):.3f}"
    )


if __name__ == "__main__":
    main()
