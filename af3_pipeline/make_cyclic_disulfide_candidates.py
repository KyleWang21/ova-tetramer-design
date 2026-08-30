#!/usr/bin/env python3
"""Create pure-OVA cyclic-disulfide candidates from a C4 interface geometry."""

from __future__ import annotations

import argparse
import csv
import pathlib


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:], []
        else:
            parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def mutate(sequence: str, replacements: dict[int, str]) -> str:
    output = list(sequence)
    for position, aa in replacements.items():
        output[position - 1] = aa
    return "".join(output)


def mutation_text(parent: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(parent, sequence), 1) if old != new
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent-manifest", type=pathlib.Path, required=True)
    ap.add_argument("--ai-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    parent_rows = list(csv.DictReader(args.parent_manifest.open(), delimiter="\t"))
    if len(parent_rows) != 1:
        raise SystemExit("parent manifest must contain exactly one sequence")
    parent = parent_rows[0]["sequence"]
    ai = read_fasta(args.ai_fasta)
    if len(parent) != 386 or len(ai) < 3:
        raise SystemExit("expected a 386-aa parent and at least three AI candidates")

    # Directed C4 edges in the AF3 clamp geometry give the following C-beta
    # distances (mean across A->B->C->D->A): 94->213 3.66 A,
    # 96->190 4.01 A, and 155->196 4.47 A.  Each chain carries both
    # residues, so one site bonds to the next chain and the other receives a
    # bond from the previous chain, closing a four-member ring.
    pair_94_213 = {94: "C", 213: "C"}
    pair_96_190 = {96: "C", 190: "C"}
    pair_155_196 = {155: "C", 196: "C"}
    double_pair = {**pair_96_190, **pair_155_196}
    specifications = [
        ("parent_C94_C213", parent, pair_94_213, "CB distance 3.66 A; includes P94C"),
        ("parent_C96_C190", parent, pair_96_190, "CB distance 4.01 A"),
        ("parent_C155_C196", parent, pair_155_196, "CB distance 4.47 A"),
        ("parent_double_ring", parent, double_pair, "two independent directed C4 disulfide pairs"),
        ("ai1_C96_C190", ai[0][1], pair_96_190, "multistate AI rank 1 plus directed disulfide"),
        ("ai1_C155_C196", ai[0][1], pair_155_196, "multistate AI rank 1 plus directed disulfide"),
        ("ai3_C96_C190", ai[2][1], pair_96_190, "multistate AI rank 3 plus directed disulfide"),
        ("ai3_C155_C196", ai[2][1], pair_155_196, "multistate AI rank 3 plus directed disulfide"),
    ]
    rows: list[dict[str, object]] = []
    with (args.out / "cyclic_disulfide_candidates.fasta").open("w") as fasta:
        for index, (name, source, replacements, rationale) in enumerate(specifications, 1):
            sequence = mutate(source, replacements)
            if any(sequence[position - 1] != parent[position - 1] for position in (14, 69, 77, 79, 191, 192, 193)):
                raise ValueError(f"{name} changed a protected P3-13R mutation")
            mutations = mutation_text(parent, sequence)
            fasta.write(f">OVA-BODY-SS-{index:02d}|{name}|mut={mutations}\n{sequence}\n")
            rows.append({
                "rank": index,
                "name": name,
                "rationale": rationale,
                "mutations_vs_p3": mutations,
                "n_mutations": len(mutations.split(",")),
                "sequence": sequence,
            })
    with (args.out / "cyclic_disulfide_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(row["rank"], row["name"], row["n_mutations"], row["mutations_vs_p3"])


if __name__ == "__main__":
    main()
