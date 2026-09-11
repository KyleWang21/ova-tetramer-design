#!/usr/bin/env python3
"""Prepare tied-AF2 targets from strict AF3 seed/sample tables.

This adapter is for promising systems that have not passed the full 25-model
acceptance screen yet.  It accepts multiple independent strict tables, keeps
only geometry-clean C4 models, and exposes known clash positions to the next
round of sequence design.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re

from Bio.PDB import MMCIFParser, PDBIO

from prepare_hit_ensemble_tied_design import PROTECTED
from prepare_iptm90_tied_design import Chains, interface_positions, write_positions


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: pathlib.Path, name: str) -> str:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header = raw[1:].split()[0]
            chunks = []
        elif header is not None:
            chunks.append(raw.strip())
    if header is not None:
        records[header] = "".join(chunks)
    if name not in records:
        raise ValueError(f"{name!r} not found in {path}; available={sorted(records)}")
    return records[name]


def parse_positions(raw: str) -> set[int]:
    return {int(value) for value in re.split(r"[\s,]+", raw.strip()) if value}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--per-model", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--system", action="append", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--pairs", default="AC,AD")
    parser.add_argument("--backbones", type=int, default=4)
    parser.add_argument("--anchor-prefix", required=True)
    parser.add_argument("--extra-design-positions", default="")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    if not (len(args.summary) == len(args.per_model) == len(args.system)):
        raise ValueError("--summary, --per-model and --system must have equal counts")
    pairs = tuple(value.strip() for value in args.pairs.split(","))
    if len(pairs) != 2 or any(len(pair) != 2 for pair in pairs):
        raise ValueError("--pairs must contain two two-chain identifiers, e.g. AC,AD")

    sequences: list[str] = []
    models: list[dict[str, str]] = []
    for summary_path, model_path, system in zip(
        args.summary, args.per_model, args.system, strict=True
    ):
        summary_rows = [row for row in read_tsv(summary_path) if row["system"] == system]
        if len(summary_rows) != 1:
            raise ValueError(f"{summary_path}: expected one row for {system}, found {len(summary_rows)}")
        sequences.append(summary_rows[0]["sequence"])
        for row in read_tsv(model_path):
            if row["system"] != system:
                continue
            if not (
                int(row["atomic_clash_count_2p4"]) == 0
                and int(row["severe_ca_overlap_count_2p5"]) == 0
                and int(row["valid_c4_network"]) == 1
                and int(row["native_c74_c121_bonds"]) == 4
                and int(row["interchain_ss_bonds"]) == 0
                and int(row["two_interfaces_each_design_ge3"]) == 1
            ):
                continue
            match = re.fullmatch(r"seed-(\d+)_sample-(\d+)", row["model"])
            if not match:
                raise ValueError(f"cannot parse seed/sample from {row['model']!r}")
            models.append({**row, "seed": match.group(1), "sample": match.group(2)})
    if len(set(sequences)) != 1:
        raise ValueError("strict tables describe different sequences")
    sequence = sequences[0]
    reference = read_fasta(args.reference_fasta, args.reference_name)
    if len(sequence) != len(reference):
        raise ValueError(f"sequence/reference length mismatch: {len(sequence)} vs {len(reference)}")
    mutations = {
        index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
        if old != new
    }
    if len(mutations) > 24:
        raise ValueError(f"anchor exceeds mutation budget: {len(mutations)}")
    surface = parse_positions(args.surface_positions.read_text())
    extra = parse_positions(args.extra_design_positions)
    if len(models) < args.backbones:
        raise ValueError(f"only {len(models)} clean models for {args.backbones} backbones")

    # Maximize independent-seed coverage first, then fill with the highest-ipTM
    # remaining samples.  This avoids choosing four nearly duplicate samples
    # when an independent seed has clean representatives.
    chosen: list[dict[str, str]] = []
    for seed in sorted({int(row["seed"]) for row in models}):
        candidates = [row for row in models if int(row["seed"]) == seed]
        chosen.append(max(candidates, key=lambda row: float(row["iptm"])))
    chosen = sorted(chosen, key=lambda row: -float(row["iptm"]))[: args.backbones]
    for row in sorted(models, key=lambda row: -float(row["iptm"])):
        if row not in chosen:
            chosen.append(row)
        if len(chosen) == args.backbones:
            break

    args.out.mkdir(parents=True)
    target_dir = args.out / "targets"
    target_dir.mkdir()
    write_positions(args.out / "surface_positions.txt", surface)
    write_positions(args.out / "known_clash_positions.txt", extra)
    records: list[dict[str, object]] = []
    with (args.out / "anchors.fasta").open("w") as fasta:
        for index, row in enumerate(chosen, 1):
            anchor = f"{args.anchor_prefix}-BB{index}"
            fasta.write(f">{anchor}\n{sequence}\n")
            structure = MMCIFParser(QUIET=True).get_structure(anchor, row["cif_path"])
            model = structure[0]
            io = PDBIO()
            io.set_structure(structure)
            stem = anchor.lower().replace("-", "_")
            target_c4 = target_dir / f"{stem}_c4.pdb"
            io.save(str(target_c4), Chains(set("ABCD")))
            record: dict[str, object] = {
                "anchor": anchor,
                "anchor_n_mutations": len(mutations),
                "source_cif": row["cif_path"],
                "source_iptm": row["iptm"],
                "source_clash": row["atomic_clash_count_2p4"],
                "pair1": pairs[0],
                "pair2": pairs[1],
                "target_c4": str(target_c4),
            }
            contact_sets: list[set[int]] = []
            for state_index, pair in enumerate(pairs, 1):
                prefix = f"{stem}_{pair.lower()}"
                target = target_dir / f"{prefix}.pdb"
                io.save(str(target), Chains(set(pair)))
                contacts = interface_positions(model, pair)
                design = (mutations | extra | (contacts & surface)) - PROTECTED
                restraint = contacts - PROTECTED
                design_path = target_dir / f"{prefix}_design.txt"
                restraint_path = target_dir / f"{prefix}_restraint.txt"
                write_positions(design_path, design)
                write_positions(restraint_path, restraint)
                contact_sets.append(contacts)
                record.update({
                    f"target{state_index}": str(target),
                    f"design{state_index}": str(design_path),
                    f"restraint{state_index}": str(restraint_path),
                    f"contact_positions{state_index}": len(contacts),
                    f"design_positions{state_index}": len(design),
                })
            joint = (
                mutations | extra | ((contact_sets[0] | contact_sets[1]) & surface)
            ) - PROTECTED
            joint_path = target_dir / f"{stem}_joint_design.txt"
            write_positions(joint_path, joint)
            record["joint_design"] = str(joint_path)
            record["joint_design_positions"] = len(joint)
            records.append(record)

    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    mutation_text = ",".join(
        f"{reference[position - 1]}{position}{sequence[position - 1]}"
        for position in sorted(mutations)
    )
    (args.out / "anchor_mutations.txt").write_text(mutation_text + "\n")
    print(
        f"prepared {len(records)} strict clean backbones; seeds/samples="
        f"{[(row['seed'], row['sample']) for row in chosen]} mutations={len(mutations)} "
        f"extra_design={sorted(extra)}"
    )


if __name__ == "__main__":
    main()
