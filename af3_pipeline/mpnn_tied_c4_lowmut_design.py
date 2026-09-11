#!/usr/bin/env python3
"""Tied ProteinMPNN sampling on a complete OVA C4 AF3 backbone.

The four chains share one sequence.  Only the union of the observed two-interface
positions and the parent mutations is designable; all remaining coordinates are
fixed.  A tunable P3-13R residue bias lets ProteinMPNN explore the AF3 geometry
without automatically changing every open position.  Output sequences are only
generation proposals and still require template-free four-chain AF3 validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
from collections import Counter

import numpy as np
from Bio.PDB import PDBParser
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.mpnn.model import aa_order


AA20 = "ARNDCQEGHILKMFPSTWYV"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if header is not None and header.split()[0] == name:
                return "".join(chunks)
            header, chunks = raw[1:], []
        elif header is not None:
            chunks.append(raw.strip())
    if header is not None and header.split()[0] == name:
        return "".join(chunks)
    raise ValueError(f"FASTA record {name!r} not found in {path}")


def positions(path: pathlib.Path) -> set[int]:
    text = path.read_text().replace(",", " ").replace("\n", " ")
    return {int(token) for token in text.split() if token}


def glyco_sites(sequence: str) -> set[int]:
    return {
        index + 1
        for index in range(len(sequence) - 2)
        if sequence[index] == "N"
        and sequence[index + 1] != "P"
        and sequence[index + 2] in "ST"
    }


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def compress_positions(chain: str, selected: set[int], length: int) -> str:
    values = sorted(selected)
    if not values:
        return ""
    chunks: list[str] = []
    start = previous = values[0]
    for value in values[1:] + [length + 2]:
        if value != previous + 1:
            chunks.append(
                f"{chain}{start}" if start == previous else f"{chain}{start}-{chain}{previous}"
            )
            start = value
        previous = value
    return ",".join(chunks)


def normalize_sampled_sequence(sequence: str, length: int) -> str:
    sequence = sequence.replace("/", "").replace(":", "")
    if len(sequence) == length:
        return sequence
    if len(sequence) == 4 * length:
        chains = [sequence[index * length:(index + 1) * length] for index in range(4)]
        if len(set(chains)) != 1:
            raise ValueError("homooligomer ProteinMPNN returned non-identical chains")
        return chains[0]
    raise ValueError(f"unexpected sampled sequence length {len(sequence)}")


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--design-positions", type=pathlib.Path, required=True)
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--temperatures", default="0.08,0.12,0.20")
    parser.add_argument("--native-biases", default="0.5,1.0,1.5,2.0")
    parser.add_argument(
        "--buried-native-bias", type=float, default=1.0,
        help="additional P3-reference logit preference at buried designable positions",
    )
    parser.add_argument("--samples-per-condition", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--min-mutations", type=int, default=12)
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--min-surface-fraction", type=float, default=0.80)
    parser.add_argument("--shortlist", type=int, default=128)
    parser.add_argument("--polar-bias-positions", type=pathlib.Path)
    parser.add_argument("--polar-bias", type=float, default=0.0)
    parser.add_argument("--polar-aas", default="DEHKNQRSTY")
    args = parser.parse_args()

    if args.samples_per_condition % args.batch_size:
        raise ValueError("samples-per-condition must be divisible by batch-size")
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa reference, found {len(reference)}")
    structure = PDBParser(QUIET=True).get_structure("c4", str(args.target_pdb))[0]
    chains = [chain for chain in structure if chain.id in "ABCD"]
    if [chain.id for chain in chains] != list("ABCD"):
        raise ValueError("target must contain chains A,B,C,D")
    parent_chains = [
        "".join(AA3[residue.resname] for residue in chain if residue.id[0] == " " and residue.resname in AA3)
        for chain in chains
    ]
    if any(len(sequence) != len(reference) for sequence in parent_chains):
        raise ValueError(f"target chain lengths are {[len(sequence) for sequence in parent_chains]}")
    if len(set(parent_chains)) != 1:
        raise ValueError("target is not a sequence-identical homotetramer")
    parent = parent_chains[0]

    requested = positions(args.design_positions)
    parent_mutations = {
        index for index, (old, new) in enumerate(zip(reference, parent), 1) if old != new
    }
    designable = (requested | parent_mutations) - NATIVE_CYS - EPITOPE
    if not designable or min(designable) < 1 or max(designable) > len(reference):
        raise ValueError("invalid design position set")
    fixed = set(range(1, len(reference) + 1)) - designable
    fix_pos = ",".join(compress_positions(chain, fixed, len(reference)) for chain in "ABCD")
    surface = positions(args.surface_positions)
    polar_positions = (
        positions(args.polar_bias_positions) if args.polar_bias_positions else set()
    )
    polar_aas = set(args.polar_aas) - {"C"}
    if not polar_aas <= set(AA20):
        raise ValueError(f"invalid polar amino-acid set: {sorted(polar_aas - set(AA20))}")
    if not polar_positions <= designable:
        raise ValueError(
            "polar-bias positions must be designable; outside mask="
            f"{sorted(polar_positions - designable)}"
        )

    if [key for key, _ in sorted(aa_order.items(), key=lambda item: item[1])] != list(AA20):
        raise RuntimeError("unexpected ProteinMPNN amino-acid alphabet")
    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, seed=args.seed,
        verbose=False, weights="original",
    )
    model.prep_inputs(
        pdb_filename=str(args.target_pdb), chain="A,B,C,D", homooligomer=True,
        fix_pos=fix_pos, rm_aa="C", verbose=True,
    )
    baseline_bias = np.asarray(model._inputs["bias"]).copy()
    bias_length = int(baseline_bias.shape[0])
    if bias_length not in {len(reference), 4 * len(reference)}:
        raise ValueError(f"unexpected ProteinMPNN bias length {bias_length}")

    args.out.mkdir(parents=True, exist_ok=True)
    temperatures = [float(value) for value in args.temperatures.split(",")]
    native_biases = [float(value) for value in args.native_biases.split(",")]
    metadata = {
        "method": "tied_c4_ProteinMPNN",
        "target_pdb": str(args.target_pdb),
        "reference": str(args.reference_fasta),
        "parent_mutations": sorted(parent_mutations),
        "requested_design_positions": sorted(requested),
        "actual_design_positions": sorted(designable),
        "fixed_positions_per_chain": len(fixed),
        "native_bias_definition": (
            "add temperature*native_bias to the P3 residue logit, giving an approximately "
            "exp(native_bias) preference after temperature scaling"
        ),
        "buried_native_bias": args.buried_native_bias,
        "polar_bias_positions": sorted(polar_positions),
        "polar_bias": args.polar_bias,
        "polar_aas": "".join(sorted(polar_aas)),
        "temperatures": temperatures,
        "native_biases": native_biases,
        "samples_per_condition": args.samples_per_condition,
        "expected_raw_samples": (
            len(temperatures) * len(native_biases) * args.samples_per_condition
        ),
    }
    (args.out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    rejection = Counter()
    reference_glyco = glyco_sites(reference)
    direct_new_glyco_asn_positions = {
        position
        for position in designable
        if position <= len(reference) - 2
        and reference[position] != "P"
        and reference[position + 1] in "ST"
        and position not in reference_glyco
    }
    metadata["direct_new_glyco_asn_positions_blocked"] = sorted(
        direct_new_glyco_asn_positions
    )
    (args.out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    condition = 0
    for temperature in temperatures:
        for native_bias in native_biases:
            condition += 1
            bias = baseline_bias.copy()
            for raw_index in range(bias_length):
                position = raw_index % len(reference) + 1
                if position in designable:
                    bias[raw_index, aa_order[reference[position - 1]]] += temperature * native_bias
                    if position not in surface:
                        bias[raw_index, aa_order[reference[position - 1]]] += (
                            temperature * args.buried_native_bias
                        )
                if position in direct_new_glyco_asn_positions:
                    bias[raw_index, aa_order["N"]] -= 1.0e4
                if position in polar_positions and args.polar_bias:
                    for amino_acid in polar_aas:
                        bias[raw_index, aa_order[amino_acid]] += temperature * args.polar_bias
            model._inputs["bias"] = bias
            for batch_index in range(args.samples_per_condition // args.batch_size):
                sampled = model.sample_parallel(
                    batch=args.batch_size, temperature=temperature, rescore=True,
                )
                for sample_index, (raw_sequence, score, seqid) in enumerate(
                    zip(sampled["seq"], sampled["score"], sampled["seqid"]), 1
                ):
                    sequence = normalize_sampled_sequence(str(raw_sequence), len(reference))
                    changed = {
                        index for index, (old, new) in enumerate(zip(reference, sequence), 1)
                        if old != new
                    }
                    nmut = len(changed)
                    surface_fraction = len(changed & surface) / nmut if nmut else 0.0
                    reason = ""
                    if not args.min_mutations <= nmut <= args.max_mutations:
                        reason = "mutation_budget"
                    elif {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
                        reason = "cysteine"
                    elif changed & EPITOPE:
                        reason = "SIINFEKL"
                    elif glyco_sites(sequence) != reference_glyco:
                        reason = "glycosylation"
                    elif surface_fraction < args.min_surface_fraction:
                        reason = "surface_fraction"
                    elif not changed <= designable:
                        reason = "mutation_outside_design_positions"
                    diagnostics.append({
                        "condition": condition,
                        "temperature": temperature,
                        "native_bias": native_bias,
                        "batch": batch_index,
                        "sample": sample_index,
                        "n_mutations": nmut,
                        "surface_mutation_fraction": surface_fraction,
                        "glyco_changed_positions": ",".join(map(str, sorted(
                            glyco_sites(sequence) ^ reference_glyco
                        ))),
                        "reason": reason or "accepted",
                        "mpnn_score": float(score),
                        "mutations": mutation_string(reference, sequence),
                        "sequence": sequence,
                    })
                    if reason:
                        rejection[reason] += 1
                        continue
                    rows.append({
                        "method": "tied_c4_ProteinMPNN",
                        "condition": condition,
                        "temperature": temperature,
                        "native_bias": native_bias,
                        "batch": batch_index,
                        "sample": sample_index,
                        "mpnn_score": float(score),
                        "mpnn_seqid_design_positions": float(seqid),
                        "n_mutations": nmut,
                        "surface_mutation_fraction": surface_fraction,
                        "hamming_to_parent": sum(a != b for a, b in zip(parent, sequence)),
                        "mutations": mutation_string(reference, sequence),
                        "sequence": sequence,
                    })

    with (args.out / "sample_diagnostics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(diagnostics[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(diagnostics)
    with (args.out / "rejection_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reason", "count"]); writer.writerows(sorted(rejection.items()))

    unique: dict[str, dict[str, object]] = {}
    for row in rows:
        sequence = str(row["sequence"])
        if sequence not in unique or float(row["mpnn_score"]) < float(unique[sequence]["mpnn_score"]):
            unique[sequence] = row
    rows = list(unique.values())
    if not rows:
        raise SystemExit(f"no valid designs; rejected={dict(rejection)}")
    score_ranks = average_ranks([float(row["mpnn_score"]) for row in rows])
    mutation_ranks = average_ranks([float(row["n_mutations"]) for row in rows])
    parent_ranks = average_ranks([float(row["hamming_to_parent"]) for row in rows])
    for row, score_rank, mutation_rank, parent_rank in zip(
        rows, score_ranks, mutation_ranks, parent_ranks
    ):
        row["generation_average_rank"] = 0.70 * score_rank + 0.15 * mutation_rank + 0.15 * parent_rank
    rows.sort(key=lambda row: (float(row["generation_average_rank"]), str(row["sequence"])))
    for index, row in enumerate(rows, 1):
        row["rank"] = index

    fields = ["rank"] + [field for field in rows[0] if field != "rank"]
    with (args.out / "mpnn_all_valid.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    shortlist: list[dict[str, object]] = []
    for row in rows:
        sequence = str(row["sequence"])
        if all(sum(a != b for a, b in zip(sequence, str(old["sequence"]))) >= 2 for old in shortlist):
            shortlist.append(row)
        if len(shortlist) >= args.shortlist:
            break
    with (args.out / "mpnn_shortlist.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(shortlist)
    with (args.out / "mpnn_shortlist.fasta").open("w") as handle:
        for index, row in enumerate(shortlist, 1):
            handle.write(
                f">OVA-C4-MPNN-{index:03d} nmut={row['n_mutations']} "
                f"score={float(row['mpnn_score']):.4f} surface={float(row['surface_mutation_fraction']):.3f}\n"
                f"{row['sequence']}\n"
            )
    print(
        f"raw_samples={len(diagnostics)} valid_unique={len(rows)} "
        f"shortlist={len(shortlist)} rejections={dict(rejection)} out={args.out}"
    )


if __name__ == "__main__":
    main()
