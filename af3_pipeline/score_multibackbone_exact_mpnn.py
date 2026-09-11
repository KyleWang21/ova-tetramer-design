#!/usr/bin/env python3
"""Score exact tied homotetramer sequences on one backbone with ProteinMPNN."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import jax
import numpy as np
from scipy.special import log_softmax
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


def read_positions(path: Path) -> set[int]:
    return {int(token) for token in path.read_text().replace(",", " ").split() if token}


def compress_positions(chain: str, selected: set[int], length: int) -> str:
    values = sorted(selected)
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


def chain_sequences(path: Path) -> list[str]:
    model = PDBParser(QUIET=True).get_structure(path.stem, str(path))[0]
    return [
        "".join(
            AA3[residue.resname]
            for residue in model[chain_id]
            if residue.id[0] == " " and residue.resname in AA3
        )
        for chain_id in "ABCD"
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--score-positions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--context", required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--orders", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument("--expected-candidates", type=int, default=512)
    args = parser.parse_args()

    with args.candidates.open() as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    if len(candidates) != args.expected_candidates:
        raise ValueError(f"expected {args.expected_candidates} candidates, found {len(candidates)}")
    sequence_length = len(candidates[0]["sequence"])
    if sequence_length != 386 or any(len(row["sequence"]) != sequence_length for row in candidates):
        raise ValueError("candidate sequences must all be 386 aa")
    score_positions = read_positions(args.score_positions)
    if not score_positions or min(score_positions) < 1 or max(score_positions) > sequence_length:
        raise ValueError("invalid score positions")
    structure_sequences = chain_sequences(args.target_pdb)
    if [len(sequence) for sequence in structure_sequences] != [sequence_length] * 4:
        raise ValueError(f"target chain lengths={[len(sequence) for sequence in structure_sequences]}")
    if len(set(structure_sequences)) != 1:
        raise ValueError("target is not a sequence-identical homotetramer")
    if set(AA20) != set(aa_order):
        raise RuntimeError("unexpected ProteinMPNN amino-acid alphabet")

    fixed = set(range(1, sequence_length + 1)) - score_positions
    fix_pos = ",".join(
        compress_positions(chain, fixed, sequence_length) for chain in "ABCD"
    )
    model = mk_mpnn_model(
        model_name="v_48_020", backbone_noise=0.0, dropout=0.0,
        seed=args.seed + args.task_index, verbose=False, weights="original",
    )
    model.prep_inputs(
        pdb_filename=str(args.target_pdb), chain="A,B,C,D", homooligomer=True,
        fix_pos=fix_pos, verbose=True,
    )
    total_length = 4 * sequence_length
    design_indices = np.asarray([
        chain * sequence_length + position - 1
        for chain in range(4) for position in sorted(score_positions)
    ], dtype=np.int32)
    fixed_indices = np.asarray(sorted(set(range(total_length)) - set(design_indices)), dtype=np.int32)
    score_mask = np.zeros(total_length, dtype=bool)
    score_mask[design_indices] = True
    chain_masks = [
        score_mask[chain * sequence_length:(chain + 1) * sequence_length]
        for chain in range(4)
    ]

    decoding_orders = []
    for order_index in range(args.orders):
        rng = np.random.default_rng(args.seed + 1009 * args.task_index + order_index)
        fixed_order = fixed_indices.copy(); rng.shuffle(fixed_order)
        design_order = design_indices.copy(); rng.shuffle(design_order)
        decoding_orders.append(np.concatenate([fixed_order, design_order]).astype(np.int32))
    decoding_orders = np.asarray(decoding_orders)

    entries = [
        (candidate_index, order_index)
        for candidate_index in range(len(candidates))
        for order_index in range(args.orders)
    ]
    if len(entries) % args.batch_size:
        raise ValueError("candidate_count * orders must be divisible by batch-size")
    per_order: list[dict[str, object]] = []
    eye = np.eye(21, dtype=np.float32)
    for batch_start in range(0, len(entries), args.batch_size):
        batch_entries = entries[batch_start:batch_start + args.batch_size]
        sequence_indices = []
        batch_orders = []
        for candidate_index, order_index in batch_entries:
            sequence = candidates[candidate_index]["sequence"] * 4
            indices = np.asarray([aa_order[aa] for aa in sequence], dtype=np.int32)
            sequence_indices.append(indices)
            batch_orders.append(decoding_orders[order_index])
        sequence_indices_array = np.asarray(sequence_indices)
        one_hot = eye[sequence_indices_array]
        key = jax.random.fold_in(jax.random.PRNGKey(args.seed + args.task_index), batch_start)
        keys = jax.random.split(key, args.batch_size)
        inputs = {name: value for name, value in model._inputs.items()}
        output = model._rescore_parallel(
            keys, inputs, one_hot, np.asarray(batch_orders, dtype=np.int32)
        )
        logits = np.asarray(output["logits"])[..., :20]
        log_probabilities = log_softmax(logits, axis=-1)
        nll = -np.take_along_axis(
            log_probabilities, sequence_indices_array[..., None], axis=-1
        )[..., 0]
        for local_index, (candidate_index, order_index) in enumerate(batch_entries):
            chain_nll = []
            for chain in range(4):
                segment = nll[
                    local_index,
                    chain * sequence_length:(chain + 1) * sequence_length,
                ]
                chain_nll.append(float(segment[chain_masks[chain]].mean()))
            per_order.append({
                "candidate_index": candidate_index,
                "source_name": candidates[candidate_index]["source_name"],
                "task_index": args.task_index,
                "context": args.context,
                "family": args.family,
                "order_index": order_index,
                "design_nll": float(nll[local_index, score_mask].mean()),
                "chain_a_nll": chain_nll[0],
                "chain_b_nll": chain_nll[1],
                "chain_c_nll": chain_nll[2],
                "chain_d_nll": chain_nll[3],
                "chain_spread": max(chain_nll) - min(chain_nll),
            })
        if batch_start % (args.batch_size * 8) == 0:
            print(f"scored {batch_start + args.batch_size}/{len(entries)} sequence-orders", flush=True)

    grouped: dict[int, list[dict[str, object]]] = {}
    for row in per_order:
        grouped.setdefault(int(row["candidate_index"]), []).append(row)
    summaries: list[dict[str, object]] = []
    for candidate_index, candidate in enumerate(candidates):
        rows = sorted(grouped[candidate_index], key=lambda row: int(row["order_index"]))
        if len(rows) != args.orders:
            raise ValueError(f"candidate {candidate_index} has {len(rows)} order scores")
        values = np.asarray([float(row["design_nll"]) for row in rows])
        chain_spreads = np.asarray([float(row["chain_spread"]) for row in rows])
        summaries.append({
            "candidate_index": candidate_index,
            "source_name": candidate["source_name"],
            "sequence": candidate["sequence"],
            "n_mutations": candidate["n_mutations"],
            "surface_mutation_fraction": candidate["surface_mutation_fraction"],
            "task_index": args.task_index,
            "context": args.context,
            "family": args.family,
            "mean_design_nll": float(values.mean()),
            "worst_order_nll": float(values.max()),
            "std_order_nll": float(values.std(ddof=0)),
            "mean_chain_spread": float(chain_spreads.mean()),
            "order_nlls": ",".join(f"{value:.8f}" for value in values),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "per_order.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(per_order[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(per_order)
    with (args.out / "candidate_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summaries[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    metadata = {
        "method": "exact_tied_c4_ProteinMPNN_conditional_NLL",
        "task_index": args.task_index,
        "context": args.context,
        "family": args.family,
        "target_pdb": str(args.target_pdb),
        "candidate_count": len(candidates),
        "orders": args.orders,
        "batch_size": args.batch_size,
        "score_position_count_per_chain": len(score_positions),
        "scored_residue_count_total": int(score_mask.sum()),
        "conditioning": "all fixed scaffold positions precede all score positions",
        "backbone_noise": 0.0,
        "dropout": 0.0,
    }
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"complete {args.context}: candidates={len(summaries)} orders={args.orders}")


if __name__ == "__main__":
    main()
