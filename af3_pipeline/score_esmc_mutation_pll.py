#!/usr/bin/env python3
"""Score designed substitutions with masked ESM-C probabilities in full OVA context."""

from __future__ import annotations

import argparse
import csv
import pathlib

import torch
from esm.models.esmc import ESMC


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    ap.add_argument("--reference", type=pathlib.Path, required=True)
    ap.add_argument("--reference-id", default="D0-P3-13R")
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--model", default="esmc_300m")
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    reference_records = read_fasta(args.reference)
    reference = next(sequence for header, sequence in reference_records
                     if header.split()[0] == args.reference_id)
    candidates = read_fasta(args.fasta)
    if not candidates or any(len(sequence) != len(reference) for _, sequence in candidates):
        raise SystemExit("candidate/reference FASTA lengths do not match")

    device = args.device if args.device == "cpu" or torch.cuda.is_available() else "cpu"
    model = ESMC.from_pretrained(args.model).to(device).eval()
    tokenizer = model.tokenizer
    rows: list[dict[str, object]] = []
    with torch.no_grad():
        for rank, (header, sequence) in enumerate(candidates, 1):
            positions = [index for index, (old, new) in enumerate(zip(reference, sequence)) if old != new]
            ids = torch.tensor(tokenizer.encode(sequence), device=device)
            batch = ids.unsqueeze(0).repeat(len(positions), 1).clone()
            for row_index, position in enumerate(positions):
                batch[row_index, position + 1] = tokenizer.mask_token_id
            logits = model(sequence_tokens=batch).sequence_logits.float()
            logp = torch.log_softmax(logits, dim=-1)
            candidate_values, reference_values = [], []
            for row_index, position in enumerate(positions):
                candidate_id = int(ids[position + 1])
                reference_id = int(tokenizer.encode(reference[position])[1])
                candidate_values.append(float(logp[row_index, position + 1, candidate_id]))
                reference_values.append(float(logp[row_index, position + 1, reference_id]))
            delta = [new - old for new, old in zip(candidate_values, reference_values)]
            rows.append({
                "rank": rank,
                "header": header,
                "n_mutations_vs_p3": len(positions),
                "mutation_pll_mean": sum(candidate_values) / len(candidate_values),
                "p3_aa_pll_same_context_mean": sum(reference_values) / len(reference_values),
                "mutation_minus_p3_pll_mean": sum(delta) / len(delta),
                "worst_mutation_delta": min(delta),
                "mutations": ",".join(
                    f"{reference[pos]}{pos + 1}{sequence[pos]}" for pos in positions
                ),
                "sequence": sequence,
            })
            print(rank, len(positions), f"mean_delta={rows[-1]['mutation_minus_p3_pll_mean']:.3f}",
                  f"worst={rows[-1]['worst_mutation_delta']:.3f}", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
