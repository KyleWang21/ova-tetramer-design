#!/usr/bin/env python3
"""Apply the frozen AF3 geometry gates to a partial seed-1 C4 screen."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import statistics
from collections import defaultdict

from screen_final20_af3_structures import reference_ca, score_one


def fasta_sequence(path: pathlib.Path, name: str) -> str:
    active = False
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if active:
                return "".join(chunks)
            active = line[1:].split()[0] == name
        elif active:
            chunks.append(line.strip())
    raise ValueError(f"FASTA record {name!r} not found")


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--reference-cif", type=pathlib.Path, required=True)
    parser.add_argument("--reference-chain", default="A")
    parser.add_argument("--shard", type=int, help="score one output shard instead of all out_s*")
    parser.add_argument("--out", type=pathlib.Path)
    args = parser.parse_args()
    out = args.out or args.experiment

    manifest = {
        row["name"]: row
        for row in csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t")
    }
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    ref_ca = reference_ca(args.reference_cif, args.reference_chain)
    grouped: dict[str, list[pathlib.Path]] = defaultdict(list)
    shard_pattern = f"out_s{args.shard}" if args.shard is not None else "out_s*"
    for summary in args.experiment.glob(
        f"{shard_pattern}/**/seed-*/*_summary_confidences.json"
    ):
        name = summary.name.split("_seed-", 1)[0]
        if name in manifest:
            grouped[name].append(summary)

    model_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for name, paths in sorted(grouped.items()):
        if len(paths) != 5:
            continue
        sequence = manifest[name]["sequence"]
        mutations = {
            index for index, (old, new) in enumerate(zip(reference, sequence), 1) if old != new
        }
        current: list[dict[str, object]] = []
        for summary_path in sorted(paths):
            confidence = json.loads(summary_path.read_text())
            stem = summary_path.name.removesuffix("_summary_confidences.json")
            cif = summary_path.with_name(stem + "_model.cif")
            metrics, _ = score_one(cif, mutations, ref_ca)
            row = {
                "system": name,
                "model": summary_path.parent.name,
                "iptm": float(confidence["iptm"]),
                "ptm": float(confidence["ptm"]),
                **metrics,
            }
            current.append(row); model_rows.append(row)
        iptm = [float(row["iptm"]) for row in current]
        design_fraction = [float(row["design_interface_fraction"]) for row in current]
        summary = {
            "system": name,
            "method": manifest[name]["description"].split(" method=", 1)[-1].split()[0],
            "n_mutations": len(mutations),
            "mean_iptm": statistics.mean(iptm),
            "median_iptm": statistics.median(iptm),
            "min_iptm": min(iptm),
            "max_iptm": max(iptm),
            "zero_atomic_clash_models": sum(int(row["atomic_clash_count_2p4"]) == 0 for row in current),
            "terminal_only_clash_models": sum(
                int(row["atomic_clash_count_2p4"]) > 0
                and int(row["terminal_atomic_clash_count_2p4"])
                == int(row["atomic_clash_count_2p4"])
                for row in current
            ),
            "nonterminal_clash_models": sum(
                int(row["nonterminal_atomic_clash_count_2p4"]) > 0 for row in current
            ),
            "design_involved_clash_models": sum(
                int(row["design_involved_atomic_clash_count_2p4"]) > 0 for row in current
            ),
            "zero_ca_overlap_models": sum(int(row["severe_ca_overlap_count_2p5"]) == 0 for row in current),
            "valid_c4_network_models": sum(int(row["valid_c4_network"]) for row in current),
            "native_ss_complete_models": sum(int(row["native_c74_c121_bonds"]) == 4 for row in current),
            "zero_interchain_ss_models": sum(int(row["interchain_ss_bonds"]) == 0 for row in current),
            "two_interfaces_design_ge3_models": sum(int(row["two_interfaces_each_design_ge3"]) for row in current),
            "mean_design_interface_fraction": statistics.mean(design_fraction),
            "max_chain_ca_rmsd_to_1ova": max(float(row["max_chain_ca_rmsd_to_1ova"]) for row in current),
            "sequence": sequence,
        }
        gates = {
            "gate_mean_iptm_ge_0p80": summary["mean_iptm"] >= 0.80,
            "gate_min_iptm_ge_0p75": summary["min_iptm"] >= 0.75,
            "gate_5_zero_atomic_clash": summary["zero_atomic_clash_models"] == 5,
            "gate_5_zero_ca_overlap": summary["zero_ca_overlap_models"] == 5,
            "gate_5_valid_c4_network": summary["valid_c4_network_models"] == 5,
            "gate_5_native_ss_complete": summary["native_ss_complete_models"] == 5,
            "gate_5_zero_interchain_ss": summary["zero_interchain_ss_models"] == 5,
            "gate_5_two_interfaces_design_ge3": summary["two_interfaces_design_ge3_models"] == 5,
            "gate_max_chain_rmsd_le_2p0": summary["max_chain_ca_rmsd_to_1ova"] <= 2.0,
        }
        summary.update({key: int(value) for key, value in gates.items()})
        failed = [key for key, value in gates.items() if not value]
        summary["seed1_strict_pass"] = int(not failed)
        summary["failed_seed1_gates"] = ";".join(failed)
        summary_rows.append(summary)

    if not summary_rows:
        raise SystemExit("no complete five-sample systems")
    summary_rows.sort(key=lambda row: (
        -int(row["seed1_strict_pass"]), -float(row["mean_iptm"]), str(row["system"])
    ))
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "seed1_strict_per_model.tsv", model_rows)
    write_tsv(out / "seed1_strict_summary.tsv", summary_rows)
    print(
        f"complete={len(summary_rows)} strict_pass={sum(int(row['seed1_strict_pass']) for row in summary_rows)}"
    )
    for row in summary_rows[:10]:
        print(
            row["system"], f"mean={float(row['mean_iptm']):.3f}",
            f"zero_clash={row['zero_atomic_clash_models']}/5",
            f"two_interfaces={row['two_interfaces_design_ge3_models']}/5",
            f"pass={row['seed1_strict_pass']}",
        )


if __name__ == "__main__":
    main()
