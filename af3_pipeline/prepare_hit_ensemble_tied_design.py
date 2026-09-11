#!/usr/bin/env python3
"""Prepare four AF3 backbone states from one multiseed C4 hit for tied-AF2 design."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, PDBIO

from prepare_iptm90_tied_design import Chains, interface_positions, write_positions


PROTECTED = {12, 31, 74, 121, 368, 383} | set(range(258, 266))


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final-candidates", type=pathlib.Path, required=True)
    parser.add_argument("--af3-per-model", type=pathlib.Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--pairs", default="AC,AD")
    parser.add_argument("--backbones", type=int, default=4)
    parser.add_argument(
        "--anchor-prefix",
        default="OVA-C4-073",
        help="Prefix for generated backbone FASTA records (default preserves E169 names).",
    )
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    final = {row["candidate"]: row for row in read_tsv(args.final_candidates)}
    if args.candidate not in final:
        raise ValueError(f"candidate {args.candidate!r} is absent from final table")
    base = final[args.candidate]
    sequence = base["sequence"]
    mutation_positions = {
        int(item[1:-1]) for item in base["mutations"].split(",") if item
    }
    surface = {
        int(value)
        for value in args.surface_positions.read_text().replace(",", " ").split()
    }
    pairs = tuple(value.strip() for value in args.pairs.split(","))
    if len(pairs) != 2 or any(len(pair) != 2 for pair in pairs):
        raise ValueError("--pairs must contain two two-chain identifiers, e.g. AC,AD")

    rows = [
        row for row in read_tsv(args.af3_per_model)
        if row["candidate"] == args.candidate
        and int(row["atomic_clash_count_2p4"]) == 0
        and int(row["severe_ca_overlap_count_2p5"]) == 0
        and int(row["valid_c4_network"]) == 1
        and int(row["native_c74_c121_bonds"]) == 4
        and int(row["interchain_ss_bonds"]) == 0
        and int(row["two_interfaces_each_design_ge3"]) == 1
        and int(row["in_main_topology_cluster"]) == 1
    ]
    if len(rows) < args.backbones:
        raise ValueError(f"only {len(rows)} eligible AF3 models for {args.backbones} backbones")
    medoids = [row for row in rows if int(row["main_topology_medoid"]) == 1]
    chosen: list[dict[str, str]] = sorted(
        medoids, key=lambda row: (-float(row["iptm"]), float(row["interface_pae"]))
    )[:1]
    chosen_seeds = {int(row["seed"]) for row in chosen}
    best_by_seed: list[dict[str, str]] = []
    for seed in sorted({int(row["seed"]) for row in rows}):
        candidates = [row for row in rows if int(row["seed"]) == seed]
        best_by_seed.append(max(candidates, key=lambda row: (
            float(row["iptm"]), float(row["min_incident_iptm"]),
            -float(row["interface_pae"]), float(row["mean_chain_plddt"]),
        )))
    for row in sorted(best_by_seed, key=lambda row: (
        -float(row["iptm"]), -float(row["min_incident_iptm"]), float(row["interface_pae"])
    )):
        seed = int(row["seed"])
        if seed in chosen_seeds:
            continue
        chosen.append(row)
        chosen_seeds.add(seed)
        if len(chosen) == args.backbones:
            break
    if len(chosen) < args.backbones:
        for row in sorted(rows, key=lambda row: -float(row["iptm"])):
            if row not in chosen:
                chosen.append(row)
            if len(chosen) == args.backbones:
                break

    args.out.mkdir(parents=True, exist_ok=True)
    write_positions(args.out / "surface_positions.txt", surface)
    target_dir = args.out / "targets"
    target_dir.mkdir(exist_ok=True)
    records: list[dict[str, object]] = []
    with (args.out / "anchors.fasta").open("w") as fasta:
        for index, row in enumerate(chosen, 1):
            anchor = f"{args.anchor_prefix}-BB{index}"
            fasta.write(f">{anchor}\n{sequence}\n")
            structure = MMCIFParser(QUIET=True).get_structure(anchor, row["cif_path"])
            model = structure[0]
            io = PDBIO(); io.set_structure(structure)
            stem = anchor.lower().replace("-", "_")
            target_c4 = target_dir / f"{stem}_c4.pdb"
            io.save(str(target_c4), Chains(set("ABCD")))
            record: dict[str, object] = {
                "anchor": anchor,
                "anchor_n_mutations": len(mutation_positions),
                "source_cif": row["cif_path"],
                "source_iptm": row["iptm"],
                "source_atomic_clashes": row["atomic_clash_count_2p4"],
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
                design = (mutation_positions | (contacts & surface)) - PROTECTED
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
            joint = (mutation_positions | ((contact_sets[0] | contact_sets[1]) & surface)) - PROTECTED
            joint_path = target_dir / f"{stem}_joint_design.txt"
            write_positions(joint_path, joint)
            record["joint_design"] = str(joint_path)
            record["joint_design_positions"] = len(joint)
            records.append(record)

    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    print(
        f"prepared {len(records)} AF3 backbone states for {args.candidate}; "
        f"seeds={[row['seed'] for row in chosen]} pairs={pairs} out={args.out}"
    )


if __name__ == "__main__":
    main()
