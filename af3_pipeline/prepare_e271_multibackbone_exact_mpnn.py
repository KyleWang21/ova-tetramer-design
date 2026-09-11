#!/usr/bin/env python3
"""Freeze the E271 exact-sequence, eight-backbone ProteinMPNN scoring input."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio.PDB import PDBParser


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def read_fasta(path: Path, wanted: str) -> str:
    records: dict[str, str] = {}
    name: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += line
    if wanted not in records:
        raise ValueError(f"missing FASTA record {wanted!r} in {path}")
    return records[wanted]


def read_positions(path: Path) -> set[int]:
    return {
        int(token)
        for token in path.read_text().replace(",", " ").split()
        if token.strip()
    }


def chain_sequences(path: Path) -> list[str]:
    model = PDBParser(QUIET=True).get_structure(path.stem, str(path))[0]
    sequences = []
    for chain_id in "ABCD":
        chain = model[chain_id]
        sequence = "".join(
            AA3[residue.resname]
            for residue in chain
            if residue.id[0] == " " and residue.resname in AA3
        )
        sequences.append(sequence)
    return sequences


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--source-ranking", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=Path, required=True)
    parser.add_argument("--min-mutations", type=int, default=12)
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--min-surface-fraction", type=float, default=0.80)
    parser.add_argument("--expected-candidates", type=int, default=512)
    parser.add_argument("--hotspots", default="99,155,157,182,184,334,336,338")
    args = parser.parse_args()

    marker = args.out / "INPUT_FROZEN"
    if marker.exists():
        print(f"{args.out} already frozen")
        return
    project = args.project.resolve()
    reference = read_fasta(args.reference_fasta, args.reference_name)
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa P3-13R reference, found {len(reference)}")
    surface = read_positions(args.surface_positions)
    native_cys = {index for index, aa in enumerate(reference, 1) if aa == "C"}
    reference_glyco = glyco_sites(reference)
    epitope_start = reference.index("SIINFEKL")

    with args.source_ranking.open() as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != args.expected_candidates:
        raise ValueError(
            f"expected {args.expected_candidates} source candidates, found {len(source_rows)}"
        )
    seen: set[str] = set()
    candidates: list[dict[str, object]] = []
    score_positions: set[int] = {
        int(value) for value in args.hotspots.split(",") if value.strip()
    }
    for source_index, row in enumerate(source_rows, 1):
        sequence = row["sequence"]
        if sequence in seen:
            raise ValueError(f"duplicate source sequence at row {source_index}")
        seen.add(sequence)
        if len(sequence) != len(reference):
            raise ValueError(f"row {source_index} length={len(sequence)}")
        changed = {
            index
            for index, (old, new) in enumerate(zip(reference, sequence), 1)
            if old != new
        }
        nmut = len(changed)
        surface_fraction = len(changed & surface) / nmut if nmut else 0.0
        if not args.min_mutations <= nmut <= args.max_mutations:
            raise ValueError(f"row {source_index} mutation count={nmut}")
        if surface_fraction < args.min_surface_fraction:
            raise ValueError(f"row {source_index} surface fraction={surface_fraction:.3f}")
        if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != native_cys:
            raise ValueError(f"row {source_index} changes cysteine pattern")
        if sequence[epitope_start:epitope_start + 8] != "SIINFEKL":
            raise ValueError(f"row {source_index} changes SIINFEKL")
        if glyco_sites(sequence) != reference_glyco:
            raise ValueError(f"row {source_index} changes glycosylation pattern")
        score_positions |= changed
        candidates.append({
            "candidate_index": source_index - 1,
            "source_ai_rank": row.get("ai_rank", source_index),
            "source_name": row["name"],
            "sequence": sequence,
            "n_mutations": nmut,
            "surface_mutation_fraction": surface_fraction,
            "mutations": mutation_string(reference, sequence),
            "generation_average_rank": row["generation_average_rank"],
            "predicted_af3_mean_iptm": row["predicted_af3_mean_iptm"],
            "predicted_zero_clash_fraction": row["predicted_zero_clash_fraction"],
            "af3_strict_surrogate_rank_score": row["af3_strict_surrogate_rank_score"],
        })

    ai11_root = project / "experiments/e221_c4_ai11_af3_protenix_crossmodel_tied_design"
    ova073_root = project / "experiments/e225_c4_073_rosetta_polar_tied_mpnn"
    backbones: list[dict[str, object]] = []
    for index in range(1, 5):
        target = ai11_root / f"targets/ai11_xm_bb{index}_af3_c4.pdb"
        mask = ai11_root / f"targets/ai11_xm_bb{index}_joint_design.txt"
        backbones.append({
            "task_index": index - 1,
            "context": f"ai11_af3_bb{index}",
            "family": "AI11_AF3",
            "target_pdb": str(target.relative_to(project)),
            "source_design_positions": str(mask.relative_to(project)),
        })
    for index in range(1, 5):
        target = ova073_root / f"targets/bb{index}_c4.pdb"
        mask = ova073_root / "weak_interface_design_positions.txt"
        backbones.append({
            "task_index": index + 3,
            "context": f"ova073_rosetta_bb{index}",
            "family": "OVA073_ROSETTA",
            "target_pdb": str(target.relative_to(project)),
            "source_design_positions": str(mask.relative_to(project)),
        })

    backbone_audit = []
    for row in backbones:
        target = project / str(row["target_pdb"])
        mask = project / str(row["source_design_positions"])
        if not target.is_file() or not mask.is_file():
            raise FileNotFoundError(f"missing backbone input: {target} or {mask}")
        sequences = chain_sequences(target)
        if any(len(sequence) != len(reference) for sequence in sequences):
            raise ValueError(f"{target}: chain lengths={[len(sequence) for sequence in sequences]}")
        if len(set(sequences)) != 1:
            raise ValueError(f"{target}: chains are not sequence-identical")
        parent_mutations = {
            index
            for index, (old, new) in enumerate(zip(reference, sequences[0]), 1)
            if old != new
        }
        design_positions = read_positions(mask)
        score_positions |= parent_mutations | design_positions
        backbone_audit.append({
            **row,
            "parent_mutations": sorted(parent_mutations),
            "source_design_positions_values": sorted(design_positions),
        })

    score_positions -= native_cys
    score_positions -= set(range(epitope_start + 1, epitope_start + 9))
    if not score_positions or min(score_positions) < 1 or max(score_positions) > len(reference):
        raise ValueError("invalid global score-position union")
    if any(
        any(old != new and index not in score_positions
            for index, (old, new) in enumerate(zip(reference, str(row["sequence"])), 1))
        for row in candidates
    ):
        raise ValueError("candidate mutation outside global score-position union")

    args.out.mkdir(parents=True, exist_ok=False)
    candidate_fields = list(candidates[0])
    with (args.out / "candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, candidate_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(candidates)
    with (args.out / "candidates.fasta").open("w") as handle:
        for row in candidates:
            handle.write(f">E271-SOURCE-{int(row['candidate_index']) + 1:03d}\n{row['sequence']}\n")
    with (args.out / "backbone_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(backbones[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(backbones)
    (args.out / "score_positions.txt").write_text(
        " ".join(map(str, sorted(score_positions))) + "\n"
    )
    metadata = {
        "method": "exact_tied_c4_ProteinMPNN_conditional_NLL",
        "candidate_count": len(candidates),
        "backbone_count": len(backbones),
        "backbone_families": ["AI11_AF3", "OVA073_ROSETTA"],
        "score_position_count": len(score_positions),
        "score_positions": sorted(score_positions),
        "conditioning": "all non-score scaffold positions decoded before score positions",
        "decoding_orders_per_candidate_per_backbone": 4,
        "source_ranking": str(args.source_ranking),
        "source_ranking_sha256": hashlib.sha256(args.source_ranking.read_bytes()).hexdigest(),
        "candidate_table_sha256": hashlib.sha256(
            (args.out / "candidates.tsv").read_bytes()
        ).hexdigest(),
        "backbone_audit": backbone_audit,
    }
    (args.out / "input_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    marker.touch()
    print(
        f"frozen E271: candidates={len(candidates)} backbones={len(backbones)} "
        f"score_positions={len(score_positions)}"
    )


if __name__ == "__main__":
    main()
