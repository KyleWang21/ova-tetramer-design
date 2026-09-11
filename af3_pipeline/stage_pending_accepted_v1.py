#!/usr/bin/env python3
"""Materialize every newly Accepted_AF3 sequence for the generic v1.1 core screen."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_sequence(path: pathlib.Path, wanted: str) -> str:
    records: dict[str, str] = {}
    name = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += line
    return records[wanted]


def candidate_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=pathlib.Path, required=True)
    parser.add_argument("--registry", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    args = parser.parse_args()
    project = args.project.resolve()
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    pending = [row for row in read_tsv(args.registry) if int(row["needs_v1_screen"]) == 1]
    args.out.mkdir(parents=True, exist_ok=True)
    staged = 0
    for registry_row in pending:
        candidate = registry_row["candidate"]
        # Historical Accepted_AF3 entries can predate the generic full25
        # directory layout.  Their missing stages are backfilled by dedicated
        # experiment runners so that completed Protenix/Rosetta/OpenDDE work is
        # reused instead of being submitted again.  Do not make one historical
        # row abort staging of all future generic promotions.
        if registry_row.get("source_kind") != "generic_25model_promotion":
            print(
                f"deferred historical backfill {candidate}: "
                f"source_kind={registry_row.get('source_kind', '')}"
            )
            continue
        token = candidate_token(candidate)
        screen = args.out / token
        if (screen / "STAGING_READY").exists():
            continue
        source = pathlib.Path(registry_row["source_experiment"])
        if not source.is_absolute():
            source = project / source
        summary_path = source / "full25" / "full25_summary.tsv"
        models_path = source / "full25" / "full25_per_model.tsv"
        if not summary_path.exists() or not models_path.exists():
            raise FileNotFoundError(f"{candidate}: missing generic full25 evidence under {source}")
        summaries = [row for row in read_tsv(summary_path) if row["system"] == candidate]
        if len(summaries) != 1 or int(summaries[0]["accepted_af3"]) != 1:
            raise ValueError(f"{candidate}: registry/full25 mismatch")
        summary = summaries[0]
        sequence = summary["sequence"]
        mutations = ",".join(
            f"{old}{index}{new}"
            for index, (old, new) in enumerate(zip(reference, sequence), 1)
            if old != new
        )
        if len(sequence) != 386 or int(summary["n_mutations"]) != len(mutations.split(",")):
            raise ValueError(f"{candidate}: invalid sequence or mutation count")
        models = [row for row in read_tsv(models_path) if row["system"] == candidate]
        if len(models) != 25:
            raise ValueError(f"{candidate}: expected 25 AF3 models, found {len(models)}")
        representative = pathlib.Path(registry_row["representative_af3_cif"])
        if not representative.is_absolute():
            representative = project / representative
        if not representative.exists():
            raise FileNotFoundError(representative)

        screen.mkdir(parents=True, exist_ok=True)
        final_row = {
            "rank": 1,
            "candidate": candidate,
            "sequence": sequence,
            "mutations": mutations,
            "passes_final": 1,
        }
        write_tsv(screen / "final_candidates.tsv", [final_row])
        normalized_summary = {
            "candidate": candidate,
            "accepted_af3_recomputed": 1,
            **{key: value for key, value in summary.items() if key not in {"system", "sequence"}},
            "sequence": sequence,
        }
        write_tsv(screen / "af3_candidate_recomputed_summary.tsv", [normalized_summary])
        normalized_models = []
        for row in models:
            normalized_models.append({
                **row,
                "candidate": candidate,
                "system": candidate,
            })
        write_tsv(screen / "af3_per_model_recomputed.tsv", normalized_models)
        write_tsv(screen / "af3_prolif_structure_manifest.tsv", [{
            "source": f"{candidate}__AF3",
            "cif_path": str(representative.resolve()),
        }])
        write_tsv(screen / "source_evidence.tsv", [{
            "candidate": candidate,
            "source_experiment": str(source),
            "source_full25_summary": str(summary_path),
            "source_full25_models": str(models_path),
            "representative_af3_cif": str(representative.resolve()),
        }])
        (screen / "STAGING_READY").touch()
        staged += 1
        print(f"staged {candidate} -> {screen}")
    print(f"pending={len(pending)} newly_staged={staged}")


if __name__ == "__main__":
    main()
