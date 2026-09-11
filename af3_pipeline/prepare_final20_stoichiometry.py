#!/usr/bin/env python3
"""Prepare C1/C2/C3/C5/C6 AF3 inputs for the final 20 OVA C4 sequences."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib


CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def source_input_from_summary(summary_path: pathlib.Path, system: str) -> pathlib.Path:
    experiment = next(parent for parent in summary_path.parents if parent.name.startswith("e"))
    out_shard = next(parent.name for parent in summary_path.parents if parent.name.startswith("out_s"))
    shard = out_shard.removeprefix("out_s")
    path = experiment / f"in_s{shard}" / f"{system}.json"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def source_input_from_model(
    row: dict[str, str], candidate: str, expected_sequence: str
) -> pathlib.Path:
    """Resolve the AF3 input for both historical and generic per-model tables."""
    if row.get("summary_path") and row.get("source_system"):
        return source_input_from_summary(
            pathlib.Path(row["summary_path"]), row["source_system"]
        )
    cif = pathlib.Path(row["cif_path"])
    if not cif.is_absolute():
        cif = pathlib.Path.cwd() / cif
    try:
        out_shard = next(parent for parent in cif.parents if parent.name.startswith("out_s"))
    except StopIteration as exc:
        raise FileNotFoundError(f"cannot infer source input from {cif}") from exc
    shard = out_shard.name.removeprefix("out_s")
    experiment = out_shard.parent
    names = [row.get("source_system", ""), row.get("system", ""), candidate]
    candidates = [experiment / f"in_s{shard}" / f"{name}.json" for name in names if name]
    candidates.extend(sorted((experiment / f"in_s{shard}").glob("*.json")))
    seen = set()
    for path in candidates:
        if path in seen or not path.exists():
            continue
        seen.add(path)
        payload = json.loads(path.read_text())
        if payload["sequences"][0]["protein"]["sequence"] == expected_sequence:
            return path
    raise FileNotFoundError(
        f"no source AF3 input for {candidate}; inferred experiment={experiment} shard={shard}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--per-model", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--stoichiometries", default="1,2,3,5,6")
    parser.add_argument("--seeds", default="1,2,3")
    parser.add_argument("--shards", type=int, default=8)
    args = parser.parse_args()

    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment already prepared: {args.out}")
    final = read_tsv(args.final)
    per_model = read_tsv(args.per_model)
    first_model: dict[str, dict[str, str]] = {}
    for row in per_model:
        first_model.setdefault(row["candidate"], row)

    stoichiometries = [int(item) for item in args.stoichiometries.split(",")]
    seeds = [int(item) for item in args.seeds.split(",")]
    args.out.mkdir(parents=True)
    for shard in range(args.shards):
        (args.out / f"in_s{shard}").mkdir()

    # Greedy balance by an n^2 proxy because AF3 cost grows steeply with token count.
    shard_load = [0] * args.shards
    jobs: list[tuple[int, int, dict[str, str]]] = []
    for base in final:
        for n in stoichiometries:
            shard = min(range(args.shards), key=lambda idx: (shard_load[idx], idx))
            shard_load[shard] += n * n
            jobs.append((shard, n, base))

    rows: list[dict[str, object]] = []
    for shard, n, base in jobs:
        candidate = base["candidate"]
        source = first_model[candidate]
        source_json = source_input_from_model(source, candidate, base["sequence"])
        source_payload = json.loads(source_json.read_text())
        protein = dict(source_payload["sequences"][0]["protein"])
        sequence = base["sequence"]
        if protein["sequence"] != sequence:
            raise ValueError(f"{candidate}: source input sequence mismatch")
        ids = CHAIN_IDS[:n]
        proteins = []
        for chain in ids:
            item = dict(protein)
            item["id"] = chain
            proteins.append({"protein": item})
        name = f"{candidate.lower().replace('-', '_')}_c{n}"
        payload = {
            "name": name,
            "modelSeeds": seeds,
            "dialect": "alphafold3",
            "version": 2,
            "sequences": proteins,
        }
        input_json = args.out / f"in_s{shard}" / f"{name}.json"
        input_json.write_text(json.dumps(payload, indent=2) + "\n")
        rows.append({
            "name": name,
            "candidate": candidate,
            "shard": shard,
            "kind": "final20_stoichiometry",
            "n_chains": n,
            "chain_ids": ids,
            "chain_length": ",".join(["386"] * n),
            "ova_ranges": ";".join(["1-386"] * n),
            "module_ranges": ";".join(["0-0"] * n),
            "architecture": "pure OVA homooligomer competition",
            "seeds": len(seeds),
            "expected_samples": len(seeds) * 5,
            "description": f"{candidate} forced C{n}; no templates; matched OVA MSA",
            "source_input_json": str(source_json),
            "input_json": str(input_json),
        })

    fields = list(rows[0])
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"prepared {len(rows)} AF3 systems across {args.shards} shards")
    print("shard_load_proxy=" + ",".join(map(str, shard_load)))


if __name__ == "__main__":
    main()
