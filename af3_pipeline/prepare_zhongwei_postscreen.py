#!/usr/bin/env python3
"""Prepare remote Zhongwei manifests for the remaining v1.1 core screens.

The source candidate directories are copied to Zhongwei first.  This helper
creates independent input copies with remote MSA paths, so the local absolute
paths embedded in historical JSON files are never used on the compute node.
OpenDDE is intentionally absent from both manifests.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib


def rewrite_json(
    source: pathlib.Path,
    destination: pathlib.Path,
    local_root: str,
    remote_root: str,
    msa_override: str | None = None,
) -> None:
    payload = json.loads(source.read_text())

    def walk(value, key: str | None = None):
        if isinstance(value, dict):
            return {item_key: walk(item, item_key) for item_key, item in value.items()}
        if isinstance(value, list):
            return [walk(item, key) for item in value]
        if isinstance(value, str):
            if key in {"unpairedMsaPath", "pairedMsaPath"} and msa_override is not None:
                return msa_override
            return value.replace(local_root, remote_root)
        return value

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(walk(payload), indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, required=True,
                        help="remote .../experiments/current_pending_v1_screen")
    parser.add_argument("--local-root", default="/root/400083/ova_p3_tetramer_design",
                        help="prefix to replace in copied JSON paths")
    args = parser.parse_args()
    root = args.root.resolve()
    remote_root = str(root.parents[1])  # .../ova_p3_tetramer_design
    candidates = []
    for directory in sorted(root.glob("ova_body_c4_*")):
        if not (directory / "final_candidates.tsv").exists():
            continue
        if (directory / "SCREEN_V1_COMPLETE").exists():
            continue
        candidates.append(directory)
    if not candidates:
        raise SystemExit("no pending candidates found")

    protenix_rows = []
    stoich_rows = []
    candidate_rows = []
    seeds = (11, 21, 31, 41, 51)
    pro_index = 0
    sto_index = 0
    for directory in candidates:
        candidate = directory.name
        candidate_rows.append({"candidate": candidate, "directory": candidate})
        # Copy/rewrite the single Protenix input.  The MSA remains in the
        # candidate's copied protenix_multiseed/msa directory.
        source_inputs = sorted((directory / "protenix_multiseed" / "input").glob("*.json"))
        if len(source_inputs) != 1:
            raise ValueError(f"{candidate}: expected one Protenix input, found {len(source_inputs)}")
        pro_input = root / candidate / "zh_protenix_input" / source_inputs[0].name
        rewrite_json(source_inputs[0], pro_input, args.local_root, remote_root)
        copied_msa = root / candidate / "protenix_multiseed" / "msa" / (
            source_inputs[0].stem.replace("_protenix", "_protenix") + ".a3m"
        )
        if not copied_msa.exists():
            msa_candidates = sorted((directory / "protenix_multiseed" / "msa").glob("*.a3m"))
            if len(msa_candidates) != 1:
                raise FileNotFoundError(f"{candidate}: copied Protenix MSA not found")
            copied_msa = root / candidate / "protenix_multiseed" / "msa" / msa_candidates[0].name
        remote_msa = str(copied_msa)
        for seed in seeds:
            protenix_rows.append({
                "shard": pro_index % 8,
                "candidate": candidate,
                "seed": seed,
                "input_json": str(pro_input.relative_to(root)),
            })
            pro_index += 1

        source_manifest = directory / "stoichiometry_af3" / "manifest.tsv"
        if not source_manifest.exists():
            raise FileNotFoundError(source_manifest)
        original_stoich = {
            row["name"]: row for row in csv.DictReader(source_manifest.open(), delimiter="\t")
        }
        candidate_stoich_rows = []
        for source in sorted((directory / "stoichiometry_af3").glob("in_s*/*.json")):
            system = source.stem
            destination = root / candidate / "zh_stoich_input" / source.name
            rewrite_json(source, destination, args.local_root, remote_root, msa_override=remote_msa)
            if system not in original_stoich:
                raise ValueError(f"{candidate}: system {system} missing from stoichiometry manifest")
            assigned_shard = sto_index % 8
            candidate_stoich = dict(original_stoich[system])
            candidate_stoich["shard"] = assigned_shard
            candidate_stoich["input_json"] = str(destination.relative_to(root))
            candidate_stoich_rows.append(candidate_stoich)
            stoich_rows.append({
                "shard": assigned_shard,
                "candidate": candidate,
                "system": system,
                "input_json": str(destination.relative_to(root)),
            })
            sto_index += 1
        write_path = root / candidate / "zh_stoich_out" / "manifest.tsv"
        write_path.parent.mkdir(parents=True, exist_ok=True)
        with write_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(candidate_stoich_rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader(); writer.writerows(candidate_stoich_rows)

    def write(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader(); writer.writerows(rows)

    write(root / "zh_candidates.tsv", candidate_rows)
    write(root / "zh_protenix_manifest.tsv", protenix_rows)
    write(root / "zh_stoich_manifest.tsv", stoich_rows)
    print(f"candidates={len(candidates)} protenix_runs={len(protenix_rows)} stoich_systems={len(stoich_rows)}")
    print("protenix_shards=" + ",".join(str(sum(int(row["shard"]) == shard for row in protenix_rows)) for shard in range(8)))
    print("stoich_shards=" + ",".join(str(sum(int(row["shard"]) == shard for row in stoich_rows)) for shard in range(8)))


if __name__ == "__main__":
    main()
