#!/usr/bin/env python3
"""Prepare OpenDDE-v1 C4 calibration controls and final-candidate inputs.

Calibration (WT, P3-13R, and process-positive designs) is retained for
diagnostic reporting. OpenDDE is not promoted to a required filter in v1.1.
"""

import argparse
import csv
import json
import pathlib
import re


def first_fasta(path: pathlib.Path) -> str:
    sequence = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if sequence:
                break
        else:
            sequence.append(line.strip())
    return "".join(sequence)


def reverse_reference(sequence: str, mutations: str) -> str:
    output = list(sequence)
    for mutation in mutations.split(","):
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mutation)
        old, position, new = match.groups()
        position = int(position)
        if output[position - 1] != new:
            raise ValueError(mutation)
        output[position - 1] = old
    return "".join(output)


def replace_query(a3m: str, sequence: str) -> str:
    """Replace only the first A3M record, preserving the frozen homolog set."""
    lines = a3m.replace("\x00", "").splitlines()
    first_header = next(index for index, line in enumerate(lines) if line.startswith(">"))
    next_header = next(
        (index for index in range(first_header + 1, len(lines)) if lines[index].startswith(">")),
        len(lines),
    )
    return "\n".join(lines[: first_header + 1] + [sequence] + lines[next_header:]) + "\n"


def write_manifest(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--wt-a3m", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-msa-dir", type=pathlib.Path, required=True)
    parser.add_argument("--remote-candidate-msa-dir", required=True)
    parser.add_argument("--remote-output-msa-dir", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seeds", default="101,102,103,104,105,106,107,108,109,110")
    args = parser.parse_args()
    rows = list(csv.DictReader(args.final.open(), delimiter="\t"))
    by_name = {row["candidate"]: row for row in rows}
    p3 = reverse_reference(rows[0]["sequence"], rows[0]["mutations"])
    wt = first_fasta(args.wt_a3m).replace("-", "")
    if len(p3) != 386 or len(wt) != 386:
        raise ValueError(f"expected 386-aa P3 and WT, got {len(p3)}, {len(wt)}")
    seeds = [int(value) for value in args.seeds.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)
    input_dir = args.out / "input"; input_dir.mkdir(exist_ok=True)
    msa_dir = args.out / "msa"; msa_dir.mkdir(exist_ok=True)
    base_candidate_msa = args.candidate_msa_dir / (
        "ova_iptm80_" + rows[0]["candidate"].lower().replace("-", "_") +
        "_protenix_1024.a3m"
    )
    base_a3m = base_candidate_msa.read_text()
    control_msa = {}
    for control_name, sequence in (("CTRL-WT-C4", wt), ("CTRL-P3-13R-C4", p3)):
        path = msa_dir / (control_name.lower().replace("-", "_") + "_1024.a3m")
        path.write_text(replace_query(base_a3m, sequence))
        control_msa[control_name] = f"{args.remote_output_msa_dir.rstrip('/')}/{path.name}"
    systems = [
        ("CTRL-WT-C4", "negative_background", wt),
        ("CTRL-P3-13R-C4", "negative_background", p3),
        ("OVA-C4-IPTM80-19K", "process_positive", by_name["OVA-C4-IPTM80-19K"]["sequence"]),
        ("OVA-C4-IPTM80-21B", "process_positive", by_name["OVA-C4-IPTM80-21B"]["sequence"]),
    ] + [(row["candidate"], "final_candidate", row["sequence"]) for row in rows]
    unique = {}
    for name, kind, sequence in systems:
        unique.setdefault(name, (kind, sequence))
    manifest = []
    for index, (name, (kind, sequence)) in enumerate(unique.items()):
        job = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_") + "_opendde_c4"
        if name in control_msa:
            msa_path = control_msa[name]
        else:
            msa_name = (
                "ova_iptm80_" + name.lower().replace("-", "_") +
                "_protenix_1024.a3m"
            )
            local_msa = args.candidate_msa_dir / msa_name
            if not local_msa.exists():
                raise FileNotFoundError(local_msa)
            msa_path = f"{args.remote_candidate_msa_dir.rstrip('/')}/{msa_name}"
        payload = [{
            "name": job,
            "modelSeeds": seeds,
            "sequences": [{"proteinChain": {
                "sequence": sequence, "count": 4, "unpairedMsaPath": msa_path,
            }}],
        }]
        path = input_dir / f"{job}.json"
        path.write_text(json.dumps(payload, indent=2) + "\n")
        manifest.append({
            "shard": index % 8, "name": name, "kind": kind,
            "n_seeds": len(seeds), "input_json": str(path.relative_to(args.out)),
        })
    write_manifest(args.out / "manifest.tsv", manifest)
    write_manifest(
        args.out / "calibration_manifest.tsv",
        [row for row in manifest if row["kind"] != "final_candidate"],
    )
    write_manifest(
        args.out / "candidate_remaining_manifest.tsv",
        [row for row in manifest if row["kind"] == "final_candidate"],
    )
    print(f"prepared {len(manifest)} OpenDDE systems; calibration controls=4")


if __name__ == "__main__":
    main()
