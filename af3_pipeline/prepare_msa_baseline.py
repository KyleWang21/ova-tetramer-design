#!/usr/bin/env python3
"""Fetch one OVA MSA and prepare the WT/P3 AF3 dimer confirmation run."""

from __future__ import annotations

import argparse
import csv
import io
import json
import pathlib
import re
import tarfile
import time

import requests


API = "https://api.colabfold.com"
USER_AGENT = "ova-tetramer-design/1.0"


def fetch_colabfold_a3m(sequence: str, cache_dir: pathlib.Path) -> str:
    """Run the public ColabFold unpaired-MSA search, with a local result cache."""
    a3m_path = cache_dir / "ova_wt_colabfold.a3m"
    if a3m_path.exists() and a3m_path.stat().st_size:
        a3m = a3m_path.read_text().replace("\x00", "")
        a3m = re.sub(r"(?<!\n)>", "\n>", a3m).lstrip()
        a3m_path.write_text(a3m)
        return a3m
    cache_dir.mkdir(parents=True, exist_ok=True)
    headers = {"User-Agent": USER_AGENT}
    query = f">101\n{sequence}\n"
    ticket = None
    for attempt in range(12):
        response = requests.post(
            f"{API}/ticket/msa",
            data={"q": query, "mode": "env"},
            headers=headers,
            timeout=(15, 90),
        )
        response.raise_for_status()
        ticket = response.json()
        status = ticket.get("status")
        print(f"submit attempt={attempt + 1} status={status}", flush=True)
        if status not in {"UNKNOWN", "RATELIMIT"}:
            break
        time.sleep(8)
    if not ticket or ticket.get("status") in {"UNKNOWN", "RATELIMIT", "ERROR", "MAINTENANCE"}:
        raise RuntimeError(f"ColabFold MSA submission failed: {ticket}")

    job_id = ticket["id"]
    status = ticket.get("status")
    while status in {"PENDING", "RUNNING", "UNKNOWN"}:
        time.sleep(8)
        response = requests.get(
            f"{API}/ticket/{job_id}", headers=headers, timeout=(15, 90)
        )
        response.raise_for_status()
        ticket = response.json()
        status = ticket.get("status")
        print(f"msa job={job_id} status={status}", flush=True)
    if status != "COMPLETE":
        raise RuntimeError(f"ColabFold MSA job failed: {ticket}")

    response = requests.get(
        f"{API}/result/download/{job_id}", headers=headers, timeout=(15, 180)
    )
    response.raise_for_status()
    parts: list[str] = []
    with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as archive:
        for name in ("uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m"):
            member = archive.extractfile(name)
            if member is not None:
                parts.append(member.read().decode())
    if not parts:
        raise RuntimeError("ColabFold result contained no A3M")
    # Each database A3M starts with >101; ensure the two files do not meet on
    # the same line if the first archive member lacks a terminal newline.
    a3m = "\n".join(part.rstrip("\n") for part in parts) + "\n"
    a3m_path.write_text(a3m)
    return a3m


def replace_first_query(a3m: str, query: str) -> str:
    """Replace only the first A3M record with the exact AF3 query sequence."""
    lines = a3m.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("invalid A3M")
    end = 1
    while end < len(lines) and not lines[end].startswith(">"):
        end += 1
    return "\n".join([lines[0], query, *lines[end:]]) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seeds", type=int, default=5)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for shard in (0, 1):
        (args.out / f"in_s{shard}").mkdir(exist_ok=True)

    wt_payload = json.loads((args.source / "in_s2/wt_dimer.json").read_text())
    p3_payload = json.loads((args.source / "in_s3/p3_dimer.json").read_text())
    wt = wt_payload["sequences"][0]["protein"]["sequence"]
    p3 = p3_payload["sequences"][0]["protein"]["sequence"]
    wt_msa = replace_first_query(
        fetch_colabfold_a3m(wt, args.out.parent / "msa_cache"), wt
    )
    p3_msa = replace_first_query(wt_msa, p3)

    systems = []
    for shard, (name, sequence, a3m) in enumerate((
        ("wt_dimer_msa", wt, wt_msa),
        ("p3_dimer_msa", p3, p3_msa),
    )):
        payload = {
            "name": name,
            "modelSeeds": list(range(1, args.seeds + 1)),
            "dialect": "alphafold3",
            "version": 2,
            "sequences": [
                {"protein": {"id": chain, "sequence": sequence,
                 "unpairedMsa": a3m, "pairedMsa": "", "templates": []}}
                for chain in "AB"
            ],
        }
        json_path = args.out / f"in_s{shard}" / f"{name}.json"
        json_path.write_text(json.dumps(payload) + "\n")
        systems.append({
            "name": name,
            "shard": shard,
            "kind": "dimer_msa_confirmation",
            "n_chains": 2,
            "chain_ids": "AB",
            "chain_length": f"{len(sequence)},{len(sequence)}",
            "ova_ranges": f"1-{len(sequence)};1-{len(sequence)}",
            "module_ranges": "0-0;0-0",
            "architecture": name,
            "seeds": args.seeds,
            "expected_samples": args.seeds * 5,
            "description": "OVA homodimer with ColabFold unpaired MSA and no templates",
            "input_json": str(json_path),
        })
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(systems[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(systems)
    print(f"prepared MSA confirmation: {args.out}", flush=True)


if __name__ == "__main__":
    main()
