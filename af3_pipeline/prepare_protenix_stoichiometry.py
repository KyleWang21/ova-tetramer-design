#!/usr/bin/env python3
"""Make Protenix count controls from an existing JSON input."""

from __future__ import annotations

import argparse
import copy
import json
import pathlib


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    parser.add_argument("--counts", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument("--fasta", type=pathlib.Path, help="replace the protein sequence")
    parser.add_argument("--msa-path", type=pathlib.Path, help="replace unpairedMsaPath")
    args = parser.parse_args()
    payloads = json.loads(args.input.read_text())
    if len(payloads) != 1 or len(payloads[0]["sequences"]) != 1:
        raise SystemExit("expected one job containing one counted proteinChain")
    if args.fasta:
        sequence = "".join(
            line.strip() for line in args.fasta.read_text().splitlines()
            if line.strip() and not line.startswith(">")
        )
        payloads[0]["sequences"][0]["proteinChain"]["sequence"] = sequence
    if args.msa_path:
        payloads[0]["sequences"][0]["proteinChain"]["unpairedMsaPath"] = str(
            args.msa_path.resolve()
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for count in args.counts:
        payload = copy.deepcopy(payloads)
        payload[0]["name"] = f"ova_nc_c{count}_protenix_control"
        payload[0]["sequences"][0]["proteinChain"]["count"] = count
        path = args.out_dir / f"ova_nc_c{count}_protenix_control.json"
        path.write_text(json.dumps(payload, indent=2) + "\n")
        print(path)


if __name__ == "__main__":
    main()
