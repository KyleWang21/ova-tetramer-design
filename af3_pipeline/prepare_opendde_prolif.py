#!/usr/bin/env python3
"""Build a ProLIF manifest from complete, eligible OpenDDE representatives."""

from __future__ import annotations

import argparse
import csv
import pathlib


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-summary", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    final = {row["candidate"] for row in read_tsv(args.final)}
    rows = []
    for row in read_tsv(args.opendde_summary):
        candidate = row["candidate"]
        if candidate not in final or row.get("opendde_evidence_available") != "1":
            continue
        cif = pathlib.Path(row["selected_cif_path"])
        if not cif.exists():
            raise FileNotFoundError(f"{candidate}: selected OpenDDE CIF is missing: {cif}")
        rows.append({
            "source": f"{candidate}__OpenDDE",
            "candidate": candidate,
            "model_system": "OpenDDE-v1",
            "cif_path": str(cif),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["source", "candidate", "model_system", "cif_path"]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} eligible OpenDDE representatives to {args.out}")


if __name__ == "__main__":
    main()
