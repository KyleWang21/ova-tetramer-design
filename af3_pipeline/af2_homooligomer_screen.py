#!/usr/bin/env python3
"""Fast AF2-multimer homooligomer screen used only as an auxiliary AI signal."""

from __future__ import annotations

import argparse
import csv
import pathlib

from colabdesign import mk_afdesign_model


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--params", type=pathlib.Path, required=True)
    ap.add_argument("--copies", type=int, default=4)
    ap.add_argument("--recycles", type=int, default=3)
    ap.add_argument("--model", type=int, default=0)
    ap.add_argument("--limit", type=int, default=8)
    args = ap.parse_args()
    records = read_fasta(args.fasta)[:args.limit]
    if not records or any(len(sequence) != 386 for _, sequence in records):
        raise SystemExit("expected non-empty 386-aa candidate FASTA")
    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"; pdb_dir.mkdir(exist_ok=True)
    af = mk_afdesign_model(
        protocol="hallucination", use_multimer=True, use_templates=False,
        data_dir=str(args.params), num_recycles=args.recycles,
    )
    af.prep_inputs(length=386, copies=args.copies)
    rows: list[dict[str, object]] = []
    for index, (header, sequence) in enumerate(records, 1):
        af.predict(seq=sequence, num_recycles=args.recycles, models=[args.model],
                   seed=index, verbose=False)
        log = af.aux["log"]
        name = f"candidate_{index:02d}"
        af.save_current_pdb(str(pdb_dir / f"{name}.pdb"))
        row = {
            "candidate": name,
            "header": header,
            "plddt": float(log.get("plddt", 0.0)),
            "ptm": float(log.get("ptm", 0.0)),
            "iptm": float(log.get("i_ptm", 0.0)),
            "pae_scaled": float(log.get("pae", 0.0)),
            "interface_pae_scaled": float(log.get("i_pae", 0.0)),
            "sequence": sequence,
        }
        rows.append(row)
        print(index, f"pLDDT={row['plddt']:.3f}", f"pTM={row['ptm']:.3f}",
              f"ipTM={row['iptm']:.3f}", flush=True)
    with (args.out / "af2_multimer_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
