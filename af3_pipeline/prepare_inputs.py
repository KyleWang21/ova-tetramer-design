#!/usr/bin/env python3
"""Prepare small, MSA-free AF3 experiments for OVA oligomer design."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_tetramer_design import make_constructs  # noqa: E402


GCN4_PLI = "RMKQIEDKLEEILSKLYHIENELARIKKLLGER"
CC_TET = "GELAAIKQELAAIKKELAAIKWELAAIKQGAG"
LINKERS = {
    "f00": "",
    "f05": "GGGGS",
    "f10": "GGGGSGGGGS",
    "f15": "GGGGSGGGGSGGGGS",
    "f20": "GGGGSGGGGSGGGGSGGGGS",
    "h10": "EAAAKEAAAK",
    "h15": "EAAAKEAAAKEAAAK",
}


def protein(chain_id: str, sequence: str, msa_path: str | None = None) -> dict[str, object]:
    payload: dict[str, object] = {
        "id": chain_id,
        "sequence": sequence,
        "pairedMsa": "",
        "templates": [],
    }
    if msa_path:
        payload["unpairedMsaPath"] = msa_path
    else:
        payload["unpairedMsa"] = ""
    return {"protein": payload}


def portable(path: pathlib.Path) -> str:
    return str(path.resolve()).replace(
        "/vepfs-mlp2/c20250508/400083", "/root/400083", 1
    )


def parse_a3m_records(a3m: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for line in a3m.replace("\x00", "").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header, sequence = line, []
        elif header is not None:
            sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records:
        raise ValueError("A3M has no records")
    return records


def combine_domain_a3ms(
    base_a3m: str,
    query: str,
    ova_range: tuple[int, int],
    module_a3m: str | None = None,
    module_range: tuple[int, int] = (0, 0),
) -> str:
    """Combine gap-padded OVA and optional module hits into a fusion A3M."""
    records = parse_a3m_records(base_a3m)
    lo, hi = ova_range
    prefix = lo - 1
    suffix = len(query) - hi
    output = [records[0][0], query]
    for header, aligned in records[1:]:
        output.extend((header, "-" * prefix + aligned + "-" * suffix))
    if module_a3m is not None:
        mod_records = parse_a3m_records(module_a3m)
        mod_lo, mod_hi = module_range
        if module_range == (0, 0):
            raise ValueError("module MSA supplied without module range")
        for header, aligned in mod_records[1:]:
            output.extend((header, "-" * (mod_lo - 1) + aligned + "-" * (len(query) - mod_hi)))
    return "\n".join(output) + "\n"


def baseline_systems() -> list[dict[str, object]]:
    constructs, _, _ = make_constructs()
    by_id = {str(c["id"]): str(c["aa"]) for c in constructs}
    wt = by_id["D0-WT"]
    p3 = by_id["D0-P3-13R"]
    return [
        {"name": "wt_monomer", "kind": "monomer_control", "seqs": [wt],
         "ova_ranges": [(1, len(wt))], "description": "OVA WT single-chain fold control"},
        {"name": "p3_monomer", "kind": "monomer_control", "seqs": [p3],
         "ova_ranges": [(1, len(p3))], "description": "OVA P3-13R single-chain fold control"},
        {"name": "wt_dimer", "kind": "dimer_control", "seqs": [wt, wt],
         "ova_ranges": [(1, len(wt))] * 2, "description": "OVA WT two-chain control"},
        {"name": "p3_dimer", "kind": "dimer_baseline", "seqs": [p3, p3],
         "ova_ranges": [(1, len(p3))] * 2, "description": "Required P3-13R dimer baseline"},
        {"name": "p3_tetramer_nodomain", "kind": "tetramer_negative", "seqs": [p3] * 4,
         "ova_ranges": [(1, len(p3))] * 4, "description": "Four P3 chains without added tetramer domain"},
        {"name": "gcn4pli_tetramer", "kind": "tetramer_positive", "seqs": [GCN4_PLI] * 4,
         "ova_ranges": [(0, 0)] * 4, "module_ranges": [(1, len(GCN4_PLI))] * 4,
         "description": "PDB 1GCL tetrameric coiled-coil positive control"},
        {"name": "cctet_tetramer", "kind": "tetramer_positive", "seqs": [CC_TET] * 4,
         "ova_ranges": [(0, 0)] * 4, "module_ranges": [(1, len(CC_TET))] * 4,
         "description": "Experimentally defined CC-Tet positive control"},
    ]


def screen1_systems() -> list[dict[str, object]]:
    constructs, _, _ = make_constructs()
    by_id = {str(c["id"]): str(c["aa"]) for c in constructs}
    p3 = by_id["D0-P3-13R"]
    systems: list[dict[str, object]] = []

    def add(module_name: str, module: str, terminus: str, linker_name: str) -> None:
        linker = LINKERS[linker_name]
        if terminus == "c":
            chain = p3 + linker + module
            ova_range = (1, len(p3))
            module_range = (len(p3) + len(linker) + 1, len(chain))
            architecture = f"P3-{linker_name}-{module_name}"
        else:
            chain = module + linker + p3
            module_range = (1, len(module))
            ova_range = (len(module) + len(linker) + 1, len(chain))
            architecture = f"{module_name}-{linker_name}-P3"
        name = f"p3_{module_name}_{terminus}_{linker_name}"
        systems.append({
            "name": name,
            "kind": "tetramer_candidate",
            "seqs": [chain] * 4,
            "ova_ranges": [ova_range] * 4,
            "module_ranges": [module_range] * 4,
            "architecture": architecture,
            "description": f"Four-chain {architecture}; {len(chain)} aa per chain",
        })

    # GCN4-pLI: broad C-terminal linker scan plus two N-terminal orientations.
    for linker_name in ("f00", "f05", "f10", "f15", "f20", "h10", "h15"):
        add("pli", GCN4_PLI, "c", linker_name)
    for linker_name in ("f10", "f20"):
        add("pli", GCN4_PLI, "n", linker_name)
    # Independent CC-Tet sequence guards against success being pLI-specific.
    for linker_name in ("f00", "f10", "f20"):
        add("cct", CC_TET, "c", linker_name)
    for linker_name in ("f10", "f20"):
        add("cct", CC_TET, "n", linker_name)
    return systems


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("baseline", "screen1"), default="baseline")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument(
        "--ova-msa", type=pathlib.Path,
        help="Optional OVA A3M to gap-pad across fusion modules (recommended for screen1)",
    )
    parser.add_argument("--pli-msa", type=pathlib.Path, help="Optional pLI domain A3M")
    parser.add_argument("--names", help="Comma-separated system names to prepare")
    args = parser.parse_args()
    if (args.out / "manifest.tsv").exists():
        raise SystemExit(f"experiment already exists: {args.out / 'manifest.tsv'}")
    if args.seeds < 1 or args.shards < 1 or args.shards > 8:
        raise SystemExit("seeds must be positive and shards must be 1..8")

    systems = baseline_systems() if args.phase == "baseline" else screen1_systems()
    if args.names:
        wanted = set(args.names.split(","))
        available_names = {str(system["name"]) for system in systems}
        if missing := wanted - available_names:
            raise SystemExit(f"unknown systems: {sorted(missing)}")
        systems = [system for system in systems if system["name"] in wanted]
    args.out.mkdir(parents=True, exist_ok=True)
    for shard in range(args.shards):
        (args.out / f"in_s{shard}").mkdir()

    rows: list[dict[str, object]] = []
    base_a3m = args.ova_msa.read_text() if args.ova_msa else None
    pli_a3m = args.pli_msa.read_text() if args.pli_msa else None
    msa_dir = args.out / "msas"
    if base_a3m is not None:
        msa_dir.mkdir()
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for index, system in enumerate(systems):
        shard = index % args.shards
        seqs = list(system["seqs"])
        ids = chain_ids[: len(seqs)]
        msa_path = None
        if base_a3m is not None:
            if any(lo == 0 and hi == 0 for lo, hi in system["ova_ranges"]):
                raise ValueError(f"{system['name']} has no OVA range for --ova-msa")
            a3m_path = msa_dir / f"{system['name']}.a3m"
            module_a3m = pli_a3m if "_pli_" in str(system["name"]) else None
            a3m_path.write_text(combine_domain_a3ms(
                base_a3m,
                str(seqs[0]),
                system["ova_ranges"][0],
                module_a3m,
                system.get("module_ranges", [(0, 0)])[0],
            ))
            msa_path = portable(a3m_path)
        payload = {
            "name": system["name"],
            "modelSeeds": list(range(1, args.seeds + 1)),
            "dialect": "alphafold3",
            "version": 2,
            "sequences": [
                protein(cid, seq, msa_path) for cid, seq in zip(ids, seqs, strict=True)
            ],
        }
        json_path = args.out / f"in_s{shard}" / f"{system['name']}.json"
        json_path.write_text(json.dumps(payload, indent=2) + "\n")
        rows.append({
            "name": system["name"],
            "shard": shard,
            "kind": system["kind"],
            "n_chains": len(seqs),
            "chain_ids": ids,
            "chain_length": ",".join(str(len(s)) for s in seqs),
            "ova_ranges": ";".join(f"{lo}-{hi}" for lo, hi in system["ova_ranges"]),
            "module_ranges": ";".join(
                f"{lo}-{hi}" for lo, hi in system.get("module_ranges", [(0, 0)] * len(seqs))
            ),
            "architecture": system.get("architecture", system["name"]),
            "seeds": args.seeds,
            "expected_samples": args.seeds * 5,
            "description": system["description"],
            "input_json": str(json_path),
        })

    fields = list(rows[0])
    with (args.out / "manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    fasta_parts = []
    for system in systems:
        sequence = str(system["seqs"][0])
        fasta_parts.append(f">{system['name']} | {system.get('architecture', system['name'])}\n")
        fasta_parts.extend(sequence[i : i + 80] + "\n" for i in range(0, len(sequence), 80))
    (args.out / "representative_chains.fasta").write_text("".join(fasta_parts))
    print(f"prepared {len(rows)} systems across {args.shards} shards -> {args.out}")
    for row in rows:
        print(f"s{row['shard']}\t{row['name']}\tchains={row['n_chains']}\tseeds={row['seeds']}")


if __name__ == "__main__":
    main()
