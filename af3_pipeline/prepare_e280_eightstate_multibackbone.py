#!/usr/bin/env python3
"""Prepare eight diverse low-mutation anchors and eight AF3-derived dimer states."""

from __future__ import annotations

import csv
import hashlib
import pathlib
import re
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}


def fasta_sequence(path: pathlib.Path, prefix: str) -> str:
    current = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines() + [">"]:
        if raw.startswith(">"):
            if current is not None and current.startswith(prefix):
                return "".join(chunks)
            current, chunks = raw[1:], []
        elif current is not None:
            chunks.append(raw.strip())
    raise ValueError(f"missing FASTA record {prefix} in {path}")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
        if old != new
    )


def glyco_sites(sequence: str) -> tuple[int, ...]:
    return tuple(
        index + 1 for index in range(len(sequence) - 2)
        if sequence[index] == "N" and sequence[index + 1] != "P"
        and sequence[index + 2] in "ST"
    )


def best_e265(source: pathlib.Path) -> dict[str, str]:
    eligible: list[dict[str, str]] = []
    for path in sorted(source.glob("run_cons2_*/**/trajectory.tsv")):
        for row in read_tsv(path):
            if row.get("stage") != "hard" or float(row.get("soft", 0.0)) < 0.999:
                continue
            if int(row["n_mutations"]) > 20 or float(row["surface_mutation_fraction"]) < 0.8:
                continue
            row = dict(row)
            row["source_path"] = str(path.relative_to(source))
            eligible.append(row)
    if not eligible:
        raise RuntimeError("E265 has no eligible fully discrete state")
    return max(
        eligible,
        key=lambda row: (
            min(float(row["state1_iptm"]), float(row["state2_iptm"])),
            (float(row["state1_iptm"]) + float(row["state2_iptm"])) / 2,
            min(float(row["state1_plddt"]), float(row["state2_plddt"])),
            -int(row["n_mutations"]),
        ),
    )


def main() -> None:
    project = pathlib.Path(__file__).resolve().parents[1]
    out = project / "experiments/e280_c4_eightstate_multibackbone_tied_design"
    frozen = out / "INPUT_FROZEN"
    if frozen.exists():
        print(f"{out} already frozen")
        return

    reference = fasta_sequence(project / "OVA_P3-13R_四聚体候选_AA.fasta", "D0-P3-13R")
    source = project / "experiments/e268_c4_e251_softbasin_stabilization_tied_design"
    surface = set(map(int, re.findall(r"\d+", (source / "surface_positions.txt").read_text())))
    registry = read_tsv(project / "experiments/current_accepted_af3_registry/accepted_af3_registry.tsv")
    accepted = {row["candidate"]: row["sequence"] for row in registry}
    e268 = read_tsv(
        project / "experiments/e277_c4_e268_hardhit_crossbackbone_tied_design/source_selection.tsv"
    )[0]
    e265 = best_e265(project / "experiments/e265_c4_crossbackbone_consensus_tied_design")
    e271_rows = {
        int(row["e271_selection_rank"]): row
        for row in read_tsv(
            project / "experiments/e271_c4_multibackbone_exact_mpnn_score/multibackbone_exact_mpnn_shortlist.tsv"
        )
    }
    sources = [
        ("E280-A01", "Accepted_AF3:OVA-C4-073", accepted["ova_body_c4_073_ms"]),
        ("E280-A02", "Accepted_AF3:AI-MULTIHIT-11", accepted["ova_body_c4_multihit11_ms"]),
        ("E280-A03", "E268:best-hard77", e268["sequence"]),
        ("E280-A04", f"E265:{e265['source_path']}/hard{e265['iteration']}", e265["sequence"]),
    ]
    for anchor_index, rank in enumerate((1, 7, 10, 12), 5):
        row = e271_rows[rank]
        sources.append((f"E280-A{anchor_index:02d}", f"E271:rank{rank}:{row['source_name']}", row["sequence"]))

    if len({sequence for _, _, sequence in sources}) != 8:
        raise RuntimeError("E280 anchors are not eight unique sequences")
    audit_rows: list[dict[str, object]] = []
    for anchor, origin, sequence in sources:
        if len(sequence) != len(reference):
            raise RuntimeError(f"{anchor} has wrong length")
        changed = [
            index for index, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
            if old != new
        ]
        surface_fraction = sum(index in surface for index in changed) / len(changed)
        if len(changed) > 20 or surface_fraction < 0.8:
            raise RuntimeError(f"{anchor} violates mutation gates")
        if {index for index, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            raise RuntimeError(f"{anchor} violates native-only Cys")
        if sequence[257:265] != "SIINFEKL" or glyco_sites(sequence) != glyco_sites(reference):
            raise RuntimeError(f"{anchor} violates protected sequence features")
        audit_rows.append({
            "anchor": anchor,
            "origin": origin,
            "n_mutations": len(changed),
            "surface_mutation_fraction": surface_fraction,
            "mutations": mutation_string(reference, sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "sequence": sequence,
        })

    context_rows = read_tsv(source / "target_manifest.tsv")
    if len(context_rows) != 4:
        raise RuntimeError("expected four independent backbone contexts")
    state_rows: list[dict[str, object]] = []
    design_positions: set[int] = set()
    for context_index, row in enumerate(context_rows, 1):
        design_positions.update(map(int, re.findall(r"\d+", (project / row["joint_design"]).read_text())))
        for interface_index, suffix in enumerate(("1", "2"), 1):
            state_rows.append({
                "state_index": len(state_rows) + 1,
                "context_index": context_index,
                "source_cif": row["source_cif"],
                "source_iptm": row["source_iptm"],
                "interface": row[f"pair{suffix}"],
                "target_pdb": row[f"target{suffix}"],
                "restraint_positions": row[f"restraint{suffix}"],
                "default_weight": 3.0,
            })
    if len(state_rows) != 8:
        raise RuntimeError("expected eight interface states")
    protected = NATIVE_CYS | set(range(258, 266))
    if design_positions & protected:
        raise RuntimeError("union design mask overlaps protected positions")

    out.mkdir(parents=True, exist_ok=False)
    with (out / "anchors.fasta").open("w") as handle:
        for row in audit_rows:
            handle.write(f">{row['anchor']} {row['origin']}\n{row['sequence']}\n")
    with (out / "anchor_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(audit_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit_rows)
    with (out / "state_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(state_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(state_rows)
    (out / "joint_design.txt").write_text(",".join(map(str, sorted(design_positions))) + "\n")
    shutil.copy2(source / "surface_positions.txt", out / "surface_positions.txt")
    frozen.touch()
    print(
        f"prepared E280 anchors=8 states=8 design_positions={len(design_positions)} "
        f"mutation_range={min(int(row['n_mutations']) for row in audit_rows)}-"
        f"{max(int(row['n_mutations']) for row in audit_rows)}"
    )


if __name__ == "__main__":
    main()
