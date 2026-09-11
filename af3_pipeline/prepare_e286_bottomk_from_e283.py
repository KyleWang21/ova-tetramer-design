#!/usr/bin/env python3
"""Select eight exact E283 states for bottom-k robust multistate continuation."""

from __future__ import annotations

import csv
import pathlib
import re
import shutil


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
ORIGINAL_HOTSPOTS = {99, 155, 157, 182, 184, 334, 336, 338}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def fasta_sequence(path: pathlib.Path, prefix: str) -> str:
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None and name.startswith(prefix):
                return "".join(chunks)
            name, chunks = line[1:], []
        elif name is not None:
            chunks.append(line.strip())
    raise ValueError(f"missing FASTA record {prefix} in {path}")


def main() -> None:
    project = pathlib.Path(__file__).resolve().parents[1]
    source = project / "experiments/e283_c4_chainmapped_eightstate_tied_design"
    out = project / "experiments/e286_c4_chainmapped_bottomk_tied_design"
    if (out / "INPUT_FROZEN").exists():
        print(f"{out} already frozen")
        return
    if not (source / "REMOTE_ARCHIVE_RECOVERED").exists():
        raise RuntimeError("E283 complete archive has not been recovered")
    if out.exists():
        raise RuntimeError(f"refusing to overwrite partial {out}")

    state_pattern = re.compile(r"state(\d+)_iptm$")
    candidates: dict[str, dict[str, object]] = {}
    result_paths = sorted(source.glob("run_chainmapped_*/**/trajectory.tsv"))
    result_paths += sorted(source.glob("run_chainmapped_*/**/af2diff_candidates.tsv"))
    for path in result_paths:
        for row in read_tsv(path):
            if path.name == "trajectory.tsv":
                if row.get("stage") != "hard" or float(row.get("hard", 0.0)) < 0.999:
                    continue
                source_stage = "hard"
                source_iteration = int(row["iteration"])
                origin_suffix = f"hard{source_iteration}"
            else:
                # af2diff_candidates.tsv is a deterministic, dropout-free
                # reevaluation of the exact argmax sequence selected from the
                # trajectory.  It is the authoritative cross-backbone score
                # for choosing continuation anchors.
                source_stage = "deterministic_final"
                source_iteration = int(row.get("selected_iteration", 0))
                origin_suffix = (
                    f"final_from_{row.get('selected_stage', 'unknown')}"
                    f"{source_iteration}"
                )
            if int(row["n_mutations"]) > 20 or float(row["surface_mutation_fraction"]) < 0.80:
                continue
            sequence = row["sequence"]
            state_columns = sorted(
                (int(match.group(1)), field)
                for field in row
                if (match := state_pattern.fullmatch(field))
            )
            if len(state_columns) != 8:
                raise RuntimeError(f"expected eight state ipTM columns in {path}")
            iptms = [float(row[field]) for _, field in state_columns]
            plddts = [float(row[f"state{index}_plddt"]) for index, _ in state_columns]
            record: dict[str, object] = {
                "origin": f"{path.relative_to(source)}/{origin_suffix}",
                "source_stage": source_stage,
                "source_iteration": source_iteration,
                "n_mutations": int(row["n_mutations"]),
                "surface_mutation_fraction": float(row["surface_mutation_fraction"]),
                "mutations": row["mutations"],
                "source_min_iptm": min(iptms),
                "source_mean_iptm": sum(iptms) / len(iptms),
                "source_min_plddt": min(plddts),
                "sequence": sequence,
            }
            old = candidates.get(sequence)
            score = (
                float(record["source_min_iptm"]),
                float(record["source_mean_iptm"]),
                float(record["source_min_plddt"]),
                -int(record["n_mutations"]),
            )
            if old is None or score > (
                float(old["source_min_iptm"]),
                float(old["source_mean_iptm"]),
                float(old["source_min_plddt"]),
                -int(old["n_mutations"]),
            ):
                candidates[sequence] = record

    ranked = sorted(
        candidates.values(),
        key=lambda row: (
            -float(row["source_min_iptm"]),
            -float(row["source_mean_iptm"]),
            -float(row["source_min_plddt"]),
            int(row["n_mutations"]),
            str(row["sequence"]),
        ),
    )
    selected: list[dict[str, object]] = []
    for minimum_distance in (3, 2, 1):
        for row in ranked:
            if row in selected:
                continue
            if all(
                hamming(str(row["sequence"]), str(old["sequence"])) >= minimum_distance
                for old in selected
            ):
                selected.append(row)
            if len(selected) == 8:
                break
        if len(selected) == 8:
            break

    # A repeated hard endpoint across tasks is possible.  Preserve continuity
    # by filling any missing anchors with unused, already-audited E283 inputs.
    reference = fasta_sequence(project / "OVA_P3-13R_四聚体候选_AA.fasta", "D0-P3-13R")
    surface = {
        int(value)
        for value in (source / "surface_positions.txt").read_text().replace("\n", ",").split(",")
        if value.strip()
    }
    for old in read_tsv(source / "anchor_manifest.tsv"):
        if len(selected) == 8:
            break
        sequence = old["sequence"]
        if any(str(row["sequence"]) == sequence for row in selected):
            continue
        changed = [i for i, (a, b) in enumerate(zip(reference, sequence, strict=True), 1) if a != b]
        surface_fraction = sum(i in surface for i in changed) / len(changed)
        if (
            not 1 <= len(changed) <= 20
            or surface_fraction < 0.80
            or {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS
        ):
            raise RuntimeError(
                f"fallback anchor {old['anchor']} violates frozen constraints: "
                f"nmut={len(changed)} surface={surface_fraction:.3f}"
            )
        selected.append({
            "origin": f"fallback:{old['anchor']}:{old['origin']}",
            "source_stage": "audited_input_fallback",
            "source_iteration": 0,
            "n_mutations": len(changed),
            "surface_mutation_fraction": surface_fraction,
            "mutations": old["mutations"],
            "source_min_iptm": -1.0,
            "source_mean_iptm": -1.0,
            "source_min_plddt": -1.0,
            "sequence": sequence,
        })
    if len(selected) != 8 or len({str(row["sequence"]) for row in selected}) != 8:
        raise RuntimeError(f"could not select eight unique continuation anchors: {len(selected)}")

    out.mkdir(parents=True)
    shutil.copy2(source / "surface_positions.txt", out / "surface_positions.txt")
    shutil.copy2(source / "joint_design.txt", out / "joint_design.txt")
    shutil.copy2(source / "joint_design.txt", out / "joint_design_surface.txt")
    shutil.copytree(source / "targets", out / "targets")
    state_rows = read_tsv(source / "state_manifest.tsv")
    source_prefix = "experiments/e283_c4_chainmapped_eightstate_tied_design/targets/"
    target_prefix = "experiments/e286_c4_chainmapped_bottomk_tied_design/targets/"
    for row in state_rows:
        row["target_pdb"] = row["target_pdb"].replace(source_prefix, target_prefix)
        row["restraint_positions"] = row["restraint_positions"].replace(source_prefix, target_prefix)
    write_tsv(out / "state_manifest.tsv", state_rows)
    surface_design = {
        int(value)
        for value in (out / "joint_design_surface.txt").read_text().replace("\n", ",").split(",")
        if value.strip()
    }
    physically_reproduced_hotspots: set[int] | None = None
    for state_index, row in enumerate(state_rows, 1):
        # Odd states are the canonical first interface.  Retain only hotspots
        # observed in every one of the four independent first-interface maps.
        if state_index % 2 == 0:
            continue
        positions = {
            int(value)
            for value in (project / row["restraint_positions"]).read_text().replace("\n", ",").split(",")
            if value.strip()
        }
        if physically_reproduced_hotspots is None:
            physically_reproduced_hotspots = positions & ORIGINAL_HOTSPOTS
        else:
            physically_reproduced_hotspots &= positions
    if physically_reproduced_hotspots is None:
        raise RuntimeError("no canonical first-interface states found")
    hotspot_expanded = surface_design | physically_reproduced_hotspots
    (out / "joint_design_hotspot.txt").write_text(
        ",".join(map(str, sorted(hotspot_expanded))) + "\n"
    )
    (out / "hotspot_expansion_audit.tsv").write_text(
        "original_hotspots\tall_four_first_interface_hotspots\tnewly_enabled_hotspots\n"
        f"{','.join(map(str, sorted(ORIGINAL_HOTSPOTS)))}\t"
        f"{','.join(map(str, sorted(physically_reproduced_hotspots)))}\t"
        f"{','.join(map(str, sorted(physically_reproduced_hotspots - surface_design)))}\n"
    )
    output_rows: list[dict[str, object]] = []
    with (out / "anchors.fasta").open("w") as handle:
        for index, row in enumerate(selected, 1):
            sequence = str(row["sequence"])
            if len(sequence) != 386 or {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
                raise RuntimeError(f"invalid selected anchor {index}")
            output = {"anchor": f"E286-A{index:02d}", **row}
            output_rows.append(output)
            handle.write(f">E286-A{index:02d} {row['origin']}\n{sequence}\n")
    write_tsv(out / "anchor_manifest.tsv", output_rows)
    (out / "INPUT_FROZEN").touch()
    print(
        f"prepared E286 anchors=8 deterministic_final_sources="
        f"{sum(row['source_stage'] == 'deterministic_final' for row in output_rows)} "
        f"hard_sources={sum(row['source_stage'] == 'hard' for row in output_rows)} "
        f"best_min_iptm={max(float(row['source_min_iptm']) for row in output_rows):.3f} "
        f"hotspot_expansion={sorted(physically_reproduced_hotspots - surface_design)}"
    )


if __name__ == "__main__":
    main()
