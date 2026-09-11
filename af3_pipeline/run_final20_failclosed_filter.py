#!/usr/bin/env python3
"""Join all OVA-C4 screening stages and apply the current hard gates fail-closed.

The script is intentionally rerunnable: optional stage tables can be omitted while jobs
are running, but every missing required stage is explicitly marked MISSING and fails.
OpenDDE is retained as an optional diagnostic and is never part of the current
``Accepted_v1`` hard-gate set.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re
import statistics
from collections import defaultdict

from Bio import Align
from Bio.PDB import MMCIFParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa


MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
    "SEP": 155.0,
}
NONVDW = {"Hydrophobic", "HydrogenBond", "SaltBridge", "PiStacking", "CationPi"}


def read_tsv(path: pathlib.Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("empty output")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def row_number(row: dict[str, str], *keys: str, default: float = 0.0) -> float:
    for key in keys:
        value = row.get(key, "")
        if value not in {"", None}:
            return float(value)
    return default


def reverse_reference(sequence: str, mutations: str) -> str:
    output = list(sequence)
    for mutation in mutations.split(","):
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mutation)
        if not match:
            raise ValueError(f"invalid mutation: {mutation}")
        old, position, new = match.groups()
        position = int(position)
        if output[position - 1] != new:
            raise ValueError(f"mutation does not match candidate sequence: {mutation}")
        output[position - 1] = old
    return "".join(output)


def glyco_sites(sequence: str) -> set[int]:
    return {index + 1 for index in range(len(sequence) - 2)
            if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in "ST"}


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEP": "S",
}


def reference_rsasa(path: pathlib.Path, chain_id: str, target_sequence: str) -> dict[int, float]:
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    model = structure[0]
    # Surface-design gate is defined on an isolated monomer. 1OVA contains four
    # crystallographic protein chains, so computing SASA on the whole model would
    # incorrectly label the native crystal-contact surface as buried.
    ShrakeRupley(probe_radius=1.4, n_points=200).compute(model[chain_id], level="R")
    records = []
    for residue in model[chain_id]:
        if not is_aa(residue, standard=False) or residue.resname not in MAX_ASA:
            continue
        records.append((AA3_TO_1[residue.resname], min(1.0, float(residue.sasa) / MAX_ASA[residue.resname])))
    crystal_sequence = "".join(record[0] for record in records)
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(crystal_sequence, target_sequence)[0]
    output = {}
    for (crystal_start, crystal_end), (target_start, target_end) in zip(*alignment.aligned):
        if crystal_end - crystal_start != target_end - target_start:
            raise ValueError("unexpected unequal aligned block")
        for crystal_index, target_index in zip(
            range(crystal_start, crystal_end), range(target_start, target_end)
        ):
            output[target_index + 1] = records[crystal_index][1]
    return output


def bool_int(value: bool) -> int:
    return int(bool(value))


def finite(value: str | float | None) -> float | None:
    try:
        output = float(value)
        return output if math.isfinite(output) else None
    except (TypeError, ValueError):
        return None


def ranks(rows: list[dict[str, object]], field: str, higher: bool) -> dict[str, float]:
    values = [(str(row["candidate"]), finite(row.get(field))) for row in rows]
    values = [(candidate, value) for candidate, value in values if value is not None]
    ordered = sorted((value for _, value in values), reverse=higher)
    positions: dict[float, list[int]] = defaultdict(list)
    for index, value in enumerate(ordered, 1):
        positions[value].append(index)
    rank_by_value = {value: statistics.mean(indices) for value, indices in positions.items()}
    return {candidate: float(rank_by_value[value]) for candidate, value in values}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--af3-summary", type=pathlib.Path, required=True)
    parser.add_argument("--af3-models", type=pathlib.Path, required=True)
    parser.add_argument("--prolif-unique", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--reference-chain", default="A")
    parser.add_argument("--protenix-summary", type=pathlib.Path)
    parser.add_argument("--stoichiometry-summary", type=pathlib.Path)
    parser.add_argument("--rosetta-summary", type=pathlib.Path)
    parser.add_argument("--opendde-summary", type=pathlib.Path)
    parser.add_argument("--opendde-prolif-summary", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    final_rows = read_tsv(args.final)
    final = {row["candidate"]: row for row in final_rows}
    reference = reverse_reference(final_rows[0]["sequence"], final_rows[0]["mutations"])
    for row in final_rows[1:]:
        if reverse_reference(row["sequence"], row["mutations"]) != reference:
            raise ValueError(f"{row['candidate']}: mutation list implies a different P3-13R reference")
    rsasa = reference_rsasa(args.reference, args.reference_chain, reference)
    af3 = {row["candidate"]: row for row in read_tsv(args.af3_summary)}
    af3_models = defaultdict(list)
    for row in read_tsv(args.af3_models):
        af3_models[row["candidate"]].append(row)
    protenix = {row["candidate"]: row for row in read_tsv(args.protenix_summary)}
    stoich = {row["candidate"]: row for row in read_tsv(args.stoichiometry_summary)}
    rosetta = {row["candidate"]: row for row in read_tsv(args.rosetta_summary)}
    opendde = {row["candidate"]: row for row in read_tsv(args.opendde_summary)}
    opendde_prolif = {
        row["candidate"]: row for row in read_tsv(args.opendde_prolif_summary)
    }

    interactions = defaultdict(list)
    for event in read_tsv(args.prolif_unique):
        interactions[event["candidate"]].append(event)

    output = []
    for base in final_rows:
        candidate = base["candidate"]
        sequence = base["sequence"]
        mutations = base["mutations"].split(",")
        positions = {int(mutation[1:-1]) for mutation in mutations}
        cys_reference = {index + 1 for index, aa in enumerate(reference) if aa == "C"}
        cys_candidate = {index + 1 for index, aa in enumerate(sequence) if aa == "C"}
        sorted_positions = sorted(positions)
        surface_values = [rsasa.get(position) for position in sorted_positions]
        surface_fraction = sum(value is not None and value >= 0.20 for value in surface_values) / len(positions)
        pg_mutations = [mutation for mutation in mutations if mutation[0] in "PG" or mutation[-1] in "PG"]
        pg_core_mutations = [
            mutation for mutation in pg_mutations
            if rsasa.get(int(mutation[1:-1]), 0.0) < 0.20
        ]
        glyco_change = glyco_sites(sequence) ^ glyco_sites(reference)
        sequence_gates = {
            "gate_seq_length_386": len(sequence) == 386,
            "gate_seq_mutations_le25": len(mutations) <= 25,
            "gate_seq_native_cys_unchanged": cys_candidate == cys_reference,
            "gate_seq_siinfecl_intact": "SIINFEKL" in sequence,
            "gate_seq_standard_aa_only": set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY"),
            "gate_seq_no_glyco_motif_change": not glyco_change,
            "gate_seq_surface_fraction_ge0p80": surface_fraction >= 0.80,
            "gate_seq_no_unapproved_core_pro_gly_change": not pg_core_mutations,
        }

        af = af3.get(candidate)
        af3_pass = bool(af and int(af.get("accepted_af3_recomputed", "0")))

        events = interactions.get(candidate, [])
        any_mutations = set()
        nonvdw_mutations = set()
        chain_support = defaultdict(set)
        nonvdw_strength = defaultdict(int)
        for event in events:
            for suffix in ("a", "b"):
                mutation = event.get(f"mutation_{suffix}", "")
                if not mutation:
                    continue
                any_mutations.add(mutation)
                if event["category"] in NONVDW:
                    nonvdw_mutations.add(mutation)
                    chain_support[mutation].add(event[f"chain_{suffix}"])
                    nonvdw_strength[mutation] += 1

        representative = max(
            af3_models.get(candidate, []),
            key=lambda row: (
                int(row_number(row, "atomic_clash_count_2p4", default=1.0) == 0),
                int(row_number(row, "valid_c4_network") == 1),
                int(row_number(row, "two_interfaces_each_design_ge3") == 1),
                row_number(row, "iptm"),
                row_number(row, "min_incident_iptm", "min_chain_plddt"),
            ),
            default=None,
        )
        type_edges = []
        if representative:
            type_edges = [group.split(",") for group in representative["interface_type_edges"].split("|") if group]
        type_nonvdw = []
        for group in type_edges:
            group_mutations = set()
            for event in events:
                if event["chain_pair"] not in group or event["category"] not in NONVDW:
                    continue
                group_mutations.update(filter(None, (event.get("mutation_a", ""), event.get("mutation_b", ""))))
            type_nonvdw.append(group_mutations)
        major = sorted(nonvdw_mutations, key=lambda mutation: (-nonvdw_strength[mutation], mutation))[:4]
        prolif_gates = {
            "gate_prolif_design_any_fraction_ge0p50": len(any_mutations) / len(mutations) >= 0.50,
            "gate_prolif_design_nonvdw_fraction_ge0p25": len(nonvdw_mutations) / len(mutations) >= 0.25,
            "gate_prolif_design_nonvdw_count_ge6": len(nonvdw_mutations) >= 6,
            "gate_prolif_two_interfaces_nonvdw_ge2": sum(
                len(mutations) >= 2 for mutations in type_nonvdw
            ) >= 2,
            "gate_prolif_major_hotspots_symmetry_ge3of4": bool(major) and all(len(chain_support[item]) >= 3 for item in major),
        }

        pt = protenix.get(candidate)
        ptx_pass = bool(pt and int(pt.get("crossmodel_validated", pt.get("accepted_protenix", "0"))))
        st = stoich.get(candidate)
        stoich_pass = bool(st and int(st.get("stoichiometry_pass", "0")))
        ros = rosetta.get(candidate)
        rosetta_pass = bool(ros and int(ros.get("rosetta_pass", "0")))
        opn = opendde.get(candidate)
        opn_prolif = opendde_prolif.get(candidate)
        if opn and opn.get("opendde_evidence_available") == "1":
            opendde_status = "PASS"
        elif opn and opn.get("gate_opendde_10_seeds") == "1":
            opendde_status = "COMPLETE_NO_ELIGIBLE_REP"
        elif opn:
            opendde_status = "INCOMPLETE"
        else:
            opendde_status = "MISSING"

        sequence_pass = all(sequence_gates.values())
        prolif_pass = all(prolif_gates.values())
        required_stages = {
            "sequence": sequence_pass,
            "af3": af3_pass,
            "prolif": prolif_pass,
            "protenix": ptx_pass,
            "stoichiometry": stoich_pass,
            "rosetta": rosetta_pass,
            # OpenDDE is deliberately excluded from the required set. Its status
            # and metrics remain in the output for independent diagnostics.
        }
        failed_stages = [stage for stage, passed in required_stages.items() if not passed]
        row: dict[str, object] = {
            "candidate": candidate,
            "original_af3_rank": base["rank"],
            "n_mutations": len(mutations),
            "mutations": base["mutations"],
            "surface_mutation_fraction_rsasa_ge0p20": surface_fraction,
            "surface_mutation_missing_positions": ",".join(str(pos) for pos, value in zip(sorted_positions, surface_values) if value is None),
            "pro_gly_mutations": ",".join(pg_mutations),
            "core_pro_gly_mutations_requiring_review": ",".join(pg_core_mutations),
            "glyco_motif_changed_positions": ",".join(map(str, sorted(glyco_change))),
            **{key: bool_int(value) for key, value in sequence_gates.items()},
            "sequence_pass": bool_int(sequence_pass),
            **({f"af3_{key}": value for key, value in af.items() if key != "candidate"} if af else {}),
            "af3_status": "PASS" if af3_pass else ("FAIL" if af else "MISSING"),
            "prolif_design_any_count": len(any_mutations),
            "prolif_design_any_fraction": len(any_mutations) / len(mutations),
            "prolif_design_nonvdw_count": len(nonvdw_mutations),
            "prolif_design_nonvdw_fraction": len(nonvdw_mutations) / len(mutations),
            "prolif_interface_type_nonvdw_counts": ",".join(map(str, map(len, type_nonvdw))),
            "prolif_major_hotspots": ",".join(major),
            "prolif_major_hotspot_chain_support": ",".join(f"{item}:{len(chain_support[item])}" for item in major),
            **{key: bool_int(value) for key, value in prolif_gates.items()},
            "prolif_pass": bool_int(prolif_pass),
            **({f"protenix_{key}": value for key, value in pt.items() if key != "candidate"} if pt else {}),
            "protenix_status": "PASS" if ptx_pass else ("FAIL" if pt else "MISSING"),
            **({f"stoich_{key}": value for key, value in st.items() if key != "candidate"} if st else {}),
            "stoichiometry_status": "PASS" if stoich_pass else ("FAIL" if st else "MISSING"),
            **({f"rosetta_{key}": value for key, value in ros.items() if key != "candidate"} if ros else {}),
            "rosetta_status": "PASS" if rosetta_pass else ("FAIL" if ros else "MISSING"),
            **({f"opendde_{key}": value for key, value in opn.items() if key != "candidate"} if opn else {}),
            **({key: value for key, value in opn_prolif.items() if key != "candidate"} if opn_prolif else {}),
            "opendde_status": opendde_status,
            "all_required_stages_pass": bool_int(not failed_stages),
            "failed_required_stages": ";".join(failed_stages),
            "sequence": sequence,
        }
        output.append(row)

    # Provisional evidence ranking never overrides hard gates. It uses only currently
    # available, directly comparable evidence and is labeled accordingly.
    rank_specs = [
        ("af3_mean_iptm", True, 0.35),
        ("af3_mean_weakest_interface_iptm", True, 0.20),
        ("af3_mean_interface_pae", False, 0.10),
        ("af3_mean_design_interface_recurrence", True, 0.15),
        ("prolif_design_nonvdw_count", True, 0.10),
        ("n_mutations", False, 0.10),
    ]
    rank_maps = [(ranks(output, field, higher), weight) for field, higher, weight in rank_specs]
    for row in output:
        values = [(mapping[str(row["candidate"])], weight) for mapping, weight in rank_maps
                  if str(row["candidate"]) in mapping]
        row["provisional_available_evidence_borda"] = sum(value * weight for value, weight in values) / sum(weight for _, weight in values)
    output.sort(key=lambda row: (float(row["provisional_available_evidence_borda"]), int(row["original_af3_rank"])))
    for rank, row in enumerate(output, 1):
        row["provisional_rank_not_final"] = rank

    # Normalize heterogeneous optional-stage keys so DictWriter remains deterministic.
    fields = []
    for row in output:
        for key in row:
            if key not in fields:
                fields.append(key)
    normalized = [{key: row.get(key, "") for key in fields} for row in output]
    write_tsv(args.out / "final20_failclosed_screen.tsv", normalized)
    pass_rows = [row for row in normalized if int(row["all_required_stages_pass"])]
    accepted_path = args.out / "accepted_final.tsv"
    if pass_rows:
        write_tsv(accepted_path, pass_rows)
    else:
        # A rerun can turn a formerly accepted candidate into a fail once a
        # required stage arrives.  Never leave a stale accepted table behind.
        accepted_path.unlink(missing_ok=True)
    summary = [
        {"item": "candidates", "value": len(output)},
        {"item": "sequence_pass", "value": sum(int(row["sequence_pass"]) for row in output)},
        {"item": "af3_pass", "value": sum(row["af3_status"] == "PASS" for row in output)},
        {"item": "prolif_pass", "value": sum(int(row["prolif_pass"]) for row in output)},
        {"item": "protenix_pass", "value": sum(row["protenix_status"] == "PASS" for row in output)},
        {"item": "stoichiometry_pass", "value": sum(row["stoichiometry_status"] == "PASS" for row in output)},
        {"item": "rosetta_pass", "value": sum(row["rosetta_status"] == "PASS" for row in output)},
        {"item": "opendde_optional_evidence", "value": sum(row["opendde_status"] == "PASS" for row in output)},
        {"item": "all_required_pass", "value": len(pass_rows)},
    ]
    write_tsv(args.out / "screen_counts.tsv", summary)
    print("; ".join(f"{row['item']}={row['value']}" for row in summary))


if __name__ == "__main__":
    main()
