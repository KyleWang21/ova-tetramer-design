#!/usr/bin/env python3
"""Prepare full AF3 C4 targets for sparse, joint two-interface tied-AF2 design."""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, Select

from run_final20_failclosed_filter import reference_rsasa, reverse_reference


PROTECTED = {12, 31, 74, 121, 368, 383} | set(range(258, 266))
PAIR_CHOICES = {
    "OVA-C4-IPTM80-19K": ("AB", "AC"),
    "OVA-C4-IPTM80-21B": ("AC", "AD"),
    "OVA-C4-IPTM80-21S": ("AC", "AD"),
    "OVA-C4-IPTM80-20A": ("AB", "AC"),
}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Chains(Select):
    def __init__(self, selected: set[str]) -> None:
        self.selected = selected

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return str(chain.id) in self.selected

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " "


def interface_positions(model, pair: str, cutoff: float = 8.0) -> set[int]:  # noqa: ANN001
    left, right = model[pair[0]], model[pair[1]]
    right_atoms = [atom for residue in right if residue.id[0] == " " for atom in residue if atom.element != "H"]
    search = NeighborSearch(right_atoms)
    positions = set()
    for residue in left:
        if residue.id[0] != " ":
            continue
        for atom in residue:
            if atom.element == "H":
                continue
            partners = search.search(atom.coord, cutoff, level="A")
            if partners:
                positions.add(int(residue.id[1]))
                positions.update(int(other.get_parent().id[1]) for other in partners)
    return positions


def write_positions(path: pathlib.Path, values: set[int]) -> None:
    path.write_text(",".join(map(str, sorted(values))) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--af3-models", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    target_dir = args.out / "targets"; target_dir.mkdir(exist_ok=True)

    final_rows = read_tsv(args.final)
    final = {row["candidate"]: row for row in final_rows}
    reference_sequence = reverse_reference(final_rows[0]["sequence"], final_rows[0]["mutations"])
    surface = {pos for pos, value in reference_rsasa(args.reference, "A", reference_sequence).items() if value >= 0.20}
    write_positions(args.out / "surface_positions.txt", surface)

    models: dict[str, list[dict[str, str]]] = {name: [] for name in PAIR_CHOICES}
    for row in read_tsv(args.af3_models):
        if row["candidate"] in models:
            models[row["candidate"]].append(row)

    manifest = []
    with (args.out / "anchors.fasta").open("w") as fasta:
        for name, pairs in PAIR_CHOICES.items():
            fasta.write(f">{name}\n{final[name]['sequence']}\n")
            rows = models[name]
            main_cluster = max(
                {row["topology_cluster"] for row in rows},
                key=lambda cluster: sum(row["topology_cluster"] == cluster for row in rows),
            )
            eligible = [row for row in rows if row["topology_cluster"] == main_cluster]
            best = max(eligible, key=lambda row: (
                int(row["atomic_clash_count_2p4"]) == 0,
                int(row["valid_c4_network"]), float(row["iptm"]),
                float(row["min_incident_iptm"]), -float(row["interface_pae"]),
            ))
            structure = MMCIFParser(QUIET=True).get_structure(name, best["cif_path"])
            model = structure[0]
            mutation_positions = {int(item[1:-1]) for item in final[name]["mutations"].split(",")}
            target_c4 = target_dir / (name.lower().replace("-", "_") + "_c4.pdb")
            io = PDBIO(); io.set_structure(structure); io.save(str(target_c4), Chains(set("ABCD")))
            record: dict[str, object] = {
                "anchor": name,
                "anchor_n_mutations": final[name]["n_mutations"],
                "source_cif": best["cif_path"],
                "source_iptm": best["iptm"],
                "source_atomic_clashes": best["atomic_clash_count_2p4"],
                "pair1": pairs[0],
                "pair2": pairs[1],
                "target_c4": str(target_c4),
            }
            contact_sets: list[set[int]] = []
            for state_index, pair in enumerate(pairs, 1):
                prefix = name.lower().replace("-", "_") + f"_{pair.lower()}"
                target = target_dir / f"{prefix}.pdb"
                io = PDBIO(); io.set_structure(structure); io.save(str(target), Chains(set(pair)))
                contact_positions = interface_positions(model, pair)
                # All anchor mutations are mutable so buried legacy mutations can
                # revert; new positions are restricted to exposed interface residues.
                design = (mutation_positions | (contact_positions & surface)) - PROTECTED
                restraint = contact_positions - PROTECTED
                design_path = target_dir / f"{prefix}_design.txt"
                restraint_path = target_dir / f"{prefix}_restraint.txt"
                write_positions(design_path, design)
                write_positions(restraint_path, restraint)
                contact_sets.append(contact_positions)
                record.update({
                    f"target{state_index}": str(target),
                    f"design{state_index}": str(design_path),
                    f"restraint{state_index}": str(restraint_path),
                    f"contact_positions{state_index}": len(contact_positions),
                    f"design_positions{state_index}": len(design),
                })
            # The optimization is performed once on the full four-chain target.
            # Both observed interfaces share one mutable mask, while their target
            # contact maps remain separate loss terms.
            joint_design = (
                mutation_positions | ((contact_sets[0] | contact_sets[1]) & surface)
            ) - PROTECTED
            joint_path = target_dir / (name.lower().replace("-", "_") + "_joint_design.txt")
            write_positions(joint_path, joint_design)
            record["joint_design"] = str(joint_path)
            record["joint_design_positions"] = len(joint_design)
            manifest.append(record)

    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        fields = list(manifest[0])
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest)
    print(
        f"prepared {len(manifest)} full C4 targets with "
        f"{2 * len(manifest)} separately restrained interface maps in {args.out}"
    )


if __name__ == "__main__":
    main()
