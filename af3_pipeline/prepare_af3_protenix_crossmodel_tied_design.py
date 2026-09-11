#!/usr/bin/env python3
"""Prepare tied-AF2 states spanning independent AF3 and Protenix C4 models.

Each of four records combines one designed interface from a zero-clash AF3
AI-MULTIHIT-11 backbone with the complementary interface from the zero-clash
Protenix-v2 seed-51 prediction.  The two states share one sequence mask; AF3
remains the downstream template-free acceptance test.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, Select


PROTECTED = {12, 31, 74, 121, 368, 383} | set(range(258, 266))


def read_positions(path: pathlib.Path) -> set[int]:
    return {
        int(token)
        for token in path.read_text().replace("\n", ",").split(",")
        if token.strip()
    }


def write_positions(path: pathlib.Path, positions: set[int]) -> None:
    path.write_text(",".join(map(str, sorted(positions))) + "\n")


def fasta_records(path: pathlib.Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            fields = line[1:].split()
            name, chunks = (fields[0] if fields else None), []
        elif name is not None:
            chunks.append(line.strip())
    return records


class Chains(Select):
    def __init__(self, selected: set[str]) -> None:
        self.selected = selected

    def accept_chain(self, chain) -> bool:  # noqa: ANN001
        return str(chain.id) in self.selected

    def accept_residue(self, residue) -> bool:  # noqa: ANN001
        return residue.id[0] == " "


def validate_model(model) -> None:  # noqa: ANN001
    chain_ids = {str(chain.id) for chain in model}
    if not set("ABCD").issubset(chain_ids):
        raise ValueError(f"expected chains A-D, found {sorted(chain_ids)}")
    for chain_id in "ABCD":
        residues = [residue for residue in model[chain_id] if residue.id[0] == " "]
        if len(residues) != 386:
            raise ValueError(f"chain {chain_id} has {len(residues)} residues, expected 386")


def interface_positions(model, pair: str, cutoff: float = 8.0) -> set[int]:  # noqa: ANN001
    left, right = model[pair[0]], model[pair[1]]
    right_atoms = [
        atom
        for residue in right
        if residue.id[0] == " "
        for atom in residue
        if atom.element != "H"
    ]
    search = NeighborSearch(right_atoms)
    positions: set[int] = set()
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


def save_chains(structure, path: pathlib.Path, chain_ids: str) -> None:  # noqa: ANN001
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), Chains(set(chain_ids)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--af3-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--anchors-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--protenix-cif", type=pathlib.Path, required=True)
    parser.add_argument("--protenix-iptm", type=float, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    target_dir = args.out / "targets"
    target_dir.mkdir(exist_ok=True)
    rows = list(csv.DictReader(args.af3_manifest.open(), delimiter="\t"))
    if len(rows) != 4:
        raise ValueError(f"expected four AF3 backbone rows, found {len(rows)}")
    anchors = fasta_records(args.anchors_fasta)
    reference = fasta_records(args.reference_fasta)[args.reference_name]
    if len(reference) != 386:
        raise ValueError(f"reference length is {len(reference)}, expected 386")
    surface = read_positions(args.surface_positions)
    write_positions(args.out / "surface_positions.txt", surface)

    ptx_structure = MMCIFParser(QUIET=True).get_structure("protenix_seed51", args.protenix_cif)
    ptx_model = ptx_structure[0]
    validate_model(ptx_model)
    save_chains(ptx_structure, target_dir / "protenix_seed51_c4.pdb", "ABCD")
    manifest: list[dict[str, object]] = []

    with (args.out / "anchors.fasta").open("w") as fasta:
        for index, row in enumerate(rows, 1):
            anchor = row["anchor"]
            sequence = anchors[anchor]
            if len(sequence) != 386:
                raise ValueError(f"{anchor} length is {len(sequence)}, expected 386")
            fasta.write(f">{anchor}-XM{index}\n{sequence}\n")
            output_anchor = f"{anchor}-XM{index}"
            mutations = {
                position
                for position, (native, designed) in enumerate(zip(reference, sequence, strict=True), 1)
                if native != designed
            }
            af3_structure = MMCIFParser(QUIET=True).get_structure(anchor, row["source_cif"])
            af3_model = af3_structure[0]
            validate_model(af3_model)

            af3_pair = "AC" if index % 2 else "AD"
            ptx_pair = "AD" if af3_pair == "AC" else "AC"
            stem = f"ai11_xm_bb{index}"
            af3_target = target_dir / f"{stem}_af3_{af3_pair.lower()}.pdb"
            ptx_target = target_dir / f"{stem}_ptx51_{ptx_pair.lower()}.pdb"
            target_c4 = target_dir / f"{stem}_af3_c4.pdb"
            save_chains(af3_structure, af3_target, af3_pair)
            save_chains(ptx_structure, ptx_target, ptx_pair)
            save_chains(af3_structure, target_c4, "ABCD")

            af3_contacts = interface_positions(af3_model, af3_pair) - PROTECTED
            ptx_contacts = interface_positions(ptx_model, ptx_pair) - PROTECTED
            if not af3_contacts or not ptx_contacts:
                raise ValueError(f"empty interface for backbone {index}")
            joint_design = (mutations | ((af3_contacts | ptx_contacts) & surface)) - PROTECTED
            af3_design = (mutations | (af3_contacts & surface)) - PROTECTED
            ptx_design = (mutations | (ptx_contacts & surface)) - PROTECTED
            af3_design_path = target_dir / f"{stem}_af3_design.txt"
            ptx_design_path = target_dir / f"{stem}_ptx51_design.txt"
            af3_restraint_path = target_dir / f"{stem}_af3_restraint.txt"
            ptx_restraint_path = target_dir / f"{stem}_ptx51_restraint.txt"
            joint_path = target_dir / f"{stem}_joint_design.txt"
            write_positions(af3_design_path, af3_design)
            write_positions(ptx_design_path, ptx_design)
            write_positions(af3_restraint_path, af3_contacts)
            write_positions(ptx_restraint_path, ptx_contacts)
            write_positions(joint_path, joint_design)
            manifest.append({
                "anchor": output_anchor,
                "anchor_n_mutations": len(mutations),
                "source_cif": row["source_cif"],
                "source_iptm": min(float(row["source_iptm"]), args.protenix_iptm),
                "source_atomic_clashes": 0,
                "pair1": af3_pair,
                "pair2": ptx_pair,
                "target_c4": str(target_c4),
                "target1": str(af3_target),
                "design1": str(af3_design_path),
                "restraint1": str(af3_restraint_path),
                "contact_positions1": len(af3_contacts),
                "design_positions1": len(af3_design),
                "target2": str(ptx_target),
                "design2": str(ptx_design_path),
                "restraint2": str(ptx_restraint_path),
                "contact_positions2": len(ptx_contacts),
                "design_positions2": len(ptx_design),
                "joint_design": str(joint_path),
                "joint_design_positions": len(joint_design),
            })

    with (args.out / "target_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(
        f"prepared {len(manifest)} AF3/Protenix cross-model tied states in {args.out}; "
        f"joint masks={','.join(str(row['joint_design_positions']) for row in manifest)}"
    )


if __name__ == "__main__":
    main()
