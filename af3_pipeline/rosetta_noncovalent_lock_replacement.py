#!/usr/bin/env python3
"""Sequence-symmetric Rosetta replacement of covalent C4 locks.

The input is an AF3 C4 backbone previously stabilized by three engineered
inter-chain disulfide pairs.  All six engineered Cys are removed before design.
Rosetta's sequence-symmetric annealer then designs identical surface residues in
all four OVA copies.  Allowed identities combine ProteinMPNN, AF2-differentiable,
the parent surface sequence, and the P3-13R reference; Cys is always excluded.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import pyrosetta


LOCK_POSITIONS = {94, 96, 155, 190, 196, 213}
LOCK_ALLOWED = set("AFGHILMNPQSTVWY")


def read_fasta_record(path: pathlib.Path, prefix: str) -> str:
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None and header.startswith(prefix):
                return "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None and header.startswith(prefix):
        return "".join(chunks)
    raise ValueError(f"record {prefix!r} not found")


def read_positions(path: pathlib.Path) -> list[int]:
    return sorted({
        int(item)
        for item in path.read_text().replace("\n", ",").split(",")
        if item.strip()
    })


def read_sequences(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    if path.suffix in {".fa", ".faa", ".fasta"}:
        records: list[str] = []
        chunks: list[str] = []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if chunks:
                    records.append("".join(chunks))
                chunks = []
            else:
                chunks.append(line.strip())
        if chunks:
            records.append("".join(chunks))
        return records
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return []
    for column in ("full_sequence", "sequence"):
        if column in rows[0]:
            return [row[column] for row in rows if row.get(column)]
    raise ValueError(f"no full_sequence or sequence column in {path}")


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def allowed_vector(rosetta, letters: set[str]):
    vector = rosetta.utility.vector1_bool(20)
    for aa in letters:
        vector[int(rosetta.core.chemical.aa_from_oneletter_code(aa))] = True
    return vector


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--design-positions", type=pathlib.Path, required=True)
    parser.add_argument("--sequence-source", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--nstruct", type=int, default=24)
    parser.add_argument("--seed", type=int, default=1001)
    args = parser.parse_args()

    pyrosetta.init("-mute all -ex1 -ex2aro -use_input_sc -constant_seed")
    from pyrosetta import rosetta

    reference = read_fasta_record(args.reference_fasta, args.reference_name)
    design_positions = read_positions(args.design_positions)
    if len(reference) != 386 or not LOCK_POSITIONS.issubset(design_positions):
        raise ValueError("expected 386-aa reference and all six lock positions in design mask")

    pose = pyrosetta.pose_from_pdb(str(args.target_pdb))
    length = len(reference)
    if pose.num_chains() != 4 or any(len(pose.chain_sequence(i)) != length for i in range(1, 5)):
        raise ValueError("expected four complete 386-aa OVA chains")
    parent = pose.chain_sequence(1)

    source_sequences = [reference, parent]
    for path in args.sequence_source:
        source_sequences.extend(read_sequences(path))
    source_sequences = [sequence for sequence in source_sequences if len(sequence) == length]

    allowed_by_position: dict[int, set[str]] = {}
    for position in design_positions:
        if position in LOCK_POSITIONS:
            allowed = set(LOCK_ALLOWED)
        else:
            allowed = {sequence[position - 1] for sequence in source_sequences}
            allowed.update({reference[position - 1], parent[position - 1]})
        allowed.discard("C")
        if not allowed:
            raise ValueError(f"no non-Cys identities at position {position}")
        allowed_by_position[position] = allowed

    # Remove every engineered Cys before any energy calculation.  Ala is only the
    # starting identity; the sequence-symmetric packer can replace it immediately.
    for chain in range(4):
        for position in LOCK_POSITIONS:
            rosetta.protocols.simple_moves.MutateResidue(
                chain * length + position, "ALA"
            ).apply(pose)

    setup = rosetta.protocols.symmetry.SetupForSequenceSymmetryMover()
    for chain in range(4):
        selector = rosetta.core.select.residue_selector.ResidueIndexSelector(
            ",".join(str(chain * length + position) for position in range(1, length + 1))
        )
        setup.add_residue_selector(0, selector)
    setup.apply(pose)

    scorefxn = pyrosetta.get_fa_scorefxn()
    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"
    pdb_dir.mkdir(exist_ok=True)
    rows: list[dict[str, object]] = []

    for index in range(args.nstruct):
        work = pose.clone()
        rosetta.numeric.random.rg().set_seed(args.seed + index)
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(work)
        task.or_include_current(True)
        for pose_index in range(1, work.total_residue() + 1):
            position = (pose_index - 1) % length + 1
            if position in allowed_by_position:
                task.nonconst_residue_task(pose_index).restrict_absent_canonical_aas(
                    allowed_vector(rosetta, allowed_by_position[position])
                )
            else:
                task.nonconst_residue_task(pose_index).prevent_repacking()
        keep = rosetta.core.pack.task.operation.KeepSequenceSymmetry()
        keep.set_setting(True)
        keep.apply(work, task)
        rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task).apply(work)

        sequences = [work.chain_sequence(chain) for chain in range(1, 5)]
        if len(set(sequences)) != 1:
            raise RuntimeError("sequence-symmetric annealer produced non-identical chains")
        sequence = sequences[0]
        if any(sequence[position - 1] == "C" for position in design_positions):
            raise RuntimeError("engineered Cys escaped the Rosetta allowed-AA mask")

        total_score = float(scorefxn(work))
        separated_score = 0.0
        split = work.split_by_chain()
        for chain in range(1, len(split) + 1):
            separated_score += float(scorefxn(split[chain]))
        interface_energy = total_score - separated_score
        candidate = f"rosetta_nc_{index + 1:02d}"
        pdb_path = pdb_dir / f"{candidate}.pdb"
        work.dump_pdb(str(pdb_path))
        row = {
            "candidate": candidate,
            "seed": args.seed + index,
            "total_score": total_score,
            "separated_score": separated_score,
            "interface_energy": interface_energy,
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutation_string(reference, sequence),
            "lock_replacements": ",".join(
                f"{position}{sequence[position - 1]}" for position in sorted(LOCK_POSITIONS)
            ),
            "native_cys_positions": ",".join(
                str(position) for position, aa in enumerate(sequence, 1) if aa == "C"
            ),
            "pdb": str(pdb_path),
            "sequence": sequence,
        }
        rows.append(row)
        print(
            candidate,
            f"dG={interface_energy:.2f}",
            f"mut={row['n_mutations']}",
            row["lock_replacements"],
            flush=True,
        )

    # Deduplicate by sequence, retaining the lowest Rosetta interface energy.
    unique: dict[str, dict[str, object]] = {}
    for row in rows:
        sequence = str(row["sequence"])
        if sequence not in unique or float(row["interface_energy"]) < float(unique[sequence]["interface_energy"]):
            unique[sequence] = row
    ranked = sorted(unique.values(), key=lambda row: (float(row["interface_energy"]), float(row["total_score"])))

    for name, table in [("all", rows), ("unique_ranked", ranked)]:
        with (args.out / f"rosetta_{name}.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(table)
    with (args.out / "rosetta_unique_ranked.fasta").open("w") as handle:
        for rank, row in enumerate(ranked, 1):
            handle.write(
                f">rosetta_nc_rank{rank:02d} source={row['candidate']} "
                f"dG={row['interface_energy']:.3f} mutations={row['n_mutations']}\n"
                f"{row['sequence']}\n"
            )


if __name__ == "__main__":
    main()
