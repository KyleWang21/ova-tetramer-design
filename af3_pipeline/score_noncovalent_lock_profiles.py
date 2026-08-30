#!/usr/bin/env python3
"""Score chemically diverse non-Cys replacements of the three directed C4 locks."""

from __future__ import annotations

import argparse
import csv
import pathlib

import pyrosetta


POSITIONS = [94, 96, 155, 190, 196, 213]
SAFE_DISULFIDE_REMOVAL_ORDER = [96, 196, 213, 94, 155, 190]
PROFILES = {
    "ros_aanaaa": "AANAAA",
    "mpnn_pnafma": "PNAFMA",
    "mpnn_pnafse": "PNAFSE",
    "mpnn_pnnfma": "PNNFMA",
    "polar_nnnqqq": "NNNQQQ",
    "polar_qqqnnn": "QQQNNN",
    "polar_ttntst": "TTNTST",
    "polar_ssnsss": "SSNSSS",
    "hydrophobic_vvnava": "VVNAVA",
    "hydrophobic_aanvav": "AANVAV",
    "salt_kkkeee": "KKKEEE",
    "salt_eeekkk": "EEEKKK",
    "salt_rrreee": "RRREEE",
    "salt_eeerrr": "EEERRR",
}


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
    raise ValueError(prefix)


def read_first_fasta(path: pathlib.Path) -> str:
    chunks: list[str] = []
    started = False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if started:
                break
            started = True
        elif started:
            chunks.append(line.strip())
    if not chunks:
        raise ValueError(f"no FASTA sequence in {path}")
    return "".join(chunks)


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b
    )


def one_aa_vector(rosetta, aa: str):
    vector = rosetta.utility.vector1_bool(20)
    vector[int(rosetta.core.chemical.aa_from_oneletter_code(aa))] = True
    return vector


def three_letter(rosetta, aa: str) -> str:
    return str(rosetta.core.chemical.name_from_aa(
        rosetta.core.chemical.aa_from_oneletter_code(aa)
    ))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--parent-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--design-positions", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    pyrosetta.init("-mute all -ex1 -ex2aro -use_input_sc -constant_seed")
    from pyrosetta import rosetta

    reference = read_fasta_record(args.reference_fasta, args.reference_name)
    design_positions = sorted({
        int(item)
        for item in args.design_positions.read_text().replace("\n", ",").split(",")
        if item.strip()
    })
    base = pyrosetta.pose_from_pdb(str(args.target_pdb))
    length = len(reference)
    parent = read_first_fasta(args.parent_fasta)
    if len(parent) != length:
        raise ValueError("parent FASTA length does not match reference")
    # target-pdb must already be a noncovalent Rosetta output.  Restore its
    # non-lock interface identities to the covalent C4 parent while keeping the
    # six former lock positions as Ala, so no Cys connection is ever recreated.
    for chain in range(4):
        for position in design_positions:
            aa = "A" if position in POSITIONS else parent[position - 1]
            rosetta.protocols.simple_moves.MutateResidue(
                chain * length + position, three_letter(rosetta, aa)
            ).apply(base)
    scorefxn = pyrosetta.get_fa_scorefxn()
    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"
    pdb_dir.mkdir(exist_ok=True)
    rows: list[dict[str, object]] = []

    for index, (name, profile) in enumerate(PROFILES.items(), 1):
        work = base.clone()
        rosetta.numeric.random.rg().set_seed(3000 + index)
        target_sequence = list(parent)
        for position, aa in zip(POSITIONS, profile):
            target_sequence[position - 1] = aa
        for chain in range(4):
            for position, aa in zip(POSITIONS, profile):
                if aa != "A":
                    rosetta.protocols.simple_moves.MutateResidue(
                        chain * length + position, three_letter(rosetta, aa)
                    ).apply(work)

        # Repack all previously identified interface positions, but keep the parent
        # identity outside the six forced non-Cys replacements.
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(work)
        task.or_include_current(True)
        for pose_index in range(1, work.total_residue() + 1):
            position = (pose_index - 1) % length + 1
            if position in design_positions:
                task.nonconst_residue_task(pose_index).restrict_absent_canonical_aas(
                    one_aa_vector(rosetta, target_sequence[position - 1])
                )
            else:
                task.nonconst_residue_task(pose_index).prevent_repacking()
        rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task).apply(work)

        sequence = work.chain_sequence(1)
        if any(work.chain_sequence(chain) != sequence for chain in range(2, 5)):
            raise RuntimeError("chain sequences differ")
        total_score = float(scorefxn(work))
        split = work.split_by_chain()
        separated_score = sum(float(scorefxn(split[chain])) for chain in range(1, len(split) + 1))
        interface_energy = total_score - separated_score
        candidate = f"lock_profile_{index:02d}_{name}"
        pdb_path = pdb_dir / f"{candidate}.pdb"
        work.dump_pdb(str(pdb_path))
        rows.append({
            "candidate": candidate,
            "profile": profile,
            "profile_class": name.split("_", 1)[0],
            "total_score": total_score,
            "separated_score": separated_score,
            "interface_energy": interface_energy,
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutation_string(reference, sequence),
            "native_cys_positions": ",".join(str(i) for i, aa in enumerate(sequence, 1) if aa == "C"),
            "pdb": str(pdb_path),
            "sequence": sequence,
        })
        print(candidate, f"dG={interface_energy:.2f}", flush=True)

    ranked = sorted(rows, key=lambda row: float(row["interface_energy"]))
    with (args.out / "lock_profiles_ranked.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranked[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(ranked)
    with (args.out / "lock_profiles_ranked.fasta").open("w") as handle:
        for rank, row in enumerate(ranked, 1):
            handle.write(
                f">lock_profile_rank{rank:02d} source={row['candidate']} "
                f"class={row['profile_class']} dG={row['interface_energy']:.3f}\n"
                f"{row['sequence']}\n"
            )


if __name__ == "__main__":
    main()
