#!/usr/bin/env python3
"""Sparse tied-homotetramer AF2 design on an AF3 C4 backbone.

The four copies share one 386-aa sequence.  The target AF3 coordinates define
three differentiable restraints: intrachain fold contacts, contacts involving
the validated first-interface surface, and contacts involving the designed
second-interface surface.  A differentiable mutation budget is measured
against P3-13R; fixed positions are always reset to P3-13R and new Cys is
forbidden.  AF2 scores are candidate-generation signals only--all retained
sequences require template-free AF3 validation.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
from typing import Iterable

import jax
import jax.numpy as jnp
import numpy as np

from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict


AA_ORDER = residue_constants.restype_order
FIRST_INTERFACE = {
    93, 94, 96, 105, 110, 111, 112, 113, 114, 115, 116, 119, 123,
    124, 127, 131, 132, 144, 148, 149, 151, 152, 153, 155, 156, 192,
}
SECOND_INTERFACE = {
    88, 89, 182, 191, 278, 285, 337, 341, 342, 343, 345, 346, 347,
    348, 349, 350,
}
NATIVE_CYS = {12, 31, 74, 121, 368, 383}
PROTECTED = NATIVE_CYS | set(range(257, 265))


def fasta_sequence(path: pathlib.Path, prefix: str | None = None) -> str:
    header = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None and (prefix is None or header.startswith(prefix)):
                return "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line)
    if header is not None and (prefix is None or header.startswith(prefix)):
        return "".join(chunks)
    raise ValueError(f"FASTA record {prefix!r} not found in {path}")


def read_positions(path: pathlib.Path) -> list[int]:
    tokens = path.read_text().replace("\n", ",").split(",")
    positions = sorted({int(token) for token in tokens if token.strip()})
    if not positions:
        raise ValueError("empty design-position list")
    return positions


def compress_ranges(positions: Iterable[int]) -> str:
    values = sorted(set(positions))
    out: list[str] = []
    start = end = values[0]
    for value in values[1:]:
        if value == end + 1:
            end = value
        else:
            out.append(str(start) if start == end else f"{start}-{end}")
            start = end = value
    out.append(str(start) if start == end else f"{start}-{end}")
    return ",".join(out)


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{idx}{new}"
        for idx, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-pdb", type=pathlib.Path, required=True)
    ap.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--reference-name", default="D0-P3-13R")
    ap.add_argument("--current-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--design-positions", type=pathlib.Path, required=True)
    ap.add_argument("--params", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--chains", default="A,B,C,D")
    ap.add_argument("--init", choices=["current", "reference", "gumbel"], default="current")
    ap.add_argument("--mutation-budget", type=float, default=20.0)
    ap.add_argument("--w-budget", type=float, default=8.0)
    ap.add_argument("--w-sparse", type=float, default=1.0)
    ap.add_argument("--w-fold", type=float, default=0.75)
    ap.add_argument("--w-first", type=float, default=2.5)
    ap.add_argument("--w-second", type=float, default=2.5)
    ap.add_argument("--free-interchain", action="store_true")
    ap.add_argument("--logits-iters", type=int, default=10)
    ap.add_argument("--soft-iters", type=int, default=6)
    ap.add_argument("--hard-iters", type=int, default=3)
    ap.add_argument("--recycles", type=int, default=0)
    args = ap.parse_args()

    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    current = fasta_sequence(args.current_fasta)
    design_positions = read_positions(args.design_positions)
    if len(reference) != 386 or len(current) != 386:
        raise ValueError(f"expected 386-aa sequences, got {len(reference)} and {len(current)}")
    if min(design_positions) < 1 or max(design_positions) > len(reference):
        raise ValueError("design position outside sequence")
    if set(design_positions) & PROTECTED:
        raise ValueError("design mask overlaps native Cys or SIINFEKL")
    current_mutations = {i for i, (a, b) in enumerate(zip(reference, current), 1) if a != b}
    missing = current_mutations - set(design_positions)
    if missing:
        raise ValueError(f"design mask must include every current mutation: {sorted(missing)}")

    n_chains = len([x for x in args.chains.split(",") if x])
    if n_chains != 4:
        raise ValueError("this experiment requires exactly four tied chains")
    fixed_positions = [i for i in range(1, len(reference) + 1) if i not in design_positions]
    fixed_spec = compress_ranges(fixed_positions)
    ref_idx = np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32)
    design_idx = jnp.asarray(np.asarray(design_positions, dtype=np.int32) - 1)
    ref_design_idx = jnp.asarray(ref_idx[np.asarray(design_positions) - 1])

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "pdb").mkdir(exist_ok=True)
    trajectory: list[dict[str, object]] = []
    active_stage = "setup"

    # These masks are populated after prep_inputs, then captured by the loss callback.
    target_masks: dict[str, jax.Array] = {}

    def contact_probability(outputs, cutoff: float) -> jax.Array:
        logits = outputs["distogram"]["logits"]
        bins = jnp.append(0.0, jnp.linspace(2.3125, 21.6875, 63))
        return (jax.nn.softmax(logits, axis=-1) * (bins < cutoff)).sum(-1)

    def masked_contact_loss(probability: jax.Array, mask: jax.Array) -> jax.Array:
        weight = mask.astype(probability.dtype)
        return (-jnp.log(probability + 1e-8) * weight).sum() / (weight.sum() + 1e-8)

    def restrained_losses(inputs, outputs) -> dict[str, jax.Array]:
        probability = contact_probability(outputs, cutoff=14.0)
        seq_prob = inputs["seq"]["pseudo"][0, : len(reference), :20]
        p_reference = seq_prob[design_idx, ref_design_idx]
        expected_mutations = (1.0 - p_reference).sum()
        excess = jax.nn.relu(expected_mutations - args.mutation_budget)
        return {
            "fold_target_con": masked_contact_loss(probability, target_masks["fold"]),
            "first_target_con": masked_contact_loss(probability, target_masks["first"]),
            "second_target_con": masked_contact_loss(probability, target_masks["second"]),
            "mutation_budget": jnp.square(excess / max(args.mutation_budget, 1.0)),
            "mutation_sparse": expected_mutations / len(design_positions),
            "mutation_expected": expected_mutations,
        }

    def capture(model) -> None:
        nonlocal active_stage
        if not hasattr(model, "aux") or "losses" not in model.aux:
            return
        sequence = model.get_seq(get_best=False)[0]
        log = model.aux.get("log", {})
        losses = model.aux["losses"]
        trajectory.append({
            "seed": args.seed,
            "stage": active_stage,
            "iteration": len(trajectory) + 1,
            "loss": float(model.aux.get("loss", np.nan)),
            "plddt": float(np.asarray(model.aux["plddt"]).mean()),
            "ptm": float(model.aux.get("ptm", np.nan)),
            "iptm": float(model.aux.get("i_ptm", np.nan)),
            "interface_pae": float(losses.get("i_pae", np.nan)),
            "fold_target_con": float(losses.get("fold_target_con", np.nan)),
            "first_target_con": float(losses.get("first_target_con", np.nan)),
            "second_target_con": float(losses.get("second_target_con", np.nan)),
            "mutation_expected": float(losses.get("mutation_expected", np.nan)),
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutation_string(reference, sequence),
            "sequence": sequence,
            "logged_loss": float(log.get("loss", np.nan)),
        })

    af = mk_afdesign_model(
        protocol="fixbb",
        use_multimer=True,
        use_templates=True,
        data_dir=str(args.params),
        num_recycles=args.recycles,
        model_names=["model_1_multimer_v3"],
        best_metric="loss",
        loss_callback=restrained_losses,
        post_design_callback=capture,
    )
    af.prep_inputs(
        pdb_filename=str(args.target_pdb),
        chain=args.chains,
        homooligomer=True,
        rm_template=False,
        rm_template_seq=True,
        rm_template_sc=True,
        rm_template_ic=args.free_interchain,
        fix_pos=fixed_spec,
    )
    if af._len != len(reference) or af._args["copies"] != 4:
        raise ValueError(f"expected tied 4x386 input, got copies={af._args['copies']} L={af._len}")

    # Derive the three target-contact masks from the AF3 coordinates.
    total = n_chains * len(reference)
    ca = np.asarray(af._inputs["batch"]["all_atom_positions"][:, 1])
    ca_mask = np.asarray(af._inputs["batch"]["all_atom_mask"][:, 1]).astype(bool)
    d2 = np.square(ca[:, None] - ca[None, :]).sum(-1)
    chain_id = np.arange(total) // len(reference)
    residue_number = np.arange(total) % len(reference) + 1
    valid = ca_mask[:, None] & ca_mask[None, :]
    same_chain = chain_id[:, None] == chain_id[None, :]
    sequence_separation = np.abs(residue_number[:, None] - residue_number[None, :])
    interchain = ~same_chain
    first_site = np.isin(residue_number, sorted(FIRST_INTERFACE))
    second_site = np.isin(residue_number, sorted(SECOND_INTERFACE))
    fold_mask = valid & same_chain & (sequence_separation >= 6) & (d2 < 10.0**2)
    first_mask = valid & interchain & (first_site[:, None] | first_site[None, :]) & (d2 < 12.0**2)
    second_mask = valid & interchain & (second_site[:, None] | second_site[None, :]) & (d2 < 12.0**2)
    target_masks.update({
        "fold": jnp.asarray(fold_mask),
        "first": jnp.asarray(first_mask),
        "second": jnp.asarray(second_mask),
    })
    print(
        f"target directed contacts fold={int(fold_mask.sum())} "
        f"first={int(first_mask.sum())} second={int(second_mask.sum())}",
        flush=True,
    )
    if min(int(fold_mask.sum()), int(first_mask.sum()), int(second_mask.sum())) == 0:
        raise ValueError("one target restraint has no contacts")

    # Every fixed site and the mutation prior reference P3-13R, not the target-PDB sequence.
    af._wt_aatype = ref_idx
    af.opt["weights"].update({
        "dgram_cce": 0.60,
        "fape": 0.05,
        "rmsd": 0.0,
        "plddt": 0.10,
        "pae": 0.10,
        "i_pae": 0.50,
        "i_con": 0.20,
        "con": 0.0,
        "fold_target_con": args.w_fold,
        "first_target_con": args.w_first,
        "second_target_con": args.w_second,
        "mutation_budget": args.w_budget,
        "mutation_sparse": args.w_sparse,
        "mutation_expected": 0.0,
    })
    af.opt["i_con"].update({"num": 4, "cutoff": 14.0, "binary": True})
    af._opt = copy_dict(af.opt)

    if args.init == "current":
        af.restart(seed=args.seed, seq=current, rm_aa="C")
    elif args.init == "reference":
        af.restart(seed=args.seed, seq=reference, rm_aa="C")
    else:
        af.restart(seed=args.seed, mode="gumbel", rm_aa="C")
    af._wt_aatype = ref_idx

    if args.logits_iters:
        active_stage = "logits"
        af.design_logits(
            iters=args.logits_iters, e_soft=1.0, num_recycles=args.recycles,
            models=[0], num_models=1, sample_models=False, dropout=True,
            save_best=True, verbose=1,
        )
    if args.soft_iters:
        active_stage = "soft"
        af.design_soft(
            iters=args.soft_iters, temp=1.0, e_temp=0.05,
            num_recycles=args.recycles, models=[0], num_models=1,
            sample_models=False, dropout=True, save_best=True, verbose=1,
        )
    if args.hard_iters:
        active_stage = "hard"
        af.design_hard(
            iters=args.hard_iters, temp=0.05, num_recycles=args.recycles,
            models=[0], num_models=1, sample_models=False, dropout=False,
            save_best=True, verbose=1,
        )

    sequence = af.get_seq(get_best=True)[0]
    n_mutations = sum(a != b for a, b in zip(reference, sequence))
    if any(sequence[p - 1] == "C" for p in design_positions):
        raise RuntimeError("new Cys escaped the design bias")
    if any(sequence[p - 1] != reference[p - 1] for p in fixed_positions):
        raise RuntimeError("fixed P3-13R position changed")
    best = af._tmp["best"]["aux"]
    log = best["log"]
    candidate = f"c4diff_seed{args.seed:02d}_{args.init}_{'free' if args.free_interchain else 'fixed'}_b{int(args.mutation_budget):02d}"
    pdb_path = args.out / "pdb" / f"{candidate}.pdb"
    af.save_pdb(str(pdb_path), get_best=True)
    final = {
        "candidate": candidate,
        "seed": args.seed,
        "init": args.init,
        "free_interchain": int(args.free_interchain),
        "mutation_budget": args.mutation_budget,
        "loss": float(log.get("loss", np.nan)),
        "plddt": float(log.get("plddt", np.nan)),
        "ptm": float(log.get("ptm", np.nan)),
        "iptm": float(log.get("i_ptm", np.nan)),
        "interface_pae": float(log.get("i_pae", np.nan)),
        "fold_target_con": float(log.get("fold_target_con", np.nan)),
        "first_target_con": float(log.get("first_target_con", np.nan)),
        "second_target_con": float(log.get("second_target_con", np.nan)),
        "mutation_expected": float(log.get("mutation_expected", np.nan)),
        "n_mutations": n_mutations,
        "mutations": mutation_string(reference, sequence),
        "native_cys_positions": ",".join(str(i) for i, aa in enumerate(sequence, 1) if aa == "C"),
        "pdb": str(pdb_path),
        "sequence": sequence,
    }
    for path, rows in [
        (args.out / "trajectory.tsv", trajectory),
        (args.out / "af2diff_candidates.tsv", [final]),
    ]:
        if not rows:
            continue
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    (args.out / "af2diff_candidates.fasta").write_text(
        f">{candidate} mutations={n_mutations}\n{sequence}\n"
    )
    print(
        f"FINAL {candidate} loss={final['loss']:.3f} pLDDT={final['plddt']:.3f} "
        f"ipTM={final['iptm']:.3f} mutations={n_mutations}",
        flush=True,
    )


if __name__ == "__main__":
    main()
