#!/usr/bin/env python3
"""Differentiable AF2 fixed-backbone design of a tied OVA homodimer interface.

The target geometry is one adjacent interface extracted from an AF3 C4 model.  The
target side-chain identities are hidden, only selected surface positions are
mutable, and Cys is forbidden at every mutable position.  Native Cys outside the
design mask remain fixed to the P3-13R sequence.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
from typing import Iterable

import numpy as np
import jax
import jax.numpy as jnp

from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict


AA_ORDER = residue_constants.restype_order


def read_fasta_record(path: pathlib.Path, name_prefix: str) -> str:
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None and header.startswith(name_prefix):
                return "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line)
    if header is not None and header.startswith(name_prefix):
        return "".join(chunks)
    raise ValueError(f"FASTA record beginning with {name_prefix!r} not found in {path}")


def read_positions(path: pathlib.Path) -> list[int]:
    text = path.read_text().replace("\n", ",")
    positions = sorted({int(token) for token in text.split(",") if token.strip()})
    if not positions:
        raise ValueError("empty design-position file")
    return positions


def compress_ranges(positions: Iterable[int]) -> str:
    values = sorted(set(positions))
    ranges: list[str] = []
    start = end = values[0]
    for value in values[1:]:
        if value == end + 1:
            end = value
        else:
            ranges.append(str(start) if start == end else f"{start}-{end}")
            start = end = value
    ranges.append(str(start) if start == end else f"{start}-{end}")
    return ",".join(ranges)


def mutations(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=pathlib.Path, required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--design-positions", type=pathlib.Path, required=True)
    parser.add_argument("--params", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--chains", default="A,B")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--logits-iters", type=int, default=10)
    parser.add_argument("--soft-iters", type=int, default=5)
    parser.add_argument("--hard-iters", type=int, default=3)
    parser.add_argument("--recycles", type=int, default=0)
    parser.add_argument("--model", default="model_1_multimer_v3")
    parser.add_argument("--init", choices=["reference", "gumbel"], default="gumbel")
    parser.add_argument("--target-contact-cutoff", type=float, default=12.0)
    parser.add_argument("--predicted-contact-cutoff", type=float, default=14.0)
    parser.add_argument("--w-target-interface", type=float, default=5.0)
    parser.add_argument("--free-interface", action="store_true",
                        help="Keep intrachain templates but remove all template interchain geometry")
    args = parser.parse_args()

    reference = read_fasta_record(args.reference_fasta, args.reference_name)
    design_positions = read_positions(args.design_positions)
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa P3-13R, got {len(reference)}")
    if min(design_positions) < 1 or max(design_positions) > len(reference):
        raise ValueError("design position outside reference sequence")
    if set(design_positions) & {12, 31, 74, 121, 368, 383}:
        raise ValueError("native Cys must not be mutable")

    fixed_positions = [i for i in range(1, len(reference) + 1) if i not in design_positions]
    fixed_spec = compress_ranges(fixed_positions)
    seeds = [int(item) for item in args.seeds.split(",") if item.strip()]
    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"
    pdb_dir.mkdir(exist_ok=True)

    trajectory: list[dict[str, object]] = []
    active_seed = -1
    active_stage = "setup"

    design_mask_one = np.zeros(len(reference), dtype=np.float32)
    design_mask_one[np.asarray(design_positions) - 1] = 1.0
    n_chains = len([chain for chain in args.chains.split(",") if chain])
    design_mask_all = jnp.asarray(np.tile(design_mask_one, n_chains))
    design_indices_all = jnp.asarray(
        np.flatnonzero(np.tile(design_mask_one.astype(bool), n_chains))
    )

    def target_interface_contact_loss(inputs, outputs) -> dict[str, jax.Array]:
        """Reward only contacts present at the selected target A/B interface."""
        length = len(reference)
        coords = inputs["batch"]["all_atom_positions"][:, 1]
        coord_mask = inputs["batch"]["all_atom_mask"][:, 1]
        distance_sq = jnp.square(coords[:, None] - coords[None, :]).sum(-1)
        chain_id = jnp.arange(n_chains * length) // length
        interchain = chain_id[:, None] != chain_id[None, :]
        selected = (design_mask_all[:, None] + design_mask_all[None, :]) > 0
        valid = coord_mask[:, None] * coord_mask[None, :]
        target_contacts = (
            interchain
            & selected
            & (distance_sq < args.target_contact_cutoff**2)
            & (valid > 0)
        )
        bins = jnp.append(0.0, jnp.linspace(2.3125, 21.6875, 63))
        probability = (
            jax.nn.softmax(outputs["distogram"]["logits"], axis=-1)
            * (bins < args.predicted_contact_cutoff)
        ).sum(-1)
        mask = target_contacts.astype(probability.dtype)
        loss = (
            -jnp.log(probability + 1e-8) * mask
        ).sum() / (mask.sum() + 1e-8)
        return {"target_i_con": loss}

    def surface_interface_contact_loss(inputs, outputs) -> dict[str, jax.Array]:
        """Reward a best subset of mutable surface sites contacting the other chain."""
        length = len(reference)
        bins = jnp.append(0.0, jnp.linspace(2.3125, 21.6875, 63))
        probability = (
            jax.nn.softmax(outputs["distogram"]["logits"], axis=-1)
            * (bins < args.predicted_contact_cutoff)
        ).sum(-1)
        chain_id = jnp.arange(n_chains * length) // length
        interchain = chain_id[:, None] != chain_id[None, :]
        best_partner = jnp.where(interchain, probability, 0.0).max(-1)
        selected_probability = best_partner[design_indices_all]
        keep = max(1, int(round(len(design_positions) * n_chains * 0.40)))
        strongest = jnp.sort(selected_probability)[-keep:]
        return {"surface_i_con": -jnp.log(strongest + 1e-8).mean()}

    def capture(model) -> None:
        if not hasattr(model, "aux") or "losses" not in model.aux:
            return
        sequence = model.get_seq(get_best=False)[0]
        losses = model.aux["losses"]
        trajectory.append({
            "seed": active_seed,
            "stage": active_stage,
            "iteration": len(model._tmp.get("log", [])),
            "loss": float(model.aux.get("loss", np.nan)),
            "plddt": float(np.asarray(model.aux["plddt"]).mean()),
            "ptm": float(model.aux.get("ptm", np.nan)),
            "iptm": float(model.aux.get("i_ptm", np.nan)),
            "dgram_cce": float(losses.get("dgram_cce", np.nan)),
            "rmsd": float(losses.get("rmsd", np.nan)),
            "pae": float(losses.get("pae", np.nan)),
            "interface_pae": float(losses.get("i_pae", np.nan)),
            "interface_contact": float(losses.get("i_con", np.nan)),
            "target_interface_contact": float(losses.get("target_i_con", np.nan)),
            "surface_interface_contact": float(losses.get("surface_i_con", np.nan)),
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutations(reference, sequence),
            "sequence": sequence,
        })

    af = mk_afdesign_model(
        protocol="fixbb",
        use_multimer=True,
        use_templates=True,
        data_dir=str(args.params),
        num_recycles=args.recycles,
        model_names=[args.model],
        best_metric="loss",
        loss_callback=(surface_interface_contact_loss if args.free_interface
                       else target_interface_contact_loss),
        post_design_callback=capture,
    )
    af.prep_inputs(
        pdb_filename=str(args.target_pdb),
        chain=args.chains,
        homooligomer=True,
        rm_template=False,
        rm_template_seq=True,
        rm_template_sc=True,
        rm_template_ic=args.free_interface,
        fix_pos=fixed_spec,
    )
    if af._len != len(reference) or af._args["copies"] != n_chains:
        raise ValueError(
            f"expected tied {n_chains}x386 target, got copies={af._args['copies']} L={af._len}"
        )
    target_ca = np.asarray(af._inputs["batch"]["all_atom_positions"][:, 1])
    target_ca_mask = np.asarray(af._inputs["batch"]["all_atom_mask"][:, 1])
    target_d2 = np.square(target_ca[:, None] - target_ca[None, :]).sum(-1)
    chain_id = np.arange(n_chains * len(reference)) // len(reference)
    selected = np.tile(design_mask_one, n_chains)
    target_pair_mask = (
        (chain_id[:, None] != chain_id[None, :])
        & ((selected[:, None] + selected[None, :]) > 0)
        & (target_d2 < args.target_contact_cutoff**2)
        & ((target_ca_mask[:, None] * target_ca_mask[None, :]) > 0)
    )
    print(f"target directed interface pairs: {int(target_pair_mask.sum())}", flush=True)

    # The PDB came from an earlier covalent design.  Use only its backbone geometry;
    # fixed positions and every restart are anchored to the actual P3-13R sequence.
    af._wt_aatype = np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32)
    af.opt["weights"].update({
        "dgram_cce": 0.75,
        "fape": 0.05,
        "rmsd": 0.0,
        "plddt": 0.10,
        "pae": 0.05,
        "i_pae": 0.50,
        "i_con": 0.50 if args.free_interface else 0.0,
        "con": 0.0,
        "target_i_con": 0.0 if args.free_interface else args.w_target_interface,
        "surface_i_con": 2.0 if args.free_interface else 0.0,
    })
    af.opt["i_con"].update({"num": 4, "cutoff": 14.0, "binary": True})
    af._opt = copy_dict(af.opt)

    final_rows: list[dict[str, object]] = []
    for seed in seeds:
        active_seed = seed
        if args.init == "reference":
            af.restart(seed=seed, seq=reference, rm_aa="C")
        else:
            af.restart(seed=seed, mode="gumbel", rm_aa="C")
        # restart restores defaults, so make the P3 fixed sequence and Cys prohibition explicit.
        af._wt_aatype = np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32)

        active_stage = "logits"
        if args.logits_iters:
            af.design_logits(
                iters=args.logits_iters,
                e_soft=1.0,
                num_recycles=args.recycles,
                models=[0],
                num_models=1,
                sample_models=False,
                dropout=True,
                save_best=True,
                verbose=1,
            )
        active_stage = "soft"
        if args.soft_iters:
            af.design_soft(
                iters=args.soft_iters,
                temp=1.0,
                e_temp=0.05,
                num_recycles=args.recycles,
                models=[0],
                num_models=1,
                sample_models=False,
                dropout=True,
                save_best=True,
                verbose=1,
            )
        active_stage = "hard"
        if args.hard_iters:
            af.design_hard(
                iters=args.hard_iters,
                temp=0.05,
                num_recycles=args.recycles,
                models=[0],
                num_models=1,
                sample_models=False,
                dropout=False,
                save_best=True,
                verbose=1,
            )

        sequence = af.get_seq(get_best=True)[0]
        best = af._tmp["best"]["aux"]
        log = best["log"]
        if any(sequence[pos - 1] == "C" for pos in design_positions):
            raise RuntimeError("engineered Cys escaped the design bias")
        if any(sequence[pos - 1] != reference[pos - 1] for pos in fixed_positions):
            raise RuntimeError("fixed P3-13R position changed")
        name = f"af2diff_seed{seed:02d}"
        pdb_path = pdb_dir / f"{name}.pdb"
        af.save_pdb(str(pdb_path), get_best=True)
        final_rows.append({
            "candidate": name,
            "seed": seed,
            "loss": float(log.get("loss", np.nan)),
            "plddt": float(log.get("plddt", np.nan)),
            "ptm": float(log.get("ptm", np.nan)),
            "iptm": float(log.get("i_ptm", np.nan)),
            "dgram_cce": float(log.get("dgram_cce", np.nan)),
            "rmsd": float(log.get("rmsd", np.nan)),
            "pae": float(log.get("pae", np.nan)),
            "interface_pae": float(log.get("i_pae", np.nan)),
            "interface_contact": float(log.get("i_con", np.nan)),
            "target_interface_contact": float(log.get("target_i_con", np.nan)),
            "surface_interface_contact": float(log.get("surface_i_con", np.nan)),
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutations(reference, sequence),
            "native_cys_positions": ",".join(str(i) for i, aa in enumerate(sequence, 1) if aa == "C"),
            "pdb": str(pdb_path),
            "sequence": sequence,
        })
        print(
            f"FINAL seed={seed} loss={final_rows[-1]['loss']:.3f} "
            f"pLDDT={final_rows[-1]['plddt']:.3f} ipTM={final_rows[-1]['iptm']:.3f} "
            f"mutations={final_rows[-1]['n_mutations']}",
            flush=True,
        )

    for path, rows in [
        (args.out / "trajectory.tsv", trajectory),
        (args.out / "af2diff_candidates.tsv", final_rows),
    ]:
        if not rows:
            continue
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    with (args.out / "af2diff_candidates.fasta").open("w") as handle:
        for row in final_rows:
            handle.write(f">{row['candidate']} mutations={row['n_mutations']}\n{row['sequence']}\n")


if __name__ == "__main__":
    main()
