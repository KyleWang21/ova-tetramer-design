#!/usr/bin/env python3
"""Sparse tied-sequence AF2 design against two or more C4 interface states.

Full 4x386 AF2 backprop requires about 299 GB and does not fit an 80-GB A800.
This implementation keeps one continuous tied sequence parameter tensor and,
at every optimization step, sums gradients from the AF3-derived dimer
states before making one optimizer update.  It is a multistate approximation
to the C4 objective; template-free four-chain AF3 remains the success test.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib

import jax
import jax.numpy as jnp
import numpy as np
from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.utils import copy_dict


AA_ORDER = residue_constants.restype_order
ORDER_AA = {index: aa for aa, index in AA_ORDER.items()}
NATIVE_CYS = {12, 31, 74, 121, 368, 383}
PROTECTED = NATIVE_CYS | set(range(258, 266))


def fasta_sequence(path: pathlib.Path, prefix: str | None = None) -> str:
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if header is not None and (prefix is None or header.startswith(prefix)):
                return "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    raise ValueError(f"FASTA record {prefix!r} not found in {path}")


def read_positions(path: pathlib.Path) -> list[int]:
    values = [token.strip() for token in path.read_text().replace("\n", ",").split(",")]
    return sorted({int(value) for value in values if value})


def compress_ranges(values: list[int]) -> str:
    if not values:
        return ""
    output: list[str] = []
    start = previous = values[0]
    for value in values[1:] + [values[-1] + 2]:
        if value != previous + 1:
            output.append(str(start) if start == previous else f"{start}-{previous}")
            start = value
        previous = value
    return ",".join(output)


def mutation_string(reference: str, sequence: str) -> str:
    return ",".join(
        f"{old}{index}{new}"
        for index, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    )


def sequence_from_aux(aux: dict[str, object]) -> str:
    indices = np.asarray(aux["seq"]["hard"]).argmax(-1)[0]
    return "".join(ORDER_AA[int(index)] for index in indices)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-pdb", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--chains", action="append", required=True)
    parser.add_argument("--restraint-positions", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--state-weight", type=float, action="append", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--base-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--base-name", required=True)
    parser.add_argument("--design-positions", type=pathlib.Path, required=True)
    parser.add_argument("--surface-positions", type=pathlib.Path, required=True)
    parser.add_argument("--params", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mutation-budget", type=float, default=24.0)
    parser.add_argument("--min-surface-mutation-fraction", type=float, default=0.80)
    parser.add_argument("--w-budget", type=float, default=20.0)
    parser.add_argument("--w-sparse", type=float, default=1.5)
    parser.add_argument("--w-buried-mutation", type=float, default=12.0)
    parser.add_argument("--w-fold", type=float, default=1.0)
    parser.add_argument(
        "--worst-state-gradient-boost", type=float, default=1.0,
        help="multiply the current lower-ipTM state's sequence gradient",
    )
    parser.add_argument(
        "--worst-state-count", type=int, default=1,
        help="number of current lowest-ipTM states receiving the gradient boost",
    )
    parser.add_argument("--logits-iters", type=int, default=18)
    parser.add_argument("--soft-iters", type=int, default=10)
    parser.add_argument("--hard-iters", type=int, default=6)
    parser.add_argument("--recycles", type=int, default=0)
    args = parser.parse_args()

    state_count = len(args.target_pdb)
    if not (
        state_count >= 2
        and state_count == len(args.chains)
        and state_count == len(args.restraint_positions)
        and state_count == len(args.state_weight)
    ):
        raise ValueError("at least two target/interface states with matching weights are required")
    if not 0.0 <= args.min_surface_mutation_fraction <= 1.0:
        raise ValueError("--min-surface-mutation-fraction must lie in [0, 1]")
    if args.worst_state_gradient_boost < 1.0:
        raise ValueError("--worst-state-gradient-boost must be at least 1.0")
    if not 1 <= args.worst_state_count <= state_count:
        raise ValueError("--worst-state-count must be in [1, number of states]")
    reference = fasta_sequence(args.reference_fasta, args.reference_name)
    base = fasta_sequence(args.base_fasta, args.base_name)
    if len(reference) != 386 or len(base) != 386:
        raise ValueError(f"expected 386-aa sequences, got {len(reference)} and {len(base)}")
    design_positions = read_positions(args.design_positions)
    surface_positions = set(read_positions(args.surface_positions))
    if set(design_positions) & PROTECTED:
        raise ValueError("design mask overlaps a native Cys or SIINFEKL")
    if not design_positions or min(design_positions) < 1 or max(design_positions) > len(reference):
        raise ValueError("invalid or empty design position set")
    if any(len([value for value in chains.split(",") if value]) != 2 for chains in args.chains):
        raise ValueError("each multistate target must contain exactly two tied chains")
    fixed_positions = [pos for pos in range(1, len(reference) + 1) if pos not in design_positions]
    fixed_spec = compress_ranges(fixed_positions)
    ref_index = jnp.asarray(np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32))
    base_index = np.asarray([AA_ORDER[aa] for aa in base], dtype=np.int32)
    buried = sorted(set(design_positions) - surface_positions)
    buried_index = jnp.asarray(np.asarray(buried, dtype=np.int32) - 1)

    def contact_probability(outputs: dict[str, object], cutoff: float) -> jax.Array:
        logits = outputs["distogram"]["logits"]
        bins = jnp.append(0.0, jnp.linspace(2.3125, 21.6875, 63))
        return (jax.nn.softmax(logits, axis=-1) * (bins < cutoff)).sum(-1)

    def masked_contact_loss(probability: jax.Array, mask: jax.Array) -> jax.Array:
        weight = mask.astype(probability.dtype)
        return (-jnp.log(probability + 1e-8) * weight).sum() / (weight.sum() + 1e-8)

    def restrained_losses(inputs, outputs) -> dict[str, jax.Array]:  # noqa: ANN001
        probability = contact_probability(outputs, 14.0)
        seq_probability = inputs["seq"]["pseudo"][0, : len(reference), :20]
        p_reference = seq_probability[jnp.arange(len(reference)), ref_index]
        expected_mutations = (1.0 - p_reference).sum()
        expected_buried = (
            (1.0 - p_reference)[buried_index].sum()
            if buried else jnp.asarray(0.0, dtype=p_reference.dtype)
        )
        excess = jax.nn.relu(expected_mutations - args.mutation_budget)
        return {
            "fold_target_con": masked_contact_loss(probability, inputs["ova_fold_target_mask"]),
            "interface_target_con": masked_contact_loss(probability, inputs["ova_interface_target_mask"]),
            "mutation_budget": jnp.square(excess / max(args.mutation_budget, 1.0)),
            "mutation_sparse": expected_mutations / len(reference),
            "mutation_buried": expected_buried / max(len(buried), 1),
            "mutation_expected": expected_mutations,
        }

    af = mk_afdesign_model(
        protocol="fixbb", use_multimer=True, use_templates=True,
        data_dir=str(args.params), num_recycles=args.recycles,
        model_names=["model_1_multimer_v3"], best_metric="loss",
        loss_callback=restrained_losses,
    )
    state_inputs: list[dict[str, object]] = []
    mask_counts: list[tuple[int, int]] = []
    for target, chains, restraint_path in zip(args.target_pdb, args.chains, args.restraint_positions):
        af.prep_inputs(
            pdb_filename=str(target), chain=chains, homooligomer=True,
            rm_template=False, rm_template_seq=True, rm_template_sc=True,
            rm_template_ic=False, fix_pos=fixed_spec,
        )
        if af._len != len(reference) or af._args["copies"] != 2:
            raise ValueError(f"expected tied 2x386 state, got copies={af._args['copies']} L={af._len}")
        total = 2 * len(reference)
        ca = np.asarray(af._inputs["batch"]["all_atom_positions"][:, 1])
        ca_mask = np.asarray(af._inputs["batch"]["all_atom_mask"][:, 1]).astype(bool)
        d2 = np.square(ca[:, None] - ca[None, :]).sum(-1)
        chain_id = np.arange(total) // len(reference)
        residue_number = np.arange(total) % len(reference) + 1
        valid = ca_mask[:, None] & ca_mask[None, :]
        same_chain = chain_id[:, None] == chain_id[None, :]
        separation = np.abs(residue_number[:, None] - residue_number[None, :])
        restrained = np.isin(residue_number, read_positions(restraint_path))
        fold_mask = valid & same_chain & (separation >= 6) & (d2 < 10.0**2)
        interface_mask = (
            valid & ~same_chain & (restrained[:, None] | restrained[None, :]) & (d2 < 12.0**2)
        )
        if min(int(fold_mask.sum()), int(interface_mask.sum())) == 0:
            raise ValueError(f"target state {target} has an empty contact mask")
        af._inputs["ova_fold_target_mask"] = jnp.asarray(fold_mask)
        af._inputs["ova_interface_target_mask"] = jnp.asarray(interface_mask)
        state_inputs.append(copy_dict(af._inputs))
        mask_counts.append((int(fold_mask.sum()), int(interface_mask.sum())))

    # Each forward contributes 1/N of shared fold/sequence/standard losses;
    # interface terms retain the requested relative weights.
    af._wt_aatype = base_index
    af.opt["weights"].update({
        "dgram_cce": 0.30, "fape": 0.025, "rmsd": 0.0,
        "plddt": 0.05, "pae": 0.05, "i_pae": 0.25,
        "i_con": 0.10, "con": 0.0,
        "fold_target_con": args.w_fold / state_count,
        "interface_target_con": 1.0,
        "mutation_budget": args.w_budget / state_count,
        "mutation_sparse": args.w_sparse / state_count,
        "mutation_buried": args.w_buried_mutation / state_count,
        "mutation_expected": 0.0,
    })
    af.opt["i_con"].update({"num": 4, "cutoff": 14.0, "binary": True})
    af._opt = copy_dict(af.opt)
    af._inputs = state_inputs[0]
    af.restart(seed=args.seed, seq=base, rm_aa="C")
    shared_bias = copy_dict(af._inputs["bias"])
    for state in state_inputs:
        state["bias"] = shared_bias

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "pdb").mkdir(exist_ok=True)
    trajectory: list[dict[str, object]] = []
    best: dict[str, object] | None = None

    def optimize_stage(stage: str, iterations: int) -> None:
        nonlocal best
        for iteration in range(iterations):
            fraction = (iteration + 1) / max(iterations, 1)
            if stage == "logits":
                soft = fraction
                hard = 0.0
                temperature = 1.0
                dropout = True
            elif stage == "soft":
                soft = 1.0
                hard = 0.0
                temperature = 0.05 + 0.95 * (1.0 - fraction) ** 2
                dropout = True
            else:
                soft = 1.0
                hard = 1.0
                temperature = 0.05
                dropout = False
            af.set_opt(soft=soft, hard=hard, temp=temperature, dropout=dropout)
            state_results: list[dict[str, object]] = []
            gradients = []
            for state_index, state in enumerate(state_inputs):
                af._inputs = state
                af.opt["weights"]["interface_target_con"] = args.state_weight[state_index]
                aux = af.run(
                    num_recycles=args.recycles, models=[0], num_models=1,
                    sample_models=False, backprop=True, return_aux=True,
                )
                gradients.append(copy_dict(aux["grad"]))
                losses = aux["losses"]
                state_results.append({
                    "loss": float(aux["loss"]),
                    "iptm": float(aux.get("i_ptm", math.nan)),
                    "plddt": float(np.asarray(aux["plddt"]).mean()),
                    "interface_pae": float(losses.get("i_pae", math.nan)),
                    "fold_target_con": float(losses["fold_target_con"]),
                    "interface_target_con": float(losses["interface_target_con"]),
                    "mutation_expected": float(losses["mutation_expected"]),
                    "sequence": sequence_from_aux(aux),
                })
            state_iptms = [float(result["iptm"]) for result in state_results]
            weakest_state = int(np.nanargmin(state_iptms))
            gradient_scales = [1.0] * state_count
            boosted_states = np.argsort(np.nan_to_num(state_iptms, nan=-math.inf))[
                : args.worst_state_count
            ]
            for boosted_state in boosted_states:
                gradient_scales[int(boosted_state)] = args.worst_state_gradient_boost
            combined_gradient = jax.tree_util.tree_map(
                lambda *leaves: sum(
                    scale * leaf for scale, leaf in zip(gradient_scales, leaves, strict=True)
                ),
                *gradients,
            )
            af.aux = {"grad": combined_gradient}
            if af.opt["norm_seq_grad"]:
                af._norm_seq_grad()
            af._state, transformed = af._optimizer(af._state, af.aux["grad"], af._params)
            lr_scale = (1.0 - soft) + soft * temperature
            learning_rate = af.opt["learning_rate"] * lr_scale
            af._params = jax.tree_util.tree_map(
                lambda parameter, gradient: parameter - learning_rate * gradient,
                af._params, transformed,
            )
            af._k += 1
            combined_loss = sum(float(result["loss"]) for result in state_results)
            sequence = str(state_results[-1]["sequence"])
            mutation_positions = {
                position
                for position, (old, new) in enumerate(
                    zip(reference, sequence, strict=True), 1
                )
                if old != new
            }
            nmut = len(mutation_positions)
            surface_fraction = (
                len(mutation_positions & surface_positions) / nmut if nmut else 1.0
            )
            record: dict[str, object] = {
                "seed": args.seed, "stage": stage,
                "iteration": len(trajectory) + 1,
                "combined_loss": combined_loss,
                "n_mutations": nmut,
                "surface_mutation_fraction": surface_fraction,
                "mutations": mutation_string(reference, sequence),
                "sequence": sequence,
                "soft": soft, "hard": hard, "temperature": temperature,
                "weakest_state": weakest_state + 1,
                "worst_state_gradient_boost": args.worst_state_gradient_boost,
                "worst_state_count": args.worst_state_count,
            }
            for state_index, result in enumerate(state_results, 1):
                record.update({
                    f"state{state_index}_loss": result["loss"],
                    f"state{state_index}_iptm": result["iptm"],
                    f"state{state_index}_plddt": result["plddt"],
                    f"state{state_index}_interface_pae": result["interface_pae"],
                    f"state{state_index}_fold_target_con": result["fold_target_con"],
                    f"state{state_index}_interface_target_con": result["interface_target_con"],
                    f"state{state_index}_mutation_expected": result["mutation_expected"],
                })
            trajectory.append(record)
            eligible = (
                stage in {"soft", "hard"}
                and nmut <= int(args.mutation_budget)
                and surface_fraction >= args.min_surface_mutation_fraction
            )
            if eligible and (best is None or combined_loss < float(best["combined_loss"])):
                best = dict(record)
            iptm_text = "/".join(f"{float(result['iptm']):.3f}" for result in state_results)
            print(
                f"{stage} {iteration + 1}/{iterations} loss={combined_loss:.3f} "
                f"ipTM={iptm_text} "
                f"nmut={nmut}", flush=True,
            )

    optimize_stage("logits", args.logits_iters)
    optimize_stage("soft", args.soft_iters)
    optimize_stage("hard", args.hard_iters)
    if best is None:
        eligible = [
            row for row in trajectory
            if int(row["n_mutations"]) <= int(args.mutation_budget)
            and float(row["surface_mutation_fraction"])
            >= args.min_surface_mutation_fraction
        ]
        if not eligible:
            raise RuntimeError(
                "no discrete state met the mutation budget and surface-fraction gate"
            )
        best = min(eligible, key=lambda row: float(row["combined_loss"]))
    sequence = str(best["sequence"])
    if {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
        raise RuntimeError("native-only Cys constraint failed")
    if any(sequence[pos - 1] != base[pos - 1] for pos in fixed_positions):
        raise RuntimeError("a fixed base-sequence position changed")
    final_mutation_positions = {
        position
        for position, (old, new) in enumerate(zip(reference, sequence, strict=True), 1)
        if old != new
    }
    final_surface_fraction = (
        len(final_mutation_positions & surface_positions) / len(final_mutation_positions)
        if final_mutation_positions else 1.0
    )
    if final_surface_fraction < args.min_surface_mutation_fraction:
        raise RuntimeError("selected sequence failed the surface-mutation-fraction gate")

    final_metrics: list[dict[str, float]] = []
    for state_index, state in enumerate(state_inputs, 1):
        af._inputs = state
        aux = af.predict(
            seq=sequence, num_models=1, num_recycles=args.recycles,
            models=[0], sample_models=False, dropout=False,
            hard=True, soft=False, temp=0.05, return_aux=True, verbose=False,
        )
        final_metrics.append({
            "iptm": float(aux["i_ptm"]),
            "ptm": float(aux["ptm"]),
            "plddt": float(np.asarray(aux["plddt"]).mean()),
            "interface_pae": float(aux["losses"].get("i_pae", math.nan)),
            "fold_target_con": float(aux["losses"]["fold_target_con"]),
            "interface_target_con": float(aux["losses"]["interface_target_con"]),
        })
        af.save_pdb(str(args.out / "pdb" / f"selected_state{state_index}.pdb"), get_best=False)

    with (args.out / "trajectory.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(trajectory[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(trajectory)
    final: dict[str, object] = {
        "candidate": f"multistate_seed{args.seed:04d}_b{int(args.mutation_budget):02d}",
        "seed": args.seed,
        "state_weights": ",".join(map(str, args.state_weight)),
        "worst_state_gradient_boost": args.worst_state_gradient_boost,
        "worst_state_count": args.worst_state_count,
        "mask_counts": ";".join(f"{fold},{interface}" for fold, interface in mask_counts),
        "selected_stage": best["stage"],
        "selected_iteration": best["iteration"],
        "combined_loss": best["combined_loss"],
        "n_mutations": sum(old != new for old, new in zip(reference, sequence)),
        "surface_mutation_fraction": final_surface_fraction,
        "mutations": mutation_string(reference, sequence),
        "native_cys_positions": ",".join(str(i) for i, aa in enumerate(sequence, 1) if aa == "C"),
        "sequence": sequence,
    }
    for state_index, metrics in enumerate(final_metrics, 1):
        final.update({f"state{state_index}_{key}": value for key, value in metrics.items()})
    with (args.out / "af2diff_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(final), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerow(final)
    (args.out / "af2diff_candidates.fasta").write_text(
        f">{final['candidate']} mutations={final['n_mutations']}\n{sequence}\n"
    )
    final_iptm_text = "/".join(f"{metrics['iptm']:.3f}" for metrics in final_metrics)
    print(
        f"FINAL {final['candidate']} state_ipTM="
        f"{final_iptm_text} "
        f"mutations={final['n_mutations']}", flush=True,
    )


if __name__ == "__main__":
    main()
