#!/usr/bin/env python3
"""Template-free AF2/BindCraft-style design of a tied noncovalent OVA homodimer."""

from __future__ import annotations

import argparse
import csv
import pathlib

import jax
import jax.numpy as jnp
import numpy as np
from colabdesign import mk_afdesign_model
from colabdesign.af.alphafold.common import residue_constants


AA_ORDER = residue_constants.restype_order


def read_fasta(path: pathlib.Path, prefix: str) -> str:
    header = None; sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None and header.startswith(prefix): return "".join(sequence)
            header, sequence = line[1:], []
        elif header is not None:
            sequence.append(line.strip())
    if header is not None and header.startswith(prefix): return "".join(sequence)
    raise ValueError(f"missing FASTA record {prefix}")


def read_positions(path: pathlib.Path) -> list[int]:
    return sorted({int(value) for value in path.read_text().replace("\n", ",").split(",") if value.strip()})


def mutation_text(reference: str, sequence: str) -> str:
    return ",".join(f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--reference-name", default="D0-P3-13R")
    ap.add_argument("--design-positions", type=pathlib.Path, required=True)
    ap.add_argument("--params", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--logits-iters", type=int, default=12)
    ap.add_argument("--soft-iters", type=int, default=6)
    ap.add_argument("--hard-iters", type=int, default=3)
    ap.add_argument("--recycles", type=int, default=0)
    ap.add_argument("--model", default="model_1_multimer_v3")
    ap.add_argument("--coverage-fraction", type=float, default=0.40)
    args = ap.parse_args()

    reference = read_fasta(args.reference_fasta, args.reference_name)
    positions = read_positions(args.design_positions)
    if len(reference) != 386 or not positions:
        raise SystemExit("expected 386-aa reference and nonempty design surface")
    native_cys = {12, 31, 74, 121, 368, 383}
    if set(positions) & native_cys:
        raise SystemExit("native Cys cannot be mutable")
    fixed = np.asarray([index for index in range(386) if index + 1 not in positions], dtype=np.int32)
    selected_one = np.zeros(386, dtype=bool); selected_one[np.asarray(positions) - 1] = True
    selected_two = jnp.asarray(np.tile(selected_one, 2))

    def surface_contact_loss(inputs, outputs) -> dict[str, jax.Array]:  # noqa: ANN001
        bins = jnp.append(0.0, jnp.linspace(2.3125, 21.6875, 63))
        probability = (jax.nn.softmax(outputs["distogram"]["logits"], -1) * (bins < 14.0)).sum(-1)
        length = 386
        chain = jnp.arange(2 * length) // length
        interchain = chain[:, None] != chain[None, :]
        probability = jnp.where(interchain, probability, 0.0)
        best_partner = probability.max(-1)
        selected_probability = best_partner[selected_two]
        keep = max(1, int(round(len(positions) * 2 * args.coverage_fraction)))
        strongest = jnp.sort(selected_probability)[-keep:]
        return {"surface_i_con": -jnp.log(strongest + 1e-8).mean()}

    args.out.mkdir(parents=True, exist_ok=True)
    pdb_dir = args.out / "pdb"; pdb_dir.mkdir(exist_ok=True)
    model = mk_afdesign_model(
        protocol="hallucination", use_multimer=True, use_templates=False,
        data_dir=str(args.params), num_recycles=args.recycles,
        model_names=[args.model], best_metric="loss", loss_callback=surface_contact_loss,
    )
    model.prep_inputs(length=386, copies=2)
    model._wt_aatype = np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32)
    model.opt["fix_pos"] = fixed
    model.opt["weights"].update({
        "plddt": 0.20, "pae": 0.05, "i_pae": 0.50,
        "con": 0.05, "i_con": 0.50, "surface_i_con": 2.0,
    })
    model.opt["i_con"].update({"num": 4, "cutoff": 14.0, "binary": True})

    rows = []
    for seed in [int(value) for value in args.seeds.split(",") if value.strip()]:
        model.restart(seed=seed, mode="gumbel", rm_aa="C")
        model._wt_aatype = np.asarray([AA_ORDER[aa] for aa in reference], dtype=np.int32)
        model.opt["fix_pos"] = fixed
        model.opt["weights"].update({
            "plddt": 0.20, "pae": 0.05, "i_pae": 0.50,
            "con": 0.05, "i_con": 0.50, "surface_i_con": 2.0,
        })
        model.design_logits(args.logits_iters, e_soft=1.0, models=[0], num_models=1,
                            sample_models=False, dropout=True, save_best=True, verbose=1)
        model.design_soft(args.soft_iters, temp=1.0, e_temp=0.05, models=[0], num_models=1,
                          sample_models=False, dropout=True, save_best=True, verbose=1)
        model.design_hard(args.hard_iters, temp=0.05, models=[0], num_models=1,
                          sample_models=False, dropout=False, save_best=True, verbose=1)
        sequence = model.get_seq(get_best=True)[0]
        if any(sequence[index] != reference[index] for index in fixed):
            raise RuntimeError("fixed reference position changed")
        if set(i + 1 for i, aa in enumerate(sequence) if aa == "C") != native_cys:
            raise RuntimeError("Cys constraint failed")
        best = model._tmp["best"]["aux"]; log = best["log"]
        name = f"af2hall_dimer_seed{seed:02d}"
        pdb = pdb_dir / f"{name}.pdb"; model.save_pdb(str(pdb), get_best=True)
        rows.append({
            "candidate": name, "seed": seed, "loss": float(log.get("loss", np.nan)),
            "plddt": float(log.get("plddt", np.nan)), "ptm": float(log.get("ptm", np.nan)),
            "iptm": float(log.get("i_ptm", np.nan)), "interface_pae": float(log.get("i_pae", np.nan)),
            "surface_contact_loss": float(log.get("surface_i_con", np.nan)),
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": mutation_text(reference, sequence), "pdb": str(pdb), "sequence": sequence,
        })
        print(name, f"pLDDT={rows[-1]['plddt']:.3f}", f"ipTM={rows[-1]['iptm']:.3f}",
              f"mut={rows[-1]['n_mutations']}", flush=True)
    with (args.out / "af2hall_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    with (args.out / "af2hall_candidates.fasta").open("w") as handle:
        for row in rows: handle.write(f">{row['candidate']}|mut={row['mutations']}\n{row['sequence']}\n")


if __name__ == "__main__":
    main()
