#!/usr/bin/env python3
"""Generate clash-aware low-mutation C4 proposals from AF3 and cross-model labels.

This is an active-learning proposal generator, not a structure validator. It fits
small Hamming-kernel ensembles to completed AF3 screens and, when informative,
to Protenix/OpenDDE/geometry/ProLIF labels. Every output still requires a new
template-free AF3 run and the frozen fail-closed screen.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from dataclasses import dataclass

import numpy as np

from run_final20_failclosed_filter import glyco_sites, reference_rsasa


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: index for index, aa in enumerate(AA_ORDER)}
NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = set(range(258, 266))


def read_tsv(path: pathlib.Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: pathlib.Path, prefix: str | None = None) -> str:
    header = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None and (prefix is None or header.startswith(prefix)):
                return "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None and (prefix is None or header.startswith(prefix)):
        return "".join(chunks)
    raise ValueError(f"FASTA record {prefix!r} not found in {path}")


def average_ranks(values: np.ndarray, higher: bool = True) -> np.ndarray:
    order = np.argsort(-values if higher else values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 4 or np.std(left) == 0 or np.std(right) == 0:
        return float("nan")
    return float(np.corrcoef(average_ranks(left), average_ranks(right))[0, 1])


def hamming(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left[:, None, :] != right[None, :, :]).sum(axis=2)


@dataclass
class KernelModel:
    gamma: float
    ridge: float
    train_x: np.ndarray
    alpha: np.ndarray

    def predict(self, test_x: np.ndarray, batch: int = 4096) -> np.ndarray:
        output = []
        for start in range(0, len(test_x), batch):
            distance = hamming(test_x[start:start + batch], self.train_x)
            output.append(np.exp(-self.gamma * distance) @ self.alpha)
        return np.concatenate(output)


def fit_model(train_x: np.ndarray, train_y: np.ndarray, gamma: float, ridge: float) -> KernelModel:
    kernel = np.exp(-gamma * hamming(train_x, train_x))
    alpha = np.linalg.solve(kernel + ridge * np.eye(len(train_x)), train_y)
    return KernelModel(gamma, ridge, train_x, alpha)


def fit_ensemble(
    train_x: np.ndarray, train_y: np.ndarray, seed: int
) -> tuple[list[KernelModel], float, list[tuple[float, float, float]]]:
    rng = np.random.default_rng(seed)
    folds = np.arange(len(train_x)) % min(5, len(train_x))
    rng.shuffle(folds)
    trials = []
    for gamma in (0.06, 0.10, 0.16, 0.24, 0.36):
        for ridge in (0.03, 0.10, 0.30, 1.00):
            predicted = np.zeros(len(train_y))
            for fold in sorted(set(folds)):
                test = folds == fold
                fit = ~test
                predicted[test] = fit_model(train_x[fit], train_y[fit], gamma, ridge).predict(train_x[test])
            trials.append((spearman(train_y, predicted), gamma, ridge))
    trials.sort(key=lambda item: (-(-math.inf if math.isnan(item[0]) else item[0]), item[1], item[2]))
    chosen = trials[:3]
    return (
        [fit_model(train_x, train_y, gamma, ridge) for _, gamma, ridge in chosen],
        chosen[0][0],
        trials,
    )


def ensemble_predict(models: list[KernelModel], test_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.stack([model.predict(test_x) for model in models])
    return values.mean(axis=0), values.std(axis=0)


def mutations(reference: str, sequence: str) -> list[str]:
    return [f"{old}{index}{new}" for index, (old, new) in enumerate(zip(reference, sequence), 1) if old != new]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--reference-name", default="D0-P3-13R")
    parser.add_argument("--reference-cif", type=pathlib.Path, required=True)
    parser.add_argument("--ranking", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--final", type=pathlib.Path, required=True)
    parser.add_argument("--anchor", action="append", required=True)
    parser.add_argument("--geometry-summary", type=pathlib.Path, required=True)
    parser.add_argument("--protenix-summary", type=pathlib.Path, required=True)
    parser.add_argument("--crossmodel-prolif", type=pathlib.Path, required=True)
    parser.add_argument("--opendde-summary", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--proposals", type=int, default=120000)
    parser.add_argument("--top", type=int, default=128)
    parser.add_argument("--min-mutations", type=int, default=15)
    parser.add_argument("--max-mutations", type=int, default=24)
    parser.add_argument("--seed", type=int, default=20260831)
    args = parser.parse_args()

    reference = read_fasta(args.reference_fasta, args.reference_name)
    if len(reference) != 386:
        raise ValueError(f"expected 386-aa reference, found {len(reference)}")
    final_rows = read_tsv(args.final)
    final = {row["candidate"]: row for row in final_rows}
    missing = set(args.anchor) - set(final)
    if missing:
        raise ValueError(f"unknown anchors: {sorted(missing)}")
    anchors = [(name, final[name]["sequence"]) for name in args.anchor]

    labelled: dict[str, float] = {}
    for path in args.ranking:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            score = row.get("mean_iptm") or row.get("median_iptm")
            if len(sequence) == 386 and score:
                labelled[sequence] = max(float(score), labelled.get(sequence, -math.inf))
    if len(labelled) < 50:
        raise ValueError(f"expected at least 50 completed AF3 labels, found {len(labelled)}")

    geometry = {row["candidate"]: row for row in read_tsv(args.geometry_summary)}
    protenix = {row["candidate"]: row for row in read_tsv(args.protenix_summary)}
    prolif = {row["candidate"]: row for row in read_tsv(args.crossmodel_prolif)}
    opendde = {row["candidate"]: row for row in read_tsv(args.opendde_summary)}
    outcomes: dict[str, dict[str, float]] = {
        "af3": labelled,
        "protenix": {}, "clean_fraction": {}, "crossmodel_prolif": {}, "opendde": {},
    }
    for name, row in final.items():
        sequence = row["sequence"]
        if name in protenix:
            outcomes["protenix"][sequence] = float(protenix[name]["mean_iptm"])
        if name in geometry:
            outcomes["clean_fraction"][sequence] = (
                float(geometry[name]["zero_atomic_clash_models"]) / float(geometry[name]["n_models"])
            )
        if name in prolif:
            outcomes["crossmodel_prolif"][sequence] = float(
                prolif[name]["shared_nonvdw_position_type_count"]
            )
        if name in opendde and opendde[name].get("gate_opendde_10_seeds") == "1":
            outcomes["opendde"][sequence] = float(opendde[name]["mean_iptm"])

    top_af3 = sorted(labelled.items(), key=lambda item: item[1], reverse=True)[:32]
    positions = sorted({
        index for sequence in list(labelled) + [seq for _, seq in anchors]
        for index, (old, new) in enumerate(zip(reference, sequence)) if old != new
    })
    # Interface blocks preserve observed multi-residue epistasis during crossover.
    blocks = [
        [pos for pos in positions if pos + 1 <= 170],
        [pos for pos in positions if 170 < pos + 1 < 270],
        [pos for pos in positions if pos + 1 >= 270],
    ]
    weighted_alphabet = {}
    for pos in positions:
        weights = {reference[pos]: 0.10}
        for sequence, score in top_af3:
            weights[sequence[pos]] = weights.get(sequence[pos], 0.0) + math.exp(9 * (score - 0.72))
        for _, sequence in anchors:
            weights[sequence[pos]] = weights.get(sequence[pos], 0.0) + 6.0
        total = sum(weights.values())
        weighted_alphabet[pos] = (list(weights), np.asarray(list(weights.values())) / total)

    rsasa = reference_rsasa(args.reference_cif, "A", reference)
    rng = np.random.default_rng(args.seed)
    excluded = set(labelled)
    proposal_set: set[str] = set()
    donor_pool = [sequence for _, sequence in anchors] + [sequence for sequence, _ in top_af3]

    def repair(sequence: list[str]) -> str:
        """Project a crossover back into the frozen sequence/surface constraints."""
        changed = [pos for pos in positions if sequence[pos] != reference[pos]]
        # Any buried Pro/Gly change is disallowed outright by the v1.0 sequence gate.
        for pos in list(changed):
            if rsasa.get(pos + 1, 0.0) < 0.20 and (
                reference[pos] in "PG" or sequence[pos] in "PG"
            ):
                sequence[pos] = reference[pos]
        changed = [pos for pos in positions if sequence[pos] != reference[pos]]
        core = [pos for pos in changed if rsasa.get(pos + 1, 0.0) < 0.20]
        surface = [pos for pos in changed if rsasa.get(pos + 1, 0.0) >= 0.20]
        rng.shuffle(core)
        while core and len(surface) / max(len(surface) + len(core), 1) < 0.80:
            pos = core.pop()
            sequence[pos] = reference[pos]
            changed = [pos for pos in positions if sequence[pos] != reference[pos]]
            core = [pos for pos in changed if rsasa.get(pos + 1, 0.0) < 0.20]
            surface = [pos for pos in changed if rsasa.get(pos + 1, 0.0) >= 0.20]
        changed = [pos for pos in positions if sequence[pos] != reference[pos]]
        while len(changed) > args.max_mutations:
            candidates = sorted(changed, key=lambda pos: (rsasa.get(pos + 1, 0.0) >= 0.20, rng.random()))
            sequence[candidates[0]] = reference[candidates[0]]
            changed = [pos for pos in positions if sequence[pos] != reference[pos]]
        surface_positions = [pos for pos in positions if rsasa.get(pos + 1, 0.0) >= 0.20]
        attempts_local = 0
        while len(changed) < args.min_mutations and attempts_local < 200:
            attempts_local += 1
            pos = int(rng.choice(surface_positions))
            donor_aa = donor_pool[int(rng.integers(len(donor_pool)))][pos]
            if donor_aa == reference[pos]:
                aas, probabilities = weighted_alphabet[pos]
                donor_aa = str(rng.choice(aas, p=probabilities))
            if donor_aa != "C":
                sequence[pos] = donor_aa
            changed = [item for item in positions if sequence[item] != reference[item]]
        return "".join(sequence)

    def valid(sequence: str) -> bool:
        changed = [index + 1 for index, (old, new) in enumerate(zip(reference, sequence)) if old != new]
        if not args.min_mutations <= len(changed) <= args.max_mutations:
            return False
        if sequence in excluded or {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            return False
        if any(sequence[pos - 1] != reference[pos - 1] for pos in EPITOPE):
            return False
        if glyco_sites(sequence) != glyco_sites(reference):
            return False
        surface_fraction = sum(rsasa.get(pos, 0.0) >= 0.20 for pos in changed) / len(changed)
        if surface_fraction < 0.80:
            return False
        if any(rsasa.get(pos, 0.0) < 0.20 and (reference[pos - 1] in "PG" or sequence[pos - 1] in "PG")
               for pos in changed):
            return False
        return True

    # Exhaustive anchor block recombination guarantees obvious cross-family hybrids
    # are represented before stochastic active-learning proposals are added.
    for base_name, base in anchors:
        del base_name
        for donor_name, donor in anchors:
            del donor_name
            for mask in range(1, 1 << len(blocks)):
                sequence = list(base)
                for block_index, block in enumerate(blocks):
                    if mask & (1 << block_index):
                        for pos in block:
                            sequence[pos] = donor[pos]
                value = repair(sequence)
                if valid(value):
                    proposal_set.add(value)

    attempts = 0
    maximum_attempts = max(1000000, args.proposals * 100)
    while len(proposal_set) < args.proposals and attempts < maximum_attempts:
        attempts += 1
        sequence = list(anchors[int(rng.integers(len(anchors)))][1])
        mode = rng.random()
        if mode < 0.55:
            donor = donor_pool[int(rng.integers(len(donor_pool)))]
            chosen_blocks = rng.choice(len(blocks), size=int(rng.integers(1, 3)), replace=False)
            for block_index in chosen_blocks:
                block = blocks[int(block_index)]
                for pos in block:
                    sequence[pos] = donor[pos]
            for pos in rng.choice(positions, size=int(rng.integers(0, 4)), replace=False):
                sequence[int(pos)] = donor_pool[int(rng.integers(len(donor_pool)))][int(pos)]
        elif mode < 0.85:
            donor = donor_pool[int(rng.integers(len(donor_pool)))]
            for pos in rng.choice(positions, size=int(rng.integers(2, 9)), replace=False):
                sequence[int(pos)] = donor[int(pos)]
        else:
            for pos in rng.choice(positions, size=int(rng.integers(1, 7)), replace=False):
                aas, probabilities = weighted_alphabet[int(pos)]
                sequence[int(pos)] = str(rng.choice(aas, p=probabilities))
        value = repair(sequence)
        if valid(value):
            proposal_set.add(value)
    if len(proposal_set) < args.top:
        raise RuntimeError(f"proposal space exhausted: {len(proposal_set)} valid sequences")

    proposal_sequences = sorted(proposal_set)
    feature_positions = positions
    proposal_x = np.asarray([[AA_INDEX[sequence[pos]] for pos in feature_positions] for sequence in proposal_sequences], dtype=np.uint8)
    predictions: dict[str, tuple[np.ndarray, np.ndarray, float]] = {}
    cv_rows = []
    for offset, (name, labels) in enumerate(outcomes.items()):
        if len(labels) < 15:
            continue
        train_sequences = list(labels)
        train_x = np.asarray([[AA_INDEX[sequence[pos]] for pos in feature_positions] for sequence in train_sequences], dtype=np.uint8)
        train_y = np.asarray([labels[sequence] for sequence in train_sequences], dtype=float)
        models, cv, trials = fit_ensemble(train_x, train_y, args.seed + offset)
        mean, sd = ensemble_predict(models, proposal_x)
        predictions[name] = (mean, sd, cv)
        for rank, (correlation, gamma, ridge) in enumerate(trials, 1):
            cv_rows.append({"outcome": name, "rank": rank, "cv_spearman": correlation,
                            "gamma": gamma, "ridge": ridge, "n_labels": len(labels)})

    weights = {"af3": 0.55, "protenix": 0.20, "clean_fraction": 0.15, "crossmodel_prolif": 0.10}
    if "opendde" in predictions:
        weights["af3"] -= 0.10
        weights["opendde"] = 0.10
    usable = {
        name: weight for name, weight in weights.items()
        if name in predictions and (name == "af3" or predictions[name][2] >= 0.20)
    }
    if "af3" not in usable:
        raise RuntimeError("AF3 surrogate was not fitted")
    total_weight = sum(usable.values())
    acquisition = np.zeros(len(proposal_sequences))
    for name, weight in usable.items():
        mean, sd, _ = predictions[name]
        robust = mean - (0.20 if name == "af3" else 0.10) * sd
        acquisition += weight / total_weight * average_ranks(robust, higher=True)
    nmut = np.asarray([len(mutations(reference, sequence)) for sequence in proposal_sequences])
    acquisition = 0.97 * acquisition + 0.03 * average_ranks(nmut, higher=False)
    order = np.argsort(acquisition)

    selected = []
    for minimum_distance in (4, 3, 2, 1):
        for index in order:
            if int(index) in selected:
                continue
            sequence = proposal_sequences[int(index)]
            if all(sum(a != b for a, b in zip(sequence, proposal_sequences[old])) >= minimum_distance
                   for old in selected):
                selected.append(int(index))
            if len(selected) == args.top:
                break
        if len(selected) == args.top:
            break

    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for rank, index in enumerate(selected, 1):
        sequence = proposal_sequences[index]
        changed = mutations(reference, sequence)
        row: dict[str, object] = {
            "rank": rank, "name": f"OVA-C4-XM-{rank:03d}",
            "acquisition_average_rank": acquisition[index],
            "n_mutations": len(changed), "mutations": ",".join(changed),
            "surface_mutation_fraction": sum(rsasa.get(int(item[1:-1]), 0.0) >= 0.20 for item in changed) / len(changed),
            "nearest_anchor_distance": min(sum(a != b for a, b in zip(sequence, anchor)) for _, anchor in anchors),
            "sequence": sequence,
        }
        for outcome, (mean, sd, cv) in predictions.items():
            row[f"predicted_{outcome}"] = mean[index]
            row[f"prediction_sd_{outcome}"] = sd[index]
            row[f"cv_spearman_{outcome}"] = cv
            row[f"used_in_acquisition_{outcome}"] = int(outcome in usable)
        rows.append(row)
    fields = []
    for row in rows:
        fields.extend(field for field in row if field not in fields)
    with (args.out / "crossmodel_surrogate_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "crossmodel_surrogate_candidates.fasta").open("w") as handle:
        for row in rows:
            handle.write(f">{row['name']} nmut={row['n_mutations']} acquisition_rank={row['acquisition_average_rank']:.3f}\n{row['sequence']}\n")
    with (args.out / "surrogate_cv.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(cv_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(cv_rows)
    with (args.out / "model_usage.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["outcome", "weight_before_normalization", "cv_spearman", "used"])
        for name in weights:
            writer.writerow([name, weights[name], predictions.get(name, (None, None, ""))[2], int(name in usable)])
    print(
        f"AF3 labels={len(labelled)} proposals={len(proposal_sequences)} selected={len(rows)}; "
        f"used outcomes={','.join(usable)}; AF3 CV={predictions['af3'][2]:.3f}"
    )


if __name__ == "__main__":
    main()
