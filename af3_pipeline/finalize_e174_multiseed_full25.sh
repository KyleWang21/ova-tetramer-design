#!/usr/bin/env bash
# Build one strict 25-model summary after an E174 reseed experiment completes.
set -euo pipefail

EXPERIMENT="${1:?multiseed experiment path required}"
PROJECT=/root/400083/ova_p3_tetramer_design
SEED1="$PROJECT/experiments/e174_c4_073_multihit_ai_af3_seed1"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_BIO=/vepfs-mlp2/c20250508/400083/GIPR-baker/env/mamba/envs/mpnn_fr/bin/python
cd "$PROJECT"

deadline=$((SECONDS + 172800))
while [[ ! -f "$EXPERIMENT/RESEED_DECISION_COMPLETE" ]]; do
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for $EXPERIMENT reseed decision" >&2
    exit 2
  fi
  sleep 30
done

strict_args=()
for seed in 2 3 4 5; do
  path="$EXPERIMENT/seed${seed}_strict/seed1_strict_per_model.tsv"
  if [[ ! -s "$path" ]]; then
    echo "$EXPERIMENT stopped before a full seeds2-5 set; no 25-model summary"
    touch "$EXPERIMENT/FULL25_NOT_AVAILABLE_FAILFAST"
    exit 0
  fi
  strict_args+=(--multiseed-strict-per-model "$path")
done

if [[ ! -s "$SEED1/model_scores.tsv" ]]; then
  "$PY_AF3" af3_pipeline/score_outputs.py --experiment "$SEED1"
fi
if [[ ! -s "$SEED1/disulfide_model_scores.tsv" ]]; then
  "$PY_BIO" af3_pipeline/score_disulfide_ring.py \
    --experiment "$SEED1" --no-designed-pairs
fi
"$PY_AF3" af3_pipeline/score_outputs.py --experiment "$EXPERIMENT"
"$PY_BIO" af3_pipeline/score_disulfide_ring.py \
  --experiment "$EXPERIMENT" --no-designed-pairs
"$PY_AF3" af3_pipeline/combine_seed1_multiseed_validation.py \
  --seed1-experiment "$SEED1" \
  --seed1-strict-per-model \
    "$SEED1/seed1_strict_final_classified/seed1_strict_per_model.tsv" \
  --multiseed-experiment "$EXPERIMENT" \
  "${strict_args[@]}" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --out "$EXPERIMENT/full25_summary.tsv"
touch "$EXPERIMENT/FULL25_COMBINE_COMPLETE"
echo "$EXPERIMENT strict 25-model summary complete"
