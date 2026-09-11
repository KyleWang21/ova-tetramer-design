#!/usr/bin/env bash
# Prioritize both E178 strict hits; use local slack for clash-informed 073 reversions.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E178="$PROJECT/experiments/e178_c4_073_af3ensemble_ai_af3_seed1"
E200="$PROJECT/experiments/e200_c4_e178_04_multiseed"
E201="$PROJECT/experiments/e201_c4_e178_15_multiseed"
E168="$PROJECT/experiments/e168_c4_073_reversion_ai_af3_seed1"
LOCK="$PROJECT/experiments/VOLC_PRIORITY_BATCH.lock"
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

cleanup() {
  rm -f "$LOCK"
}
trap cleanup EXIT
touch "$LOCK"

deadline=$((SECONDS + 172800))
while true; do
  complete=$(find "$E178" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
  if [[ "$complete" -eq 80 ]]; then
    touch "$E178/SEED1_ALL_MODELS_COMPLETE"
    break
  fi
  if (( SECONDS >= deadline )); then
    echo "timed out waiting for E178: $complete/80" >&2
    exit 2
  fi
  sleep 30
done

ran_e168=0
set +e
"$PY" af3_pipeline/submit_e200_e201_priority_if_capacity.py
remote_rc=$?
set -e
if [[ "$remote_rc" -eq 0 ]]; then
  while true; do
    e200_complete=$(find "$E200" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    e201_complete=$(find "$E201" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$e200_complete" -eq 20 && "$e201_complete" -eq 20 ]] && break
    if (( SECONDS >= deadline )); then
      echo "timed out waiting for priority Volc jobs: E200=$e200_complete/20 E201=$e201_complete/20" >&2
      exit 2
    fi
    sleep 30
  done
  for shard in 0 1 2 3; do
    seed=$((shard + 2))
    for experiment in "$E200" "$E201"; do
      "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
        --experiment "$experiment" \
        --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
        --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
        --shard "$shard" --out "$experiment/seed${seed}_strict"
    done
  done
  touch "$E200/RESEED_DECISION_COMPLETE"
  touch "$E201/RESEED_DECISION_COMPLETE"
else
  if [[ "$remote_rc" -ne 10 ]]; then
    echo "priority remote submission failed rc=$remote_rc" >&2
    exit "$remote_rc"
  fi
  bash af3_pipeline/run_single_multiseed_failfast.sh \
    "$E200" "$E178/SEED1_ALL_MODELS_COMPLETE"
  bash af3_pipeline/run_single_multiseed_failfast.sh \
    "$E201" "$E178/SEED1_ALL_MODELS_COMPLETE"
  for shard in 0 1 2 3; do
    found=$(find "$E168/out_s$shard" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$found" -eq 10 ]] || bash af3_pipeline/run_worker.sh "$E168" "$shard"
  done
  ran_e168=1
fi

if [[ "$ran_e168" -eq 1 ]]; then
  for shard in 0 1 2 3; do
    "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
      --experiment "$E168" \
      --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
      --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
      --shard "$shard" --out "$E168/seed1_strict_s$shard"
  done
  touch "$E168/PRIORITY_SHARDS_0_3_SCREEN_COMPLETE"
fi
echo "E200/E201 priority multiseed complete; E168_local_screen=$ran_e168"
