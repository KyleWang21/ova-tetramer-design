#!/usr/bin/env bash
# Promote up to four strict seed-1 hits through independent AF3 seeds2--5.
# Seed2 is fail-fast; seeds3+4 run together (<=8 GPUs), then seed5.
set -euo pipefail

SOURCE="${1:?source seed1 experiment required}"
MULTISEED="${2:?multiseed experiment required}"
JOB_BASE="${3:?Volc job-number base required}"
TOP="${4:-4}"
PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/alphafold3/repo/.venv/bin/python
SURFACE="$PROJECT/experiments/e155_c4_iptm90_crossmodel_design/surface_positions.txt"
LOCK="$PROJECT/experiments/VOLC_PRIORITY_BATCH.lock"
cd "$PROJECT"

update_registry() {
  "$PY" af3_pipeline/update_accepted_af3_registry.py \
    --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
  pending=$(awk -F '\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="needs_v1_screen") need=i
        if ($i=="source_kind") kind=i
      }
      next
    }
    $kind=="generic_25model_promotion" {n += $need}
    END {print n+0}
  ' "$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv")
  if [[ "$pending" -gt 0 ]]; then
    touch "$PROJECT/experiments/current_accepted_af3_registry/NEW_ACCEPTED_AF3_PENDING_V1"
  fi
  bash af3_pipeline/refresh_current_accepted_delivery.sh
}

if [[ -f "$MULTISEED/MULTISEED_VALIDATION_COMPLETE" ]]; then
  echo "$MULTISEED already complete"
  update_registry
  exit 0
fi

"$PY" af3_pipeline/prepare_seed1_hits_multiseed.py \
  --source "$SOURCE" --out "$MULTISEED" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R --surface-positions "$SURFACE" --top "$TOP"

if [[ -f "$MULTISEED/NO_SEED1_HITS" ]]; then
  mkdir -p "$MULTISEED/full25"
  "$PY" af3_pipeline/summarize_seed1_hits_multiseed.py \
    --source "$SOURCE" --multiseed "$MULTISEED" --out "$MULTISEED/full25"
  touch "$MULTISEED/MULTISEED_VALIDATION_COMPLETE"
  update_registry
  exit 0
fi

while [[ -e "$LOCK" ]]; do sleep 30; done
printf '%s\t%s\n' "$$" "$MULTISEED" > "$LOCK"
cleanup() {
  if [[ -f "$LOCK" ]] && grep -q "^$$" "$LOCK"; then rm -f "$LOCK"; fi
}
trap cleanup EXIT

join_csv() {
  local IFS=,
  echo "$*"
}

submit_phase() {
  local label="$1" expected="$2"
  shift 2
  local shards=("$@")
  [[ "$expected" -gt 0 ]] || return 0
  local marker="VOLC_${label}_SHARDS_SUBMITTED"
  local record="volc_${label}_jobs.tsv"
  local spec
  spec=$(join_csv "${shards[@]}")
  while [[ ! -f "$MULTISEED/$marker" ]]; do
    "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
      --experiment "$MULTISEED" --job-number-start "$JOB_BASE" \
      --shards "$spec" --expected "$expected" --record-file "$record" \
      --marker "$marker" --priority
    [[ -f "$MULTISEED/$marker" ]] && break
    sleep 30
  done
  for shard in "${shards[@]}"; do
    while true; do
      if [[ -f "$MULTISEED/WORKER_FAILED_s$shard" ]]; then
        echo "AF3 worker failure marker detected for $MULTISEED shard $shard" >&2
        exit 1
      fi
      complete=$(find "$MULTISEED/out_s$shard" \
        -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
      [[ "$complete" -eq 5 ]] && break
      sleep 30
    done
    if [[ ! -s "$MULTISEED/strict_s$shard/seed1_strict_summary.tsv" ]]; then
      "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
        --experiment "$MULTISEED" \
        --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
        --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
        --shard "$shard" --out "$MULTISEED/strict_s$shard"
    fi
  done
}

strict_pass() {
  local shard="$1"
  awk -F '\t' '
    NR==1 {for (i=1; i<=NF; i++) if ($i=="seed1_strict_pass") c=i; next}
    NR==2 {print $c}
  ' "$MULTISEED/strict_s$shard/seed1_strict_summary.tsv"
}

mapfile -t base_shards < <(awk -F '\t' 'NR>1 {print $4}' "$MULTISEED/selected_candidates.tsv")
submit_phase seed2 "${#base_shards[@]}" "${base_shards[@]}"

seed2_survivors=()
for base in "${base_shards[@]}"; do
  [[ "$(strict_pass "$base")" == 1 ]] && seed2_survivors+=("$base")
done

seed34_shards=()
for base in "${seed2_survivors[@]}"; do
  seed34_shards+=("$((base + 1))" "$((base + 2))")
done
submit_phase seed34 "${#seed34_shards[@]}" "${seed34_shards[@]}"

seed34_survivors=()
for base in "${seed2_survivors[@]}"; do
  if [[ "$(strict_pass "$((base + 1))")" == 1 && \
        "$(strict_pass "$((base + 2))")" == 1 ]]; then
    seed34_survivors+=("$base")
  fi
done

seed5_shards=()
for base in "${seed34_survivors[@]}"; do seed5_shards+=("$((base + 3))"); done
submit_phase seed5 "${#seed5_shards[@]}" "${seed5_shards[@]}"

"$PY" af3_pipeline/summarize_seed1_hits_multiseed.py \
  --source "$SOURCE" --multiseed "$MULTISEED" --out "$MULTISEED/full25"
touch "$MULTISEED/MULTISEED_VALIDATION_COMPLETE"
update_registry
echo "$MULTISEED seeds2-5 fail-fast validation complete"
