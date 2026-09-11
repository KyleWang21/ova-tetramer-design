#!/usr/bin/env bash
# Promote strict seed1 hits after the first Top4 while keeping <=8 concurrent GPUs.
set -euo pipefail

SOURCE="${1:?source seed1 experiment required}"
MULTISEED="${2:?multiseed experiment required}"
JOB_BASE="${3:?Volc job-number base required}"
SKIP="${4:-4}"
TOP="${5:-5}"
PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/alphafold3/repo/.venv/bin/python
SURFACE="$PROJECT/experiments/e155_c4_iptm90_crossmodel_design/surface_positions.txt"
LOCK="$PROJECT/experiments/VOLC_PRIORITY_BATCH.lock"
cd "$PROJECT"

update_registry() {
  "$PY" af3_pipeline/update_accepted_af3_registry.py \
    --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
  pending=$(awk -F '\t' '
    NR==1 {for(i=1;i<=NF;i++){if($i=="needs_v1_screen")need=i;if($i=="source_kind")kind=i};next}
    $kind=="generic_25model_promotion" {n += $need}
    END {print n+0}
  ' "$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv")
  [[ "$pending" -gt 0 ]] && touch "$PROJECT/experiments/current_accepted_af3_registry/NEW_ACCEPTED_AF3_PENDING_V1"
  bash af3_pipeline/refresh_current_accepted_delivery.sh
}

if [[ -f "$MULTISEED/MULTISEED_VALIDATION_COMPLETE" ]]; then
  update_registry
  exit 0
fi

if [[ ! -f "$MULTISEED/CANDIDATES_FROZEN" && ! -f "$MULTISEED/NO_SEED1_HITS" ]]; then
  "$PY" af3_pipeline/prepare_seed1_hits_multiseed.py \
    --source "$SOURCE" --out "$MULTISEED" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --surface-positions "$SURFACE" \
    --skip "$SKIP" --top "$TOP"
fi

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

strict_pass() {
  local shard="$1"
  awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="seed1_strict_pass")c=i;next} NR==2{print $c}' \
    "$MULTISEED/strict_s$shard/seed1_strict_summary.tsv"
}

submit_chunk() {
  local label="$1"
  shift
  local shards=("$@")
  [[ "${#shards[@]}" -le 8 ]] || { echo "chunk exceeds 8 GPUs" >&2; exit 2; }
  [[ "${#shards[@]}" -gt 0 ]] || return 0
  local csv marker record
  csv=$(IFS=,; echo "${shards[*]}")
  marker="VOLC_${label}_SHARDS_SUBMITTED"
  record="volc_${label}_jobs.tsv"
  while [[ ! -f "$MULTISEED/$marker" ]]; do
    submitted=()
    if [[ -s "$MULTISEED/$record" ]]; then
      mapfile -t submitted < <(awk -F '\t' '
        NR==1 {for(i=1;i<=NF;i++){if($i=="shard")h=i;if($i=="mode")m=i;if($i=="status")s=i};next}
        $m=="submit" && $s=="ok" {print $h}
      ' "$MULTISEED/$record")
    fi
    missing=()
    for shard in "${shards[@]}"; do
      present=0
      for old in "${submitted[@]}"; do [[ "$old" == "$shard" ]] && present=1; done
      [[ "$present" -eq 1 ]] || missing+=("$shard")
    done
    if [[ "${#missing[@]}" -eq 0 ]]; then
      touch "$MULTISEED/$marker"
      break
    fi
    csv=$(IFS=,; echo "${missing[*]}")
    "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
      --experiment "$MULTISEED" --job-number-start "$JOB_BASE" \
      --shards "$csv" --expected "${#shards[@]}" \
      --record-file "$record" --marker "$marker" --priority
    [[ -f "$MULTISEED/$marker" ]] || sleep 30
  done
  for shard in "${shards[@]}"; do
    while true; do
      [[ ! -f "$MULTISEED/WORKER_FAILED_s$shard" ]] || {
        echo "AF3 worker failure shard=$shard" >&2; exit 1;
      }
      complete=$(find "$MULTISEED/out_s$shard" \
        -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
      [[ "$complete" -eq 5 ]] && break
      sleep 30
    done
    "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
      --experiment "$MULTISEED" \
      --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
      --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
      --shard "$shard" --out "$MULTISEED/strict_s$shard"
  done
}

mapfile -t bases < <(awk -F '\t' 'NR>1{print $4}' "$MULTISEED/selected_candidates.tsv")
submit_chunk seed2 "${bases[@]}"

seed34=()
for base in "${bases[@]}"; do
  [[ "$(strict_pass "$base")" == 1 ]] && seed34+=("$((base+1))" "$((base+2))")
done
chunk_index=0
for ((start=0; start<${#seed34[@]}; start+=8)); do
  submit_chunk "seed34_chunk$chunk_index" "${seed34[@]:start:8}"
  chunk_index=$((chunk_index+1))
done

seed5=()
for base in "${bases[@]}"; do
  if [[ "$(strict_pass "$base")" == 1 && \
        "$(strict_pass "$((base+1))")" == 1 && \
        "$(strict_pass "$((base+2))")" == 1 ]]; then
    seed5+=("$((base+3))")
  fi
done
submit_chunk seed5 "${seed5[@]}"

"$PY" af3_pipeline/summarize_seed1_hits_multiseed.py \
  --source "$SOURCE" --multiseed "$MULTISEED" --out "$MULTISEED/full25"
touch "$MULTISEED/MULTISEED_VALIDATION_COMPLETE"
update_registry
echo "$MULTISEED remaining strict hits validation complete"
