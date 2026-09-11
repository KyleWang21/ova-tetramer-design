#!/usr/bin/env bash
# OpenDDE-v1 evidence run. Use only after installing the official preview release
# and downloading opendde.pt under OPENDDE_ROOT_DIR.
set -uo pipefail

export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

EXP="${1:?usage: run_opendde_c4_shard.sh EXP SHARD}"
SHARD="${2:?usage: run_opendde_c4_shard.sh EXP SHARD}"
OPENDDE_BIN="${OPENDDE_BIN:?set OPENDDE_BIN to the installed official CLI}"
OUT="$EXP/out"
MANIFEST="${MANIFEST:-$EXP/manifest.tsv}"
mkdir -p "$OUT" "$EXP/logs"
exec >> "$EXP/logs/opendde_s${SHARD}.log" 2>&1

while IFS=$'\t' read -r row_shard name kind n_seeds input_json; do
  input_json="${input_json%$'\r'}"
  [ "$row_shard" = "shard" ] && continue
  [ "$row_shard" = "$SHARD" ] || continue
  [[ "$input_json" = /* ]] || input_json="$EXP/$input_json"
  echo "$(date -u '+%F %T UTC') OpenDDE name=$name kind=$kind seeds=$n_seeds"
  job_name="$(${OPENDDE_BIN%/opendde}/python - "$input_json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))[0]["name"])
PY
)"
  CUDA_VISIBLE_DEVICES=0 "$OPENDDE_BIN" pred \
    -i "$input_json" -o "$OUT" -n opendde_v1 \
    --use_msa true --use_template false --use_rna_msa false \
    --trimul_kernel cuequivariance --triatt_kernel cuequivariance \
    --dtype bf16 --sample 1 --step 200 --cycle 10
  rc=$?
  produced=$(find "$OUT" -path "*/$job_name/seed_*/predictions/*_summary_confidence_sample_0.json" | wc -l)
  if [ "$produced" -lt "$n_seeds" ]; then
    echo "$(date -u '+%F %T UTC') name=$name rc=$rc produced=$produced expected=$n_seeds fail_closed=1" >&2
    exit 90
  fi
  echo "$(date -u '+%F %T UTC') name=$name rc=$rc produced=$produced"
  [ "$rc" -eq 0 ] || exit "$rc"
done < "$MANIFEST"
