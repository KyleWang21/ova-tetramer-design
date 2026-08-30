#!/bin/bash
# Run one OVA AF3 shard on a local or Volcengine A100 worker.
set -uo pipefail

ROOT=/root/400083
PROJECT="$ROOT/ova_p3_tetramer_design"
AF3="$ROOT/alphafold3"
EXP="${1:?usage: run_worker.sh EXPERIMENT_DIR SHARD}"
SHARD="${2:?usage: run_worker.sh EXPERIMENT_DIR SHARD}"
IN="$EXP/in_s$SHARD"
OUT="$EXP/out_s$SHARD"

mkdir -p "$OUT" "$EXP/logs" "$PROJECT/experiments/_jaxcache"
exec >> "$EXP/logs/af3_s${SHARD}.log" 2>&1
echo "$(date -u '+%F %T UTC') start shard=$SHARD host=$(hostname) input=$IN"
nvidia-smi -L

if [ ! -e /vepfs-mlp2/c20250508/400083 ]; then
  mkdir -p /vepfs-mlp2/c20250508
  ln -s /root/400083 /vepfs-mlp2/c20250508/400083
fi

count=$(find "$IN" -maxdepth 1 -name '*.json' | wc -l)
if [ "$count" -eq 0 ]; then
  echo "$(date -u '+%F %T UTC') no inputs for shard=$SHARD"
  exit 0
fi

PY="$AF3/repo/.venv/bin/python"
XLA_PYTHON_CLIENT_PREALLOCATE=false \
env -u LD_LIBRARY_PATH CUDA_VISIBLE_DEVICES=0 \
  "$PY" "$AF3/repo/run_alphafold.py" \
  --input_dir="$IN" --output_dir="$OUT" \
  --model_dir="$AF3/models" --norun_data_pipeline \
  --noresolve_msa_overlaps \
  --jax_compilation_cache_dir="$PROJECT/experiments/_jaxcache"
rc=$?

# Full atom-pair confidence JSONs are large and are not needed for the initial
# oligomer gate. Keep summary confidences, CIFs, terms of use, and input JSONs.
find "$OUT" -name '*_confidences.json' ! -name '*_summary_confidences.json' -delete 2>/dev/null
echo "$(date -u '+%F %T UTC') finish shard=$SHARD rc=$rc"
exit "$rc"
