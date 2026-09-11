#!/usr/bin/env bash
# Score the same 512 exact sequences on one of eight independent C4 backbones.
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="${OVA_MPNN_EXPERIMENT:-$PROJECT_ROOT/experiments/e271_c4_multibackbone_exact_mpnn_score}"
if (( TASK_INDEX < 0 || TASK_INDEX >= 8 )); then
  echo "invalid TASK_INDEX=$TASK_INDEX" >&2
  exit 2
fi
row=$((TASK_INDEX + 2))
IFS=$'\t' read -r manifest_index context family target_pdb source_mask \
  < <(sed -n "${row}p" "$EXPERIMENT/backbone_manifest.tsv")
[[ "$manifest_index" == "$TASK_INDEX" ]] || {
  echo "manifest/task mismatch: $manifest_index != $TASK_INDEX" >&2
  exit 2
}
out="$EXPERIMENT/task_${TASK_INDEX}_${context}"
mkdir -p "$out" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$out/candidate_scores.tsv" ]]; then
  echo "E271 task $TASK_INDEX already complete: $out"
  exit 0
fi
export PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"
"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/score_multibackbone_exact_mpnn.py" \
  --target-pdb "$PROJECT_ROOT/$target_pdb" \
  --candidates "$EXPERIMENT/candidates.tsv" \
  --score-positions "$EXPERIMENT/score_positions.txt" \
  --out "$out" --task-index "$TASK_INDEX" --context "$context" --family "$family" \
  --orders 4 --batch-size 32 --seed 20260902 --expected-candidates 512
echo "E271 task $TASK_INDEX complete: $context"
