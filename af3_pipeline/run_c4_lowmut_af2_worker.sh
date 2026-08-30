#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="$PROJECT_ROOT/experiments/e109_tied_c4_lowmut_af2"

INITS=(current current current current reference reference reference gumbel current current current current reference reference reference gumbel)
BUDGETS=(16 18 20 22 16 20 20 24 16 18 20 24 18 18 22 24)
FREE_IC=(0 0 0 0 0 0 1 0 1 1 1 1 0 1 1 0)
W_BUDGET=(14 12 10 8 14 10 10 7 14 12 10 7 12 12 8 7)
W_SPARSE=(2.0 1.5 1.2 1.0 2.0 1.2 1.2 0.8 2.0 1.5 1.2 0.8 1.5 1.5 1.0 0.8)

if (( TASK_INDEX < 0 || TASK_INDEX >= ${#INITS[@]} )); then
  echo "invalid task index $TASK_INDEX" >&2
  exit 2
fi

OUT="$EXPERIMENT/run_$(printf '%02d' "$TASK_INDEX")"
mkdir -p "$OUT" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$OUT/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"

FREE_ARG=()
if [[ "${FREE_IC[$TASK_INDEX]}" == 1 ]]; then
  FREE_ARG=(--free-interchain)
fi

exec "$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$PROJECT_ROOT/experiments/e106_final_noncovalent_c4/interface_hotspots/c4_ova_only_backbone.pdb" \
  --reference-fasta "$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta" \
  --reference-name D0-P3-13R \
  --current-fasta "$PROJECT_ROOT/experiments/e106_final_noncovalent_c4/OVA-NC-C4-D2-92-AI01.fasta" \
  --design-positions "$EXPERIMENT/design_positions.txt" \
  --params "$BINDCRAFT_ROOT/params" \
  --out "$OUT" --seed "$TASK_INDEX" \
  --init "${INITS[$TASK_INDEX]}" \
  --mutation-budget "${BUDGETS[$TASK_INDEX]}" \
  --w-budget "${W_BUDGET[$TASK_INDEX]}" \
  --w-sparse "${W_SPARSE[$TASK_INDEX]}" \
  --w-fold 0.75 --w-first 2.5 --w-second 2.5 \
  --logits-iters 10 --soft-iters 6 --hard-iters 3 --recycles 0 \
  "${FREE_ARG[@]}"
