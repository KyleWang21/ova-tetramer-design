#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="$PROJECT_ROOT/experiments/e109_tied_c4_lowmut_af2"
TOTAL_BUDGET=(16 18 20 22 18 20 22 24 16 18 20 24 16 20 22 24)
SECOND_FREE=(0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0)

if (( TASK_INDEX < 0 || TASK_INDEX >= ${#TOTAL_BUDGET[@]} )); then
  echo "invalid task index $TASK_INDEX" >&2
  exit 2
fi

RUN="$EXPERIMENT/run_$(printf '%02d' "$TASK_INDEX")"
STAGE1="$RUN/stage1_first_interface"
OUT="$RUN/stage2_second_interface_warm"
REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
CURRENT="$PROJECT_ROOT/experiments/e106_final_noncovalent_c4/OVA-NC-C4-D2-92-AI01.fasta"
[[ -s "$STAGE1/af2diff_candidates.fasta" ]] || { echo "missing stage1 for task $TASK_INDEX" >&2; exit 3; }
mkdir -p "$OUT" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_warm_$TASK_INDEX"
[[ ! -s "$OUT/af2diff_candidates.tsv" ]] || { echo "warm task $TASK_INDEX already complete"; exit 0; }

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_warm_$TASK_INDEX"
FREE_ARG=(); [[ "${SECOND_FREE[$TASK_INDEX]}" == 1 ]] && FREE_ARG=(--free-interchain)

# Fixed positions come from the sparse first-interface result, while mutable
# second-interface positions start from the validated 42-mutation C4 sequence.
# The optimizer can therefore remove redundant second-interface mutations
# without having to discover that interface from an all-P3 initialization.
"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$EXPERIMENT/target_second_AC.pdb" --chains A,C \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" \
  --init current --init-fasta "$CURRENT" \
  --design-positions "$EXPERIMENT/second_positions.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$OUT" --seed "$TASK_INDEX" \
  --mutation-budget "${TOTAL_BUDGET[$TASK_INDEX]}" \
  --w-budget 12 --w-sparse 1.5 --w-fold 0.75 --w-first 0 --w-second 3.0 \
  --logits-iters 10 --soft-iters 6 --hard-iters 3 --recycles 0 \
  "${FREE_ARG[@]}"

echo "warm second-interface task $TASK_INDEX complete"
