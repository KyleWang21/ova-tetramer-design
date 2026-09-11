#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="$PROJECT_ROOT/experiments/e109_tied_c4_lowmut_af2"

# Each GPU runs one complete low-mutation dimer-of-dimers design trajectory:
# first interface on AF3 chains A/B, then second interface on A/C while holding
# the first-stage sequence fixed.  The mutation loss in stage 2 is still total
# mutations versus P3-13R, not just newly mutable positions.
INITS=(current current current current reference reference reference gumbel current current current current reference reference reference gumbel)
FIRST_BUDGET=(8 10 12 14 8 10 12 14 8 10 12 14 8 10 12 14)
TOTAL_BUDGET=(16 18 20 22 18 20 22 24 16 18 20 24 16 20 22 24)
FIRST_FREE=(0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0)
SECOND_FREE=(0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0)

if (( TASK_INDEX < 0 || TASK_INDEX >= ${#INITS[@]} )); then
  echo "invalid task index $TASK_INDEX" >&2
  exit 2
fi

OUT="$EXPERIMENT/run_$(printf '%02d' "$TASK_INDEX")"
STAGE1="$OUT/stage1_first_interface"
STAGE2="$OUT/stage2_second_interface"
mkdir -p "$OUT" "$STAGE1" "$STAGE2" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$STAGE2/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"

REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
CURRENT="$PROJECT_ROOT/experiments/e106_final_noncovalent_c4/OVA-NC-C4-D2-92-AI01.fasta"
INIT_ARG=(--init current --init-fasta "$CURRENT")
if [[ "${INITS[$TASK_INDEX]}" == reference ]]; then
  INIT_ARG=(--init reference --init-fasta "$REFERENCE")
elif [[ "${INITS[$TASK_INDEX]}" == gumbel ]]; then
  INIT_ARG=(--init gumbel --init-fasta "$REFERENCE")
fi
FIRST_FREE_ARG=(); SECOND_FREE_ARG=()
[[ "${FIRST_FREE[$TASK_INDEX]}" == 1 ]] && FIRST_FREE_ARG=(--free-interchain)
[[ "${SECOND_FREE[$TASK_INDEX]}" == 1 ]] && SECOND_FREE_ARG=(--free-interchain)

"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$EXPERIMENT/target_first_AB.pdb" --chains A,B \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$REFERENCE" "${INIT_ARG[@]}" \
  --design-positions "$EXPERIMENT/first_positions.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE1" --seed "$TASK_INDEX" \
  --mutation-budget "${FIRST_BUDGET[$TASK_INDEX]}" \
  --w-budget 12 --w-sparse 1.5 --w-fold 0.75 --w-first 3.0 --w-second 0 \
  --logits-iters 10 --soft-iters 6 --hard-iters 3 --recycles 0 \
  "${FIRST_FREE_ARG[@]}"

"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$EXPERIMENT/target_second_AC.pdb" --chains A,C \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" \
  --init current --init-fasta "$STAGE1/af2diff_candidates.fasta" \
  --design-positions "$EXPERIMENT/second_positions.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE2" --seed "$TASK_INDEX" \
  --mutation-budget "${TOTAL_BUDGET[$TASK_INDEX]}" \
  --w-budget 12 --w-sparse 1.5 --w-fold 0.75 --w-first 0 --w-second 3.0 \
  --logits-iters 10 --soft-iters 6 --hard-iters 3 --recycles 0 \
  "${SECOND_FREE_ARG[@]}"

echo "pipeline task $TASK_INDEX complete"
