#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
SOURCE_EXP="$PROJECT_ROOT/experiments/e115_c4_iptm80_design"
EXPERIMENT="$PROJECT_ROOT/experiments/e119_c4_iptm80_hit_refine"

BASE=(H H H H F F F F)
STATE1=(AB AC AB AC AB AC AB AC)
STATE2=(AC AB AD AD AC AB AD AD)
BUDGET=(20 20 22 22 20 20 22 22)
W1=(6 6 7 7 6 6 7 7)
W2=(7 7 6 6 7 7 6 6)
FREE1=(0 0 0 1 0 0 0 1)
FREE2=(0 0 1 1 0 0 1 1)

if (( TASK_INDEX < 0 || TASK_INDEX >= ${#BASE[@]} )); then
  echo "invalid task index $TASK_INDEX" >&2
  exit 2
fi

RUN="$EXPERIMENT/run_$(printf '%02d' "$TASK_INDEX")"
S1="${STATE1[$TASK_INDEX]}"
S2="${STATE2[$TASK_INDEX]}"
STAGE1="$RUN/stage1_${S1}"
STAGE2="$RUN/stage2_${S2}"
mkdir -p "$STAGE1" "$STAGE2" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$STAGE2/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi

case "${BASE[$TASK_INDEX]}" in
  H) BASE_FASTA="$EXPERIMENT/hit_h.fasta" ;;
  F) BASE_FASTA="$EXPERIMENT/hit_f.fasta" ;;
  *) exit 3 ;;
esac

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"
REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
FREE1_ARG=(); FREE2_ARG=()
[[ "${FREE1[$TASK_INDEX]}" == 1 ]] && FREE1_ARG=(--free-interchain)
[[ "${FREE2[$TASK_INDEX]}" == 1 ]] && FREE2_ARG=(--free-interchain)

echo "task=$TASK_INDEX base=${BASE[$TASK_INDEX]} states=$S1->$S2 budget=${BUDGET[$TASK_INDEX]}"

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S1}.pdb" --chains "${S1:0:1},${S1:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$BASE_FASTA" --init current --init-fasta "$BASE_FASTA" \
  --design-positions "$SOURCE_EXP/design_${S1}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S1}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE1" --seed "$((310 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 16 --w-sparse 1.4 --w-fold 1.0 --w-first "${W1[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0 \
  "${FREE1_ARG[@]}"

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S2}.pdb" --chains "${S2:0:1},${S2:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" --init current --init-fasta "$BASE_FASTA" \
  --design-positions "$SOURCE_EXP/design_${S2}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S2}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE2" --seed "$((410 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 16 --w-sparse 1.4 --w-fold 1.0 --w-first "${W2[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0 \
  "${FREE2_ARG[@]}"

echo "task $TASK_INDEX complete"
