#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
SOURCE_EXP="$PROJECT_ROOT/experiments/e115_c4_iptm80_design"
EXPERIMENT="$PROJECT_ROOT/experiments/e127_c4_iptm80_expand_design"
FINAL_FASTA="$PROJECT_ROOT/experiments/e123_final_iptm80_c4/final_candidates.fasta"

BASE=(19A 19A 19A 19B 19B 19C 19C 19C)
STATE1=(AB AC AB AB AC AB AC AC)
STATE2=(AC AB AD AC AD AC AB AD)
BUDGET=(19 20 21 19 21 20 20 22)
W1=(7 7 7 7 7 7 7 7)
W2=(7 7 6 7 6 7 7 6)
FREE1=(0 0 0 0 0 0 0 0)
FREE2=(0 0 1 0 1 0 0 1)

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

BASE_NAME="OVA-C4-IPTM80-${BASE[$TASK_INDEX]}"
REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
FREE1_ARG=(); FREE2_ARG=()
[[ "${FREE1[$TASK_INDEX]}" == 1 ]] && FREE1_ARG=(--free-interchain)
[[ "${FREE2[$TASK_INDEX]}" == 1 ]] && FREE2_ARG=(--free-interchain)

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"
echo "task=$TASK_INDEX base=$BASE_NAME states=$S1->$S2 budget=${BUDGET[$TASK_INDEX]}"

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S1}.pdb" --chains "${S1:0:1},${S1:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$FINAL_FASTA" --base-name "$BASE_NAME" \
  --init current --init-fasta "$FINAL_FASTA" --init-name "$BASE_NAME" \
  --design-positions "$SOURCE_EXP/design_${S1}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S1}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE1" --seed "$((510 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 18 --w-sparse 1.5 --w-fold 1.0 --w-first "${W1[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0 \
  "${FREE1_ARG[@]}"

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S2}.pdb" --chains "${S2:0:1},${S2:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" \
  --init current --init-fasta "$STAGE1/af2diff_candidates.fasta" \
  --design-positions "$SOURCE_EXP/design_${S2}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S2}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE2" --seed "$((610 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 18 --w-sparse 1.5 --w-fold 1.0 --w-first "${W2[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0 \
  "${FREE2_ARG[@]}"

echo "task $TASK_INDEX complete"
