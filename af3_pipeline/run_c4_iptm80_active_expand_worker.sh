#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
SOURCE_EXP="$PROJECT_ROOT/experiments/e115_c4_iptm80_design"
EXPERIMENT="${OVA_EXPAND_EXPERIMENT:-$PROJECT_ROOT/experiments/e137_c4_iptm80_active_expand_design}"
SEED_OFFSET="${OVA_EXPAND_SEED_OFFSET:-0}"
ANCHORS="$EXPERIMENT/anchors.fasta"

BASE=(MS01 MS02 MS03 MS04 MS05 MS06 MS07 MS08)
STATE1=(AB AC AB AC AB AC AB AC)
STATE2=(AC AB AD AD AC AB AD AD)
BUDGET=(20 20 22 22 21 21 22 22)
W1=(8 8 7 7 8 8 7 7)
W2=(8 8 7 7 8 8 7 7)
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

BASE_NAME="OVA-IPTM80-${BASE[$TASK_INDEX]}"
REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
FREE2_ARG=()
[[ "${FREE2[$TASK_INDEX]}" == 1 ]] && FREE2_ARG=(--free-interchain)
export PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"
echo "task=$TASK_INDEX base=$BASE_NAME states=$S1->$S2 budget=${BUDGET[$TASK_INDEX]}"

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S1}.pdb" --chains "${S1:0:1},${S1:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$ANCHORS" --base-name "$BASE_NAME" \
  --init current --init-fasta "$ANCHORS" --init-name "$BASE_NAME" \
  --design-positions "$SOURCE_EXP/design_${S1}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S1}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE1" --seed "$((710 + SEED_OFFSET + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 18 --w-sparse 1.5 --w-fold 1.0 --w-first "${W1[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0

"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE_EXP/target_${S2}.pdb" --chains "${S2:0:1},${S2:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" \
  --init current --init-fasta "$STAGE1/af2diff_candidates.fasta" \
  --design-positions "$SOURCE_EXP/design_${S2}.txt" \
  --first-restraint-positions "$SOURCE_EXP/restraint_${S2}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE2" --seed "$((810 + SEED_OFFSET + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 18 --w-sparse 1.5 --w-fold 1.0 --w-first "${W2[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 15 --soft-iters 8 --hard-iters 5 --recycles 0 \
  "${FREE2_ARG[@]}"

echo "task $TASK_INDEX complete"
