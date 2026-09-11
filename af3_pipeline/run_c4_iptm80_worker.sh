#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="$PROJECT_ROOT/experiments/e115_c4_iptm80_design"

# The two optimized pair states are deliberately varied across the 16 GPUs.
# AB and AC are the two recurrent strong interfaces in the final 17-B AF3
# backbone; AD is the weaker diagonal interface and is included in half of the
# trajectories as a direct route to raising global ipTM.
BASE=(B A B A B A B A B A B A B A B A)
STATE1=(AB AB AC AC AB AB AC AC AB AC AB AC AB AC AB AC)
STATE2=(AC AC AB AB AD AD AD AD AC AB AD AD AC AB AD AD)
INIT1=(OLD B OLD A A OLD B OLD OLD A B OLD A B OLD A)
INIT2=(B A A B OLD OLD B A A B OLD B B A A OLD)
BUDGET=(20 21 22 23 24 22 23 24 21 22 23 24 24 23 22 24)
W1=(4 5 4 5 4 5 4 5 6 6 5 5 6 6 6 6)
W2=(5 4 5 4 6 6 5 5 5 5 7 7 6 6 7 7)
FREE1=(0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1)
FREE2=(0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1)

if (( TASK_INDEX < 0 || TASK_INDEX >= ${#BASE[@]} )); then
  echo "invalid task index $TASK_INDEX" >&2
  exit 2
fi

RUN="$EXPERIMENT/run_$(printf '%02d' "$TASK_INDEX")"
STAGE1="$RUN/stage1_${STATE1[$TASK_INDEX]}"
STAGE2="$RUN/stage2_${STATE2[$TASK_INDEX]}"
mkdir -p "$STAGE1" "$STAGE2" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$STAGE2/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi

export PYTHONUNBUFFERED=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"

REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
SEQ_B="$PROJECT_ROOT/experiments/e114_final_lowmut_c4_preferred/OVA-LM-C4-17-B.fasta"
SEQ_A="$PROJECT_ROOT/experiments/e113_final_lowmut_c4/OVA-LM-C4-17-A.fasta"
SEQ_OLD="$PROJECT_ROOT/experiments/e106_final_noncovalent_c4/OVA-NC-C4-D2-92-AI01.fasta"

sequence_path() {
  case "$1" in
    A) echo "$SEQ_A" ;;
    B) echo "$SEQ_B" ;;
    OLD) echo "$SEQ_OLD" ;;
    *) echo "unknown sequence alias $1" >&2; exit 3 ;;
  esac
}

BASE_FASTA="$(sequence_path "${BASE[$TASK_INDEX]}")"
INIT1_FASTA="$(sequence_path "${INIT1[$TASK_INDEX]}")"
INIT2_FASTA="$(sequence_path "${INIT2[$TASK_INDEX]}")"
S1="${STATE1[$TASK_INDEX]}"
S2="${STATE2[$TASK_INDEX]}"
FREE1_ARG=(); FREE2_ARG=()
[[ "${FREE1[$TASK_INDEX]}" == 1 ]] && FREE1_ARG=(--free-interchain)
[[ "${FREE2[$TASK_INDEX]}" == 1 ]] && FREE2_ARG=(--free-interchain)

echo "task=$TASK_INDEX base=${BASE[$TASK_INDEX]} states=$S1->$S2 budget=${BUDGET[$TASK_INDEX]} init=${INIT1[$TASK_INDEX]}/${INIT2[$TASK_INDEX]}"

"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$EXPERIMENT/target_${S1}.pdb" --chains "${S1:0:1},${S1:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$BASE_FASTA" --init current --init-fasta "$INIT1_FASTA" \
  --design-positions "$EXPERIMENT/design_${S1}.txt" \
  --first-restraint-positions "$EXPERIMENT/restraint_${S1}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE1" --seed "$((100 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 14 --w-sparse 1.2 --w-fold 0.9 --w-first "${W1[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 12 --soft-iters 7 --hard-iters 4 --recycles 0 \
  "${FREE1_ARG[@]}"

"$BINDCRAFT_ROOT/env/bin/python" \
  "$PROJECT_ROOT/af3_pipeline/af2_tied_c4_lowmut_design.py" \
  --target-pdb "$EXPERIMENT/target_${S2}.pdb" --chains "${S2:0:1},${S2:1:1}" \
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R \
  --base-fasta "$STAGE1/af2diff_candidates.fasta" --init current --init-fasta "$INIT2_FASTA" \
  --design-positions "$EXPERIMENT/design_${S2}.txt" \
  --first-restraint-positions "$EXPERIMENT/restraint_${S2}.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$STAGE2" --seed "$((200 + TASK_INDEX))" \
  --mutation-budget "${BUDGET[$TASK_INDEX]}" \
  --w-budget 14 --w-sparse 1.2 --w-fold 0.9 --w-first "${W2[$TASK_INDEX]}" --w-second 0 \
  --logits-iters 12 --soft-iters 7 --hard-iters 4 --recycles 0 \
  "${FREE2_ARG[@]}"

echo "task $TASK_INDEX complete"
