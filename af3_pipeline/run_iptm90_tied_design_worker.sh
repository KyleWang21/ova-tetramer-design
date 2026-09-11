#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="${OVA_IPTM90_EXPERIMENT:-$PROJECT_ROOT/experiments/e155_c4_iptm90_crossmodel_design}"
REFERENCE="$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta"
MANIFEST="$EXPERIMENT/target_manifest.tsv"
cd "$PROJECT_ROOT"

if (( TASK_INDEX < 0 || TASK_INDEX >= 8 )); then
  echo "invalid TASK_INDEX=$TASK_INDEX" >&2
  exit 2
fi
anchor_row=$((TASK_INDEX / 2 + 2))
replicate=$((TASK_INDEX % 2))
IFS=$'\t' read -r anchor anchor_n source_cif source_iptm source_clash pair1 pair2 target_c4 target1 design1 restraint1 contacts1 ndesign1 target2 design2 restraint2 contacts2 ndesign2 joint_design joint_ndesign \
  < <(sed -n "${anchor_row}p" "$MANIFEST")
if (( replicate == 0 )); then
  w1="${OVA_STATE1_WEIGHT_R0:-7}"
  w2="${OVA_STATE2_WEIGHT_R0:-5}"
else
  w1="${OVA_STATE1_WEIGHT_R1:-5}"
  w2="${OVA_STATE2_WEIGHT_R1:-7}"
fi

run_prefix="${OVA_IPTM90_RUN_PREFIX:-run}"
run="$EXPERIMENT/${run_prefix}_$(printf '%02d' "$TASK_INDEX")"
joint="$run/joint_${pair1}_${pair2}_r${replicate}"
mkdir -p "$joint" "$EXPERIMENT/logs" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$joint/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi
export PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"

common=(
  --reference-fasta "$REFERENCE" --reference-name D0-P3-13R
  --surface-positions "$EXPERIMENT/surface_positions.txt"
  --params "$BINDCRAFT_ROOT/params"
  --mutation-budget "${OVA_MUTATION_BUDGET:-24}"
  --min-surface-mutation-fraction "${OVA_MIN_SURFACE_MUTATION_FRACTION:-0.80}"
  --w-budget "${OVA_W_BUDGET:-20}"
  --w-sparse "${OVA_W_SPARSE:-1.5}"
  --w-buried-mutation "${OVA_W_BURIED_MUTATION:-12}"
  --w-fold 1.0 \
  --worst-state-gradient-boost "${OVA_WORST_STATE_GRADIENT_BOOST:-1.0}" \
  --logits-iters "${OVA_LOGITS_ITERS:-18}" \
  --soft-iters "${OVA_SOFT_ITERS:-10}" \
  --hard-iters "${OVA_HARD_ITERS:-6}" --recycles "${OVA_RECYCLES:-0}"
)

echo "task=$TASK_INDEX anchor=$anchor shared_tied_sequence=1 interfaces=$pair1+$pair2 weights=$w1,$w2 source_iptm=$source_iptm"
"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_multistate_tied_lowmut_design.py" \
  --target-pdb "$target1" --chains "${pair1:0:1},${pair1:1:1}" \
  --restraint-positions "$restraint1" --state-weight "$w1" \
  --target-pdb "$target2" --chains "${pair2:0:1},${pair2:1:1}" \
  --restraint-positions "$restraint2" --state-weight "$w2" \
  --base-fasta "$EXPERIMENT/anchors.fasta" --base-name "$anchor" \
  --design-positions "$joint_design" \
  --out "$joint" --seed "$((1201 + TASK_INDEX + ${OVA_SEED_OFFSET:-0}))" \
  "${common[@]}"

echo "task $TASK_INDEX complete"
