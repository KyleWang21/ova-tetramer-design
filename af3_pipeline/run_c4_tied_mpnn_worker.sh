#!/usr/bin/env bash
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
SOURCE="$PROJECT_ROOT/experiments/e155_c4_iptm90_crossmodel_design"
EXPERIMENT="${OVA_MPNN_EXPERIMENT:-$PROJECT_ROOT/experiments/e166_c4_tied_proteinmpnn}"

anchors=(19k 19k 21b 21b 21s 21s 20a 20a)
if (( TASK_INDEX < 0 || TASK_INDEX >= ${#anchors[@]} )); then
  echo "invalid TASK_INDEX=$TASK_INDEX" >&2
  exit 2
fi
anchor="${anchors[$TASK_INDEX]}"
replicate=$((TASK_INDEX % 2))
prefix="${OVA_MPNN_RUN_PREFIX:-run}"
out="$EXPERIMENT/${prefix}_${anchor}_r${replicate}"
mkdir -p "$out" "$EXPERIMENT/logs"
if [[ -s "$out/mpnn_shortlist.tsv" ]]; then
  echo "ProteinMPNN task $TASK_INDEX already complete: $out"
  exit 0
fi

export XLA_PYTHON_CLIENT_PREALLOCATE=false
"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/mpnn_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE/targets/ova_c4_iptm80_${anchor}_c4.pdb" \
  --reference-fasta "$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta" \
  --reference-name D0-P3-13R \
  --design-positions "$SOURCE/targets/ova_c4_iptm80_${anchor}_joint_design.txt" \
  --surface-positions "$SOURCE/surface_positions.txt" \
  --out "$out" \
  --temperatures "${OVA_MPNN_TEMPERATURES:-0.08,0.12,0.20}" \
  --native-biases "${OVA_MPNN_NATIVE_BIASES:-0.5,1.0,1.5,2.0}" \
  --buried-native-bias "${OVA_MPNN_BURIED_NATIVE_BIAS:-1.0}" \
  --samples-per-condition "${OVA_MPNN_SAMPLES_PER_CONDITION:-128}" \
  --batch-size "${OVA_MPNN_BATCH_SIZE:-32}" \
  --seed "$((20260901 + TASK_INDEX + ${OVA_MPNN_SEED_OFFSET:-0}))" \
  --min-mutations 12 --max-mutations 24 --min-surface-fraction 0.80 --shortlist 128

echo "ProteinMPNN task $TASK_INDEX complete: $out"
