#!/usr/bin/env bash
# Polar-network-biased tied ProteinMPNN on four OVA-C4-073 Rosetta backbones.
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
SOURCE="${OVA_MPNN_SOURCE:-$PROJECT_ROOT/experiments/e225_c4_073_rosetta_polar_tied_mpnn}"
EXPERIMENT="${OVA_MPNN_EXPERIMENT:-$SOURCE}"

if (( TASK_INDEX < 0 || TASK_INDEX >= 8 )); then
  echo "invalid TASK_INDEX=$TASK_INDEX" >&2
  exit 2
fi
backbones=(bb1 bb1 bb2 bb2 bb3 bb3 bb4 bb4)
polar_biases=(1.0 1.8 1.0 1.8 1.0 1.8 1.0 1.8)
backbone="${backbones[$TASK_INDEX]}"
polar_bias="${polar_biases[$TASK_INDEX]}"
replicate=$((TASK_INDEX % 2))
prefix="${OVA_MPNN_RUN_PREFIX:-run_073polar1}"
out="$EXPERIMENT/${prefix}_${backbone}_p${polar_bias/./p}_r${replicate}"
mkdir -p "$out" "$EXPERIMENT/logs"
if [[ -s "$out/mpnn_shortlist.tsv" ]]; then
  echo "Rosetta-polar ProteinMPNN task $TASK_INDEX already complete: $out"
  exit 0
fi

export XLA_PYTHON_CLIENT_PREALLOCATE=false
"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/mpnn_tied_c4_lowmut_design.py" \
  --target-pdb "$SOURCE/targets/${backbone}_c4.pdb" \
  --reference-fasta "$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta" \
  --reference-name D0-P3-13R \
  --design-positions "$SOURCE/weak_interface_design_positions.txt" \
  --surface-positions "$SOURCE/surface_positions.txt" \
  --polar-bias-positions "$SOURCE/weak_interface_polar_positions.txt" \
  --polar-bias "$polar_bias" --polar-aas DEHKNQRSTY \
  --out "$out" \
  --temperatures "${OVA_MPNN_TEMPERATURES:-0.08,0.12}" \
  --native-biases "${OVA_MPNN_NATIVE_BIASES:-1.5,2.0,2.5}" \
  --buried-native-bias "${OVA_MPNN_BURIED_NATIVE_BIAS:-12.0}" \
  --samples-per-condition "${OVA_MPNN_SAMPLES_PER_CONDITION:-512}" \
  --batch-size "${OVA_MPNN_BATCH_SIZE:-32}" \
  --seed "$((20260901 + TASK_INDEX + ${OVA_MPNN_SEED_OFFSET:-68000}))" \
  --min-mutations 12 --max-mutations 22 \
  --min-surface-fraction 0.80 --shortlist 128

echo "Rosetta-polar ProteinMPNN task $TASK_INDEX complete: $out"
