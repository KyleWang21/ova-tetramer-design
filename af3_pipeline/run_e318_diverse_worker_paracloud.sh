#!/usr/bin/env bash
# E318 second-generation eight-state tied-AF2 worker on Zhongwei.
# Tasks 0-3 reinforce the two weakest states; tasks 4-7 additionally open
# the cross-backbone hotspot mask and weight the second interface more.
#SBATCH --partition=vip_gpu_a800_scwb671
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

set -euo pipefail
TASK_INDEX="${TASK_INDEX:?submit with TASK_INDEX=0..7}"
ROOT=/data/home/scwb671/run/wky
export OVA_TETRA_ROOT="$ROOT/ova_p3_tetramer_design"
export BINDCRAFT_MOD="$ROOT/antibody-design/pipeline_stage6_maskcap_bindcraft"
export OVA_E318_EXPERIMENT="${OVA_E318_EXPERIMENT:-$OVA_TETRA_ROOT/experiments/e318_c4_e311_diverse_tied_design}"
cd "$OVA_TETRA_ROOT"
if (( TASK_INDEX < 0 || TASK_INDEX >= 8 )); then
  echo "invalid TASK_INDEX=$TASK_INDEX" >&2
  exit 2
fi

row=$((TASK_INDEX + 2))
anchor=$(sed -n "${row}p" "$OVA_E318_EXPERIMENT/anchor_manifest.tsv" | cut -f1)
if (( TASK_INDEX < 4 )); then
  bottom_k=2
  boost=2.5
  design_mask="$OVA_E318_EXPERIMENT/joint_design_surface.txt"
  mask_mode=surface
else
  bottom_k=3
  boost=2.0
  design_mask="$OVA_E318_EXPERIMENT/joint_design_hotspot.txt"
  mask_mode=hotspot_expanded
fi
run="$OVA_E318_EXPERIMENT/run_e318_$(printf '%02d' "$TASK_INDEX")"
joint="$run/joint_E318_k${bottom_k}"
mkdir -p "$joint" "$OVA_E318_EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$joint/af2diff_candidates.tsv" ]]; then
  echo "task $TASK_INDEX already complete"
  exit 0
fi
export PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$OVA_E318_EXPERIMENT/jaxcache_$TASK_INDEX"

state_args=()
while IFS=$'\t' read -r state_index context_index source_cif source_iptm interface target restraint weight; do
  [[ "$state_index" == "state_index" ]] && continue
  # Odd tasks deliberately increase the AD/second-interface target weight.
  if (( TASK_INDEX % 2 == 1 )) && [[ "$interface" == AD ]]; then
    state_weight=4.2
  elif (( TASK_INDEX % 2 == 1 )) && [[ "$interface" == AC ]]; then
    state_weight=2.6
  else
    state_weight=3.0
  fi
  state_args+=(--target-pdb "$target" --chains "${interface:0:1},${interface:1:1}" \
    --restraint-positions "$restraint" --state-weight "$state_weight")
done < "$OVA_E318_EXPERIMENT/state_manifest.tsv"
if (( ${#state_args[@]} != 64 )); then
  echo "expected 8 state argument groups" >&2
  exit 2
fi

echo "task=$TASK_INDEX anchor=$anchor bottom_k=$bottom_k boost=$boost mask=$mask_mode"
"$BINDCRAFT_MOD/env/bin/python" "$OVA_TETRA_ROOT/af3_pipeline/af2_multistate_tied_lowmut_design.py" \
  "${state_args[@]}" \
  --reference-fasta "$OVA_TETRA_ROOT/OVA_P3-13R_四聚体候选_AA.fasta" \
  --reference-name D0-P3-13R --base-fasta "$OVA_E318_EXPERIMENT/anchors.fasta" \
  --base-name "$anchor" --design-positions "$design_mask" \
  --surface-positions "$OVA_E318_EXPERIMENT/surface_positions.txt" \
  --params "$BINDCRAFT_MOD/params" --out "$joint" \
  --seed "$((192001 + TASK_INDEX))" --mutation-budget 24 \
  --min-surface-mutation-fraction 0.80 --w-budget 22 --w-sparse 1.0 \
  --w-buried-mutation 10 --w-fold 1.5 \
  --worst-state-gradient-boost "$boost" --worst-state-count "$bottom_k" \
  --logits-iters 14 --soft-iters 8 --hard-iters 4 --recycles 1
echo "task $TASK_INDEX complete"
