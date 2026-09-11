#!/usr/bin/env bash
# One tied sequence optimized simultaneously against four backbones / eight interfaces.
set -euo pipefail

TASK_INDEX="${1:-${TASK_INDEX:?TASK_INDEX is required}}"
PROJECT_ROOT="${OVA_TETRA_ROOT:-/root/400083/ova_p3_tetramer_design}"
BINDCRAFT_ROOT="${BINDCRAFT_MOD:-/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft}"
EXPERIMENT="${OVA_E280_EXPERIMENT:-$PROJECT_ROOT/experiments/e280_c4_eightstate_multibackbone_tied_design}"
cd "$PROJECT_ROOT"
if (( TASK_INDEX < 0 || TASK_INDEX >= 8 )); then echo "invalid TASK_INDEX=$TASK_INDEX" >&2; exit 2; fi

row=$((TASK_INDEX + 2))
anchor=$(sed -n "${row}p" "$EXPERIMENT/anchor_manifest.tsv" | cut -f1)
run="$EXPERIMENT/run_eightstate_$(printf '%02d' "$TASK_INDEX")"
joint="$run/joint_AC_AD_r$((TASK_INDEX % 2))"
mkdir -p "$joint" "$EXPERIMENT/jaxcache_$TASK_INDEX"
if [[ -s "$joint/af2diff_candidates.tsv" ]]; then echo "task $TASK_INDEX already complete"; exit 0; fi
export PYTHONUNBUFFERED=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export JAX_COMPILATION_CACHE_DIR="$EXPERIMENT/jaxcache_$TASK_INDEX"

state_args=()
while IFS=$'\t' read -r state_index context_index source_cif source_iptm interface target restraint weight; do
  [[ "$state_index" == "state_index" ]] && continue
  state_args+=(--target-pdb "$target" --chains "${interface:0:1},${interface:1:1}" \
    --restraint-positions "$restraint" --state-weight "$weight")
done < "$EXPERIMENT/state_manifest.tsv"
if (( ${#state_args[@]} != 64 )); then echo "expected 8 complete state argument groups" >&2; exit 2; fi

echo "task=$TASK_INDEX anchor=$anchor shared_tied_sequence=1 contexts=4 interface_states=8"
"$BINDCRAFT_ROOT/env/bin/python" "$PROJECT_ROOT/af3_pipeline/af2_multistate_tied_lowmut_design.py" \
  "${state_args[@]}" \
  --reference-fasta "$PROJECT_ROOT/OVA_P3-13R_四聚体候选_AA.fasta" \
  --reference-name D0-P3-13R --base-fasta "$EXPERIMENT/anchors.fasta" --base-name "$anchor" \
  --design-positions "$EXPERIMENT/joint_design.txt" \
  --surface-positions "$EXPERIMENT/surface_positions.txt" \
  --params "$BINDCRAFT_ROOT/params" --out "$joint" \
  --seed "$((171201 + TASK_INDEX))" --mutation-budget 20 \
  --min-surface-mutation-fraction 0.80 --w-budget 34 --w-sparse 2.2 \
  --w-buried-mutation 16 --w-fold 1.0 --worst-state-gradient-boost 2.5 \
  --logits-iters 16 --soft-iters 8 --hard-iters 4 --recycles 1
echo "task $TASK_INDEX complete"
