#!/usr/bin/env bash
# Capacity-gated AF3 chain for the frozen v7/AI11-v1/v8 batches with 25-model promotion.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/alphafold3/repo/.venv/bin/python
cd "$PROJECT"

while [[ ! -f experiments/e195_c4_e177_02_multiseed/RESEED_DECISION_COMPLETE || \
        ! -f experiments/e196_c4_e177_23_multiseed/RESEED_DECISION_COMPLETE ]]; do
  sleep 30
done

# The newest cross-validated epistasis/Hamming ensemble has the highest current
# expected information gain and therefore runs before the older frozen waves.
while [[ ! -f experiments/e252_c4_epistasis_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e254_c4_epistasis_tied_joint_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e257_c4_recovered_trajectory_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e236_c4_crossmodel_active_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e244_c4_crossensemble_active_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e248_c4_fullcondition_mpnn_joint_ai_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e262_c4_e243_crossbackbone_continuation_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e264_c4_multihit11_stoichiometry_v1/STOICHIOMETRY_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e266_c4_crossbackbone_consensus_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e269_c4_e251_softbasin_stabilization_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e272_c4_multibackbone_exact_mpnn_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e275_c4_e261_softbasin_stabilization_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e284_c4_chainmapped_eightstate_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e287_c4_chainmapped_bottomk_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done
while [[ ! -f experiments/e278_c4_e268_hardhit_crossbackbone_af3_seed1/SEED1_STRICT_SCREEN_COMPLETE ]]; do
  sleep 30
done

run_one() {
  local experiment="$1" seed1_base="$2" multiseed="$3" multiseed_base="$4"
  local expected_models
  expected_models=$(( $(tail -n +2 "$experiment/manifest.tsv" | wc -l) * 5 ))
  while [[ ! -f "$experiment/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
    "$PY" af3_pipeline/submit_af3_volc_if_capacity.py \
      --experiment "$experiment" --job-number-start "$seed1_base"
    [[ -f "$experiment/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
    sleep 30
  done
  while true; do
    complete=$(find "$experiment" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
    [[ "$complete" -eq "$expected_models" ]] && break
    sleep 30
  done
  "$PY" af3_pipeline/summarize_iptm90_seed1_strict.py \
    --experiment "$experiment" \
    --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
    --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
    --out "$experiment/seed1_strict_final"
  touch "$experiment/SEED1_STRICT_RESULTS_COMPLETE"
  bash af3_pipeline/run_seed1_hits_multiseed_volc.sh \
    "$experiment" "$multiseed" "$multiseed_base" 4
  touch "$experiment/SEED1_STRICT_SCREEN_COMPLETE"
}

run_one experiments/e194_c4_v7_ai_af3_seed1 626 \
  experiments/e238_e194_top4_af3_multiseed 1
run_one experiments/e199_c4_ai11v1_ai_af3_seed1 659 \
  experiments/e239_e199_top4_af3_multiseed 17
run_one experiments/e197_c4_v8_ai_af3_seed1 642 \
  experiments/e240_e197_top4_af3_multiseed 33
echo "E194/E199/E197 AF3 and promoted seeds2-5 chain complete"
