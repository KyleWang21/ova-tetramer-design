#!/usr/bin/env bash
# Combine 25-model AF3 evidence for every E168 strict-seed1 hit that survives
# the independent seed2 fail-fast gate.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
SOURCE="$PROJECT/experiments/e168_c4_073_reversion_ai_af3_seed1"
cd "$PROJECT"

experiments=(
  e206_c4_e168_rev002_multiseed
  e207_c4_e168_rev010_multiseed
  e208_c4_e168_rev003_multiseed
  e209_c4_e168_rev006_multiseed
  e211_c4_e168_rev001_multiseed
)
seed1_systems=(ova_body_c4_01 ova_body_c4_08 ova_body_c4_03 ova_body_c4_10 ova_body_c4_04)
seed1_tables=(
  seed1_strict_s0/seed1_strict_per_model.tsv
  seed1_strict_s0/seed1_strict_per_model.tsv
  seed1_strict_s2/seed1_strict_per_model.tsv
  seed1_strict_s2/seed1_strict_per_model.tsv
  seed1_strict_s3/seed1_strict_per_model.tsv
)

for index in "${!experiments[@]}"; do
  experiment="$PROJECT/experiments/${experiments[$index]}"
  while [[ ! -f "$experiment/RESEED_DECISION_COMPLETE" ]]; do
    sleep 30
  done
  candidate=$(awk -F '\t' 'NR==2 {print $1}' "$experiment/manifest.tsv")
  if [[ -f "$experiment/seed5_strict/seed1_strict_per_model.tsv" ]]; then
    "$PY" af3_pipeline/combine_strict_per_model_c4.py \
      --seed1-per-model "$SOURCE/${seed1_tables[$index]}" \
      --seed1-system "${seed1_systems[$index]}" \
      --multiseed-experiment "$experiment" --candidate-system "$candidate" \
      --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
      --reference-name D0-P3-13R --reference-cif references/1OVA.cif \
      --out "$experiment/full25_strict"
  else
    echo "${experiments[$index]} stopped at seed2; no 25-model combination"
  fi
done
touch "$SOURCE/MULTISEED_FOLLOWUPS_FINALIZED"
echo "all E168 multiseed follow-ups finalized"
