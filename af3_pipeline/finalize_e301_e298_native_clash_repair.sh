#!/usr/bin/env bash
# Freeze R340-E341 local repairs after E298 fail-fast validation.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e298_e297_top4_af3_multiseed"
WORK="$PROJECT/experiments/e301_c4_e298_nativeclash_repair_generation"
OUT="$PROJECT/experiments/e301_c4_e298_nativeclash_repair_af3_seed1"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$OUT/CANDIDATES_FROZEN" ]]; then echo "E301 already frozen"; exit 0; fi
while [[ ! -f "$SOURCE/MULTISEED_VALIDATION_COMPLETE" ]]; do sleep 30; done
mkdir -p "$WORK"
"$PY" af3_pipeline/prepare_e301_e298_native_clash_repair.py \
  --source "$SOURCE" --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --out "$WORK/generation_library.tsv"
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$WORK/generation_library.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$OUT" --top 16 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$OUT/CANDIDATES_FROZEN"
echo "E301 E298 native-clash repair candidates frozen"
