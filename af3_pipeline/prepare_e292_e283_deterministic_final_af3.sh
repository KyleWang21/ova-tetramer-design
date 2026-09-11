#!/usr/bin/env bash
# Freeze exact deterministic E283 endpoints for a dedicated AF3 quota.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
OUT="$PROJECT/experiments/e292_c4_e283_deterministic_final_af3_seed1"
POOL="$PROJECT/experiments/e283_c4_chainmapped_eightstate_tied_design/candidate_pool_deterministic_final"
cd "$PROJECT"
if [[ -f "$OUT/CANDIDATES_FROZEN" ]]; then echo "E292 already frozen"; exit 0; fi
"$PY" af3_pipeline/prepare_e292_e283_deterministic_final_af3.py
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$POOL/deterministic_final_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$OUT" --top 8 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$OUT/CANDIDATES_FROZEN"
echo "E292 deterministic E283 endpoints frozen"
