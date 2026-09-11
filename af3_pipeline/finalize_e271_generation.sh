#!/usr/bin/env bash
# Aggregate E271 exact scores and freeze E272 for template-free AF3.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
E271="$PROJECT/experiments/e271_c4_multibackbone_exact_mpnn_score"
E272="$PROJECT/experiments/e272_c4_multibackbone_exact_mpnn_af3_seed1"
SOURCE="$PROJECT/experiments/e246_c4_ai11_crossmodel_tied_mpnn_full/candidate_pool_joint_e247_novel/hamming_ai_ranking_current_weighted.tsv"
PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$E272/CANDIDATES_FROZEN" ]]; then echo "E272 already frozen"; exit 0; fi
while [[ ! -f "$E271/REMOTE_ARCHIVE_RECOVERED" ]]; do sleep 30; done
"$PY" af3_pipeline/rank_e271_multibackbone_exact_mpnn.py \
  --project "$PROJECT" --experiment "$E271" --source-ranking "$SOURCE" --out "$E271" \
  --top 32 --minimum-distance 3 \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --exclude-glob 'deliverables/**/*.fasta'
"$PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$E271/multibackbone_exact_mpnn_shortlist.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$E272" --top 32 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$E272/CANDIDATES_FROZEN"
echo "E271 multi-backbone exact MPNN consensus frozen as E272"
