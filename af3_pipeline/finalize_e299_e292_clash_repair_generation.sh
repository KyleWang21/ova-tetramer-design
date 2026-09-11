#!/usr/bin/env bash
# Freeze direct repairs of the E292 high-AF3/low-clash deterministic endpoint.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
SOURCE="$PROJECT/experiments/e292_c4_e283_deterministic_final_af3_seed1"
OUT="$PROJECT/experiments/e299_c4_e292_clashrepair_af3_seed1"
WORK="$PROJECT/experiments/e299_c4_e292_clashrepair_generation"
AF3_PY=/root/400083/alphafold3/repo/.venv/bin/python
BIO_PY=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"
if [[ -f "$OUT/CANDIDATES_FROZEN" ]]; then echo "E299 already frozen"; exit 0; fi
while [[ ! -f "$SOURCE/SEED1_STRICT_RESULTS_COMPLETE" ]]; do sleep 30; done

mkdir -p "$WORK"
"$AF3_PY" af3_pipeline/audit_af3_interchain_clashes.py \
  --experiment "$SOURCE" --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --out "$WORK/clash_audit"
"$BIO_PY" af3_pipeline/prepare_e289_clash_repair_library.py \
  --experiment "$SOURCE" \
  --strict-summary "$SOURCE/seed1_strict_final/seed1_strict_summary.tsv" \
  --clash-recurrence "$WORK/clash_audit/clash_residue_pair_recurrence.tsv" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --accepted-fasta deliverables/current_accepted_af3_ranked_20260903/all_sequences.fasta \
  --surface-positions experiments/e283_c4_chainmapped_eightstate_tied_design/surface_positions.txt \
  --exclude-glob 'experiments/e*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/manifest.tsv' \
  --exclude-glob 'experiments/e*/*/*/manifest.tsv' \
  --out "$WORK/generation_library.tsv" --max-parents 2 --max-edits 2 \
  --max-mutations 20 --min-surface-fraction 0.80 --name-prefix E299-REPAIR

# This is a complete six-member mechanistic matrix, not a large library: test
# every unique repair instead of allowing a surrogate rank to discard an
# epistatic combination before AF3.
"$BIO_PY" af3_pipeline/prepare_ova_body_af3.py \
  --fasta "$WORK/generation_library.fasta" \
  --ova-msa experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m \
  --out "$OUT" --top 8 --seeds 1 --seed-start 1 --seed-offset 0 \
  --shards 8 --n-chains 4
touch "$OUT/CANDIDATES_FROZEN"
echo "E299 E292 AF3 clash-directed repair candidates frozen"
