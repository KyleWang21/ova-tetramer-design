#!/usr/bin/env bash
# Score the ten-seed OVA-C4-073 OpenDDE-v1 run after Zhongwei recovery.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
EXP="$PROJECT/experiments/e214_c4_073_opendde_v1"
AF3_EXP="$PROJECT/experiments/e165_c4_iptm90_073_multiseed"
OUT="$EXP/summary"
PY=/root/400083/alphafold3/repo/.venv/bin/python
DOCKQ="$PROJECT/experiments/e153_dockq/env/bin/DockQ"
cd "$PROJECT"

while [[ ! -f "$EXP/REMOTE_ARCHIVE_RECOVERED" ]]; do
  sleep 30
done
count=$(find "$EXP/out" -path '*/seed_*/predictions/*_sample_0.cif' | wc -l)
[[ "$count" -eq 10 ]] || { echo "expected 10 recovered OpenDDE CIFs, found $count" >&2; exit 2; }

export PYTHONPATH=/root/400083/alphafold3/repo/src${PYTHONPATH:+:$PYTHONPATH}
"$PY" af3_pipeline/summarize_opendde_c4.py \
  --manifest "$EXP/manifest.tsv" --out-root "$EXP/out" \
  --final "$AF3_EXP/final25_screen_inputs/final_candidates.tsv" \
  --af3-models "$AF3_EXP/final25_screen/af3_per_model_recomputed.tsv" \
  --reference references/1OVA.cif --dockq-bin "$DOCKQ" --out "$OUT"
touch "$EXP/OPEN_DDE_SUMMARY_COMPLETE"
echo "E214 ten-seed OpenDDE scoring complete"
