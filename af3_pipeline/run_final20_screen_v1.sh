#!/usr/bin/env bash
# Idempotent collector for the OVA noncovalent-C4 screening plan v1.1.
# Remote AF3/Protenix/OpenDDE jobs are submitted separately; this driver scores
# every available result. OpenDDE, when present, is retained as optional
# diagnostic evidence and is not a required stage.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
FINAL=experiments/e144_final20_iptm80_c4/final_candidates.tsv
AF3_MODELS=experiments/e151_final20_screen_v1/geometry/af3_per_model_recomputed.tsv
AF3_SUMMARY=experiments/e151_final20_screen_v1/geometry/af3_candidate_recomputed_summary.tsv
STOICH=experiments/e149_final20_stoichiometry_v1
PTX31=experiments/e148_final20_protenix_v2_all_seed31
PTXNEW=experiments/e150_final20_protenix_multiseed_v1
SCREEN=experiments/e151_final20_screen_v1
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_PROLIF=experiments/e124_prolif_interface/env/bin/python

# Exact geometry is always recomputed from coordinates; no model-native clash flag
# is accepted as a substitute for the frozen 2.4-A definition.
"$PY_AF3" af3_pipeline/screen_final20_af3_structures.py \
  --final "$FINAL" --per-model experiments/e144_final20_iptm80_c4/per_model_scores.tsv \
  --reference references/1OVA.cif --out "$SCREEN/geometry"

stoich_args=()
stoich_count=$("$PY_AF3" - "$STOICH" <<'PY'
import csv, pathlib, sys
root = pathlib.Path(sys.argv[1])
rows = csv.DictReader((root / "manifest.tsv").open(), delimiter="\t")
print(sum(len(list((root / f"out_s{row['shard']}" / row["name"]).glob(
    "seed-*/*_summary_confidences.json"
))) for row in rows))
PY
)
if [ "$stoich_count" -eq 1500 ]; then
  "$PY_AF3" af3_pipeline/score_outputs.py --experiment "$STOICH"
  "$PY_AF3" af3_pipeline/summarize_final20_stoichiometry.py \
    --manifest "$STOICH/manifest.tsv" --stoich-models "$STOICH/model_scores.tsv" \
    --c4-models "$AF3_MODELS" --reference references/1OVA.cif --out "$SCREEN/stoichiometry"
  stoich_args=(--stoichiometry-summary "$SCREEN/stoichiometry/stoichiometry_candidate_summary.tsv")
else
  echo "stoichiometry incomplete: $stoich_count/1500; stage remains MISSING"
fi

protenix_args=()
ptx_count=$(find "$PTX31/out" "$PTX31/out_volc" "$PTXNEW/out" \
  -name '*_summary_confidence_sample_0.json' 2>/dev/null | wc -l)
if [ "$ptx_count" -eq 100 ]; then
  mkdir -p "$SCREEN/protenix" "$SCREEN/prolif_crossmodel"
  "$PY_AF3" af3_pipeline/summarize_final20_protenix_multiseed.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-root "$PTX31/out" --protenix-root "$PTX31/out_volc" \
    --protenix-root "$PTXNEW/out" --reference references/1OVA.cif --out "$SCREEN/protenix"
  "$PY_AF3" af3_pipeline/prepare_crossmodel_prolif.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-summary "$SCREEN/protenix/protenix_candidate_summary.tsv" \
    --out "$SCREEN/prolif_crossmodel/structure_manifest.tsv"
  "$PY_PROLIF" af3_pipeline/analyze_prolif_interfaces.py \
    --structure-manifest "$SCREEN/prolif_crossmodel/structure_manifest.tsv" \
    --final-candidates "$FINAL" --out "$SCREEN/prolif_crossmodel/results"
  "$PY_AF3" af3_pipeline/summarize_crossmodel_prolif.py \
    --unique "$SCREEN/prolif_crossmodel/results/unique_residue_interactions.tsv" \
    --out "$SCREEN/prolif_crossmodel/shared_nonvdw.tsv"
  "$PY_AF3" af3_pipeline/summarize_final20_protenix_multiseed.py \
    --final "$FINAL" --af3-models "$AF3_MODELS" \
    --protenix-root "$PTX31/out" --protenix-root "$PTX31/out_volc" \
    --protenix-root "$PTXNEW/out" --reference references/1OVA.cif \
    --shared-nonvdw "$SCREEN/prolif_crossmodel/shared_nonvdw.tsv" --out "$SCREEN/protenix"
  protenix_args=(--protenix-summary "$SCREEN/protenix/protenix_candidate_summary.tsv")
else
  echo "Protenix-v2 incomplete: $ptx_count/100; stage remains MISSING"
fi

rosetta_args=()
if [ -f "$SCREEN/rosetta_manifest.tsv" ]; then
  "$PY_AF3" af3_pipeline/collect_rosetta_c4_summaries.py \
    --manifest "$SCREEN/rosetta_manifest.tsv" --root "$SCREEN/rosetta_diagnostic" \
    --out "$SCREEN/rosetta_diagnostic_summary.tsv"
  rosetta_args=(--rosetta-summary "$SCREEN/rosetta_diagnostic_summary.tsv")
fi

opendde_args=()
if [ -f "$SCREEN/opendde/opendde_candidate_summary.tsv" ]; then
  opendde_args=(--opendde-summary "$SCREEN/opendde/opendde_candidate_summary.tsv")
fi
opendde_prolif_args=()
if [ -f "$SCREEN/opendde/opendde_prolif_summary.tsv" ]; then
  opendde_prolif_args=(--opendde-prolif-summary "$SCREEN/opendde/opendde_prolif_summary.tsv")
fi

"$PY_PROLIF" af3_pipeline/run_final20_failclosed_filter.py \
  --final "$FINAL" --af3-summary "$AF3_SUMMARY" --af3-models "$AF3_MODELS" \
  --prolif-unique experiments/e145_final20_prolif/results/unique_residue_interactions.tsv \
  --reference references/1OVA.cif "${protenix_args[@]}" "${stoich_args[@]}" \
  "${rosetta_args[@]}" "${opendde_args[@]}" "${opendde_prolif_args[@]}" --out "$SCREEN/unified"

report_args=()
if [ -f "$SCREEN/opendde/opendde_calibration_summary.tsv" ]; then
  report_args=(--opendde-calibration "$SCREEN/opendde/opendde_calibration_summary.tsv")
fi
"$PY_AF3" af3_pipeline/make_final20_screen_report.py \
  --unified "$SCREEN/unified/final20_failclosed_screen.tsv" \
  "${report_args[@]}" --out "$SCREEN/report"
