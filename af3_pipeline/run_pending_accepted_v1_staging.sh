#!/usr/bin/env bash
# Prepare all applicable v1.1 core stages and run the CPU AF3/ProLIF pass.
# OpenDDE is no longer prepared or submitted here; existing results remain
# available as optional diagnostics.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
ROOT="$PROJECT/experiments/current_pending_v1_screen"
REGISTRY="$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
PY_PROLIF="$PROJECT/experiments/e124_prolif_interface/env/bin/python"
BASE_MSA="$PROJECT/experiments/e202_c4_ai11v2_ai_af3_seed1/msas/ova_body_c4_01.a3m"
cd "$PROJECT"

"$PY_AF3" af3_pipeline/stage_pending_accepted_v1.py \
  --project "$PROJECT" --registry "$REGISTRY" --out "$ROOT" \
  --reference-fasta OVA_P3-13R_四聚体候选_AA.fasta \
  --reference-name D0-P3-13R

for screen in "$ROOT"/*; do
  [[ -f "$screen/STAGING_READY" ]] || continue
  relative=${screen#"$PROJECT/"}
  if [[ ! -f "$screen/INPUTS_PREPARED" ]]; then
    "$PY_AF3" af3_pipeline/prepare_protenix_candidates.py \
      --validation "$screen/final_candidates.tsv" --base-a3m "$BASE_MSA" \
      --out "$screen/protenix_base" --passing-only
    "$PY_AF3" af3_pipeline/prepare_final20_protenix_multiseed.py \
      --final "$screen/final_candidates.tsv" --input-dir "$screen/protenix_base/input" \
      --out "$screen/protenix_multiseed" --seeds 11,21,31,41,51 --shards 5 \
      --expected-candidates 1 \
      --msa-path-prefix "$PROJECT/$relative/protenix_multiseed/msa"
    "$PY_AF3" af3_pipeline/prepare_final20_stoichiometry.py \
      --final "$screen/final_candidates.tsv" \
      --per-model "$screen/af3_per_model_recomputed.tsv" \
      --out "$screen/stoichiometry_af3" --stoichiometries 1,2,3,5,6 \
      --seeds 1,2,3 --shards 8
    "$PY_AF3" af3_pipeline/prepare_rosetta_c4_manifest.py \
      --models "$screen/af3_per_model_recomputed.tsv" \
      --out "$screen/rosetta_manifest.tsv"
    touch "$screen/INPUTS_PREPARED"
  fi
  if [[ ! -f "$screen/INITIAL_AF3_PROLIF_COMPLETE" ]]; then
    "$PY_PROLIF" af3_pipeline/analyze_prolif_interfaces.py \
      --structure-manifest "$screen/af3_prolif_structure_manifest.tsv" \
      --final-candidates "$screen/final_candidates.tsv" \
      --out "$screen/prolif_af3"
    "$PY_FILTER" af3_pipeline/run_final20_failclosed_filter.py \
      --final "$screen/final_candidates.tsv" \
      --af3-summary "$screen/af3_candidate_recomputed_summary.tsv" \
      --af3-models "$screen/af3_per_model_recomputed.tsv" \
      --prolif-unique "$screen/prolif_af3/unique_residue_interactions.tsv" \
      --reference references/1OVA.cif --out "$screen/unified_initial"
    touch "$screen/INITIAL_AF3_PROLIF_COMPLETE"
  fi
done

"$PY_AF3" af3_pipeline/update_accepted_af3_registry.py \
  --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
"$PY_AF3" af3_pipeline/allocate_pending_v1_job_bases.py \
  --root "$ROOT" --start 2000 --block 32 \
  --reserved 2200-2263 --reserved 3000-3063
