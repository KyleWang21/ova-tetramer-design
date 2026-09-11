#!/usr/bin/env bash
# Complete the Fire-resident core of every staged generic Accepted_AF3 v1.1
# screen. All GPU stages run first while the caller owns the priority lock;
# the lock is then released before CPU-only Rosetta. OpenDDE is no longer a
# required or automatically submitted stage. Every GPU submit remains gated
# on free>10.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
ROOT="$PROJECT/experiments/current_pending_v1_screen"
LOCK="$PROJECT/experiments/VOLC_PRIORITY_BATCH.lock"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_FILTER=/root/400083/antibody-design/pipeline_stage6_maskcap_bindcraft/env/bin/python
cd "$PROJECT"

successful_jobs() {
  local path="$1"
  [[ -s "$path" ]] || { echo 0; return; }
  awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++){if($i=="mode")m=i;if($i=="status")s=i;if($i=="job_id")j=i}next} $m=="submit"&&$s=="ok"&&length($j){n++} END{print n+0}' "$path"
}

for screen in "$ROOT"/*; do
  [[ -f "$screen/INPUTS_PREPARED" && -f "$screen/v1_job_base.tsv" ]] || continue
  [[ ! -f "$screen/SCREEN_V1_COMPLETE" ]] || continue
  base=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="job_base")c=i;next} NR==2{print $c}' "$screen/v1_job_base.tsv")
  candidate=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="candidate")c=i;next} NR==2{print $c}' "$screen/final_candidates.tsv")

  stoich="$screen/stoichiometry_af3"
  if [[ ! -f "$screen/STOICHIOMETRY_COMPLETE" ]]; then
    while [[ ! -f "$stoich/VOLC_AF3_SHARDS_SUBMITTED" ]]; do
      "$PY_AF3" af3_pipeline/submit_af3_volc_if_capacity.py \
        --experiment "$stoich" --job-number-start "$base" --priority
      [[ -f "$stoich/VOLC_AF3_SHARDS_SUBMITTED" ]] && break
      sleep 30
    done
    while true; do
      complete=$(find "$stoich" -path '*/seed-*/*_summary_confidences.json' 2>/dev/null | wc -l || true)
      [[ "$complete" -eq 75 ]] && break
      sleep 30
    done
    "$PY_AF3" af3_pipeline/score_outputs.py --experiment "$stoich"
    "$PY_AF3" af3_pipeline/summarize_final20_stoichiometry.py \
      --manifest "$stoich/manifest.tsv" --stoich-models "$stoich/model_scores.tsv" \
      --c4-models "$screen/af3_per_model_recomputed.tsv" \
      --reference references/1OVA.cif --out "$screen/stoichiometry_summary"
    touch "$screen/STOICHIOMETRY_COMPLETE"
  fi

  protenix="$screen/protenix_multiseed"
  if [[ ! -f "$screen/PROTENIX_COMPLETE" ]]; then
    expected=$(awk -F '\t' 'NR>1{seen[$1]=1} END{for(x in seen)n++;print n+0}' "$protenix/volc_manifest.tsv")
    while [[ "$(successful_jobs "$protenix/volc_jobs.tsv")" -lt "$expected" ]]; do
      "$PY_AF3" af3_pipeline/submit_volc_protenix_multiseed.py \
        --experiment "$protenix" --job-number-start "$((base + 8))" --submit
      [[ "$(successful_jobs "$protenix/volc_jobs.tsv")" -ge "$expected" ]] && break
      sleep 30
    done
    while true; do
      complete=$(find "$protenix/out" -path '*/seed_*/*_summary_confidence_sample_0.json' 2>/dev/null | wc -l || true)
      [[ "$complete" -eq 5 ]] && break
      sleep 30
    done
    export PYTHONPATH=/root/400083/alphafold3/repo/src${PYTHONPATH:+:$PYTHONPATH}
    "$PY_AF3" af3_pipeline/summarize_final20_protenix_multiseed.py \
      --final "$screen/final_candidates.tsv" --af3-models "$screen/af3_per_model_recomputed.tsv" \
      --protenix-root "$protenix/out" --reference references/1OVA.cif \
      --out "$screen/protenix_summary_preprolif"
    "$PY_AF3" af3_pipeline/prepare_crossmodel_prolif.py \
      --final "$screen/final_candidates.tsv" --af3-models "$screen/af3_per_model_recomputed.tsv" \
      --protenix-summary "$screen/protenix_summary_preprolif/protenix_candidate_summary.tsv" \
      --out "$screen/prolif_crossmodel/structure_manifest.tsv"
    "$PY_FILTER" af3_pipeline/analyze_prolif_interfaces.py \
      --structure-manifest "$screen/prolif_crossmodel/structure_manifest.tsv" \
      --final-candidates "$screen/final_candidates.tsv" --out "$screen/prolif_crossmodel"
    "$PY_AF3" af3_pipeline/summarize_crossmodel_prolif.py \
      --unique "$screen/prolif_crossmodel/unique_residue_interactions.tsv" \
      --out "$screen/prolif_crossmodel/crossmodel_nonvdw.tsv"
    "$PY_AF3" af3_pipeline/summarize_final20_protenix_multiseed.py \
      --final "$screen/final_candidates.tsv" --af3-models "$screen/af3_per_model_recomputed.tsv" \
      --protenix-root "$protenix/out" --reference references/1OVA.cif \
      --shared-nonvdw "$screen/prolif_crossmodel/crossmodel_nonvdw.tsv" \
      --out "$screen/protenix_summary"
    touch "$screen/PROTENIX_COMPLETE"
  fi

done

# No Fire GPU work remains in this invocation.  Do not make AF3 design wait
# for the much slower CPU FastRelax stage.  The parent watcher's ownership-
# checked cleanup remains safe if another AF3 process acquires the lock next.
if [[ -f "$LOCK" ]] && awk -F '\t' 'NR==1 && $2=="generic_pending_v1"{ok=1} END{exit !ok}' "$LOCK"; then
  rm -f "$LOCK"
fi

for screen in "$ROOT"/*; do
  [[ -f "$screen/INPUTS_PREPARED" && -f "$screen/v1_job_base.tsv" ]] || continue
  [[ ! -f "$screen/SCREEN_V1_COMPLETE" ]] || continue
  candidate=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="candidate")c=i;next} NR==2{print $c}' "$screen/final_candidates.tsv")
  if [[ ! -f "$screen/ROSETTA_COMPLETE" ]]; then
    manifest_values=$(awk -F '\t' 'NR==2 {
      if (NF >= 8) printf "%s|%s|%s|%s|%s|%s|%s|%s", $1,$2,$3,$4,$5,$6,$7,$8;
      else printf "%s|%s|%s|%s|%s||||%s", $1,$2,$3,$4,$4,$5;
    }' "$screen/rosetta_manifest.tsv")
    IFS='|' read -r ros_candidate source_model iptm effective_edges design_edges background_edges design_selection cif_path <<<"$manifest_values"
    OVA_ROSETTA_CHUNKS=32 OVA_ROSETTA_TOTAL=200 \
      bash af3_pipeline/run_rosetta_c4_parallel.sh \
      "$cif_path" "$ros_candidate" "$effective_edges" "$design_edges" "$source_model" \
      "$screen/rosetta_chunks" "$screen/rosetta"
    touch "$screen/ROSETTA_COMPLETE"
  fi

  touch "$screen/V1_CORE_COMPLETE"
  echo "$candidate Fire v1.1 core stages complete; OpenDDE is optional diagnostic"
done

"$PY_AF3" af3_pipeline/update_accepted_af3_registry.py \
  --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
