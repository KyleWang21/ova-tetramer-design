#!/usr/bin/env bash
# Historical OpenDDE/DockQ/ProLIF collector for the frozen 20-candidate screen.
# It is not part of the current v1.1 pipeline and is retained only for audit replay.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
EXP=${EXP:-experiments/e154_final20_opendde_msa_v1}
SCREEN=${SCREEN:-experiments/e151_final20_screen_v1}
FINAL=experiments/e144_final20_iptm80_c4/final_candidates.tsv
AF3_MODELS="$SCREEN/geometry/af3_per_model_recomputed.tsv"
PY_AF3=/root/400083/alphafold3/repo/.venv/bin/python
PY_PROLIF=experiments/e124_prolif_interface/env/bin/python
DOCKQ=experiments/e153_dockq/env/bin/DockQ

# Validate each manifest system and each frozen seed, not only a global file count.
"$PY_AF3" - "$EXP" <<'PY'
import csv, json, pathlib, re, sys

root = pathlib.Path(sys.argv[1])
manifest = list(csv.DictReader((root / "manifest.tsv").open(), delimiter="\t"))
if len(manifest) != 22:
    raise SystemExit(f"fail_closed: expected 22 OpenDDE systems, found {len(manifest)}")
total = 0
for item in manifest:
    payload = json.loads((root / item["input_json"]).read_text())[0]
    job = payload["name"]
    expected = set(map(int, payload["modelSeeds"]))
    found = {}
    for summary in (root / "out").glob(
        f"**/{job}/seed_*/predictions/*_summary_confidence_sample_0.json"
    ):
        match = re.search(r"seed_(\d+)", str(summary))
        if not match:
            continue
        seed = int(match.group(1))
        cif = summary.with_name(summary.name.replace(
            "_summary_confidence_sample_0.json", "_sample_0.cif"
        ))
        if not cif.exists():
            raise SystemExit(f"fail_closed: missing CIF for {summary}")
        if seed in found and found[seed] != summary:
            raise SystemExit(f"fail_closed: duplicate output for {job} seed {seed}")
        found[seed] = summary
    if set(found) != expected:
        raise SystemExit(
            f"fail_closed: {item['name']} expected seeds {sorted(expected)}, found {sorted(found)}"
        )
    total += len(found)
if total != 220:
    raise SystemExit(f"fail_closed: expected 220 models, found {total}")
print(f"validated OpenDDE completeness: {len(manifest)} systems, {total} models")
PY

mkdir -p "$SCREEN/opendde"
"$PY_AF3" af3_pipeline/summarize_opendde_c4.py \
  --manifest "$EXP/manifest.tsv" --out-root "$EXP/out" --final "$FINAL" \
  --af3-models "$AF3_MODELS" --reference references/1OVA.cif \
  --dockq-bin "$DOCKQ" --out "$SCREEN/opendde"

"$PY_AF3" af3_pipeline/prepare_opendde_prolif.py \
  --final "$FINAL" --opendde-summary "$SCREEN/opendde/opendde_candidate_summary.tsv" \
  --out "$SCREEN/opendde/prolif_structure_manifest.tsv"

representatives=$(($(wc -l < "$SCREEN/opendde/prolif_structure_manifest.tsv") - 1))
opendde_unique=()
if [ "$representatives" -gt 0 ]; then
  "$PY_PROLIF" af3_pipeline/analyze_prolif_interfaces.py \
    --structure-manifest "$SCREEN/opendde/prolif_structure_manifest.tsv" \
    --final-candidates "$FINAL" --out "$SCREEN/opendde/prolif_results"
  opendde_unique=(--opendde-unique "$SCREEN/opendde/prolif_results/unique_residue_interactions.tsv")
fi

"$PY_AF3" af3_pipeline/summarize_opendde_prolif.py \
  --final "$FINAL" --opendde-summary "$SCREEN/opendde/opendde_candidate_summary.tsv" \
  --opendde-per-seed "$SCREEN/opendde/opendde_per_seed_recomputed.tsv" \
  "${opendde_unique[@]}" \
  --crossmodel-unique "$SCREEN/prolif_crossmodel/results/unique_residue_interactions.tsv" \
  --out "$SCREEN/opendde/opendde_prolif_summary.tsv"

bash af3_pipeline/run_final20_screen_v1.sh
echo "finalized OpenDDE/DockQ/ProLIF and regenerated $SCREEN/report"
