#!/usr/bin/env bash
# Atomically refresh the ranked Accepted_AF3 delivery only when its content changes.
set -euo pipefail

PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/alphafold3/repo/.venv/bin/python
REGISTRY="$PROJECT/experiments/current_accepted_af3_registry/accepted_af3_registry.tsv"
DATE_TAG=$(TZ=Asia/Shanghai date +%Y%m%d)
NAME="current_accepted_af3_ranked_${DATE_TAG}"
DEST="$PROJECT/deliverables/$NAME"
ARCHIVE="$PROJECT/deliverables/$NAME.tgz"
BACKUP_ROOT="$PROJECT/deliverables/archive"
cd "$PROJECT"

"$PY" af3_pipeline/update_accepted_af3_registry.py \
  --project "$PROJECT" --out "$PROJECT/experiments/current_accepted_af3_registry"
[[ -s "$REGISTRY" ]] || { echo "Accepted_AF3 registry is empty" >&2; exit 2; }

TMP_ROOT=$(mktemp -d "$PROJECT/deliverables/.accepted-refresh.XXXXXX")
cleanup() {
  if [[ "$TMP_ROOT" == "$PROJECT/deliverables/.accepted-refresh."* && -d "$TMP_ROOT" ]]; then
    rm -rf -- "$TMP_ROOT"
  fi
}
trap cleanup EXIT
PAYLOAD="$TMP_ROOT/$NAME"
"$PY" af3_pipeline/export_ranked_accepted_af3.py --registry "$REGISTRY" --out "$PAYLOAD"

expected=$(tail -n +2 "$REGISTRY" | wc -l)
table_rows=$(tail -n +2 "$PAYLOAD/unified_candidates.tsv" | wc -l)
fasta_rows=$(grep -c '^>' "$PAYLOAD/all_sequences.fasta")
cif_rows=$(find "$PAYLOAD/cif" -maxdepth 1 -type f -name 'RANK*.cif' | wc -l)
[[ "$expected" -gt 0 && "$table_rows" -eq "$expected" && \
   "$fasta_rows" -eq "$expected" && "$cif_rows" -eq "$expected" ]] || {
  echo "delivery audit failed: registry=$expected table=$table_rows fasta=$fasta_rows cif=$cif_rows" >&2
  exit 3
}
printf 'registry_rows\t%s\ntable_rows\t%s\nfasta_records\t%s\ncif_files\t%s\n' \
  "$expected" "$table_rows" "$fasta_rows" "$cif_rows" > "$PAYLOAD/delivery_audit.tsv"
tar -C "$TMP_ROOT" -czf "$TMP_ROOT/$NAME.tgz" "$NAME"

if [[ -s "$DEST/unified_candidates.tsv" ]] && \
   cmp -s "$DEST/unified_candidates.tsv" "$PAYLOAD/unified_candidates.tsv" && \
   [[ -s "$DEST/SUMMARY.md" && -s "$DEST/delivery_audit.tsv" && -s "$ARCHIVE" ]]; then
  echo "delivery unchanged: $ARCHIVE"
  exit 0
fi

timestamp=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$BACKUP_ROOT"
if [[ -e "$DEST" || -e "$ARCHIVE" ]]; then
  backup="$BACKUP_ROOT/${NAME}_${timestamp}"
  mkdir -p "$backup"
  [[ -e "$DEST" ]] && mv "$DEST" "$backup/"
  [[ -e "$ARCHIVE" ]] && mv "$ARCHIVE" "$backup/"
  echo "previous delivery preserved under $backup"
fi
mv "$PAYLOAD" "$DEST"
mv "$TMP_ROOT/$NAME.tgz" "$ARCHIVE"
echo "refreshed $ARCHIVE with $expected ranked Accepted_AF3 candidates"
