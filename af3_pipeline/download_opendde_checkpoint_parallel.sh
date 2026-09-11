#!/usr/bin/env bash
# Parallel, resumable downloader for the released OpenDDE-v1 checkpoint.
# The frozen size and SHA-256 are from OpenDDE 1.1.0 model_manifest.json.
set -euo pipefail

ROOT="${OPENDDE_ROOT_DIR:?set OPENDDE_ROOT_DIR}"
URL="${OPENDDE_CHECKPOINT_URL:-https://huggingface.co/aurekaresearch/OpenDDE/resolve/eddd563ce96571f784012edd8f045181c8f8627d/opendde.pt}"
SIZE=2625249069
SHA256=7b826620390afad877ee2babc6a4d0df81b94d3a0be030959853d6a7da0807cc
N_CHUNKS="${N_CHUNKS:-64}"
PARALLEL="${PARALLEL:-32}"
CHUNK_DIR="$ROOT/checkpoint/opendde.pt.parts"
TARGET="$ROOT/checkpoint/opendde.pt"
mkdir -p "$CHUNK_DIR" "$(dirname "$TARGET")"

download_one() {
  index="$1"
  chunk_size=$(( (SIZE + N_CHUNKS - 1) / N_CHUNKS ))
  start=$(( index * chunk_size ))
  end=$(( start + chunk_size - 1 ))
  [ "$end" -lt "$SIZE" ] || end=$(( SIZE - 1 ))
  expected=$(( end - start + 1 ))
  output=$(printf '%s/%03d.part' "$CHUNK_DIR" "$index")
  if [ -f "$output" ] && [ "$(stat -c %s "$output")" -eq "$expected" ]; then
    echo "chunk=$index already complete"
    return
  fi
  tmp="$output.tmp"
  attempts=0
  while true; do
    have=$([ -f "$tmp" ] && stat -c %s "$tmp" || echo 0)
    if [ "$have" -gt "$expected" ]; then
      mv "$tmp" "$tmp.oversize"
      have=0
    fi
    [ "$have" -lt "$expected" ] || break
    attempts=$((attempts + 1))
    [ "$attempts" -le 200 ] || { echo "chunk=$index retry limit" >&2; return 4; }
    resume_start=$((start + have))
    segment="$tmp.segment"
    rm -f "$segment"
    # Keep partial bytes even when the proxy closes a long transfer; the next
    # attempt asks only for the remaining absolute byte range.
    curl --http1.1 --silent --show-error -fL --connect-timeout 30 --max-time 300 \
      --range "$resume_start-$end" --output "$segment" "$URL" || true
    if [ -s "$segment" ]; then
      cat "$segment" >> "$tmp"
      rm -f "$segment"
    else
      sleep 2
    fi
  done
  actual=$(stat -c %s "$tmp")
  if [ "$actual" -ne "$expected" ]; then
    echo "chunk=$index expected=$expected actual=$actual" >&2
    return 2
  fi
  mv "$tmp" "$output"
  echo "chunk=$index bytes=$actual complete"
}
export -f download_one
export ROOT URL SIZE SHA256 N_CHUNKS CHUNK_DIR TARGET
seq 0 $((N_CHUNKS - 1)) | xargs -n1 -P "$PARALLEL" bash -c 'download_one "$1"' _

parts=$(find "$CHUNK_DIR" -name '*.part' | wc -l)
[ "$parts" -eq "$N_CHUNKS" ] || { echo "only $parts/$N_CHUNKS parts" >&2; exit 3; }
find "$CHUNK_DIR" -name '*.part' -print0 | sort -z | xargs -0 cat > "$TARGET.tmp"
[ "$(stat -c %s "$TARGET.tmp")" -eq "$SIZE" ]
printf '%s  %s\n' "$SHA256" "$TARGET.tmp" | sha256sum -c -
mv "$TARGET.tmp" "$TARGET"
echo "checkpoint ready: $TARGET"
