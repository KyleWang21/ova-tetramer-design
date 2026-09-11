#!/usr/bin/env bash
# Run one Rosetta C4 screen as independent, contiguous seed chunks.
set -euo pipefail

SOURCE="${1:?usage: run_rosetta_c4_parallel.sh CIF CANDIDATE ALL_EDGES [DESIGN_EDGES] SOURCE_MODEL CHUNKS_ROOT OUT}"
CANDIDATE="${2:?usage: run_rosetta_c4_parallel.sh CIF CANDIDATE ALL_EDGES [DESIGN_EDGES] SOURCE_MODEL CHUNKS_ROOT OUT}"
EDGES="${3:?usage: run_rosetta_c4_parallel.sh CIF CANDIDATE ALL_EDGES [DESIGN_EDGES] SOURCE_MODEL CHUNKS_ROOT OUT}"
if (( $# == 6 )); then
  # Legacy invocation: chemistry was evaluated on every effective edge.
  DESIGN_EDGES="$EDGES"
  SOURCE_MODEL="${4:?}"
  CHUNKS_ROOT="${5:?}"
  OUT="${6:?}"
elif (( $# == 7 )); then
  DESIGN_EDGES="${4:?}"
  SOURCE_MODEL="${5:?}"
  CHUNKS_ROOT="${6:?}"
  OUT="${7:?}"
else
  echo "usage: run_rosetta_c4_parallel.sh CIF CANDIDATE ALL_EDGES [DESIGN_EDGES] SOURCE_MODEL CHUNKS_ROOT OUT" >&2
  exit 2
fi
PROJECT=/root/400083/ova_p3_tetramer_design
PY=/root/400083/antibody-design/pipeline_stage1_interface_energetics/env/bin/python
# PyRosetta imports NumPy/BLAS and can otherwise create one large thread pool
# per process. This launcher already runs independent Relax workers in
# parallel, so nested threading causes severe oversubscription (64 workers x
# approximately 127 threads each was observed on the dev instance). Keep each
# worker single-threaded and use process-level parallelism instead.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export BLIS_NUM_THREADS="${BLIS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
NCHUNKS="${OVA_ROSETTA_CHUNKS:-8}"
TOTAL="${OVA_ROSETTA_TOTAL:-200}"
SEED_START=7001

if (( NCHUNKS < 1 || NCHUNKS > TOTAL )); then
  echo "invalid NCHUNKS=$NCHUNKS for TOTAL=$TOTAL" >&2
  exit 2
fi
BASE=$((TOTAL / NCHUNKS))
REMAINDER=$((TOTAL % NCHUNKS))

cd "$PROJECT"
mkdir -p "$CHUNKS_ROOT" "$OUT"
pids=()
for chunk in $(seq 0 $((NCHUNKS - 1))); do
  chunk_name=$(printf 'chunk%02d' "$chunk")
  chunk_out="$CHUNKS_ROOT/$chunk_name"
  chunk_extra=0
  if (( chunk < REMAINDER )); then
    chunk_extra=1
  fi
  chunk_count=$((BASE + chunk_extra))
  preceding_extra=$chunk
  if (( preceding_extra > REMAINDER )); then
    preceding_extra=$REMAINDER
  fi
  chunk_seed=$((SEED_START + chunk * BASE + preceding_extra))
  mkdir -p "$chunk_out"
  (
    "$PY" af3_pipeline/rosetta_c4_relax_screen.py \
      --in "$SOURCE" --candidate "$CANDIDATE" --effective-edges "$EDGES" \
      --design-interface-edges "$DESIGN_EDGES" \
      --source-model "$SOURCE_MODEL" --nrelax "$chunk_count" \
      --seed-start "$chunk_seed" --out "$chunk_out"
  ) >"$chunk_out/worker.log" 2>&1 &
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    failed=1
  fi
done
if (( failed != 0 )); then
  echo "one or more Rosetta chunks failed" >&2
  exit 1
fi

"$PY" af3_pipeline/combine_rosetta_c4_chunks.py \
  --chunks-root "$CHUNKS_ROOT" --candidate "$CANDIDATE" \
  --source-model "$SOURCE_MODEL" --expected "$TOTAL" --expected-chunks "$NCHUNKS" \
  --seed-start "$SEED_START" --out "$OUT"
