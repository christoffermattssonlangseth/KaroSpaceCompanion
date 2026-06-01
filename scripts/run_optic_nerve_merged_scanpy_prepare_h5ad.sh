#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/moldiassd/optic_nerve_merged.scanpy.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/moldiassd/optic_nerve_merged.scanpy.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type CellCharter_10 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_10,CellCharter_12,CellCharter_15,CellCharter_20,CellCharter_25,CellCharter_30,CellCharter_6,CellCharter_8,leiden,leiden_0_2,leiden_0_4,leiden_0_6,leiden_0_8,leiden_1_0,leiden_1_5,leiden_2_0,leiden_2_5,leiden_3_0,leiden_3_5,leiden_4_0 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
