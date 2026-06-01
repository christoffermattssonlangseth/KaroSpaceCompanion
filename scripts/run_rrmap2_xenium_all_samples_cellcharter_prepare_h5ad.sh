#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.cellcharter.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.cellcharter.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby kmeans_split_id \
  --composition-cell-type CellCharter_10 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_10,CellCharter_15,CellCharter_20,CellCharter_25,CellCharter_30,CellCharter_35,CellCharter_40,CellCharter_45,CellCharter_5,CellCharter_50,leiden_0.5,leiden_1,leiden_1.5,leiden_2,leiden_2.5,leiden_3,leiden_3.5,leiden_4.0 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
