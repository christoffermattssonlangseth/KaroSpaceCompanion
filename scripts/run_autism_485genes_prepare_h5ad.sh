#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing2/autism/autism_concatenated_filtered_sparse_485genes.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing2/autism/autism_concatenated_filtered_sparse_485genes.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby source_file \
  --composition-cell-type tangram_cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns tangram_cell_type,anatomical_region \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
