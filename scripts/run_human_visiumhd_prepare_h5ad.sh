#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/Human_VisiumHD_compressed_v1.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/Human_VisiumHD_compressed_v1.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby Sample_id \
  --composition-cell-type Final_celltype \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns Final_celltype,Subclass_id,Class_id,Clinical_Brain_Diagnosis_Short,Region,Disease_region,module_label \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
