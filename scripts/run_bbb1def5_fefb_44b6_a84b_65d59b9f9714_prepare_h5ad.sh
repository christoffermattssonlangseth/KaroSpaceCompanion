#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/bbb1def5-fefb-44b6-a84b-65d59b9f9714.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/bbb1def5-fefb-44b6-a84b-65d59b9f9714.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby Sample \
  --composition-cell-type cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cell_type,cluster,IR_VDJ_1_c_call,IR_VDJ_2_c_call,IR_VDJ_1_d_call,IR_VDJ_2_d_call,IR_VJ_1_j_call,IR_VJ_2_j_call \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
