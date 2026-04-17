#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/GSE282203_combined.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/GSE282203_combined.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample \
  --composition-cell-type leiden_1.0 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns leiden_1.0,leiden_0.5,gsm_id,capture_area,genotype,age,sex,zeitgeber_time,tissue \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
