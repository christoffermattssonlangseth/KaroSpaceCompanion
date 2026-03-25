#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing2/RRmap/data/RRmap_metadata_fixed_update.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing2/RRmap/data/rrmap.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type anno_L2 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns anno_L3,anno_L2,anno_L1,leiden_3.5,compartment_mana,stage,condition \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
