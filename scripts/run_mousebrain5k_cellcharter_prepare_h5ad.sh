#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/mouseBrain5k/mouseBrain5k_cellcharter.h5ad}"
OUTPUT_H5AD="${2:-/tmp/mouseBrain5k_cellcharter.companion.ready.h5ad}"

cargo run --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --k 8 \
  --groupby sample_id \
  --composition-cell-type CellCharter_10 \
  --normalize-from layer:counts \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_5,CellCharter_15,CellCharter_20,leiden_0.1,leiden_0.5,leiden_1.0 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

