#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/baloMS_indep_clust_balo_MANA_balo_annot.h5ad}"
OUTPUT_H5AD="${2:-/tmp/baloMS_companion.ready.h5ad}"

cargo run --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --k 8 \
  --groupby sample_id \
  --composition-cell-type leiden_2_names_sub \
  --normalize-from layer:counts \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns leiden_2_names_sub,leiden_0.5 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
