#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/baloMS_indep_clust_balo_MANA_SC.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/baloMS_indep_clust_balo_MANA_SC.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type niche \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns niche,segmentation_method,run,leiden_0.5,leiden_1,leiden_1.5,leiden_2,gmm_mana_10 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
