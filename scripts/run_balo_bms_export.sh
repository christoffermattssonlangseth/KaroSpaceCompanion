#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/baloMS_indep_clust_balo_MANA_balo_annot.h5ad}"
OUTPUT_H5AD="${2:-/tmp/baloMS_full.fast.enriched.h5ad}"
OUTPUT_JSON="${3:-/tmp/baloMS_full.fast.viewer.json}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type leiden_2_names_sub \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --viewer-json "$OUTPUT_JSON" \
  --initial-color leiden_2_names_sub \
  --viewer-analytics-columns leiden_2_names_sub \
  --skip-viewer-interaction-markers \
  --viewer-genes MBP,GFAP,CLDN5,SLC1A3,SNAP25 \
  --viewer-cluster-de-method t-test

echo
echo "Outputs:"
echo "  H5AD: $OUTPUT_H5AD"
echo "  JSON: $OUTPUT_JSON"
ls -lh "$OUTPUT_H5AD" "$OUTPUT_JSON"
