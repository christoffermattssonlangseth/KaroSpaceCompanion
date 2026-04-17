#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/talbot_xenium_tumor_annotated_updated_cellcharter.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/talbot_xenium_tumor_annotated_updated_cellcharter.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type cytetype_annotation_leiden_4 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_20,CellCharter_15,condition,genotype,cytetype_annotation_leiden_4,cytetype_cellState_leiden_4,leiden_4\
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
