#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing2/ST_BRICHOS/data/ST_BRICHOS_region_subcluster.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing2/ST_BRICHOS/data/ST_BRICHOS_region_subcluster.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type leiden \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns leiden,treatment,region_annotation,re_annotation_regions,leiden_0.75 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
