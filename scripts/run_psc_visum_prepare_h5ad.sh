#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/cellxgene/wrangled/psc-visum.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/cellxgene/wrangled/psc-visum.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cell_type,author_cell_type,donor_id,disease,sex,sample_preservation_method,sample_derivation_process \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
