#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/cellxgene/cellxgene_visium_merged/hindlimb_10.1038_s41586-023-06806-x.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/cellxgene/cellxgene_visium_merged/hindlimb_10.1038_s41586-023-06806-x.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby source_file \
  --composition-cell-type cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cell_type,clusters,donor_id,development_stage,sex \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
