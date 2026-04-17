#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/cellxgene/cellxgene_visium_merged/lung_10.1016_j.cell.2022.11.005.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/cellxgene/cellxgene_visium_merged/lung_10.1016_j.cell.2022.11.005.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby source_file \
  --composition-cell-type cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cell_type,annotation,clusters,leiden,regions,development_stage,donor_id,sex \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
