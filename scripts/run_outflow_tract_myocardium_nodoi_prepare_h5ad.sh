#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/cellxgene/cellxgene_visium_merged/outflow_tract_myocardium_nodoi.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/cellxgene/cellxgene_visium_merged/outflow_tract_myocardium_nodoi.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby source_file \
  --composition-cell-type cell_type \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns annotation,cell_type,cell_type_ontology_term_id \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
