#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/LABEL_mouse_all_concat_shared_genes_log1p_non3d.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/LABEL_mouse_all_concat_shared_genes_log1p_non3d.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby source_file \
  --composition-cell-type clusters \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns clusters,organ_standardized,Subregion,sample_id,Tissue_ID \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
