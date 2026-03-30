#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/GSE248904_All_Samples_HVG.ready.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/GSE248904_All_Samples_HVG.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby Sample \
  --composition-cell-type clusters \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns clusters,Treatment,Organ_Full_Name,Subregion \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
