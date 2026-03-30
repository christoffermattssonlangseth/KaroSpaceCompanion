#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/Nanostring_CosMX_v1.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/Nanostring_CosMX_v1.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby Sample.FOV \
  --composition-cell-type new_identity \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns new_identity,Subclass_label,Pathology,module_label,Max_10x_Mg_cluster_mapped \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
