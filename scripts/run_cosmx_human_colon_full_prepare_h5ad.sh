#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/notebooks/data/processed/cosmx-human-colon/cosmx-human-colon-full.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/notebooks/data/processed/cosmx-human-colon/cosmx-human-colon-full.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby fov \
  --composition-cell-type cluster \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cluster \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
