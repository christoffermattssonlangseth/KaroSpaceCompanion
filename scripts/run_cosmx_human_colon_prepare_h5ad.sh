#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/notebooks/data/processed/cosmx-human-colon/cosmx-human-colon.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/notebooks/data/processed/cosmx-human-colon/cosmx-human-colon.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
