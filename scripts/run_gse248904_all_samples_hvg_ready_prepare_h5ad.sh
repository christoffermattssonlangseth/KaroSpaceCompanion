#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/christoffer/work/karolinska/development/KaroSpaceDataWrangling/data/gse248904/GSE248904_All_Samples_HVG.ready.h5ad}"
OUTPUT_H5AD="${2:-/Users/christoffer/work/karolinska/development/KaroSpaceDataWrangling/data/gse248904/GSE248904_All_Samples_HVG.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby Sample \
  --composition-cell-type clusters \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns 'clusters,Treatment,Organ_Full_Name,Subregion,Organ,Other Annotation' \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
