#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/muBrainRelease_seurat.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/muBrainRelease_seurat.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type RNA_nbclust_clusters_long \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns RNA_nbclust_clusters_long,RNA_nbclust_clusters,spatialClusteringAssignments,tissue,qcFlagsCellComplex,spatialClusteringAssignments \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
