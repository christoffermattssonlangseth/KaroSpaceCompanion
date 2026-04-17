#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/processing/RRmap_with_refined_anno_and_myeloid_lineage_and_UCell_signature_scoring_scanpy_scores.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/processing/RRmap_with_refined_anno_and_myeloid_lineage_and_UCell_signature_scoring_scanpy_scores.companion.ready.h5ad}"

SCORE_COLUMNS="Surveillance_score,Neuroprotection_score,Phagocytosis_score,Inflammation_score,Cytokine_score,Antigen_presentation_score,IFN_score,Proliferation_score,ECM remodelling_score,Myelin debris_score,Anti-inflammatory_score"
python scripts/bin_obs_numeric_to_tertiles.py "$INPUT" --columns "$SCORE_COLUMNS" --suffix "_tertile"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type anno_L2 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns anno_L3,anno_L2,anno_L1,anno_hanna,anno_hanna_new,myeloid_lineage,UCell_supercluster_myeloid,compartment_mana,compartment_anno,Class \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"
