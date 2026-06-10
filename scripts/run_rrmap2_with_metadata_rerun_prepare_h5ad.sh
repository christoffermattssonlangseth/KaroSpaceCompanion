#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.cellcharter.companion.ready.with_metadata.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.cellcharter.companion.ready.with_metadata.rerun.h5ad}"

# Re-run of the companion pipeline on the already-processed + metadata-enriched
# file (~1.38M cells, 54 samples). Same clusterings as the recorded
# uns/karospace_companion config, PLUS the sample-level metadata columns
# (stage, condition, region, sex, model) added as analytics columns so they get
# cluster-DE (differential expression between those biological groups). X is
# log-normalized (raw counts in layers["counts"]) so we skip recomputing the
# normalized layer. No uns/spatial => no DAPI/HE images. At >200k cells the
# neighbor-stats z-score pass auto-off.
#
# NOTE: prepare copies the input verbatim then refreshes only the
# uns/karospace_companion fields, so the existing uns/karospace_polygons and
# uns/karospace_category_labels and all metadata columns are preserved.
cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby kmeans_split_id \
  --composition-cell-type CellCharter_10 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_5,CellCharter_10,CellCharter_15,CellCharter_20,CellCharter_25,CellCharter_30,CellCharter_35,CellCharter_40,CellCharter_45,CellCharter_50,leiden_0.5,leiden_1,leiden_1.5,leiden_2,leiden_2.5,leiden_3,leiden_3.5,leiden_4.0,stage,condition,region,sex,model \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# Confirm metadata extras carried through, and report on images.
echo
echo "Preserved extras + images in output:"
python3 - "$OUTPUT_H5AD" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    for key in ["uns/karospace_polygons", "uns/karospace_category_labels"]:
        print(f"  {key}: {'present' if key in f else 'MISSING'}")
    base = "uns/spatial"
    if base not in f:
        print("  uns/spatial images: (none) — this dataset has no images")
        sys.exit(0)
    found = False
    for sample in f[base]:
        img = f"{base}/{sample}/images"
        if img in f:
            for name in f[img]:
                d = f[f"{img}/{name}"]
                print(f"  {sample}/{name}: shape={d.shape} dtype={d.dtype}")
                found = True
    print("  images OK" if found else "  uns/spatial images: (none)")
PY