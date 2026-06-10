#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.annotated.filtered.processed.cellcharter.h5ad}"
OUTPUT_H5AD="${2:-/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/RRMAP2_xenium_all_samples.annotated.filtered.processed.cellcharter.companion.ready.h5ad}"

# Large Xenium 5k dataset (~1.42M cells, 54 samples). X is already
# log-normalized (raw counts in layers["counts"]), so we skip recomputing the
# normalized layer and let analytics use X directly. No uns/spatial group =>
# no DAPI/HE images. At >200k cells the neighbor-stats z-score pass auto-off.
# Graph is built per kmeans_split_id tile (mirrors the prior RRMAP2 config).
cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby kmeans_split_id \
  --composition-cell-type CellCharter_10 \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns CellCharter_5,CellCharter_10,CellCharter_15,CellCharter_20,CellCharter_25,CellCharter_30,CellCharter_35,CellCharter_40,CellCharter_45,CellCharter_50,leiden_0.5,leiden_1,leiden_1.5,leiden_2,leiden_2.5,leiden_3,leiden_3.5,leiden_4.0 \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# This dataset has no uns/spatial images. Report honestly.
echo
echo "Spatial images in output:"
python3 - "$OUTPUT_H5AD" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    base = "uns/spatial"
    if base not in f:
        print("  (none) — this dataset has no uns/spatial images")
        sys.exit(0)
    found = False
    for sample in f[base]:
        img = f"{base}/{sample}/images"
        if img in f:
            for name in f[img]:
                d = f[f"{img}/{name}"]
                print(f"  {sample}/{name}: shape={d.shape} dtype={d.dtype}")
                found = True
    print("  OK" if found else "  (none) — no images found under uns/spatial")
PY