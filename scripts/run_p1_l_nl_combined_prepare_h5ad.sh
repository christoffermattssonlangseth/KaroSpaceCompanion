#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/xenium-P1-annotated/P1_L_NL_combined.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/xenium-P1-annotated/P1_L_NL_combined.companion.ready.h5ad}"

# This dataset is a Seurat export: X holds log-normalized values and raw counts
# live in layers["counts"], so we normalize from the counts layer. There is no
# uns/spatial group, so this dataset carries no DAPI/HE images.
cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type predicted.celltype \
  --normalize-from layer:counts \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns predicted.celltype,seurat_clusters,niches,condition \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# This dataset has no uns/spatial images (Seurat export). Report honestly.
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