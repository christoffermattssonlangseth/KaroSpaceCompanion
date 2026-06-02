#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/xenium-dapi-he/breast_cancer_2samples.companion.annotated.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/xenium-dapi-he/breast_cancer_2samples.companion.annotated.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type leiden \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns leiden,leiden_label \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# Confirm the KaroSpace spatial images survived the copy (prepare does a full
# fs::copy of the input, so uns/spatial/.../images is preserved verbatim).
echo
echo "Spatial images in output:"
python3 - "$OUTPUT_H5AD" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    base = "uns/spatial"
    if base not in f:
        print("  WARNING: uns/spatial missing from output!"); sys.exit(1)
    found = False
    for sample in f[base]:
        img = f"{base}/{sample}/images"
        if img in f:
            for name in f[img]:
                d = f[f"{img}/{name}"]
                print(f"  {sample}/{name}: shape={d.shape} dtype={d.dtype}")
                found = True
    print("  OK" if found else "  WARNING: no images found under uns/spatial!")
PY