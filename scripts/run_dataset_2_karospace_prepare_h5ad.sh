#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/dataset_2_karospace.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/dataset_2_karospace.companion.ready.h5ad}"

# Mouse-brain neuroinflammation Xenium dataset (~744k cells; conditions
# sal/lps/poly x 3 replicates = 9 physical sections). X is log-normalized and
# raw counts live in layers["counts"], so we normalize from the counts layer.
# No uns/spatial group => no DAPI/HE images. At >200k cells the neighbor-stats
# z-score permutation pass auto-disables.
#
# groupby     : sample_repl (9 physical tissue sections)
# composition : celltype_major (27 major cell types)
# analytics   : cell-type annotations + spatial GMM domains + leiden
#               resolutions + condition/replicate metadata.
cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_repl \
  --composition-cell-type celltype_major \
  --normalize-from layer:counts \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns celltype_major,cytetype_annotation_leiden_0.2,cytetype_leiden_0.2_astrosub,astro_subclass,leiden_0.2,leiden_0.5,leiden_1,gmm_CC_15,gmm_CC_20,gmm_CC_25,gmm_CC_30,gmm_CC_35,sample_id,replicate \
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