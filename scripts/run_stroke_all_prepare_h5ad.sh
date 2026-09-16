#!/usr/bin/env bash
set -euo pipefail

# ===========================================================================
# stroke_all.h5ad: ~3.66M cells x 500 genes, 49 samples, raw counts in X,
# and a `hires` tissue image per sample under uns/spatial (preserved on copy).
#
# DISK: input is ~7.8 GB. `prepare` copies input->output (>=7.8 GB) and the
# cluster->categorical step needs another ~7.8 GB temp => budget ~16 GB free on
# the OUTPUT volume. Pass a path with more space as $2 if needed.
#
# `cluster` is stored numeric (int64, values 1/2/3). The companion rejects
# numeric columns as analytics/cluster columns, so we convert it to categorical
# (same approach used for banksy) into a temp input that `prepare` consumes.
# ===========================================================================

INPUT="${1:-/Users/chrislangseth/Downloads/stroke_all.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/stroke_all.companion.ready.h5ad}"

PREP_INPUT="$(dirname "$OUTPUT_H5AD")/$(basename "${INPUT%.h5ad}").cluster_cat.h5ad"

echo "Preparing curated input (cluster -> categorical) at: $PREP_INPUT"
cp -f "$INPUT" "$PREP_INPUT"

python3 - "$PREP_INPUT" <<'PY'
import sys, h5py, numpy as np
dst = sys.argv[1]
with h5py.File(dst, "r+") as f:
    obs = f["obs"]
    col = obs["cluster"]
    if isinstance(col, h5py.Group):
        print("  cluster already categorical; nothing to do")
    else:
        vals = col[:].astype(np.int64)
        cats = sorted(set(vals.tolist()))
        code_map = {c: i for i, c in enumerate(cats)}
        codes = np.array([code_map[v] for v in vals], dtype=np.int8)
        del obs["cluster"]
        g = obs.create_group("cluster")
        g.attrs["encoding-type"] = "categorical"
        g.attrs["encoding-version"] = "0.2.0"
        g.attrs["ordered"] = np.bool_(False)
        str_dt = h5py.string_dtype(encoding="utf-8")
        cds = g.create_dataset("categories",
                               data=np.array([str(c) for c in cats], dtype=object),
                               dtype=str_dt)
        cds.attrs["encoding-type"] = "string-array"
        cds.attrs["encoding-version"] = "0.2.0"
        cod = g.create_dataset("codes", data=codes)
        cod.attrs["encoding-type"] = "array"
        cod.attrs["encoding-version"] = "0.2.0"
        print(f"  converted cluster -> categorical ({len(cats)} categories: {[str(c) for c in cats]})")
PY

# X is raw counts => normalize from X (default) to build layers["normalized"].
# groupby per sample (49 tissue sections). >200k cells => neighbor z-score off.
cargo run --release --offline -- prepare "$PREP_INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample \
  --composition-cell-type cluster \
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

# Clean up the temporary curated input.
rm -f "$PREP_INPUT"

# Confirm the per-sample tissue images survived the copy.
echo
echo "Spatial images in output:"
python3 - "$OUTPUT_H5AD" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    base = "uns/spatial"
    if base not in f:
        print("  (none) — no uns/spatial group")
        sys.exit(0)
    n_with_img = 0
    examples = []
    for sample in f[base]:
        img = f"{base}/{sample}/images"
        if img in f:
            n_with_img += 1
            if len(examples) < 3:
                for name in f[img]:
                    d = f[f"{img}/{name}"]
                    examples.append(f"{sample}/{name} shape={d.shape} dtype={d.dtype}")
    print(f"  {n_with_img} sample(s) carry images. e.g.:")
    for e in examples:
        print(f"    {e}")
PY