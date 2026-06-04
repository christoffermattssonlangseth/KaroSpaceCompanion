#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-/Users/chrislangseth/Downloads/deitaled_wBanky-lam08res03_pID.h5ad}"
OUTPUT_H5AD="${2:-/Users/chrislangseth/Downloads/deitaled_wBanky-lam08res03_pID.companion.ready.h5ad}"

# ---------------------------------------------------------------------------
# Preprocess: `banksy` is stored as a numeric int64 column (values 0-9). The
# companion treats numeric obs columns as continuous and rejects them as
# analytics/cluster columns, so we convert it to a categorical column first.
#
# We edit a COPY with h5py directly (not anndata): this file's
# uns/cytetype_jobDetails/auth_token is `null`-encoded and trips
# anndata.read_h5ad, so we avoid a full AnnData load and touch only obs/banksy.
# ---------------------------------------------------------------------------
PREP_INPUT="${TMPDIR:-/tmp}/$(basename "${INPUT%.h5ad}").banksy_cat.h5ad"

echo "Preparing curated input (banksy -> categorical)..."
cp -f "$INPUT" "$PREP_INPUT"

python3 - "$PREP_INPUT" <<'PY'
import sys, h5py, numpy as np
dst = sys.argv[1]
with h5py.File(dst, "r+") as f:
    obs = f["obs"]
    col = obs["banksy"]
    if isinstance(col, h5py.Group):
        print("  banksy already categorical; nothing to do")
    else:
        vals = col[:].astype(np.int64)
        cats = sorted(set(vals.tolist()))
        code_map = {c: i for i, c in enumerate(cats)}
        codes = np.array([code_map[v] for v in vals], dtype=np.int8)
        del obs["banksy"]
        g = obs.create_group("banksy")
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
        print(f"  converted banksy -> categorical ({len(cats)} categories)")
PY

# ---------------------------------------------------------------------------
# Prepare. ~212k cells, 26 sections. X is log-normalized; raw counts live in
# layers["counts"], so we normalize from the counts layer. No uns/spatial group
# => no DAPI/HE images. At >200k cells the neighbor-stats z-score pass auto-off.
#
# Column choices (edits to the previous detailedAnnot viewer config):
#   - groupby            : sample_id   -> sample_NewNAME
#   - composition        : cell_annotation_detailed
#   - analytics/colors   : cell_class, cell_annotation_detailed, banksy
#     (removed: condition, status, cell_subclass, cluster_cellcharter_*, leiden_*)
#   - patient_id stays a categorical obs column (selectable as metadata/color).
# ---------------------------------------------------------------------------
cargo run --release --offline -- prepare "$PREP_INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_NewNAME \
  --composition-cell-type cell_annotation_detailed \
  --normalize-from layer:counts \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns cell_class,cell_annotation_detailed,banksy \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# Clean up the temporary curated input.
rm -f "$PREP_INPUT"

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