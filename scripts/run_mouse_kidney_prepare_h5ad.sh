#!/usr/bin/env bash
set -euo pipefail

# ===========================================================================
# mouse_kidney.h5ad: ~276.9k cells x 56,748 genes, single sample
# (`mouse_kidney`), raw counts in X, with a pre-computed `lognorm` layer and an
# H&E `he_hires` tissue image under uns/spatial (preserved on copy).
#
# Cluster/domain columns (leiden, cellcharter_domains, cellcharter_k*) are
# already stored categorical, so no cluster->categorical conversion is needed.
#
# MEMORY: this is a big single section (276.9k cells x 56,748 genes). X alone is
# ~4-5 GB sparse; rebuilding a full `normalized` layer doubles the peak and gets
# the process OOM-killed (Jetsam) on a 16 GB machine when swap is tight. The file
# already ships a `lognorm` layer, so we pass --skip-normalized-layer: peak stays
# ~4.3 GB and the viewer still has normalized expression from `lognorm`.
# groupby is sample_id (one section). >200k cells => neighbor z-score off.
#
# Analytics columns are kept lean (leiden + cellcharter_domains) for the same
# reason; add the other cellcharter_k* columns back if running with more RAM.
# ===========================================================================

INPUT="${1:-/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpaceDataWrangling/data/processed/illumina-novaseqx-spatial/mouse_kidney.h5ad}"
OUTPUT_H5AD="${2:-$(dirname "$INPUT")/mouse_kidney.companion.ready.h5ad}"

cargo run --release --offline -- prepare "$INPUT" \
  --output "$OUTPUT_H5AD" \
  --delaunay \
  --groupby sample_id \
  --composition-cell-type leiden \
  --skip-normalized-layer \
  --skip-aggregation \
  --overwrite-derived \
  --persist-analytics-in-h5ad \
  --viewer-analytics-columns leiden,cellcharter_domains \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test

echo
echo "Output:"
echo "  H5AD: $OUTPUT_H5AD"
ls -lh "$OUTPUT_H5AD"

# Confirm the H&E tissue image survived the copy.
echo
echo "Spatial images in output:"
python3 - "$OUTPUT_H5AD" <<'PY'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:
    base = "uns/spatial"
    if base not in f:
        print("  (none) — no uns/spatial group")
        sys.exit(0)
    for sample in f[base]:
        img = f"{base}/{sample}/images"
        if img in f:
            for name in f[img]:
                d = f[f"{img}/{name}"]
                print(f"    {sample}/{name} shape={d.shape} dtype={d.dtype}")
PY
