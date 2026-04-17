#!/usr/bin/env python3
"""
Generate a KaroSpaceCompanion runner shell script for an H5AD dataset.

Usage:
    python3 scripts/generate_runner.py DATASET.h5ad [--repo-root .] [--groupby COL] ...
"""

import argparse
import os
import re
import sys
from pathlib import Path

GROUPBY_PRIORITY = [
    "Sample.FOV", "source_file", "Sample_id", "sample_id", "sample",
    "Sample", "slide", "orig.ident", "FOV",
]

COMPOSITION_PRIORITY = [
    "Final_celltype", "new_identity", "majorCluster_final", "tangram_cell_type",
    "anno_L2", "niches", "areas", "CellTypes", "cell_type_E", "Subclass_label",
    "Class_id", "cell_type", "CellCharter_10", "RNA_snn_res.2", "seurat_clusters",
    "leiden", "louvain", "clusters", "cluster",
]

ANALYTICS_EXCLUDE_RE = re.compile(
    r"^(x|y|X|Y)$"
    r"|^(x_|y_|X_|Y_)"
    r"|(_x|_y)$"
    r"|^n_"
    r"|(_count|_counts)$"
    r"|^(pct_|log1p_|total_|percent_|__)",
    re.IGNORECASE,
)


def pick_column(cols, priority):
    col_set = set(cols)
    for c in priority:
        if c in col_set:
            return c
    return None


def is_good_analytics_col(col, obs, groupby, composition,
                           max_categories=200, max_unique_ratio=0.3):
    if col in (groupby, composition):
        return False
    if ANALYTICS_EXCLUDE_RE.search(col):
        return False
    series = obs[col]
    dtype_str = str(series.dtype)
    if dtype_str in ("float32", "float64"):
        return False
    try:
        n_unique = series.nunique()
    except Exception:
        return False
    n_total = len(series)
    if n_unique < 2 or n_unique > max_categories:
        return False
    if n_total > 0 and n_unique / n_total > max_unique_ratio:
        return False
    return True


def looks_normalized(adata):
    """Return True if X appears to be log-normalised (skip raw normalisation step)."""
    try:
        import numpy as np
        import scipy.sparse as sp

        X = adata.X
        n = min(200, X.shape[0])
        sample = X[:n].toarray() if sp.issparse(X) else np.asarray(X[:n])
        max_val = float(sample.max())
        pos = sample[sample > 0]
        mean_pos = float(pos.mean()) if pos.size else 0.0
        # Looks normalised: small values, not integer-like
        if max_val < 20 and mean_pos < 5:
            return True
        frac_int = float((sample == sample.astype(int)).mean())
        if frac_int > 0.9 and max_val > 10:
            return False   # raw counts
        return True
    except Exception:
        return True  # conservative: assume normalised


def script_stem(h5ad_path):
    stem = Path(h5ad_path).stem
    stem = re.sub(r"\.(h5ad|ready|companion)$", "", stem, flags=re.IGNORECASE)
    return re.sub(r"[^a-zA-Z0-9]+", "_", stem).lower().strip("_")


class _H5adProxy:
    """Minimal stand-in for AnnData backed object, built from raw h5py."""

    def __init__(self, path):
        import h5py
        import numpy as np

        self._f = h5py.File(path, "r")
        obs_grp = self._f["obs"]

        # obs columns (skip _index)
        self._obs_cols = [k for k in obs_grp.keys() if k != "_index"]

        # n_obs from X indptr (CSR) or first obs column length
        if "X" in self._f and "indptr" in self._f["X"]:
            self.n_obs = len(self._f["X"]["indptr"]) - 1
        else:
            first = self._obs_cols[0]
            v = obs_grp[first]
            self.n_obs = (len(v["codes"]) if isinstance(v, h5py.Group) and "codes" in v
                          else len(v))

        # n_vars
        var_grp = self._f["var"]
        var_cols = [k for k in var_grp.keys() if k != "_index"]
        if var_cols:
            v = var_grp[var_cols[0]]
            self.n_vars = (len(v["codes"]) if isinstance(v, h5py.Group) and "codes" in v
                           else len(v))
        else:
            self.n_vars = 0

        # Build obs DataFrame proxy (dict-like for nunique / dtype checks)
        self._obs_grp = obs_grp
        self.obs = self   # reuse self as the "obs" object

        # layers / uns presence
        self.layers = list(self._f["layers"].keys()) if "layers" in self._f else []
        self.uns_keys = list(self._f["uns"].keys()) if "uns" in self._f else []

    @property
    def columns(self):
        return self._obs_cols

    def __getitem__(self, key):
        return _H5ColProxy(self._obs_grp[key])

    def nunique_col(self, key):
        return self[key].nunique()

    def dtype_str(self, key):
        return self[key].dtype_str

    def close(self):
        self._f.close()


class _H5ColProxy:
    """Lightweight proxy for a single obs column stored in h5py."""

    def __init__(self, ds_or_grp):
        import h5py
        self._v = ds_or_grp
        self._is_cat = isinstance(ds_or_grp, h5py.Group) and "categories" in ds_or_grp

    @property
    def dtype(self):
        if self._is_cat:
            return type("_DT", (), {"__str__": lambda s: "category"})()
        ds = self._v
        return type("_DT", (), {"__str__": lambda s: str(ds.dtype)})()

    @property
    def dtype_str(self):
        return "category" if self._is_cat else str(self._v.dtype)

    def nunique(self):
        if self._is_cat:
            return len(self._v["categories"])
        import numpy as np
        return int(np.unique(self._v[:]).shape[0])

    def __len__(self):
        if self._is_cat:
            return len(self._v["codes"])
        return len(self._v)


def _read_dataset(path):
    """Try anndata first; fall back to h5py proxy on codec errors."""
    # anndata path
    try:
        import anndata as ad
        adata = ad.read_h5ad(path, backed="r")
        return adata, False   # (obj, is_proxy)
    except ImportError:
        sys.exit("ERROR: anndata not installed — run: pip install anndata")
    except Exception as exc:
        print(f"  anndata failed ({exc.__class__.__name__}: {exc}); falling back to h5py.",
              file=sys.stderr)

    return _H5adProxy(path), True


def _obs_cols(obj, is_proxy):
    if is_proxy:
        return obj.obs.columns
    return list(obj.obs.columns)


def _is_good_analytics_col_generic(col, obj, is_proxy, groupby, composition):
    if col in (groupby, composition):
        return False
    if ANALYTICS_EXCLUDE_RE.search(col):
        return False
    if is_proxy:
        proxy = obj.obs[col]
        dtype_s = proxy.dtype_str
        if dtype_s in ("float32", "float64"):
            return False
        try:
            n_unique = proxy.nunique()
        except Exception:
            return False
    else:
        return is_good_analytics_col(col, obj.obs, groupby, composition)
    n_total = obj.n_obs
    if n_unique < 2 or n_unique > 200:
        return False
    if n_total > 0 and n_unique / n_total > 0.3:
        return False
    return True


def generate(args):
    print(f"Reading {args.dataset_path} …", file=sys.stderr)
    obj, is_proxy = _read_dataset(args.dataset_path)
    obs_cols = list(_obs_cols(obj, is_proxy))

    print(f"  {obj.n_obs:,} cells × {obj.n_vars:,} genes", file=sys.stderr)
    print(f"  obs columns: {obs_cols}", file=sys.stderr)

    # ── groupby ──────────────────────────────────────────────────────────────
    groupby = args.groupby or pick_column(obs_cols, GROUPBY_PRIORITY)
    if groupby is None:
        groupby = obs_cols[0] if obs_cols else "sample"
        print(f"  WARNING: no standard groupby column found, defaulting to '{groupby}'",
              file=sys.stderr)

    # ── composition ───────────────────────────────────────────────────────────
    composition = args.composition or pick_column(obs_cols, COMPOSITION_PRIORITY)
    if composition is None:
        composition = obs_cols[1] if len(obs_cols) > 1 else "cell_type"
        print(f"  WARNING: no standard composition column found, defaulting to '{composition}'",
              file=sys.stderr)

    # ── analytics columns ─────────────────────────────────────────────────────
    if args.analytics:
        analytics = [c.strip() for c in args.analytics.split(",")]
    else:
        candidates = [
            c for c in obs_cols
            if _is_good_analytics_col_generic(c, obj, is_proxy, groupby, composition)
        ]
        prio_set = {c: i for i, c in enumerate(COMPOSITION_PRIORITY)}
        candidates.sort(key=lambda c: prio_set.get(c, len(COMPOSITION_PRIORITY)))
        analytics = [composition] + [c for c in candidates if c != composition]
        analytics = analytics[:8]

    # ── normalisation ─────────────────────────────────────────────────────────
    if args.skip_normalized_layer:
        skip_norm = True
    elif is_proxy:
        # Check uns for log1p key, or sample X data
        skip_norm = "log1p" in obj.uns_keys
        if not skip_norm:
            skip_norm = looks_normalized(obj._f)  # pass h5py file
    else:
        skip_norm = looks_normalized(obj)

    # ── paths ─────────────────────────────────────────────────────────────────
    repo_root = args.repo_root or os.getcwd()
    stem = script_stem(args.dataset_path)
    script_path = args.script_path or os.path.join(repo_root, "scripts",
                                                    f"run_{stem}_prepare_h5ad.sh")
    input_path = args.dataset_path
    output_h5ad = (args.output_h5ad
                   or str(Path(input_path).parent / f"{Path(input_path).stem}.companion.ready.h5ad"))

    # ── build script ──────────────────────────────────────────────────────────
    analytics_str = ",".join(analytics)

    body_args = [
        '  --output "$OUTPUT_H5AD" \\',
        '  --delaunay \\',
        f"  --groupby {groupby} \\",
        f"  --composition-cell-type {composition} \\",
    ]
    if skip_norm:
        body_args.append("  --skip-normalized-layer \\")
    body_args += [
        "  --skip-aggregation \\",
        "  --overwrite-derived \\",
        "  --persist-analytics-in-h5ad \\",
        f"  --viewer-analytics-columns {analytics_str} \\",
        "  --skip-viewer-interaction-markers \\",
    ]
    if args.viewer_neighbor_permutations:
        body_args.append(f"  --viewer-neighbor-permutations {args.viewer_neighbor_permutations} \\")
    body_args.append("  --viewer-cluster-de-method t-test")

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f'INPUT="${{1:-{input_path}}}"',
        f'OUTPUT_H5AD="${{2:-{output_h5ad}}}"',
        "",
        'cargo run --release --offline -- prepare "$INPUT" \\',
        *body_args,
        "",
        "echo",
        'echo "Output:"',
        'echo "  H5AD: $OUTPUT_H5AD"',
        'ls -lh "$OUTPUT_H5AD"',
        "",
    ]

    script_content = "\n".join(lines)
    Path(script_path).parent.mkdir(parents=True, exist_ok=True)
    Path(script_path).write_text(script_content)
    os.chmod(script_path, 0o755)

    # ── validate syntax ───────────────────────────────────────────────────────
    rc = os.system(f"bash -n {script_path}")  # noqa: S605
    if rc != 0:
        print(f"WARNING: bash -n reported a syntax error in {script_path}", file=sys.stderr)

    print(f"\nWrote: {script_path}", file=sys.stderr)
    print("\n--- Chosen defaults ---")
    print(f"  groupby:              {groupby}")
    print(f"  composition:          {composition}")
    print(f"  analytics:            {analytics_str}")
    print(f"  skip-normalized-layer:{skip_norm}")
    print(f"  output h5ad:          {output_h5ad}")
    print(f"  script:               {script_path}")


def main():
    p = argparse.ArgumentParser(
        description="Generate a KaroSpaceCompanion runner shell script for an H5AD dataset.",
    )
    p.add_argument("dataset_path", help="Path to the .h5ad file")
    p.add_argument("--repo-root", default=None,
                   help="Repo root directory (default: cwd)")
    p.add_argument("--groupby", default=None,
                   help="Override the groupby column")
    p.add_argument("--composition", default=None,
                   help="Override the composition-cell-type column")
    p.add_argument("--analytics", default=None,
                   help="Override analytics columns (comma-separated)")
    p.add_argument("--skip-normalized-layer", action="store_true", default=False,
                   help="Force --skip-normalized-layer in the generated script")
    p.add_argument("--viewer-neighbor-permutations", type=int, default=None,
                   help="Add --viewer-neighbor-permutations N")
    p.add_argument("--script-path", default=None,
                   help="Override output script file path")
    p.add_argument("--output-h5ad", default=None,
                   help="Override output .h5ad path baked into the script")
    args = p.parse_args()
    generate(args)


if __name__ == "__main__":
    main()
