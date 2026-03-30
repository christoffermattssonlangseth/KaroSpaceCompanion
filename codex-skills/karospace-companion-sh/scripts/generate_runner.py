#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import h5py

GROUPBY_PRIORITY = [
    "Sample.FOV",
    "source_file",
    "Sample_id",
    "sample_id",
    "sample",
    "Sample",
    "slide",
    "orig.ident",
    "FOV",
]

COMPOSITION_PRIORITY = [
    "Final_celltype",
    "new_identity",
    "majorCluster_final",
    "tangram_cell_type",
    "anno_L2",
    "niches",
    "areas",
    "CellTypes",
    "cell_type_E",
    "Subclass_label",
    "Class_id",
    "cell_type",
    "CellCharter_10",
    "RNA_snn_res.2",
    "seurat_clusters",
]

ANALYTICS_PRIORITY = [
    "Final_celltype",
    "new_identity",
    "majorCluster_final",
    "niches",
    "areas",
    "Subclass_id",
    "Subclass_label",
    "Class_id",
    "cell_type",
    "CellTypes",
    "cell_type_E",
    "anno_L3",
    "anno_L2",
    "anno_L1",
    "Clinical_Brain_Diagnosis_Short",
    "Clinical_Brain_Diagnosis",
    "Disease_region",
    "Region",
    "Tissue_Region",
    "Tissue_Type",
    "Pathology",
    "condition",
    "stage",
    "module_label",
    "patient_id",
    "Sample",
    "Sample_id",
    "sample_id",
    "Subject_ID",
    "Max_10x_Mg_cluster_mapped",
    "RNA_snn_res.2",
    "seurat_clusters",
    "seurat_clusters_Class_id",
    "CellCharter_5",
    "CellCharter_10",
    "CellCharter_15",
    "CellCharter_20",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a KaroSpaceCompanion dataset runner shell script."
    )
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--script-path", type=Path)
    parser.add_argument("--output-h5ad", type=Path)
    parser.add_argument("--groupby")
    parser.add_argument("--composition")
    parser.add_argument("--analytics")
    parser.add_argument("--viewer-neighbor-permutations", type=int)
    return parser.parse_args()


def decode_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def read_obs_column(obs: h5py.Group, name: str) -> list[str]:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = [decode_value(value) for value in obj["categories"][:]]
        codes = obj["codes"][:]
        values = []
        for code in codes:
            if code < 0:
                values.append("NA")
            else:
                values.append(categories[int(code)])
        return values
    values = obj[:]
    return [decode_value(value) for value in values]


def is_categorical_obs_column(obs: h5py.Group, name: str) -> bool:
    obj = obs[name]
    if isinstance(obj, h5py.Group):
        return "categories" in obj and "codes" in obj
    dtype = obj.dtype
    return dtype.kind in {"S", "O", "U", "b"}


def approx_integer_ratio(values) -> float:
    if len(values) == 0:
        return 1.0
    near_int = 0
    for value in values:
        if abs(float(value) - round(float(value))) < 1e-6:
            near_int += 1
    return near_int / len(values)


def normalize_slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()
    return slug or "dataset"


def default_output_h5ad(dataset: Path) -> Path:
    return dataset.with_name(f"{dataset.stem}.companion.ready.h5ad")


def choose_groupby(obs_columns: list[str]) -> str | None:
    for name in GROUPBY_PRIORITY:
        if name in obs_columns:
            return name
    return None


def choose_composition(
    obs_columns: list[str], uniqueness: dict[str, int], categorical: set[str]
) -> str | None:
    for name in COMPOSITION_PRIORITY:
        if (
            name in categorical
            and name in obs_columns
            and 2 <= uniqueness.get(name, 0) <= 200
        ):
            return name
    fallback = [
        name
        for name, n_unique in uniqueness.items()
        if name in categorical and name != "_index" and 2 <= n_unique <= 80
    ]
    return sorted(fallback, key=lambda name: (uniqueness[name], name))[0] if fallback else None


def choose_analytics(
    obs_columns: list[str],
    uniqueness: dict[str, int],
    categorical: set[str],
    groupby: str | None,
    composition: str | None,
) -> list[str]:
    selected: list[str] = []
    seen: set[str] = set()

    def maybe_add(name: str) -> None:
        if (
            name in seen
            or name == groupby
            or name not in obs_columns
            or name not in categorical
        ):
            return
        n_unique = uniqueness.get(name, 0)
        if not (2 <= n_unique <= 160):
            return
        selected.append(name)
        seen.add(name)

    if composition:
        maybe_add(composition)
    for name in ANALYTICS_PRIORITY:
        maybe_add(name)
        if len(selected) >= 7:
            return selected

    extras = [
        name
        for name, n_unique in uniqueness.items()
        if name in categorical
        and name not in seen
        and name != groupby
        and 2 <= n_unique <= 40
    ]
    for name in sorted(extras, key=lambda item: (uniqueness[item], item)):
        maybe_add(name)
        if len(selected) >= 7:
            break
    return selected


def build_script_content(
    dataset: Path,
    output_h5ad: Path,
    groupby: str | None,
    composition: str | None,
    analytics: list[str],
    skip_normalized_layer: bool,
    normalize_from_counts: bool,
    viewer_neighbor_permutations: int | None,
) -> str:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f'INPUT="${{1:-{dataset}}}"',
        f'OUTPUT_H5AD="${{2:-{output_h5ad}}}"',
        "",
        'cargo run --release --offline -- prepare "$INPUT" \\',
        '  --output "$OUTPUT_H5AD" \\',
        "  --delaunay \\",
    ]
    if groupby:
        lines.append(f"  --groupby {groupby} \\")
    if composition:
        lines.append(f"  --composition-cell-type {composition} \\")
    if normalize_from_counts:
        lines.append("  --normalize-from layer:counts \\")
    if skip_normalized_layer:
        lines.append("  --skip-normalized-layer \\")
    lines.extend(
        [
            "  --skip-aggregation \\",
            "  --overwrite-derived \\",
            "  --persist-analytics-in-h5ad \\",
        ]
    )
    if analytics:
        lines.append(f"  --viewer-analytics-columns {','.join(analytics)} \\")
    if viewer_neighbor_permutations is not None:
        lines.append(
            f"  --viewer-neighbor-permutations {viewer_neighbor_permutations} \\"
        )
    lines.extend(
        [
            "  --skip-viewer-interaction-markers \\",
            "  --viewer-cluster-de-method t-test",
            "",
            'echo',
            'echo "Output:"',
            'echo "  H5AD: $OUTPUT_H5AD"',
            'ls -lh "$OUTPUT_H5AD"',
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    dataset = args.dataset.expanduser().resolve()
    repo_root = args.repo_root.expanduser().resolve()
    if not dataset.exists():
        raise SystemExit(f"Dataset not found: {dataset}")
    if not (repo_root / "scripts").is_dir():
        raise SystemExit(f"Expected repo root with scripts/: {repo_root}")

    with h5py.File(dataset, "r") as handle:
        obs = handle["obs"]
        obs_columns = sorted(obs.keys())
        uniqueness: dict[str, int] = {}
        categorical: set[str] = set()
        for name in obs_columns:
            try:
                values = read_obs_column(obs, name)
            except Exception:
                continue
            uniqueness[name] = len(set(values))
            if is_categorical_obs_column(obs, name):
                categorical.add(name)

        groupby = args.groupby or choose_groupby(obs_columns)
        composition = args.composition or choose_composition(
            obs_columns, uniqueness, categorical
        )
        analytics = (
            [item.strip() for item in args.analytics.split(",") if item.strip()]
            if args.analytics
            else choose_analytics(obs_columns, uniqueness, categorical, groupby, composition)
        )

        x_values = handle["X"]["data"][: min(256, handle["X"]["data"].shape[0])]
        x_integer_like = approx_integer_ratio(x_values) >= 0.95
        shape = tuple(int(v) for v in handle["X"].attrs.get("shape", (0, 0)))
        matrix_size = shape[0] * shape[1]
        has_counts_layer = "counts" in handle["layers"]
        skip_normalized_layer = (not x_integer_like) or matrix_size > 50_000_000
        normalize_from_counts = has_counts_layer and not skip_normalized_layer

    output_h5ad = (
        args.output_h5ad.expanduser().resolve()
        if args.output_h5ad
        else default_output_h5ad(dataset)
    )
    script_path = (
        args.script_path.expanduser().resolve()
        if args.script_path
        else repo_root
        / "scripts"
        / f"run_{normalize_slug(dataset.stem)}_prepare_h5ad.sh"
    )

    content = build_script_content(
        dataset=dataset,
        output_h5ad=output_h5ad,
        groupby=groupby,
        composition=composition,
        analytics=analytics,
        skip_normalized_layer=skip_normalized_layer,
        normalize_from_counts=normalize_from_counts,
        viewer_neighbor_permutations=args.viewer_neighbor_permutations,
    )
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(content, encoding="utf-8")

    summary = {
        "dataset": str(dataset),
        "script_path": str(script_path),
        "output_h5ad": str(output_h5ad),
        "groupby": groupby,
        "composition": composition,
        "analytics": analytics,
        "skip_normalized_layer": skip_normalized_layer,
        "normalize_from_counts": normalize_from_counts,
        "viewer_neighbor_permutations": args.viewer_neighbor_permutations,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
