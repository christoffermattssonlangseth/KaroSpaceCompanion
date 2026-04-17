#!/usr/bin/env python3
"""Create categorical tertile columns in obs for numeric H5AD columns."""

from __future__ import annotations

import argparse
import re
from typing import Iterable

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add categorical tertile columns for numeric obs columns in an H5AD file. "
            "Each output column is stored as obs/<name>/codes + obs/<name>/categories."
        )
    )
    parser.add_argument("h5ad", help="Path to input H5AD (modified in place).")
    parser.add_argument(
        "--columns",
        required=True,
        help="Comma-separated obs column names to bin.",
    )
    parser.add_argument(
        "--suffix",
        default="_tertile",
        help="Suffix for generated columns (default: _tertile).",
    )
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def sanitize_name(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]+", "_", name).strip("_")


def read_numeric_obs_column(obs: h5py.Group, name: str) -> np.ndarray:
    if name not in obs:
        raise KeyError(f"obs/{name} is missing")
    obj = obs[name]
    if isinstance(obj, h5py.Group):
        raise TypeError(f"obs/{name} is a group (expected numeric dataset)")
    arr = obj[()]
    if arr.ndim != 1:
        raise TypeError(f"obs/{name} is not 1D")
    if arr.dtype.kind not in {"f", "i", "u", "b"}:
        raise TypeError(f"obs/{name} is not numeric (dtype={arr.dtype})")
    return np.asarray(arr, dtype=np.float64)


def bin_tertiles(values: np.ndarray) -> tuple[np.ndarray, list[str], tuple[float, float] | None]:
    codes = np.full(values.shape[0], -1, dtype=np.int16)
    finite = np.isfinite(values)
    finite_values = values[finite]
    if finite_values.size == 0:
        return codes, ["low", "mid", "high"], None

    q1, q2 = np.quantile(finite_values, [1.0 / 3.0, 2.0 / 3.0])
    if np.isclose(q1, q2):
        min_v = float(np.min(finite_values))
        max_v = float(np.max(finite_values))
        if np.isclose(min_v, max_v):
            codes[finite] = 0
            return codes, ["constant"], (q1, q2)
        # Fallback to median split when tertiles collapse.
        median = float(np.quantile(finite_values, 0.5))
        finite_codes = np.where(finite_values <= median, 0, 1).astype(np.int16)
        codes[finite] = finite_codes
        return codes, ["low", "high"], (median, median)

    finite_codes = np.empty_like(finite_values, dtype=np.int16)
    finite_codes[finite_values <= q1] = 0
    finite_codes[(finite_values > q1) & (finite_values <= q2)] = 1
    finite_codes[finite_values > q2] = 2
    codes[finite] = finite_codes
    return codes, ["low", "mid", "high"], (float(q1), float(q2))


def write_categorical_column(
    obs: h5py.Group, target: str, codes: np.ndarray, categories: Iterable[str]
) -> None:
    if target in obs:
        del obs[target]
    group = obs.create_group(target)
    group.create_dataset("codes", data=codes, dtype=np.int16)
    str_dtype = h5py.string_dtype(encoding="utf-8")
    group.create_dataset(
        "categories", data=np.asarray(list(categories), dtype=object), dtype=str_dtype
    )


def main() -> None:
    args = parse_args()
    columns = parse_csv(args.columns)
    if not columns:
        raise ValueError("No columns were provided via --columns")

    with h5py.File(args.h5ad, "r+") as f:
        obs = f["obs"]
        for name in columns:
            values = read_numeric_obs_column(obs, name)
            target = f"{sanitize_name(name)}{args.suffix}"
            codes, categories, cutoffs = bin_tertiles(values)
            write_categorical_column(obs, target, codes, categories)
            missing_count = int(np.count_nonzero(codes == -1))
            if cutoffs is None:
                print(
                    f"{name} -> {target}: categories={categories} "
                    f"(all values missing, missing={missing_count})"
                )
            else:
                print(
                    f"{name} -> {target}: categories={categories} "
                    f"cutoffs={cutoffs} missing={missing_count}"
                )


if __name__ == "__main__":
    main()
