---
name: karospace-companion-sh
description: Create or update KaroSpaceCompanion dataset shell scripts for new `.h5ad` inputs. Use when the user asks for a new companion `.sh`, a dataset-specific `run_*_prepare_h5ad.sh`, or a reusable wrapper for a KaroSpaceCompanion dataset.
---

# KaroSpace Companion Sh

Use this skill when working in the KaroSpaceCompanion repo and the user wants a new dataset runner script.

## Workflow

1. Run `scripts/generate_runner.py DATASET_PATH --repo-root REPO_ROOT`.
2. Read the generated script and confirm the chosen `groupby`, composition column, analytics columns, and normalization strategy are sensible.
3. If the dataset or user context calls for a change, rerun with overrides:
   - `--groupby COL`
   - `--composition COL`
   - `--analytics col1,col2,...`
   - `--viewer-neighbor-permutations N`
   - `--script-path PATH`
   - `--output-h5ad PATH`
4. Set the script executable and validate with `bash -n`.
5. Do not run the full dataset unless the user explicitly asks for that run.

## Default Heuristics

- Prefer section boundaries in this order: `Sample.FOV`, `source_file`, `Sample_id`, `sample_id`, `sample`, `Sample`, `slide`, `orig.ident`, `FOV`.
- Prefer composition columns in this order: `Final_celltype`, `new_identity`, `majorCluster_final`, `tangram_cell_type`, `anno_L2`, `niches`, `areas`, `CellTypes`, `cell_type_E`, `Subclass_label`, `Class_id`, `cell_type`, `CellCharter_10`, `RNA_snn_res.2`, `seurat_clusters`.
- Keep analytics columns interpretable and bounded. Favor 4-8 actual categorical columns and exclude the `groupby` column.
- If `X` already looks transformed or the matrix is large, prefer `--skip-normalized-layer`.
- Always default to `--skip-aggregation --overwrite-derived --persist-analytics-in-h5ad --skip-viewer-interaction-markers --viewer-cluster-de-method t-test`.
- Only add `--viewer-neighbor-permutations` when the user explicitly wants forced z-score permutations on datasets that would otherwise auto-disable them.

## Bundled Script

- `scripts/generate_runner.py` inspects the H5AD and writes the shell script.
- It prints the chosen defaults after writing, so review them before finishing.
