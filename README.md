# KaroSpaceCompanion

Rust CLI for preparing KaroSpace-ready `.h5ad` files.

Current entrypoint:

```bash
cargo run -- prepare input.h5ad \
  --output enriched.h5ad \
  --radius 50 \
  --groupby sample \
  --composition-cell-type cell_type
```

Viewer-ready export:

```bash
cargo run -- prepare input.h5ad \
  --output enriched.h5ad \
  --radius 50 \
  --groupby sample \
  --composition-cell-type cell_type \
  --viewer-json viewer_precompute.json \
  --initial-color cell_type \
  --viewer-genes EPCAM,KRT19,COL1A1 \
  --viewer-analytics-columns cell_type \
  --skip-viewer-interaction-markers \
  --viewer-cluster-de-method t-test \
  --viewer-interaction-method t-test
```

For the Balo test dataset in `~/Downloads`, you can avoid shell copy/paste issues entirely with:

```bash
bash scripts/run_balo_bms_export.sh
```

The v1 `prepare` command:

- copies the input `.h5ad` to a new output path
- writes `obsm["spatial"]`
- writes `obsp["spatial_connectivities"]` and `obsp["spatial_distances"]`
- writes `layers["normalized"]` unless skipped
- writes `obsm["X_karo_agg"]` unless skipped
- optionally writes `obsm["X_karo_comp"]`
- optionally writes `obsm["X_karo_nmf"]` and `varm["X_karo_nmf_loadings"]`
- records provenance in `uns["karospace_companion"]`
- can persist downstream analytics directly into `uns["karospace_companion"]`, including:
  - `marker_genes_json`
  - `cluster_de_json`
  - `neighbor_stats_json`
  - `interaction_markers_json`
  - `gene_correlations_json`
  - `spatial_variable_genes_json`
  - `cluster_gene_means_json`
- optionally writes a KaroSpace viewer sidecar JSON with:
  - per-section coordinates, obs indices, colors, UMAP coordinates, gene vectors, and edges
  - `colors_meta`, `genes_meta`, `gene_encodings`, and categorical `metadata_filters`
  - `marker_genes`, `cluster_de`, `neighbor_stats`, `interaction_markers`
  - `gene_correlations`, `spatial_variable_genes`, `cluster_gene_means`

Notes:

- `--persist-analytics-in-h5ad` computes and stores the downstream analytics in the `.h5ad` even if you do not write a viewer JSON.
- `--viewer-json` uses `layers["normalized"]` for gene vectors and analytics when present, otherwise `X`.
- Omitting `--viewer-genes` exports all genes into the viewer JSON. That is complete, but it can become very large on real datasets.
- Dense section arrays can be packed to base64 for coordinates, colors, UMAP, and edges. Sparse gene sections can be packed with `ib64`/`vb64`.
- `--viewer-analytics-columns` limits the expensive precomputations (`cluster_de`, `neighbor_stats`, `interaction_markers`, `marker_genes`, `cluster_gene_means`) to the categorical columns you name. This is the practical way to run large datasets.
- `--skip-viewer-interaction-markers` disables the contact-vs-noncontact DE block entirely. For large datasets this is often the right default.
- `--viewer-cluster-de-method` and `--viewer-interaction-method` support `t-test` and `wilcoxon`. `t-test` is the fast path for large datasets.
- `--viewer-interaction-gene-limit` restricts the gene universe scanned for interaction markers to the top variable genes, which keeps large exports tractable.
