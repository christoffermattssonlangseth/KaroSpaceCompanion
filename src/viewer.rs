use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use base64::engine::general_purpose::STANDARD;
use base64::Engine;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use serde_json::{json, Map, Value};

use crate::anndata::{read_viewer_h5ad, CsrMatrix, ObsColumnData, ViewerAnnDataRead};
use crate::progress::ProgressReporter;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneEncodingMode {
    Auto,
    Dense,
    Sparse,
}

impl GeneEncodingMode {
    pub fn parse(value: &str) -> Result<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "dense" => Ok(Self::Dense),
            "sparse" => Ok(Self::Sparse),
            other => bail!("unsupported viewer gene encoding '{other}'"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ViewerPrecomputeConfig {
    pub output_json: Option<PathBuf>,
    pub initial_color: Option<String>,
    pub genes: Option<Vec<String>>,
    pub analytics_columns: Option<Vec<String>>,
    pub include_interaction_markers: bool,
    pub gene_encoding: GeneEncodingMode,
    pub gene_sparse_zero_threshold: f32,
    pub pack_arrays: bool,
    pub pack_arrays_min_len: usize,
    pub gene_sparse_pack_min_nnz: usize,
    pub gene_correlation_top_n: usize,
    pub cluster_means_n_genes: usize,
    pub spatial_variable_genes_n: usize,
    pub marker_genes_top_n: usize,
    pub cluster_de_method: DeMethod,
    pub cluster_de_top_n: usize,
    pub cluster_de_min_cells: usize,
    pub neighbor_stats_seed: u64,
    pub interaction_markers_method: DeMethod,
    pub interaction_markers_gene_limit: usize,
    pub interaction_markers_top_targets: usize,
    pub interaction_markers_top_genes: usize,
    pub interaction_markers_min_cells: usize,
    pub interaction_markers_min_neighbors: usize,
}

#[derive(Debug, Clone)]
pub struct CompanionAnalyticsBundle {
    pub analytics_columns: Vec<String>,
    pub marker_genes: Value,
    pub cluster_de: Value,
    pub neighbor_stats: Value,
    pub interaction_markers: Value,
    pub gene_correlations: Value,
    pub spatial_variable_genes: Value,
    pub cluster_gene_means: Value,
}

impl CompanionAnalyticsBundle {
    pub fn to_json_scalars(&self) -> Result<Vec<(String, String)>> {
        Ok(vec![
            (
                "marker_genes_json".to_string(),
                serde_json::to_string(&self.marker_genes).context("serializing marker_genes")?,
            ),
            (
                "cluster_de_json".to_string(),
                serde_json::to_string(&self.cluster_de).context("serializing cluster_de")?,
            ),
            (
                "neighbor_stats_json".to_string(),
                serde_json::to_string(&self.neighbor_stats)
                    .context("serializing neighbor_stats")?,
            ),
            (
                "interaction_markers_json".to_string(),
                serde_json::to_string(&self.interaction_markers)
                    .context("serializing interaction_markers")?,
            ),
            (
                "gene_correlations_json".to_string(),
                serde_json::to_string(&self.gene_correlations)
                    .context("serializing gene_correlations")?,
            ),
            (
                "spatial_variable_genes_json".to_string(),
                serde_json::to_string(&self.spatial_variable_genes)
                    .context("serializing spatial_variable_genes")?,
            ),
            (
                "cluster_gene_means_json".to_string(),
                serde_json::to_string(&self.cluster_gene_means)
                    .context("serializing cluster_gene_means")?,
            ),
        ])
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeMethod {
    TTest,
    Wilcoxon,
}

impl DeMethod {
    pub fn parse(value: &str) -> Result<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "t-test" | "ttest" | "t_test" => Ok(Self::TTest),
            "wilcoxon" | "rank-sum" | "ranksum" => Ok(Self::Wilcoxon),
            other => bail!("unsupported DE method '{other}'"),
        }
    }
}

#[derive(Clone)]
struct GroupSummary {
    n_cells: usize,
    sums: Vec<f64>,
    sums_sq: Vec<f64>,
    nnz: Vec<usize>,
}

impl GroupSummary {
    fn zeros(n_genes: usize) -> Self {
        Self {
            n_cells: 0,
            sums: vec![0.0; n_genes],
            sums_sq: vec![0.0; n_genes],
            nnz: vec![0; n_genes],
        }
    }
}

struct CategoricalColumnSummary {
    group_summaries: Vec<GroupSummary>,
    total: GroupSummary,
}

struct NeighborStatsContext {
    codes: Vec<Option<usize>>,
    counts: Vec<Vec<f64>>,
    zscores: Option<Vec<Vec<f64>>>,
    n_cells: Vec<usize>,
}

struct DifferentialExpressionResult {
    genes: Vec<String>,
    logfoldchanges: Vec<Option<f64>>,
    pvals_adj: Vec<Option<f64>>,
    scores: Vec<Option<f64>>,
    pct_a: Vec<Option<f64>>,
    pct_b: Vec<Option<f64>>,
}

struct GeneExportMeta {
    index: usize,
    name: String,
    encoding: GeneEncodingMode,
    vmin: f64,
    vmax: f64,
}

pub fn export_viewer_precompute_json(
    h5ad_path: &Path,
    groupby: Option<&str>,
    config: &ViewerPrecomputeConfig,
) -> Result<PathBuf> {
    let reporter = ProgressReporter::new("viewer");
    let stage = reporter.stage("Loading prepared h5ad");
    let adata = read_viewer_h5ad(h5ad_path)?;
    stage.finish(format!(
        "cells={}, genes={}, expr_source={}",
        adata.obs_names.len(),
        adata.var_names.len(),
        adata.expression_source
    ));
    let analytics = compute_companion_analytics_from_loaded(&adata, groupby, config, &reporter)?;
    export_viewer_precompute_json_with_analytics(h5ad_path, groupby, config, &analytics)
}

pub fn export_viewer_precompute_json_with_analytics(
    h5ad_path: &Path,
    groupby: Option<&str>,
    config: &ViewerPrecomputeConfig,
    analytics: &CompanionAnalyticsBundle,
) -> Result<PathBuf> {
    let reporter = ProgressReporter::new("viewer");
    let stage = reporter.stage("Loading prepared h5ad");
    let adata = read_viewer_h5ad(h5ad_path)?;
    stage.finish(format!(
        "cells={}, genes={}, expr_source={}",
        adata.obs_names.len(),
        adata.var_names.len(),
        adata.expression_source
    ));
    let output_json = config
        .output_json
        .as_ref()
        .context("viewer JSON output path was not provided")?;
    let payload = build_payload(&adata, groupby, config, analytics, &reporter)?;
    if let Some(parent) = output_json.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("creating {:?}", parent))?;
        }
    }
    let stage = reporter.stage("Writing viewer JSON");
    let encoded = serde_json::to_vec(&payload).context("serializing viewer precompute JSON")?;
    fs::write(output_json, encoded)
        .with_context(|| format!("writing {:?}", output_json))?;
    stage.finish(output_json.display().to_string());
    Ok(output_json.clone())
}

pub fn compute_companion_analytics(
    h5ad_path: &Path,
    groupby: Option<&str>,
    config: &ViewerPrecomputeConfig,
) -> Result<CompanionAnalyticsBundle> {
    let reporter = ProgressReporter::new("analytics");
    let stage = reporter.stage("Loading prepared h5ad");
    let adata = read_viewer_h5ad(h5ad_path)?;
    stage.finish(format!(
        "cells={}, genes={}, expr_source={}",
        adata.obs_names.len(),
        adata.var_names.len(),
        adata.expression_source
    ));
    compute_companion_analytics_from_loaded(&adata, groupby, config, &reporter)
}

fn build_payload(
    adata: &ViewerAnnDataRead,
    groupby: Option<&str>,
    config: &ViewerPrecomputeConfig,
    analytics: &CompanionAnalyticsBundle,
    reporter: &ProgressReporter,
) -> Result<Value> {
    let stage = reporter.stage("Preparing core payload");
    let group_values = resolve_group_values(adata, groupby);
    let sections = build_sections(&group_values);
    let initial_color = resolve_initial_color(adata, groupby, config.initial_color.as_deref())?;

    let available_colors = adata.obs_order.clone();
    let metadata_columns: Vec<String> = adata
        .obs_order
        .iter()
        .filter(|name| Some(name.as_str()) != groupby)
        .filter(|name| {
            adata
                .obs_columns
                .get(*name)
                .and_then(ObsColumnData::categories)
                .is_some()
        })
        .cloned()
        .collect();
    let selected_genes = resolve_selected_genes(adata, config.genes.as_deref())?;
    let gene_export_meta = build_gene_export_meta(adata, &selected_genes, config);

    let colors_meta = build_colors_meta(adata)?;
    let metadata_filters = build_metadata_filters(adata, &metadata_columns);
    let sections_json = build_sections_json(
        adata,
        &sections,
        &available_colors,
        &metadata_columns,
        &gene_export_meta,
        config,
    )?;
    stage.finish(format!(
        "sections={}, colors={}, genes={}",
        sections.len(),
        available_colors.len(),
        gene_export_meta.len()
    ));

    let genes_meta = gene_export_meta
        .iter()
        .map(|gene| {
            (
                gene.name.clone(),
                json!({
                    "vmin": gene.vmin,
                    "vmax": gene.vmax,
                }),
            )
        })
        .collect::<Map<String, Value>>();
    let gene_encodings = gene_export_meta
        .iter()
        .map(|gene| {
            let mode = match gene.encoding {
                GeneEncodingMode::Auto => "auto",
                GeneEncodingMode::Dense => "dense",
                GeneEncodingMode::Sparse => "sparse",
            };
            (gene.name.clone(), Value::String(mode.to_string()))
        })
        .collect::<Map<String, Value>>();

    let umap_bounds = adata.umap.as_ref().map(compute_umap_bounds).unwrap_or(Value::Null);

    Ok(json!({
        "initial_color": initial_color,
        "groupby": groupby.unwrap_or("all"),
        "colors_meta": Value::Object(colors_meta),
        "genes_meta": Value::Object(genes_meta),
        "gene_encodings": Value::Object(gene_encodings),
        "metadata_filters": Value::Object(metadata_filters),
        "n_sections": sections.len(),
        "total_cells": adata.obs_names.len(),
        "sections": sections_json,
        "available_colors": available_colors,
        "available_genes": gene_export_meta.iter().map(|gene| gene.name.clone()).collect::<Vec<_>>(),
        "marker_genes": analytics.marker_genes.clone(),
        "cluster_de": analytics.cluster_de.clone(),
        "has_umap": adata.umap.is_some(),
        "umap_bounds": umap_bounds,
        "has_neighbors": adata.graph.is_some(),
        "neighbors_key": adata.graph_key,
        "neighbor_stats": analytics.neighbor_stats.clone(),
        "interaction_markers": analytics.interaction_markers.clone(),
        "gene_correlations": analytics.gene_correlations.clone(),
        "spatial_variable_genes": analytics.spatial_variable_genes.clone(),
        "cluster_gene_means": analytics.cluster_gene_means.clone(),
    }))
}

fn compute_companion_analytics_from_loaded(
    adata: &ViewerAnnDataRead,
    groupby: Option<&str>,
    config: &ViewerPrecomputeConfig,
    reporter: &ProgressReporter,
) -> Result<CompanionAnalyticsBundle> {
    let categorical_columns = collect_categorical_columns(adata, groupby);
    let analytics_columns =
        resolve_analytics_columns(&categorical_columns, config.analytics_columns.as_deref())?;
    let analytics_column_names = analytics_columns
        .iter()
        .map(|(name, _, _)| name.clone())
        .collect::<Vec<_>>();
    let selected_genes = resolve_selected_genes(adata, config.genes.as_deref())?;
    let top_variable_gene_indices =
        select_top_variable_gene_indices(&adata.expression, config.cluster_means_n_genes);
    let cluster_gene_means = if !analytics_columns.is_empty() && !top_variable_gene_indices.is_empty()
    {
        compute_cluster_gene_means(
            adata,
            &analytics_columns,
            &top_variable_gene_indices,
            Some(reporter),
        )
    } else {
        Value::Null
    };

    let gene_correlations = if config.gene_correlation_top_n > 0 && !selected_genes.is_empty() {
        compute_gene_correlations(
            &adata.expression,
            &selected_genes,
            &adata.var_names,
            config.gene_correlation_top_n,
            Some(reporter),
        )
    } else {
        Value::Object(Map::new())
    };

    let spatial_variable_genes = if config.spatial_variable_genes_n > 0 {
        if let Some(graph) = &adata.graph {
            let svg_genes = select_top_variable_gene_indices(
                &adata.expression,
                config.spatial_variable_genes_n,
            );
            compute_spatial_variable_genes(
                &adata.expression,
                graph,
                &adata.var_names,
                &svg_genes,
                Some(reporter),
            )
        } else {
            Value::Array(Vec::new())
        }
    } else {
        Value::Array(Vec::new())
    };

    let marker_genes = compute_marker_genes(
        &adata.expression,
        &adata.var_names,
        &analytics_columns,
        config.marker_genes_top_n,
        Some(reporter),
    );
    let cluster_de = compute_cluster_de(
        &adata.expression,
        &adata.var_names,
        &analytics_columns,
        config.cluster_de_top_n,
        config.cluster_de_min_cells,
        config.cluster_de_method,
        Some(reporter),
    );

    let (neighbor_stats, neighbor_contexts) = if let Some(graph) = &adata.graph {
        compute_neighbor_stats(
            graph,
            &analytics_columns,
            auto_neighbor_permutations(adata.obs_names.len()),
            config.neighbor_stats_seed,
            Some(reporter),
        )
    } else {
        (Value::Object(Map::new()), HashMap::new())
    };

    let interaction_gene_indices = if config.include_interaction_markers {
        if config.interaction_markers_gene_limit == 0 {
            (0..adata.expression.ncols()).collect::<Vec<_>>()
        } else {
            select_top_variable_gene_indices(
                &adata.expression,
                config.interaction_markers_gene_limit,
            )
        }
    } else {
        Vec::new()
    };
    let interaction_markers = if config.include_interaction_markers && adata.graph.is_some() {
        compute_interaction_markers(
            &adata.expression,
            &adata.var_names,
            &analytics_columns,
            &neighbor_contexts,
            graph_neighbors_from_csr(adata.graph.as_ref().expect("checked graph")),
            &interaction_gene_indices,
            config.interaction_markers_top_targets,
            config.interaction_markers_top_genes,
            config.interaction_markers_min_cells,
            config.interaction_markers_min_neighbors,
            config.interaction_markers_method,
            Some(reporter),
        )
    } else {
        Value::Object(Map::new())
    };

    Ok(CompanionAnalyticsBundle {
        analytics_columns: analytics_column_names,
        marker_genes,
        cluster_de,
        neighbor_stats,
        interaction_markers,
        gene_correlations,
        spatial_variable_genes,
        cluster_gene_means,
    })
}

fn resolve_group_values(adata: &ViewerAnnDataRead, groupby: Option<&str>) -> Vec<String> {
    if let Some(column) = groupby {
        if let Some(values) = adata.obs_columns.get(column) {
            return values.as_strings();
        }
    }
    vec!["all".to_string(); adata.obs_names.len()]
}

fn build_sections(group_values: &[String]) -> Vec<(String, Vec<usize>)> {
    let mut by_section: HashMap<String, Vec<usize>> = HashMap::new();
    let mut order = Vec::new();
    for (index, value) in group_values.iter().enumerate() {
        if !by_section.contains_key(value) {
            order.push(value.clone());
            by_section.insert(value.clone(), Vec::new());
        }
        by_section
            .get_mut(value)
            .expect("section inserted")
            .push(index);
    }
    order
        .into_iter()
        .map(|section| {
            let indices = by_section.remove(&section).unwrap_or_default();
            (section, indices)
        })
        .collect()
}

fn resolve_initial_color(
    adata: &ViewerAnnDataRead,
    groupby: Option<&str>,
    requested: Option<&str>,
) -> Result<String> {
    if let Some(color) = requested {
        if adata.obs_columns.contains_key(color) {
            return Ok(color.to_string());
        }
        bail!("requested initial color '{color}' not found in obs");
    }
    if let Some((name, _)) = adata
        .obs_order
        .iter()
        .filter(|name| Some(name.as_str()) != groupby)
        .find_map(|name| adata.obs_columns.get(name).map(|col| (name, col)))
    {
        if !adata
            .obs_columns
            .get(name)
            .map(|col| col.is_continuous())
            .unwrap_or(false)
        {
            return Ok(name.clone());
        }
    }
    adata
        .obs_order
        .iter()
        .find(|name| Some(name.as_str()) != groupby)
        .cloned()
        .context("no obs columns available for initial color")
}

fn resolve_selected_genes(
    adata: &ViewerAnnDataRead,
    requested: Option<&[String]>,
) -> Result<Vec<usize>> {
    match requested {
        Some(names) => {
            let lookup: HashMap<&str, usize> = adata
                .var_names
                .iter()
                .enumerate()
                .map(|(index, name)| (name.as_str(), index))
                .collect();
            let mut indices = Vec::new();
            for name in names {
                if let Some(index) = lookup.get(name.as_str()) {
                    indices.push(*index);
                }
            }
            if indices.is_empty() {
                bail!("none of the requested viewer genes were found in var_names");
            }
            Ok(indices)
        }
        None => Ok((0..adata.var_names.len()).collect()),
    }
}

fn build_gene_export_meta(
    adata: &ViewerAnnDataRead,
    gene_indices: &[usize],
    config: &ViewerPrecomputeConfig,
) -> Vec<GeneExportMeta> {
    gene_indices
        .iter()
        .map(|&index| {
            let column = adata.expression.column(index);
            let mut finite_min = f64::INFINITY;
            let mut finite_max = f64::NEG_INFINITY;
            let mut nonzero = 0usize;
            for value in column.iter().copied() {
                if value.is_finite() {
                    finite_min = finite_min.min(value as f64);
                    finite_max = finite_max.max(value as f64);
                    if value != 0.0 {
                        nonzero += 1;
                    }
                }
            }
            if !finite_min.is_finite() {
                finite_min = 0.0;
                finite_max = 0.0;
            }
            let zero_frac = if column.is_empty() {
                1.0
            } else {
                1.0 - (nonzero as f32 / column.len() as f32)
            };
            let encoding = match config.gene_encoding {
                GeneEncodingMode::Auto => {
                    if zero_frac >= config.gene_sparse_zero_threshold {
                        GeneEncodingMode::Sparse
                    } else {
                        GeneEncodingMode::Dense
                    }
                }
                mode => mode,
            };
            GeneExportMeta {
                index,
                name: adata.var_names[index].clone(),
                encoding,
                vmin: finite_min,
                vmax: finite_max,
            }
        })
        .collect()
}

fn build_colors_meta(adata: &ViewerAnnDataRead) -> Result<Map<String, Value>> {
    let mut meta = Map::new();
    for name in &adata.obs_order {
        let column = adata
            .obs_columns
            .get(name)
            .with_context(|| format!("missing cached obs column '{name}'"))?;
        let value = match column {
            ObsColumnData::Categorical { categories, .. } => json!({
                "is_continuous": false,
                "categories": categories,
                "vmin": Value::Null,
                "vmax": Value::Null,
            }),
            ObsColumnData::Numeric { values } => {
                let finite: Vec<f64> = values.iter().filter_map(|value| *value).collect();
                let vmin = finite.iter().copied().reduce(f64::min).unwrap_or(0.0);
                let vmax = finite.iter().copied().reduce(f64::max).unwrap_or(0.0);
                json!({
                    "is_continuous": true,
                    "categories": Value::Null,
                    "vmin": vmin,
                    "vmax": vmax,
                })
            }
        };
        meta.insert(name.clone(), value);
    }
    Ok(meta)
}

fn build_metadata_filters(
    adata: &ViewerAnnDataRead,
    metadata_columns: &[String],
) -> Map<String, Value> {
    let mut filters = Map::new();
    for column_name in metadata_columns {
        if let Some(column) = adata.obs_columns.get(column_name) {
            let mut values = column
                .as_strings()
                .into_iter()
                .filter(|value| !value.is_empty())
                .collect::<Vec<_>>();
            values.sort();
            values.dedup();
            filters.insert(
                column_name.clone(),
                Value::Array(values.into_iter().map(Value::String).collect()),
            );
        }
    }
    filters
}

fn build_sections_json(
    adata: &ViewerAnnDataRead,
    sections: &[(String, Vec<usize>)],
    available_colors: &[String],
    metadata_columns: &[String],
    genes: &[GeneExportMeta],
    config: &ViewerPrecomputeConfig,
) -> Result<Vec<Value>> {
    let graph_neighbors = adata.graph.as_ref().map(graph_neighbors_from_csr);
    let mut out = Vec::new();
    for (section_id, indices) in sections {
        let metadata = build_section_metadata(adata, metadata_columns, indices);
        let n_cells = indices.len();
        let coords_x = indices
            .iter()
            .map(|index| adata.coordinates[[*index, 0]] as f32)
            .collect::<Vec<_>>();
        let coords_y = indices
            .iter()
            .map(|index| adata.coordinates[[*index, 1]] as f32)
            .collect::<Vec<_>>();
        let bounds = json!({
            "xmin": coords_x.iter().copied().reduce(f32::min).unwrap_or(0.0),
            "xmax": coords_x.iter().copied().reduce(f32::max).unwrap_or(0.0),
            "ymin": coords_y.iter().copied().reduce(f32::min).unwrap_or(0.0),
            "ymax": coords_y.iter().copied().reduce(f32::max).unwrap_or(0.0),
        });

        let mut section = Map::new();
        section.insert("id".to_string(), Value::String(section_id.clone()));
        section.insert("metadata".to_string(), Value::Object(metadata));
        section.insert("n_cells".to_string(), json!(n_cells));
        section.insert("rotation_deg".to_string(), json!(0.0));
        section.insert("bounds".to_string(), bounds);

        if config.pack_arrays && n_cells >= config.pack_arrays_min_len {
            section.insert("x".to_string(), Value::Null);
            section.insert("y".to_string(), Value::Null);
            section.insert("xb64".to_string(), Value::String(encode_f32_b64(&coords_x)));
            section.insert("yb64".to_string(), Value::String(encode_f32_b64(&coords_y)));
            section.insert(
                "obs_idx".to_string(),
                Value::Null,
            );
            section.insert(
                "obs_idxb64".to_string(),
                Value::String(encode_u32_b64(
                    &indices.iter().map(|index| *index as u32).collect::<Vec<_>>(),
                )),
            );
        } else {
            section.insert("x".to_string(), f32_vec_to_json(&coords_x));
            section.insert("y".to_string(), f32_vec_to_json(&coords_y));
            section.insert(
                "xb64".to_string(),
                Value::Null,
            );
            section.insert(
                "yb64".to_string(),
                Value::Null,
            );
            section.insert(
                "obs_idx".to_string(),
                Value::Array(indices.iter().map(|index| json!(*index)).collect()),
            );
            section.insert("obs_idxb64".to_string(), Value::Null);
        }

        if let Some(umap) = &adata.umap {
            let umap_x = indices
                .iter()
                .map(|index| umap[[*index, 0]] as f32)
                .collect::<Vec<_>>();
            let umap_y = indices
                .iter()
                .map(|index| umap[[*index, 1]] as f32)
                .collect::<Vec<_>>();
            if config.pack_arrays && n_cells >= config.pack_arrays_min_len {
                section.insert("umap_x".to_string(), Value::Null);
                section.insert("umap_y".to_string(), Value::Null);
                section.insert("umap_xb64".to_string(), Value::String(encode_f32_b64(&umap_x)));
                section.insert("umap_yb64".to_string(), Value::String(encode_f32_b64(&umap_y)));
            } else {
                section.insert("umap_x".to_string(), f32_vec_to_json(&umap_x));
                section.insert("umap_y".to_string(), f32_vec_to_json(&umap_y));
                section.insert("umap_xb64".to_string(), Value::Null);
                section.insert("umap_yb64".to_string(), Value::Null);
            }
        } else {
            section.insert("umap_x".to_string(), Value::Null);
            section.insert("umap_y".to_string(), Value::Null);
            section.insert("umap_xb64".to_string(), Value::Null);
            section.insert("umap_yb64".to_string(), Value::Null);
        }

        let (colors, colors_b64) = build_section_colors(
            adata,
            available_colors,
            indices,
            config.pack_arrays && n_cells >= config.pack_arrays_min_len,
        )?;
        section.insert("colors".to_string(), Value::Object(colors));
        section.insert("colors_b64".to_string(), Value::Object(colors_b64));

        let (genes_dense, genes_sparse) = build_section_genes(adata, genes, indices, config);
        section.insert("genes".to_string(), Value::Object(genes_dense));
        section.insert("genes_sparse".to_string(), Value::Object(genes_sparse));

        let (edges, edges_b64) = build_section_edges(graph_neighbors.as_deref(), indices, config);
        section.insert("edges".to_string(), edges);
        section.insert("edges_b64".to_string(), edges_b64);

        out.push(Value::Object(section));
    }
    Ok(out)
}

fn build_section_metadata(
    adata: &ViewerAnnDataRead,
    metadata_columns: &[String],
    indices: &[usize],
) -> Map<String, Value> {
    let mut meta = Map::new();
    for column_name in metadata_columns {
        if let Some(column) = adata.obs_columns.get(column_name) {
            let values = column.as_strings();
            let mut unique = indices
                .iter()
                .filter_map(|index| {
                    let value = values.get(*index)?.clone();
                    if value.is_empty() {
                        None
                    } else {
                        Some(value)
                    }
                })
                .collect::<Vec<_>>();
            unique.sort();
            unique.dedup();
            let value = match unique.len() {
                0 => String::new(),
                1 => unique[0].clone(),
                _ => "mixed".to_string(),
            };
            meta.insert(column_name.clone(), Value::String(value));
        }
    }
    meta
}

fn build_section_colors(
    adata: &ViewerAnnDataRead,
    color_names: &[String],
    indices: &[usize],
    pack: bool,
) -> Result<(Map<String, Value>, Map<String, Value>)> {
    let mut colors = Map::new();
    let mut colors_b64 = Map::new();
    for name in color_names {
        let column = adata
            .obs_columns
            .get(name)
            .with_context(|| format!("missing cached color column '{name}'"))?;
        let values = if column.is_continuous() {
            column
                .numeric_values_f32()
                .expect("checked numeric column")
        } else {
            column
                .categorical_codes_f32()
                .expect("checked categorical column")
        };
        let section_values = indices.iter().map(|index| values[*index]).collect::<Vec<_>>();
        if pack {
            colors_b64.insert(name.clone(), Value::String(encode_f32_b64(&section_values)));
        } else {
            colors.insert(name.clone(), f32_vec_to_json(&section_values));
        }
    }
    Ok((colors, colors_b64))
}

fn build_section_genes(
    adata: &ViewerAnnDataRead,
    genes: &[GeneExportMeta],
    indices: &[usize],
    config: &ViewerPrecomputeConfig,
) -> (Map<String, Value>, Map<String, Value>) {
    let mut dense = Map::new();
    let mut sparse = Map::new();
    for gene in genes {
        let values = indices
            .iter()
            .map(|index| adata.expression[[*index, gene.index]])
            .collect::<Vec<_>>();
        match gene.encoding {
            GeneEncodingMode::Dense => {
                dense.insert(gene.name.clone(), f32_vec_to_json(&values));
            }
            GeneEncodingMode::Sparse | GeneEncodingMode::Auto => {
                sparse.insert(
                    gene.name.clone(),
                    serialize_sparse_gene_section(&values, config.gene_sparse_pack_min_nnz),
                );
            }
        }
    }
    (dense, sparse)
}

fn serialize_sparse_gene_section(values: &[f32], pack_min_nnz: usize) -> Value {
    let mut indices = Vec::new();
    let mut nonzero = Vec::new();
    let mut nan = Vec::new();
    for (index, value) in values.iter().copied().enumerate() {
        if value.is_nan() {
            nan.push(index as u32);
        } else if value != 0.0 {
            indices.push(index as u32);
            nonzero.push(value);
        }
    }
    let mut payload = Map::new();
    if indices.len() >= pack_min_nnz {
        payload.insert("ib64".to_string(), Value::String(encode_u32_b64(&indices)));
        payload.insert("vb64".to_string(), Value::String(encode_f32_b64(&nonzero)));
    } else {
        payload.insert(
            "i".to_string(),
            Value::Array(indices.into_iter().map(|index| json!(index)).collect()),
        );
        payload.insert("v".to_string(), f32_vec_to_json(&nonzero));
    }
    if !nan.is_empty() {
        payload.insert(
            "nan".to_string(),
            Value::Array(nan.into_iter().map(|index| json!(index)).collect()),
        );
    }
    Value::Object(payload)
}

fn build_section_edges(
    neighbors: Option<&[Vec<(usize, f32)>]>,
    indices: &[usize],
    config: &ViewerPrecomputeConfig,
) -> (Value, Value) {
    let Some(neighbors) = neighbors else {
        return (Value::Array(Vec::new()), Value::Null);
    };

    let local_lookup: HashMap<usize, usize> = indices
        .iter()
        .enumerate()
        .map(|(local, global)| (*global, local))
        .collect();
    let mut edge_pairs = Vec::<(u32, u32)>::new();
    for (local_i, global_i) in indices.iter().enumerate() {
        for &(global_j, _) in &neighbors[*global_i] {
            let Some(&local_j) = local_lookup.get(&global_j) else {
                continue;
            };
            if local_j <= local_i {
                continue;
            }
            edge_pairs.push((local_i as u32, local_j as u32));
        }
    }
    if config.pack_arrays && indices.len() >= config.pack_arrays_min_len {
        let mut flat = Vec::with_capacity(edge_pairs.len() * 2);
        for (i, j) in edge_pairs {
            flat.push(i);
            flat.push(j);
        }
        (Value::Null, Value::String(encode_u32_b64(&flat)))
    } else {
        (
            Value::Array(
                edge_pairs
                    .into_iter()
                    .map(|(i, j)| Value::Array(vec![json!(i), json!(j)]))
                    .collect(),
            ),
            Value::Null,
        )
    }
}

fn collect_categorical_columns(
    adata: &ViewerAnnDataRead,
    groupby: Option<&str>,
) -> Vec<(String, Vec<String>, Vec<String>)> {
    adata
        .obs_order
        .iter()
        .filter(|name| Some(name.as_str()) != groupby)
        .filter_map(|name| {
            let column = adata.obs_columns.get(name)?;
            let categories = column.categories()?.to_vec();
            Some((name.clone(), column.as_strings(), categories))
        })
        .collect()
}

fn resolve_analytics_columns(
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    requested: Option<&[String]>,
) -> Result<Vec<(String, Vec<String>, Vec<String>)>> {
    if let Some(names) = requested {
        let lookup: HashMap<&str, &(String, Vec<String>, Vec<String>)> = categorical_columns
            .iter()
            .map(|column| (column.0.as_str(), column))
            .collect();
        let mut selected = Vec::with_capacity(names.len());
        for name in names {
            let column = lookup
                .get(name.as_str())
                .with_context(|| format!("viewer analytics column '{name}' is missing or not categorical"))?;
            selected.push((*column).clone());
        }
        return Ok(selected);
    }
    Ok(categorical_columns.to_vec())
}

fn select_top_variable_gene_indices(
    matrix: &ndarray::Array2<f32>,
    n: usize,
) -> Vec<usize> {
    if n == 0 {
        return Vec::new();
    }
    let mut scored = Vec::with_capacity(matrix.ncols());
    for gene_idx in 0..matrix.ncols() {
        let column = matrix.column(gene_idx);
        if column.is_empty() {
            continue;
        }
        let mean = column.iter().copied().map(f64::from).sum::<f64>() / column.len() as f64;
        let mean_sq = column
            .iter()
            .copied()
            .map(|value| {
                let v = f64::from(value);
                v * v
            })
            .sum::<f64>()
            / column.len() as f64;
        let variance = mean_sq - mean * mean;
        scored.push((gene_idx, variance));
    }
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    scored.into_iter().take(n).map(|(index, _)| index).collect()
}

fn compute_cluster_gene_means(
    adata: &ViewerAnnDataRead,
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    gene_indices: &[usize],
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing cluster gene means"));
    let genes = gene_indices
        .iter()
        .map(|index| adata.var_names[*index].clone())
        .collect::<Vec<_>>();
    let background = gene_indices
        .iter()
        .map(|index| {
            let column = adata.expression.column(*index);
            if column.is_empty() {
                0.0
            } else {
                column.iter().copied().map(f64::from).sum::<f64>() / column.len() as f64
            }
        })
        .collect::<Vec<_>>();

    let mut columns_json = Map::new();
    for (column_idx, (column_name, values, categories)) in categorical_columns.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("columns", column_idx + 1, categorical_columns.len());
            stage.note(format!(
                "column {}/{}: {} ({} categories)",
                column_idx + 1,
                categorical_columns.len(),
                column_name,
                categories.len()
            ));
        }
        let mut means = Map::new();
        for (category_idx, category) in categories.iter().enumerate() {
            if let Some(stage) = stage.as_mut() {
                stage.progress("categories", category_idx + 1, categories.len());
            }
            let cell_indices = values
                .iter()
                .enumerate()
                .filter_map(|(idx, value)| if value == category { Some(idx) } else { None })
                .collect::<Vec<_>>();
            let gene_means = gene_indices
                .iter()
                .map(|gene_idx| mean_for_indices(&adata.expression, *gene_idx, &cell_indices))
                .collect::<Vec<_>>();
            means.insert(category.clone(), numeric_vec_to_json(&gene_means));
        }
        columns_json.insert(
            column_name.clone(),
            json!({
                "categories": categories,
                "means": Value::Object(means),
                "background": numeric_vec_to_json(&background),
            }),
        );
    }

    let result = json!({
        "genes": genes,
        "columns": Value::Object(columns_json),
    });
    if let Some(stage) = stage {
        stage.finish(format!(
            "columns={}, genes={}",
            categorical_columns.len(),
            gene_indices.len()
        ));
    }
    result
}

fn compute_gene_correlations(
    matrix: &ndarray::Array2<f32>,
    gene_indices: &[usize],
    gene_names: &[String],
    top_n: usize,
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing gene correlations"));
    let mut centered = Vec::with_capacity(gene_indices.len());
    let mut norms = Vec::with_capacity(gene_indices.len());
    for &gene_idx in gene_indices {
        let column = matrix.column(gene_idx);
        let mean = if column.is_empty() {
            0.0
        } else {
            column.iter().copied().map(f64::from).sum::<f64>() / column.len() as f64
        };
        let values = column
            .iter()
            .copied()
            .map(|value| f64::from(value) - mean)
            .collect::<Vec<_>>();
        let norm = values.iter().map(|value| value * value).sum::<f64>().sqrt();
        centered.push(values);
        norms.push(norm);
    }

    let mut out = Map::new();
    for (i, &gene_idx) in gene_indices.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("genes", i + 1, gene_indices.len());
        }
        if top_n == 0 || norms[i] == 0.0 {
            out.insert(gene_names[gene_idx].clone(), Value::Array(Vec::new()));
            continue;
        }
        let mut scores = Vec::new();
        for (j, &other_gene_idx) in gene_indices.iter().enumerate() {
            if i == j || norms[j] == 0.0 {
                continue;
            }
            let dot = centered[i]
                .iter()
                .zip(centered[j].iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let r = dot / (norms[i] * norms[j]);
            if r > 0.0 && r.is_finite() {
                scores.push((other_gene_idx, r));
            }
        }
        scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        let payload = scores
            .into_iter()
            .take(top_n)
            .map(|(other_gene_idx, r)| {
                json!({
                    "gene": gene_names[other_gene_idx],
                    "r": r,
                })
            })
            .collect::<Vec<_>>();
        out.insert(gene_names[gene_idx].clone(), Value::Array(payload));
    }
    if let Some(stage) = stage {
        stage.finish(format!("genes={}, top_n={}", gene_indices.len(), top_n));
    }
    Value::Object(out)
}

fn compute_spatial_variable_genes(
    matrix: &ndarray::Array2<f32>,
    graph: &CsrMatrix,
    gene_names: &[String],
    gene_indices: &[usize],
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing spatial variable genes"));
    let neighbors = graph_neighbors_from_csr(graph);
    let row_sums = neighbors
        .iter()
        .map(|row| row.iter().map(|(_, weight)| f64::from(*weight)).sum::<f64>())
        .collect::<Vec<_>>();
    let s0 = row_sums
        .iter()
        .map(|sum| if *sum > 0.0 { 1.0 } else { 0.0 })
        .sum::<f64>();
    if s0 == 0.0 {
        return Value::Array(Vec::new());
    }

    let n = matrix.nrows() as f64;
    let mut scored = Vec::new();
    for (index, &gene_idx) in gene_indices.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("genes", index + 1, gene_indices.len());
        }
        let column = matrix.column(gene_idx);
        let mean = column.iter().copied().map(f64::from).sum::<f64>() / column.len() as f64;
        let z = column
            .iter()
            .copied()
            .map(|value| f64::from(value) - mean)
            .collect::<Vec<_>>();
        let denom = z.iter().map(|value| value * value).sum::<f64>();
        if denom == 0.0 {
            continue;
        }
        let mut wz = vec![0.0f64; z.len()];
        for (row_idx, row) in neighbors.iter().enumerate() {
            let row_sum = row_sums[row_idx];
            if row_sum <= 0.0 {
                continue;
            }
            let mut acc = 0.0;
            for &(neighbor_idx, weight) in row {
                acc += (f64::from(weight) / row_sum) * z[neighbor_idx];
            }
            wz[row_idx] = acc;
        }
        let numerator = n
            * z.iter()
                .zip(wz.iter())
                .map(|(lhs, rhs)| lhs * rhs)
                .sum::<f64>();
        let i_value = (numerator / (s0 * denom)).clamp(-1.0, 1.0);
        scored.push((gene_idx, i_value));
    }
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    let result = Value::Array(
        scored
            .into_iter()
            .map(|(gene_idx, i_value)| {
                json!({
                    "gene": gene_names[gene_idx],
                    "I": i_value,
                })
            })
            .collect(),
    );
    if let Some(stage) = stage {
        stage.finish(format!("genes={}", gene_indices.len()));
    }
    result
}

fn compute_marker_genes(
    matrix: &ndarray::Array2<f32>,
    gene_names: &[String],
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    top_n: usize,
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing marker genes"));
    let mut out = Map::new();
    let all_gene_indices = (0..matrix.ncols()).collect::<Vec<_>>();
    for (column_idx, (column_name, values, categories)) in categorical_columns.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("columns", column_idx + 1, categorical_columns.len());
            stage.note(format!(
                "column {}/{}: {} ({} categories)",
                column_idx + 1,
                categorical_columns.len(),
                column_name,
                categories.len()
            ));
        }
        let summary = build_categorical_column_summary(matrix, values, categories);
        let mut column_out = Map::new();
        for (category_idx, category) in categories.iter().enumerate() {
            if let Some(stage) = stage.as_mut() {
                stage.progress("categories", category_idx + 1, categories.len());
            }
            let group_a = &summary.group_summaries[category_idx];
            let group_b = subtract_group_summary(&summary.total, group_a);
            if group_a.n_cells == 0 || group_b.n_cells == 0 {
                column_out.insert(category.clone(), Value::Array(Vec::new()));
                continue;
            }
            let de = differential_expression_from_summaries(
                gene_names,
                &all_gene_indices,
                group_a,
                &group_b,
                top_n,
            );
            column_out.insert(
                category.clone(),
                Value::Array(de.genes.into_iter().map(Value::String).collect()),
            );
        }
        out.insert(column_name.clone(), Value::Object(column_out));
    }
    if let Some(stage) = stage {
        stage.finish(format!("columns={}", categorical_columns.len()));
    }
    Value::Object(out)
}

fn compute_cluster_de(
    matrix: &ndarray::Array2<f32>,
    gene_names: &[String],
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    top_n: usize,
    min_cells: usize,
    method: DeMethod,
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing cluster DE"));
    let mut out = Map::new();
    let all_gene_indices = (0..matrix.ncols()).collect::<Vec<_>>();
    for (column_idx, (column_name, values, categories)) in categorical_columns.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("columns", column_idx + 1, categorical_columns.len());
            stage.note(format!(
                "column {}/{}: {} ({} categories)",
                column_idx + 1,
                categorical_columns.len(),
                column_name,
                categories.len()
            ));
        }
        let summary = if matches!(method, DeMethod::TTest) {
            Some(build_categorical_column_summary(matrix, values, categories))
        } else {
            None
        };
        let mut group_out = Map::new();
        for (source_idx, source) in categories.iter().enumerate() {
            if let Some(stage) = stage.as_mut() {
                stage.progress("source groups", source_idx + 1, categories.len());
            }
            let source_indices = values
                .iter()
                .enumerate()
                .filter_map(|(idx, value)| if value == source { Some(idx) } else { None })
                .collect::<Vec<_>>();
            let mut source_out = Map::new();
            for (reference_idx, reference) in categories.iter().enumerate() {
                if reference == source {
                    continue;
                }
                let reference_indices = values
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, value)| if value == reference { Some(idx) } else { None })
                    .collect::<Vec<_>>();
                if source_indices.len() < min_cells || reference_indices.len() < min_cells {
                    source_out.insert(
                        reference.clone(),
                        json!({
                            "available": false,
                            "reason": "insufficient_cells",
                            "genes": [],
                            "logfoldchanges": [],
                            "pvals_adj": [],
                            "scores": [],
                            "pct_source": [],
                            "pct_reference": [],
                            "n_source": source_indices.len(),
                            "n_reference": reference_indices.len(),
                            "min_cells_required": min_cells,
                        }),
                    );
                    continue;
                }
                let de = match method {
                    DeMethod::TTest => {
                        let summary = summary.as_ref().expect("summary available");
                        let group_a = &summary.group_summaries[source_idx];
                        let group_b = &summary.group_summaries[reference_idx];
                        differential_expression_from_summaries(
                            gene_names,
                            &all_gene_indices,
                            group_a,
                            group_b,
                            top_n,
                        )
                    }
                    DeMethod::Wilcoxon => differential_expression(
                        matrix,
                        gene_names,
                        &source_indices,
                        &reference_indices,
                        top_n,
                        method,
                    ),
                };
                source_out.insert(
                    reference.clone(),
                    json!({
                        "available": true,
                        "genes": de.genes,
                        "logfoldchanges": option_vec_to_json(&de.logfoldchanges),
                        "pvals_adj": option_vec_to_json(&de.pvals_adj),
                        "scores": option_vec_to_json(&de.scores),
                        "pct_source": option_vec_to_json(&de.pct_a),
                        "pct_reference": option_vec_to_json(&de.pct_b),
                        "n_source": source_indices.len(),
                        "n_reference": reference_indices.len(),
                    }),
                );
            }
            if !source_out.is_empty() {
                group_out.insert(source.clone(), Value::Object(source_out));
            }
        }
        if !group_out.is_empty() {
            out.insert(column_name.clone(), Value::Object(group_out));
        }
    }
    if let Some(stage) = stage {
        stage.finish(format!("columns={}, method={:?}", categorical_columns.len(), method));
    }
    Value::Object(out)
}

fn compute_neighbor_stats(
    graph: &CsrMatrix,
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    permutations: usize,
    seed: u64,
    reporter: Option<&ProgressReporter>,
) -> (Value, HashMap<String, NeighborStatsContext>) {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing neighbor stats"));
    let neighbors = graph_neighbors_from_csr(graph);
    let mut payload = Map::new();
    let mut contexts = HashMap::new();
    for (column_idx, (column_name, values, categories)) in categorical_columns.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("columns", column_idx + 1, categorical_columns.len());
            stage.note(format!(
                "column {}/{}: {} ({} categories)",
                column_idx + 1,
                categorical_columns.len(),
                column_name,
                categories.len()
            ));
        }
        let codes = categorical_codes(values, categories);
        let (counts, n_cells, mean_degree) = directed_neighbor_counts(&neighbors, &codes, categories.len());
        let zscores = if permutations > 0 {
            Some(permutation_neighbor_zscores(
                &neighbors,
                &codes,
                categories.len(),
                permutations,
                seed,
            ))
        } else {
            None
        };
        let mut entry = Map::new();
        entry.insert(
            "categories".to_string(),
            Value::Array(categories.iter().cloned().map(Value::String).collect()),
        );
        entry.insert("counts".to_string(), matrix_to_json(&counts));
        entry.insert(
            "n_cells".to_string(),
            Value::Array(n_cells.iter().map(|value| json!(*value)).collect()),
        );
        entry.insert(
            "mean_degree".to_string(),
            Value::Array(mean_degree.iter().map(|value| json!(*value)).collect()),
        );
        if let Some(zscore) = &zscores {
            entry.insert("perm_n".to_string(), json!(permutations));
            entry.insert("zscore".to_string(), matrix_to_json(zscore));
        }
        payload.insert(column_name.clone(), Value::Object(entry));
        contexts.insert(
            column_name.clone(),
            NeighborStatsContext {
                codes,
                counts,
                zscores,
                n_cells,
            },
        );
    }
    if let Some(stage) = stage {
        stage.finish(format!(
            "columns={}, permutations={}",
            categorical_columns.len(),
            permutations
        ));
    }
    (Value::Object(payload), contexts)
}

fn compute_interaction_markers(
    matrix: &ndarray::Array2<f32>,
    gene_names: &[String],
    categorical_columns: &[(String, Vec<String>, Vec<String>)],
    contexts: &HashMap<String, NeighborStatsContext>,
    neighbors: Vec<Vec<(usize, f32)>>,
    gene_indices: &[usize],
    top_targets: usize,
    top_genes: usize,
    min_cells: usize,
    min_neighbors: usize,
    method: DeMethod,
    reporter: Option<&ProgressReporter>,
) -> Value {
    let mut stage = reporter.map(|reporter| reporter.stage("Computing interaction markers"));
    let mut payload = Map::new();
    for (column_idx, (column_name, values, categories)) in categorical_columns.iter().enumerate() {
        if let Some(stage) = stage.as_mut() {
            stage.progress("columns", column_idx + 1, categorical_columns.len());
            stage.note(format!(
                "column {}/{}: {} ({} categories)",
                column_idx + 1,
                categorical_columns.len(),
                column_name,
                categories.len()
            ));
        }
        let Some(ctx) = contexts.get(column_name) else {
            continue;
        };
        let mut by_source = Map::new();
        for (source_idx, source) in categories.iter().enumerate() {
            if let Some(stage) = stage.as_mut() {
                stage.progress("source groups", source_idx + 1, categories.len());
            }
            if *ctx.n_cells.get(source_idx).unwrap_or(&0) == 0 {
                continue;
            }
            let source_cells = values
                .iter()
                .enumerate()
                .filter_map(|(idx, value)| if value == source { Some(idx) } else { None })
                .collect::<Vec<_>>();
            if source_cells.is_empty() {
                continue;
            }

            let mut candidate_targets = categories
                .iter()
                .enumerate()
                .filter_map(|(target_idx, _)| {
                    if target_idx == source_idx {
                        None
                    } else {
                        let count = ctx
                            .counts
                            .get(source_idx)
                            .and_then(|row| row.get(target_idx))
                            .copied()
                            .unwrap_or(0.0);
                        if count > 0.0 {
                            Some(target_idx)
                        } else {
                            None
                        }
                    }
                })
                .collect::<Vec<_>>();

            candidate_targets.sort_by(|lhs, rhs| {
                let left_z = ctx
                    .zscores
                    .as_ref()
                    .and_then(|matrix| matrix.get(source_idx).and_then(|row| row.get(*lhs)))
                    .copied()
                    .unwrap_or(f64::NEG_INFINITY);
                let right_z = ctx
                    .zscores
                    .as_ref()
                    .and_then(|matrix| matrix.get(source_idx).and_then(|row| row.get(*rhs)))
                    .copied()
                    .unwrap_or(f64::NEG_INFINITY);
                right_z
                    .partial_cmp(&left_z)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| {
                        let right_count = ctx.counts[source_idx][*rhs];
                        let left_count = ctx.counts[source_idx][*lhs];
                        right_count
                            .partial_cmp(&left_count)
                            .unwrap_or(Ordering::Equal)
                    })
            });
            candidate_targets.truncate(top_targets);

            let mut source_payload = Map::new();
            for target_idx in candidate_targets {
                let target_name = &categories[target_idx];
                let target_neighbor_counts = source_cells
                    .iter()
                    .map(|cell_idx| {
                        neighbors[*cell_idx]
                            .iter()
                            .filter(|(neighbor_idx, _)| {
                                ctx.codes.get(*neighbor_idx).and_then(|code| *code) == Some(target_idx)
                            })
                            .count()
                    })
                    .collect::<Vec<_>>();

                let mut pos = Vec::new();
                let mut neg = Vec::new();
                for (offset, cell_idx) in source_cells.iter().enumerate() {
                    let count = target_neighbor_counts[offset];
                    if count >= min_neighbors {
                        pos.push(*cell_idx);
                    } else if count == 0 {
                        neg.push(*cell_idx);
                    }
                }

                let pct_contact = if pos.len() + neg.len() == 0 {
                    0.0
                } else {
                    100.0 * pos.len() as f64 / (pos.len() + neg.len()) as f64
                };
                let mean_contact = mean_usize(
                    &source_cells
                        .iter()
                        .enumerate()
                        .filter_map(|(offset, _)| {
                            if target_neighbor_counts[offset] >= min_neighbors {
                                Some(target_neighbor_counts[offset])
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>(),
                );
                let mean_non_contact = mean_usize(
                    &source_cells
                        .iter()
                        .enumerate()
                        .filter_map(|(offset, _)| {
                            if target_neighbor_counts[offset] == 0 {
                                Some(target_neighbor_counts[offset])
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>(),
                );
                let target_edge_count = ctx.counts[source_idx][target_idx];
                let target_zscore = ctx
                    .zscores
                    .as_ref()
                    .and_then(|matrix| matrix.get(source_idx).and_then(|row| row.get(target_idx)))
                    .copied();

                if pos.len() < min_cells || neg.len() < min_cells {
                    source_payload.insert(
                        target_name.clone(),
                        json!({
                            "available": false,
                            "reason": "insufficient_cells",
                            "genes": [],
                            "logfoldchanges": [],
                            "pvals_adj": [],
                            "n_contact": pos.len(),
                            "n_non_contact": neg.len(),
                            "min_cells_required": min_cells,
                            "pct_contact": pct_contact,
                            "mean_target_neighbors_contact": mean_contact,
                            "mean_target_neighbors_non_contact": mean_non_contact,
                            "target_edge_count": target_edge_count,
                            "target_zscore": target_zscore,
                        }),
                    );
                    continue;
                }

                let de = differential_expression_subset(
                    matrix,
                    gene_names,
                    gene_indices,
                    &pos,
                    &neg,
                    top_genes,
                    method,
                );
                source_payload.insert(
                    target_name.clone(),
                    json!({
                        "available": true,
                        "genes": de.genes,
                        "logfoldchanges": option_vec_to_json(&de.logfoldchanges),
                        "pvals_adj": option_vec_to_json(&de.pvals_adj),
                        "n_contact": pos.len(),
                        "n_non_contact": neg.len(),
                        "pct_contact": pct_contact,
                        "mean_target_neighbors_contact": mean_contact,
                        "mean_target_neighbors_non_contact": mean_non_contact,
                        "target_edge_count": target_edge_count,
                        "target_zscore": target_zscore,
                    }),
                );
            }

            if !source_payload.is_empty() {
                by_source.insert(source.clone(), Value::Object(source_payload));
            }
        }
        if !by_source.is_empty() {
            payload.insert(column_name.clone(), Value::Object(by_source));
        }
    }
    if let Some(stage) = stage {
        stage.finish(format!(
            "columns={}, method={:?}, gene_limit={}",
            categorical_columns.len(),
            method,
            gene_indices.len()
        ));
    }
    Value::Object(payload)
}

fn build_categorical_column_summary(
    matrix: &ndarray::Array2<f32>,
    values: &[String],
    categories: &[String],
) -> CategoricalColumnSummary {
    let n_genes = matrix.ncols();
    let codes = categorical_codes(values, categories);
    let mut group_summaries = (0..categories.len())
        .map(|_| GroupSummary::zeros(n_genes))
        .collect::<Vec<_>>();
    let mut total = GroupSummary::zeros(n_genes);

    for (row_idx, code) in codes.into_iter().enumerate() {
        let Some(code) = code else {
            continue;
        };
        let row = matrix.row(row_idx);
        let row_iter: Box<dyn Iterator<Item = (usize, f32)>> = if let Some(slice) = row.as_slice() {
            Box::new(slice.iter().copied().enumerate())
        } else {
            Box::new(row.iter().copied().enumerate())
        };
        group_summaries[code].n_cells += 1;
        total.n_cells += 1;
        for (gene_idx, value) in row_iter {
            let value = f64::from(value);
            group_summaries[code].sums[gene_idx] += value;
            group_summaries[code].sums_sq[gene_idx] += value * value;
            total.sums[gene_idx] += value;
            total.sums_sq[gene_idx] += value * value;
            if value > 0.0 {
                group_summaries[code].nnz[gene_idx] += 1;
                total.nnz[gene_idx] += 1;
            }
        }
    }

    CategoricalColumnSummary {
        group_summaries,
        total,
    }
}

fn subtract_group_summary(total: &GroupSummary, subset: &GroupSummary) -> GroupSummary {
    GroupSummary {
        n_cells: total.n_cells.saturating_sub(subset.n_cells),
        sums: total
            .sums
            .iter()
            .zip(subset.sums.iter())
            .map(|(lhs, rhs)| lhs - rhs)
            .collect(),
        sums_sq: total
            .sums_sq
            .iter()
            .zip(subset.sums_sq.iter())
            .map(|(lhs, rhs)| lhs - rhs)
            .collect(),
        nnz: total
            .nnz
            .iter()
            .zip(subset.nnz.iter())
            .map(|(lhs, rhs)| lhs.saturating_sub(*rhs))
            .collect(),
    }
}

fn differential_expression_from_summaries(
    gene_names: &[String],
    gene_indices: &[usize],
    group_a: &GroupSummary,
    group_b: &GroupSummary,
    top_n: usize,
) -> DifferentialExpressionResult {
    let mut stats = Vec::with_capacity(gene_indices.len());
    for &gene_idx in gene_indices {
        let mean_a = group_a.sums[gene_idx] / group_a.n_cells.max(1) as f64;
        let mean_b = group_b.sums[gene_idx] / group_b.n_cells.max(1) as f64;
        let var_a = group_a.sums_sq[gene_idx] / group_a.n_cells.max(1) as f64 - mean_a * mean_a;
        let var_b = group_b.sums_sq[gene_idx] / group_b.n_cells.max(1) as f64 - mean_b * mean_b;
        let score = welch_t_score(
            mean_a,
            mean_b,
            var_a.max(0.0),
            var_b.max(0.0),
            group_a.n_cells,
            group_b.n_cells,
        );
        let p_value = two_sided_p_from_z(score);
        let logfc = ((mean_a + 1e-9) / (mean_b + 1e-9)).log2();
        let pct_a = group_a.nnz[gene_idx] as f64 / group_a.n_cells.max(1) as f64;
        let pct_b = group_b.nnz[gene_idx] as f64 / group_b.n_cells.max(1) as f64;
        stats.push((gene_idx, score, p_value, logfc, pct_a, pct_b));
    }

    finalize_de_stats(gene_names, stats, top_n)
}

fn differential_expression(
    matrix: &ndarray::Array2<f32>,
    gene_names: &[String],
    group_a: &[usize],
    group_b: &[usize],
    top_n: usize,
    method: DeMethod,
) -> DifferentialExpressionResult {
    let (mean_a, var_a, pct_a) = summarize_group(matrix, group_a);
    let (mean_b, var_b, pct_b) = summarize_group(matrix, group_b);

    let mut stats = Vec::with_capacity(matrix.ncols());
    for gene_idx in 0..matrix.ncols() {
        let score = match method {
            DeMethod::TTest => welch_t_score(
                mean_a[gene_idx],
                mean_b[gene_idx],
                var_a[gene_idx],
                var_b[gene_idx],
                group_a.len(),
                group_b.len(),
            ),
            DeMethod::Wilcoxon => wilcoxon_rank_sum_score(matrix, gene_idx, group_a, group_b),
        };
        let p_value = two_sided_p_from_z(score);
        let logfc = ((mean_a[gene_idx] + 1e-9) / (mean_b[gene_idx] + 1e-9)).log2();
        stats.push((gene_idx, score, p_value, logfc, pct_a[gene_idx], pct_b[gene_idx]));
    }

    finalize_de_stats(gene_names, stats, top_n)
}

fn differential_expression_subset(
    matrix: &ndarray::Array2<f32>,
    gene_names: &[String],
    gene_indices: &[usize],
    group_a: &[usize],
    group_b: &[usize],
    top_n: usize,
    method: DeMethod,
) -> DifferentialExpressionResult {
    match method {
        DeMethod::TTest => {
            let (mean_a, var_a, pct_a) = summarize_group_subset(matrix, group_a, gene_indices);
            let (mean_b, var_b, pct_b) = summarize_group_subset(matrix, group_b, gene_indices);
            let mut stats = Vec::with_capacity(gene_indices.len());
            for (offset, gene_idx) in gene_indices.iter().copied().enumerate() {
                let score = welch_t_score(
                    mean_a[offset],
                    mean_b[offset],
                    var_a[offset],
                    var_b[offset],
                    group_a.len(),
                    group_b.len(),
                );
                let p_value = two_sided_p_from_z(score);
                let logfc = ((mean_a[offset] + 1e-9) / (mean_b[offset] + 1e-9)).log2();
                stats.push((gene_idx, score, p_value, logfc, pct_a[offset], pct_b[offset]));
            }
            finalize_de_stats(gene_names, stats, top_n)
        }
        DeMethod::Wilcoxon => {
            let mut stats = Vec::with_capacity(gene_indices.len());
            for &gene_idx in gene_indices {
                let score = wilcoxon_rank_sum_score(matrix, gene_idx, group_a, group_b);
                let p_value = two_sided_p_from_z(score);
                let (mean_a, var_a, pct_a) = summarize_group_subset(matrix, group_a, &[gene_idx]);
                let (mean_b, var_b, pct_b) = summarize_group_subset(matrix, group_b, &[gene_idx]);
                let _ = (var_a, var_b);
                let logfc = ((mean_a[0] + 1e-9) / (mean_b[0] + 1e-9)).log2();
                stats.push((gene_idx, score, p_value, logfc, pct_a[0], pct_b[0]));
            }
            finalize_de_stats(gene_names, stats, top_n)
        }
    }
}

fn finalize_de_stats(
    gene_names: &[String],
    stats: Vec<(usize, f64, f64, f64, f64, f64)>,
    top_n: usize,
) -> DifferentialExpressionResult {
    let adjusted = benjamini_hochberg(
        &stats
            .iter()
            .map(|(_, _, p_value, _, _, _)| *p_value)
            .collect::<Vec<_>>(),
    );
    let mut scored = stats
        .into_iter()
        .enumerate()
        .map(|(idx, (gene_idx, score, _p, logfc, pct_a, pct_b))| {
            (
                gene_idx,
                score,
                adjusted[idx],
                logfc,
                pct_a,
                pct_b,
            )
        })
        .collect::<Vec<_>>();
    scored.sort_by(|lhs, rhs| rhs.1.partial_cmp(&lhs.1).unwrap_or(Ordering::Equal));

    let top = scored
        .into_iter()
        .filter(|(_, score, _, _, _, _)| score.is_finite())
        .take(top_n)
        .collect::<Vec<_>>();
    DifferentialExpressionResult {
        genes: top
            .iter()
            .map(|(gene_idx, _, _, _, _, _)| gene_names[*gene_idx].clone())
            .collect(),
        logfoldchanges: top
            .iter()
            .map(|(_, _, _, value, _, _)| finite_or_none(*value))
            .collect(),
        pvals_adj: top
            .iter()
            .map(|(_, _, value, _, _, _)| finite_or_none(*value))
            .collect(),
        scores: top
            .iter()
            .map(|(_, value, _, _, _, _)| finite_or_none(*value))
            .collect(),
        pct_a: top
            .iter()
            .map(|(_, _, _, _, value, _)| finite_or_none(*value))
            .collect(),
        pct_b: top
            .iter()
            .map(|(_, _, _, _, _, value)| finite_or_none(*value))
            .collect(),
    }
}

fn summarize_group_subset(
    matrix: &ndarray::Array2<f32>,
    indices: &[usize],
    gene_indices: &[usize],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut sums = vec![0.0f64; gene_indices.len()];
    let mut sums_sq = vec![0.0f64; gene_indices.len()];
    let mut nnz = vec![0usize; gene_indices.len()];
    for row_idx in indices {
        for (offset, gene_idx) in gene_indices.iter().copied().enumerate() {
            let value = f64::from(matrix[[*row_idx, gene_idx]]);
            sums[offset] += value;
            sums_sq[offset] += value * value;
            if value > 0.0 {
                nnz[offset] += 1;
            }
        }
    }
    let n = indices.len().max(1) as f64;
    let means = sums.iter().map(|sum| *sum / n).collect::<Vec<_>>();
    let vars = sums_sq
        .iter()
        .zip(means.iter())
        .map(|(sum_sq, mean)| (sum_sq / n) - (mean * mean))
        .collect::<Vec<_>>();
    let pct = nnz
        .iter()
        .map(|count| *count as f64 / indices.len().max(1) as f64)
        .collect::<Vec<_>>();
    (means, vars, pct)
}

fn summarize_group(
    matrix: &ndarray::Array2<f32>,
    indices: &[usize],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut sums = vec![0.0f64; matrix.ncols()];
    let mut sums_sq = vec![0.0f64; matrix.ncols()];
    let mut nnz = vec![0usize; matrix.ncols()];
    for row_idx in indices {
        for gene_idx in 0..matrix.ncols() {
            let value = f64::from(matrix[[*row_idx, gene_idx]]);
            sums[gene_idx] += value;
            sums_sq[gene_idx] += value * value;
            if value > 0.0 {
                nnz[gene_idx] += 1;
            }
        }
    }
    let n = indices.len().max(1) as f64;
    let means = sums.iter().map(|sum| *sum / n).collect::<Vec<_>>();
    let vars = sums_sq
        .iter()
        .zip(means.iter())
        .map(|(sum_sq, mean)| (sum_sq / n) - (mean * mean))
        .collect::<Vec<_>>();
    let pct = nnz
        .iter()
        .map(|count| *count as f64 / indices.len().max(1) as f64)
        .collect::<Vec<_>>();
    (means, vars, pct)
}

fn welch_t_score(
    mean_a: f64,
    mean_b: f64,
    var_a: f64,
    var_b: f64,
    n_a: usize,
    n_b: usize,
) -> f64 {
    let denom = (var_a / n_a.max(1) as f64 + var_b / n_b.max(1) as f64 + 1e-12).sqrt();
    (mean_a - mean_b) / denom
}

fn wilcoxon_rank_sum_score(
    matrix: &ndarray::Array2<f32>,
    gene_idx: usize,
    group_a: &[usize],
    group_b: &[usize],
) -> f64 {
    if group_a.is_empty() || group_b.is_empty() {
        return 0.0;
    }
    let mut values = Vec::with_capacity(group_a.len() + group_b.len());
    for idx in group_a {
        values.push((f64::from(matrix[[*idx, gene_idx]]), true));
    }
    for idx in group_b {
        values.push((f64::from(matrix[[*idx, gene_idx]]), false));
    }
    values.sort_by(|lhs, rhs| lhs.0.partial_cmp(&rhs.0).unwrap_or(Ordering::Equal));

    let mut rank_sum_a = 0.0;
    let mut start = 0usize;
    while start < values.len() {
        let value = values[start].0;
        let mut end = start + 1;
        while end < values.len() && values[end].0 == value {
            end += 1;
        }
        let avg_rank = (start + 1 + end) as f64 / 2.0;
        let a_count = values[start..end].iter().filter(|(_, is_a)| *is_a).count();
        rank_sum_a += avg_rank * a_count as f64;
        start = end;
    }

    let n1 = group_a.len() as f64;
    let n2 = group_b.len() as f64;
    let u1 = rank_sum_a - (n1 * (n1 + 1.0) / 2.0);
    let mean = n1 * n2 / 2.0;
    let variance = (n1 * n2 * (n1 + n2 + 1.0) / 12.0).max(1e-12);
    (u1 - mean) / variance.sqrt()
}

fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    if p_values.is_empty() {
        return Vec::new();
    }
    let mut indexed = p_values
        .iter()
        .copied()
        .enumerate()
        .collect::<Vec<_>>();
    indexed.sort_by(|lhs, rhs| lhs.1.partial_cmp(&rhs.1).unwrap_or(Ordering::Equal));
    let n = p_values.len() as f64;
    let mut adjusted = vec![1.0; p_values.len()];
    let mut prev = 1.0;
    for (rank, (idx, p_value)) in indexed.into_iter().enumerate().rev() {
        let raw = (p_value * n / (rank as f64 + 1.0)).min(1.0);
        let value = raw.min(prev);
        adjusted[idx] = value;
        prev = value;
    }
    adjusted
}

fn two_sided_p_from_z(z: f64) -> f64 {
    erfc(z.abs() / std::f64::consts::SQRT_2).clamp(0.0, 1.0)
}

fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0
        - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * (-x * x).exp());
    sign * y
}

fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

fn auto_neighbor_permutations(n_cells: usize) -> usize {
    if n_cells >= 200_000 {
        0
    } else {
        20
    }
}

fn categorical_codes(values: &[String], categories: &[String]) -> Vec<Option<usize>> {
    let lookup: HashMap<&str, usize> = categories
        .iter()
        .enumerate()
        .map(|(index, value)| (value.as_str(), index))
        .collect();
    values
        .iter()
        .map(|value| {
            if value.is_empty() {
                None
            } else {
                lookup.get(value.as_str()).copied()
            }
        })
        .collect()
}

fn directed_neighbor_counts(
    neighbors: &[Vec<(usize, f32)>],
    codes: &[Option<usize>],
    n_categories: usize,
) -> (Vec<Vec<f64>>, Vec<usize>, Vec<f64>) {
    let mut counts = vec![vec![0.0f64; n_categories]; n_categories];
    let mut n_cells = vec![0usize; n_categories];
    let mut degree_sums = vec![0.0f64; n_categories];
    for (cell_idx, source_code) in codes.iter().enumerate() {
        let Some(source_code) = source_code else {
            continue;
        };
        n_cells[*source_code] += 1;
        let mut degree = 0.0;
        for &(neighbor_idx, weight) in &neighbors[cell_idx] {
            let Some(target_code) = codes.get(neighbor_idx).and_then(|code| *code) else {
                continue;
            };
            counts[*source_code][target_code] += f64::from(weight);
            degree += f64::from(weight);
        }
        degree_sums[*source_code] += degree;
    }
    let mean_degree = degree_sums
        .iter()
        .enumerate()
        .map(|(idx, sum)| {
            if n_cells[idx] == 0 {
                0.0
            } else {
                *sum / n_cells[idx] as f64
            }
        })
        .collect::<Vec<_>>();
    (counts, n_cells, mean_degree)
}

fn permutation_neighbor_zscores(
    neighbors: &[Vec<(usize, f32)>],
    codes: &[Option<usize>],
    n_categories: usize,
    permutations: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let observed = directed_neighbor_counts(neighbors, codes, n_categories).0;
    let valid_positions = codes
        .iter()
        .enumerate()
        .filter_map(|(idx, code)| code.map(|_| idx))
        .collect::<Vec<_>>();
    let base_labels = valid_positions
        .iter()
        .map(|idx| codes[*idx].expect("valid positions"))
        .collect::<Vec<_>>();
    let mut mean = vec![vec![0.0f64; n_categories]; n_categories];
    let mut m2 = vec![vec![0.0f64; n_categories]; n_categories];
    let mut rng = StdRng::seed_from_u64(seed);

    for perm in 0..permutations {
        let mut shuffled = base_labels.clone();
        shuffled.shuffle(&mut rng);
        let mut perm_codes = vec![None; codes.len()];
        for (offset, global_idx) in valid_positions.iter().enumerate() {
            perm_codes[*global_idx] = Some(shuffled[offset]);
        }
        let counts = directed_neighbor_counts(neighbors, &perm_codes, n_categories).0;
        for i in 0..n_categories {
            for j in 0..n_categories {
                let delta = counts[i][j] - mean[i][j];
                mean[i][j] += delta / (perm + 1) as f64;
                m2[i][j] += delta * (counts[i][j] - mean[i][j]);
            }
        }
    }

    let mut z = vec![vec![0.0f64; n_categories]; n_categories];
    for i in 0..n_categories {
        for j in 0..n_categories {
            let variance = if permutations > 1 {
                m2[i][j] / (permutations - 1) as f64
            } else {
                0.0
            };
            let std = variance.sqrt();
            if std > 0.0 {
                z[i][j] = (observed[i][j] - mean[i][j]) / std;
            }
        }
    }
    z
}

fn graph_neighbors_from_csr(graph: &CsrMatrix) -> Vec<Vec<(usize, f32)>> {
    let mut neighbors = vec![Vec::new(); graph.nrows];
    for row in 0..graph.nrows {
        let start = graph.indptr[row].max(0) as usize;
        let end = graph.indptr[row + 1].max(0) as usize;
        for idx in start..end {
            let col = graph.indices[idx].max(0) as usize;
            neighbors[row].push((col, graph.data[idx]));
        }
    }
    neighbors
}

fn compute_umap_bounds(umap: &ndarray::Array2<f64>) -> Value {
    let xs = umap.column(0);
    let ys = umap.column(1);
    json!({
        "xmin": xs.iter().copied().reduce(f64::min).unwrap_or(0.0),
        "xmax": xs.iter().copied().reduce(f64::max).unwrap_or(0.0),
        "ymin": ys.iter().copied().reduce(f64::min).unwrap_or(0.0),
        "ymax": ys.iter().copied().reduce(f64::max).unwrap_or(0.0),
    })
}

fn mean_for_indices(matrix: &ndarray::Array2<f32>, gene_idx: usize, indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    indices
        .iter()
        .map(|row_idx| f64::from(matrix[[*row_idx, gene_idx]]))
        .sum::<f64>()
        / indices.len() as f64
}

fn mean_usize(values: &[usize]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<usize>() as f64 / values.len() as f64
    }
}

fn finite_or_none(value: f64) -> Option<f64> {
    if value.is_finite() {
        Some(value)
    } else {
        None
    }
}

fn encode_f32_b64(values: &[f32]) -> String {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    STANDARD.encode(bytes)
}

fn encode_u32_b64(values: &[u32]) -> String {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    STANDARD.encode(bytes)
}

fn f32_vec_to_json(values: &[f32]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|value| {
                if value.is_finite() {
                    json!(*value)
                } else {
                    Value::Null
                }
            })
            .collect(),
    )
}

fn numeric_vec_to_json(values: &[f64]) -> Value {
    Value::Array(values.iter().map(|value| json!(*value)).collect())
}

fn option_vec_to_json(values: &[Option<f64>]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|value| value.map_or(Value::Null, |v| json!(v)))
            .collect(),
    )
}

fn matrix_to_json(values: &[Vec<f64>]) -> Value {
    Value::Array(values.iter().map(|row| numeric_vec_to_json(row)).collect())
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process;
    use std::str::FromStr;
    use std::time::{SystemTime, UNIX_EPOCH};

    use hdf5::types::VarLenUnicode;
    use hdf5_metno as hdf5;
    use ndarray::arr2;

    use super::{export_viewer_precompute_json, DeMethod, GeneEncodingMode, ViewerPrecomputeConfig};
    use serde_json::Value;

    #[test]
    fn viewer_precompute_json_contains_core_karospace_keys() {
        with_temp_paths("viewer-json", |input, output, json_path| {
            create_test_h5ad(input);
            let config = ViewerPrecomputeConfig {
                output_json: Some(json_path.to_path_buf()),
                initial_color: Some("cell_type".to_string()),
                genes: None,
                analytics_columns: Some(vec!["cell_type".to_string()]),
                include_interaction_markers: true,
                gene_encoding: GeneEncodingMode::Auto,
                gene_sparse_zero_threshold: 0.8,
                pack_arrays: false,
                pack_arrays_min_len: 1024,
                gene_sparse_pack_min_nnz: 256,
                gene_correlation_top_n: 3,
                cluster_means_n_genes: 3,
                spatial_variable_genes_n: 3,
                marker_genes_top_n: 3,
                cluster_de_method: DeMethod::TTest,
                cluster_de_top_n: 3,
                cluster_de_min_cells: 1,
                neighbor_stats_seed: 0,
                interaction_markers_method: DeMethod::TTest,
                interaction_markers_gene_limit: 3,
                interaction_markers_top_targets: 2,
                interaction_markers_top_genes: 3,
                interaction_markers_min_cells: 1,
                interaction_markers_min_neighbors: 1,
            };

            fs::copy(input, output).unwrap();
            export_viewer_precompute_json(output, Some("sample"), &config).unwrap();

            let payload: Value = serde_json::from_slice(&fs::read(json_path).unwrap()).unwrap();
            assert!(payload.get("sections").unwrap().is_array());
            assert!(payload.get("available_colors").unwrap().is_array());
            assert!(payload.get("available_genes").unwrap().is_array());
            assert!(payload.get("colors_meta").unwrap().is_object());
            assert!(payload.get("marker_genes").unwrap().is_object());
            assert!(payload.get("cluster_de").unwrap().is_object());
            assert!(payload.get("neighbor_stats").unwrap().is_object());
            assert!(payload.get("cluster_gene_means").is_some());
            assert!(payload["available_colors"]
                .as_array()
                .unwrap()
                .iter()
                .any(|value| value.as_str() == Some("score")));
            assert!(payload["metadata_filters"].get("cell_type").is_some());
            assert!(payload["metadata_filters"].get("score").is_none());
        });
    }

    fn create_test_h5ad(path: &Path) {
        let file = hdf5::File::create(path).unwrap();

        let obs = file.create_group("obs").unwrap();
        write_string_dataset(&obs, "_index", &["cell1", "cell2", "cell3", "cell4"]);
        write_string_dataset(&obs, "sample", &["A", "A", "B", "B"]);
        write_string_dataset(&obs, "cell_type", &["t1", "t2", "t1", "t2"]);
        obs.new_dataset_builder()
            .with_data(&[1.0f32, 2.0, 3.0, 4.0])
            .create("score")
            .unwrap();

        let var = file.create_group("var").unwrap();
        write_string_dataset(&var, "_index", &["gene1", "gene2", "gene3"]);

        file.new_dataset_builder()
            .with_data(&arr2(&[
                [10.0f32, 0.0, 2.0],
                [0.0f32, 6.0, 1.0],
                [4.0f32, 0.0, 0.0],
                [0.0f32, 5.0, 3.0],
            ]))
            .create("X")
            .unwrap();

        let layers = file.create_group("layers").unwrap();
        layers
            .new_dataset_builder()
            .with_data(&arr2(&[
                [3.0f32, 0.0, 1.0],
                [0.0f32, 2.0, 1.0],
                [2.0f32, 0.0, 0.0],
                [0.0f32, 3.0, 1.0],
            ]))
            .create("normalized")
            .unwrap();

        let obsm = file.create_group("obsm").unwrap();
        obsm.new_dataset_builder()
            .with_data(&arr2(&[
                [0.0f32, 0.0],
                [1.0f32, 0.0],
                [10.0f32, 10.0],
                [11.0f32, 10.0],
            ]))
            .create("spatial")
            .unwrap();
        obsm.new_dataset_builder()
            .with_data(&arr2(&[
                [0.1f32, 0.2],
                [0.2f32, 0.3],
                [0.8f32, 0.7],
                [0.9f32, 0.8],
            ]))
            .create("X_umap")
            .unwrap();

        let obsp = file.create_group("obsp").unwrap();
        set_dict_attrs(&obsp);
        write_csr(
            &obsp,
            "spatial_connectivities",
            &[1.0, 1.0, 1.0, 1.0],
            &[1, 0, 3, 2],
            &[0, 1, 2, 3, 4],
            4,
            4,
        );
    }

    fn write_string_dataset(parent: &hdf5::Group, name: &str, values: &[&str]) {
        let encoded: Vec<VarLenUnicode> = values
            .iter()
            .map(|value| VarLenUnicode::from_str(value).unwrap())
            .collect();
        parent
            .new_dataset_builder()
            .with_data(&encoded)
            .create(name)
            .unwrap();
    }

    fn set_dict_attrs(group: &hdf5::Group) {
        let enc = VarLenUnicode::from_str("dict").unwrap();
        let ver = VarLenUnicode::from_str("0.1.0").unwrap();
        let attr = group.new_attr::<VarLenUnicode>().shape(()).create("encoding-type").unwrap();
        attr.write_scalar(&enc).unwrap();
        let attr = group.new_attr::<VarLenUnicode>().shape(()).create("encoding-version").unwrap();
        attr.write_scalar(&ver).unwrap();
    }

    fn write_csr(
        parent: &hdf5::Group,
        name: &str,
        data: &[f32],
        indices: &[i32],
        indptr: &[i32],
        nrows: usize,
        ncols: usize,
    ) {
        let group = parent.create_group(name).unwrap();
        let enc = VarLenUnicode::from_str("csr_matrix").unwrap();
        let ver = VarLenUnicode::from_str("0.1.0").unwrap();
        let attr = group.new_attr::<VarLenUnicode>().shape(()).create("encoding-type").unwrap();
        attr.write_scalar(&enc).unwrap();
        let attr = group.new_attr::<VarLenUnicode>().shape(()).create("encoding-version").unwrap();
        attr.write_scalar(&ver).unwrap();
        group
            .new_attr_builder()
            .with_data(&[nrows as i64, ncols as i64])
            .create("shape")
            .unwrap();
        group.new_dataset_builder().with_data(data).create("data").unwrap();
        group.new_dataset_builder().with_data(indices).create("indices").unwrap();
        group.new_dataset_builder().with_data(indptr).create("indptr").unwrap();
    }

    fn with_temp_paths(name: &str, f: impl FnOnce(&Path, &Path, &Path)) {
        let input = temp_path(&format!("{name}-input"), "h5ad");
        let output = temp_path(&format!("{name}-output"), "h5ad");
        let json = temp_path(&format!("{name}-viewer"), "json");
        f(&input, &output, &json);
        let _ = fs::remove_file(&input);
        let _ = fs::remove_file(&output);
        let _ = fs::remove_file(&json);
    }

    fn temp_path(name: &str, ext: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "karospace-companion-{name}-{}-{nanos}.{ext}",
            process::id()
        ))
    }
}
