pub mod anndata;
pub mod features;
pub mod graph;
pub mod nmf;
pub mod progress;
pub mod viewer;

use std::path::PathBuf;

use anndata::{
    ensure_copy_of_input, read_h5ad, write_augmented_outputs, AnnDataRead, AugmentedOutputs,
    DenseMatrixUpdate, GraphMatrixUpdate, MatrixSource, MetadataScalars, StringArrayUpdate,
};
use anyhow::{bail, Context, Result};
use features::{aggregate_neighbors, compute_composition_matrix, normalize_counts};
use graph::{build_spatial_graph, GraphMode};
use nmf::{run_nmf, NmfConfig};
use progress::ProgressReporter;
use viewer::{
    compute_companion_analytics, export_viewer_precompute_json_with_analytics,
    ViewerPrecomputeConfig,
};

pub const SCHEMA_VERSION: &str = "0.1.0";

#[derive(Debug, Clone)]
pub struct PrepareConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub graph_mode: GraphMode,
    pub groupby: Option<String>,
    pub composition_cell_type: Option<String>,
    pub normalize_from: MatrixSource,
    pub skip_normalized_layer: bool,
    pub target_sum: f32,
    pub log1p: bool,
    pub skip_aggregation: bool,
    pub aggregate_from: MatrixSource,
    pub nmf_components: Option<usize>,
    pub nmf_max_iter: usize,
    pub nmf_seed: u64,
    pub overwrite_derived: bool,
    pub viewer_precompute: Option<ViewerPrecomputeConfig>,
}

#[derive(Debug, Clone)]
pub struct PrepareSummary {
    pub n_cells: usize,
    pub n_genes: usize,
    pub n_edges: usize,
    pub output: PathBuf,
    pub viewer_json: Option<PathBuf>,
}

pub fn prepare(config: &PrepareConfig) -> Result<PrepareSummary> {
    let reporter = ProgressReporter::new("prepare");
    if config.input == config.output {
        bail!("input and output paths must be different");
    }

    if !config.skip_aggregation {
        match config.aggregate_from {
            MatrixSource::Normalized => {
                if config.skip_normalized_layer {
                    bail!("aggregate source 'normalized' requires the normalized layer to be computed");
                }
            }
            MatrixSource::Obsm(_) => {}
            MatrixSource::X | MatrixSource::Layer(_) => {
                bail!("aggregate source must be 'normalized' or 'obsm:<key>'");
            }
        }
    }

    let mut obs_cols: Vec<String> = Vec::new();
    if let Some(groupby) = &config.groupby {
        obs_cols.push(groupby.clone());
    }
    if let Some(cell_type) = &config.composition_cell_type {
        if !obs_cols.iter().any(|col| col == cell_type) {
            obs_cols.push(cell_type.clone());
        }
    }

    let mut obsm_keys: Vec<String> = Vec::new();
    if let MatrixSource::Obsm(key) = &config.aggregate_from {
        obsm_keys.push(key.clone());
    }

    let needs_expression = !config.skip_normalized_layer || config.nmf_components.is_some();
    let stage = reporter.stage("Loading input h5ad");
    let adata = read_h5ad(
        &config.input,
        &obs_cols.iter().map(String::as_str).collect::<Vec<_>>(),
        &obsm_keys.iter().map(String::as_str).collect::<Vec<_>>(),
        needs_expression,
        match &config.normalize_from {
            MatrixSource::X => None,
            MatrixSource::Layer(name) => Some(name.as_str()),
            MatrixSource::Normalized | MatrixSource::Obsm(_) => {
                bail!("normalize-from must be 'X' or 'layer:<name>'");
            }
        },
    )?;
    stage.finish(format!(
        "cells={}, genes={}, obs_cols={}, embeddings={}",
        adata.obs_names.len(),
        adata.var_names.len(),
        obs_cols.len(),
        obsm_keys.len()
    ));

    let group_values = resolve_group_values(&adata, config.groupby.as_deref())?;
    let stage = reporter.stage("Building spatial graph");
    let graph = build_spatial_graph(&adata.coordinates, &group_values, config.graph_mode.clone())?;
    stage.finish(format!(
        "groups={}, undirected_edges={}",
        unique_count(&group_values),
        graph.n_undirected_edges
    ));

    let normalized = if config.skip_normalized_layer {
        None
    } else {
        let stage = reporter.stage("Normalizing expression");
        let expr = adata
            .expression
            .as_ref()
            .context("expression matrix is required to compute a normalized layer")?;
        let normalized = normalize_counts(expr, config.target_sum, config.log1p);
        stage.finish(format!(
            "shape={}x{}, source={}",
            normalized.nrows(),
            normalized.ncols(),
            config.normalize_from.as_cli_value()
        ));
        Some(normalized)
    };

    let aggregate_matrix = if config.skip_aggregation {
        None
    } else {
        Some(match &config.aggregate_from {
            MatrixSource::Normalized => normalized
                .as_ref()
                .context("normalized matrix missing for aggregation")?
                .clone(),
            MatrixSource::Obsm(key) => adata
                .embeddings
                .get(key)
                .with_context(|| format!("requested aggregation embedding '{key}' was not loaded"))?
                .mapv(|value| value as f32),
            MatrixSource::X | MatrixSource::Layer(_) => unreachable!("validated above"),
        })
    };

    let aggregated = aggregate_matrix.as_ref().map(|matrix| {
        let stage = reporter.stage("Aggregating neighbor features");
        let aggregated = aggregate_neighbors(&graph.neighbors, matrix);
        stage.finish(format!(
            "shape={}x{}, source={}",
            aggregated.nrows(),
            aggregated.ncols(),
            config.aggregate_from.as_cli_value()
        ));
        aggregated
    });

    let composition = if let Some(cell_type_col) = &config.composition_cell_type {
        let stage = reporter.stage(format!(
            "Computing neighborhood composition ({})",
            cell_type_col
        ));
        let values = adata
            .obs
            .get(cell_type_col)
            .with_context(|| format!("obs column '{cell_type_col}' not found"))?;
        let composition = compute_composition_matrix(&graph.neighbors, values);
        stage.finish(format!(
            "shape={}x{}, categories={}",
            composition.0.nrows(),
            composition.0.ncols(),
            composition.1.len()
        ));
        Some(composition)
    } else {
        None
    };

    let nmf_result = if let Some(n_components) = config.nmf_components {
        let stage = reporter.stage("Running NMF");
        let matrix = normalized
            .as_ref()
            .context("normalized layer is required before running NMF")?;
        let nmf_config = NmfConfig {
            n_components,
            max_iter: config.nmf_max_iter,
            seed: config.nmf_seed,
            tol: 1e-4,
            epsilon: 1e-12,
        };
        let result = run_nmf(matrix, &nmf_config)?;
        stage.finish(format!(
            "components={}, iterations={}, final_error={:.4}",
            result.w.ncols(),
            result.n_iter,
            result.final_error
        ));
        Some(result)
    } else {
        None
    };

    let stage = reporter.stage("Copying input file");
    ensure_copy_of_input(&config.input, &config.output)?;
    stage.finish(config.output.display().to_string());

    let mut dense_updates = Vec::new();
    dense_updates.push(DenseMatrixUpdate {
        group_path: "obsm".to_string(),
        name: "spatial".to_string(),
        data: adata.coordinates.mapv(|value| value as f32),
        overwrite: true,
        require_flag: false,
        encoding_type: "array".to_string(),
        encoding_version: "0.2.0".to_string(),
    });

    if let Some(matrix) = normalized {
        dense_updates.push(DenseMatrixUpdate {
            group_path: "layers".to_string(),
            name: "normalized".to_string(),
            data: matrix,
            overwrite: config.overwrite_derived,
            require_flag: true,
            encoding_type: "array".to_string(),
            encoding_version: "0.2.0".to_string(),
        });
    }

    if let Some(matrix) = aggregated {
        dense_updates.push(DenseMatrixUpdate {
            group_path: "obsm".to_string(),
            name: "X_karo_agg".to_string(),
            data: matrix,
            overwrite: config.overwrite_derived,
            require_flag: true,
            encoding_type: "array".to_string(),
            encoding_version: "0.2.0".to_string(),
        });
    }

    let mut string_arrays = Vec::new();
    if let Some((matrix, categories)) = composition {
        dense_updates.push(DenseMatrixUpdate {
            group_path: "obsm".to_string(),
            name: "X_karo_comp".to_string(),
            data: matrix,
            overwrite: config.overwrite_derived,
            require_flag: true,
            encoding_type: "array".to_string(),
            encoding_version: "0.2.0".to_string(),
        });
        string_arrays.push(StringArrayUpdate {
            path: "uns/karospace_companion/X_karo_comp_categories".to_string(),
            values: categories,
            overwrite: config.overwrite_derived,
            require_flag: true,
        });
    }

    if let Some(result) = &nmf_result {
        dense_updates.push(DenseMatrixUpdate {
            group_path: "obsm".to_string(),
            name: "X_karo_nmf".to_string(),
            data: result.w.clone(),
            overwrite: config.overwrite_derived,
            require_flag: true,
            encoding_type: "array".to_string(),
            encoding_version: "0.2.0".to_string(),
        });
        dense_updates.push(DenseMatrixUpdate {
            group_path: "varm".to_string(),
            name: "X_karo_nmf_loadings".to_string(),
            data: result.h.t().to_owned(),
            overwrite: config.overwrite_derived,
            require_flag: true,
            encoding_type: "array".to_string(),
            encoding_version: "0.2.0".to_string(),
        });
    }

    let graph_updates = vec![
        GraphMatrixUpdate {
            name: "spatial_connectivities".to_string(),
            matrix: graph.connectivities,
            overwrite: config.overwrite_derived,
            require_flag: true,
        },
        GraphMatrixUpdate {
            name: "spatial_distances".to_string(),
            matrix: graph.distances,
            overwrite: config.overwrite_derived,
            require_flag: true,
        },
    ];

    let mut metadata = MetadataScalars::default();
    metadata
        .strings
        .insert("schema_version".to_string(), SCHEMA_VERSION.to_string());
    metadata
        .strings
        .insert("command".to_string(), "prepare".to_string());
    metadata.strings.insert(
        "graph_mode".to_string(),
        config.graph_mode.label().to_string(),
    );
    metadata.strings.insert(
        "groupby".to_string(),
        config.groupby.clone().unwrap_or_default(),
    );
    metadata.strings.insert(
        "composition_cell_type".to_string(),
        config.composition_cell_type.clone().unwrap_or_default(),
    );
    metadata.strings.insert(
        "normalize_from".to_string(),
        config.normalize_from.as_cli_value(),
    );
    metadata.strings.insert(
        "aggregate_from".to_string(),
        config.aggregate_from.as_cli_value(),
    );
    metadata
        .numbers
        .insert("target_sum".to_string(), config.target_sum as f64);
    metadata
        .numbers
        .insert("log1p".to_string(), if config.log1p { 1.0 } else { 0.0 });
    metadata
        .numbers
        .insert("n_cells".to_string(), adata.obs_names.len() as f64);
    metadata
        .numbers
        .insert("n_genes".to_string(), adata.var_names.len() as f64);
    metadata
        .numbers
        .insert("n_edges".to_string(), graph.n_undirected_edges as f64);
    match &config.graph_mode {
        GraphMode::Radius(radius) => {
            metadata.numbers.insert("radius".to_string(), *radius);
        }
        GraphMode::Knn(k) => {
            metadata.numbers.insert("k".to_string(), *k as f64);
        }
        GraphMode::Delaunay {
            remove_long_links_percentile,
        } => {
            metadata.numbers.insert(
                "remove_long_links_percentile".to_string(),
                *remove_long_links_percentile,
            );
        }
    }
    if let Some(result) = &nmf_result {
        metadata
            .numbers
            .insert("nmf_components".to_string(), result.w.ncols() as f64);
        metadata
            .numbers
            .insert("nmf_iterations".to_string(), result.n_iter as f64);
        metadata
            .numbers
            .insert("nmf_final_error".to_string(), result.final_error as f64);
    }

    let outputs = AugmentedOutputs {
        dense_updates,
        graph_updates,
        metadata,
        string_arrays,
    };
    let stage = reporter.stage("Writing derived outputs");
    write_augmented_outputs(&config.output, outputs)?;
    stage.finish(config.output.display().to_string());

    let companion_analytics = if let Some(viewer_config) = &config.viewer_precompute {
        let stage = reporter.stage("Computing companion analytics");
        let analytics =
            compute_companion_analytics(&config.output, config.groupby.as_deref(), viewer_config)?;
        stage.finish(format!(
            "columns={}, interaction_markers={}",
            analytics.analytics_columns.len(),
            if viewer_config.include_interaction_markers {
                "enabled"
            } else {
                "disabled"
            }
        ));

        let mut metadata = MetadataScalars::default();
        metadata.strings.insert(
            "analytics_storage".to_string(),
            "json-string-v1".to_string(),
        );
        metadata.strings.insert(
            "analytics_cluster_de_method".to_string(),
            format!("{:?}", viewer_config.cluster_de_method).to_ascii_lowercase(),
        );
        metadata.strings.insert(
            "analytics_interaction_method".to_string(),
            format!("{:?}", viewer_config.interaction_markers_method).to_ascii_lowercase(),
        );
        metadata.strings.insert(
            "analytics_interaction_markers_enabled".to_string(),
            if viewer_config.include_interaction_markers {
                "true".to_string()
            } else {
                "false".to_string()
            },
        );
        for (key, value) in analytics.to_json_scalars()? {
            metadata.strings.insert(key, value);
        }
        let string_arrays = if analytics.analytics_columns.is_empty() {
            Vec::new()
        } else {
            vec![StringArrayUpdate {
                path: "uns/karospace_companion/analytics_columns".to_string(),
                values: analytics.analytics_columns.clone(),
                overwrite: true,
                require_flag: true,
            }]
        };
        let stage = reporter.stage("Persisting analytics into h5ad");
        write_augmented_outputs(
            &config.output,
            AugmentedOutputs {
                dense_updates: Vec::new(),
                graph_updates: Vec::new(),
                metadata,
                string_arrays,
            },
        )?;
        stage.finish("uns/karospace_companion".to_string());
        Some(analytics)
    } else {
        None
    };

    let viewer_json = if let (Some(viewer_config), Some(companion_analytics)) =
        (&config.viewer_precompute, companion_analytics.as_ref())
    {
        if viewer_config.output_json.is_some() {
            let stage = reporter.stage("Exporting viewer precompute");
            let viewer_json = export_viewer_precompute_json_with_analytics(
                &config.output,
                config.groupby.as_deref(),
                viewer_config,
                companion_analytics,
            )?;
            stage.finish(viewer_json.display().to_string());
            Some(viewer_json)
        } else {
            None
        }
    } else {
        None
    };

    Ok(PrepareSummary {
        n_cells: adata.obs_names.len(),
        n_genes: adata.var_names.len(),
        n_edges: graph.n_undirected_edges,
        output: config.output.clone(),
        viewer_json,
    })
}

fn unique_count(values: &[String]) -> usize {
    let mut seen = std::collections::HashSet::new();
    for value in values {
        seen.insert(value);
    }
    seen.len()
}

fn resolve_group_values(adata: &AnnDataRead, groupby: Option<&str>) -> Result<Vec<String>> {
    if let Some(column) = groupby {
        let values = adata
            .obs
            .get(column)
            .with_context(|| format!("obs column '{column}' not found"))?;
        Ok(values.clone())
    } else {
        Ok(vec!["all".to_string(); adata.obs_names.len()])
    }
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

    use super::{prepare, PrepareConfig};
    use crate::anndata::MatrixSource;
    use crate::graph::GraphMode;
    use crate::viewer::{DeMethod, GeneEncodingMode, ViewerPrecomputeConfig};

    #[test]
    fn prepare_writes_expected_karospace_outputs() {
        with_temp_paths("prepare-basic", |input, output| {
            create_test_h5ad(input);
            let config = PrepareConfig {
                input: input.to_path_buf(),
                output: output.to_path_buf(),
                graph_mode: GraphMode::Radius(1.5),
                groupby: Some("sample".to_string()),
                composition_cell_type: Some("cell_type".to_string()),
                normalize_from: MatrixSource::X,
                skip_normalized_layer: false,
                target_sum: 10_000.0,
                log1p: true,
                skip_aggregation: false,
                aggregate_from: MatrixSource::Normalized,
                nmf_components: Some(2),
                nmf_max_iter: 25,
                nmf_seed: 7,
                overwrite_derived: true,
                viewer_precompute: None,
            };

            let summary = prepare(&config).unwrap();
            assert_eq!(summary.n_cells, 4);
            assert_eq!(summary.n_genes, 3);
            assert_eq!(summary.n_edges, 2);

            let file = hdf5::File::open(output).unwrap();
            assert!(file.group("layers").unwrap().link_exists("normalized"));
            assert!(file.group("obsm").unwrap().link_exists("spatial"));
            assert!(file.group("obsm").unwrap().link_exists("X_karo_agg"));
            assert!(file.group("obsm").unwrap().link_exists("X_karo_comp"));
            assert!(file.group("obsm").unwrap().link_exists("X_karo_nmf"));
            assert!(file
                .group("varm")
                .unwrap()
                .link_exists("X_karo_nmf_loadings"));
            assert!(file
                .group("obsp")
                .unwrap()
                .link_exists("spatial_connectivities"));
            assert!(file.group("obsp").unwrap().link_exists("spatial_distances"));
            assert!(file
                .group("uns")
                .unwrap()
                .link_exists("karospace_companion"));
        });
    }

    #[test]
    fn prepare_persists_companion_analytics_into_uns() {
        with_temp_paths("prepare-analytics", |input, output| {
            create_test_h5ad(input);
            let config = PrepareConfig {
                input: input.to_path_buf(),
                output: output.to_path_buf(),
                graph_mode: GraphMode::Radius(1.5),
                groupby: Some("sample".to_string()),
                composition_cell_type: Some("cell_type".to_string()),
                normalize_from: MatrixSource::X,
                skip_normalized_layer: true,
                target_sum: 10_000.0,
                log1p: true,
                skip_aggregation: true,
                aggregate_from: MatrixSource::Normalized,
                nmf_components: None,
                nmf_max_iter: 25,
                nmf_seed: 7,
                overwrite_derived: true,
                viewer_precompute: Some(ViewerPrecomputeConfig {
                    output_json: None,
                    initial_color: Some("cell_type".to_string()),
                    genes: Some(vec!["gene1".to_string(), "gene2".to_string()]),
                    analytics_columns: Some(vec!["cell_type".to_string()]),
                    include_interaction_markers: false,
                    gene_encoding: GeneEncodingMode::Auto,
                    gene_sparse_zero_threshold: 0.8,
                    pack_arrays: false,
                    pack_arrays_min_len: 1024,
                    gene_sparse_pack_min_nnz: 256,
                    gene_correlation_top_n: 2,
                    gene_correlation_n_genes: 2,
                    cluster_means_n_genes: 2,
                    spatial_variable_genes_n: 2,
                    marker_genes_top_n: 2,
                    cluster_de_method: DeMethod::TTest,
                    cluster_de_top_n: 2,
                    cluster_de_min_cells: 1,
                    neighbor_stats_permutations: Some(4),
                    neighbor_stats_seed: 0,
                    interaction_markers_method: DeMethod::TTest,
                    interaction_markers_gene_limit: 2,
                    interaction_markers_top_targets: 2,
                    interaction_markers_top_genes: 2,
                    interaction_markers_min_cells: 1,
                    interaction_markers_min_neighbors: 1,
                }),
            };

            let summary = prepare(&config).unwrap();
            assert_eq!(summary.viewer_json, None);

            let file = hdf5::File::open(output).unwrap();
            let group = file.group("uns/karospace_companion").unwrap();
            assert!(group.link_exists("analytics_storage"));
            assert!(group.link_exists("analytics_columns"));
            assert!(group.link_exists("marker_genes_json"));
            assert!(group.link_exists("cluster_de_json"));
            assert!(group.link_exists("neighbor_stats_json"));
            assert!(group.link_exists("interaction_markers_json"));
            assert!(group.link_exists("gene_correlations_json"));
            assert!(group.link_exists("spatial_variable_genes_json"));
            assert!(group.link_exists("cluster_gene_means_json"));
        });
    }

    #[test]
    fn prepare_reads_csc_integer_counts_layer() {
        with_temp_paths("prepare-csc-counts", |input, output| {
            create_test_h5ad(input);
            let config = PrepareConfig {
                input: input.to_path_buf(),
                output: output.to_path_buf(),
                graph_mode: GraphMode::Radius(1.5),
                groupby: Some("sample".to_string()),
                composition_cell_type: Some("cell_type".to_string()),
                normalize_from: MatrixSource::Layer("counts".to_string()),
                skip_normalized_layer: false,
                target_sum: 10_000.0,
                log1p: true,
                skip_aggregation: true,
                aggregate_from: MatrixSource::Normalized,
                nmf_components: None,
                nmf_max_iter: 25,
                nmf_seed: 7,
                overwrite_derived: true,
                viewer_precompute: None,
            };

            let summary = prepare(&config).unwrap();
            assert_eq!(summary.n_cells, 4);

            let file = hdf5::File::open(output).unwrap();
            assert!(file.group("layers").unwrap().link_exists("normalized"));
        });
    }

    fn create_test_h5ad(path: &Path) {
        let file = hdf5::File::create(path).unwrap();

        let obs = file.create_group("obs").unwrap();
        write_string_dataset(&obs, "_index", &["cell1", "cell2", "cell3", "cell4"]);
        write_string_dataset(&obs, "sample", &["A", "A", "B", "B"]);
        write_string_dataset(&obs, "cell_type", &["t1", "t2", "t1", "t2"]);
        write_categorical_i64_dataset(&obs, "cluster_id", &[0, 1, 0, 1], &[10, 20]);

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
        write_csc_u32(
            &layers,
            "counts",
            &[3, 2, 2, 3, 1, 1, 1],
            &[0, 2, 1, 3, 0, 1, 3],
            &[0, 2, 4, 7],
            4,
            3,
        );

        let obsm = file.create_group("obsm").unwrap();
        obsm.new_dataset_builder()
            .with_data(&arr2(&[
                [0.0f32, 0.0],
                [1.0f32, 0.0],
                [10.0f32, 10.0],
                [11.0f32, 10.0],
            ]))
            .create("X_spatial")
            .unwrap();
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

    fn write_csc_u32(
        parent: &hdf5::Group,
        name: &str,
        data: &[u32],
        indices: &[i32],
        indptr: &[i32],
        nrows: usize,
        ncols: usize,
    ) {
        let group = parent.create_group(name).unwrap();
        let enc = VarLenUnicode::from_str("csc_matrix").unwrap();
        let ver = VarLenUnicode::from_str("0.1.0").unwrap();
        let attr = group
            .new_attr::<VarLenUnicode>()
            .shape(())
            .create("encoding-type")
            .unwrap();
        attr.write_scalar(&enc).unwrap();
        let attr = group
            .new_attr::<VarLenUnicode>()
            .shape(())
            .create("encoding-version")
            .unwrap();
        attr.write_scalar(&ver).unwrap();
        group
            .new_attr_builder()
            .with_data(&[nrows as i64, ncols as i64])
            .create("shape")
            .unwrap();
        group
            .new_dataset_builder()
            .with_data(data)
            .create("data")
            .unwrap();
        group
            .new_dataset_builder()
            .with_data(indices)
            .create("indices")
            .unwrap();
        group
            .new_dataset_builder()
            .with_data(indptr)
            .create("indptr")
            .unwrap();
    }

    fn write_categorical_i64_dataset(
        parent: &hdf5::Group,
        name: &str,
        codes: &[i64],
        categories: &[i64],
    ) {
        let group = parent.create_group(name).unwrap();
        group
            .new_dataset_builder()
            .with_data(codes)
            .create("codes")
            .unwrap();
        group
            .new_dataset_builder()
            .with_data(categories)
            .create("categories")
            .unwrap();
    }

    fn with_temp_paths(name: &str, f: impl FnOnce(&Path, &Path)) {
        let input = temp_path(&format!("{name}-input"));
        let output = temp_path(&format!("{name}-output"));
        f(&input, &output);
        let _ = fs::remove_file(&input);
        let _ = fs::remove_file(&output);
    }

    fn temp_path(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "karospace-companion-{name}-{}-{nanos}.h5ad",
            process::id()
        ))
    }
}
