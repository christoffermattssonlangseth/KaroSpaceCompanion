use std::path::PathBuf;

use anyhow::Result;
use clap::{Args, Parser, Subcommand};

use karospace_companion::anndata::MatrixSource;
use karospace_companion::graph::GraphMode;
use karospace_companion::viewer::{DeMethod, GeneEncodingMode, ViewerPrecomputeConfig};
use karospace_companion::{prepare, PrepareConfig};

#[derive(Parser)]
#[command(name = "karospace-companion")]
#[command(about = "Prepare KaroSpace-ready h5ad files in Rust")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    Prepare(PrepareArgs),
}

#[derive(Args)]
struct PrepareArgs {
    input: PathBuf,

    #[arg(short, long)]
    output: PathBuf,

    #[arg(long, conflicts_with = "k")]
    radius: Option<f64>,

    #[arg(long, conflicts_with = "radius")]
    k: Option<usize>,

    #[arg(long)]
    groupby: Option<String>,

    #[arg(long = "composition-cell-type")]
    composition_cell_type: Option<String>,

    #[arg(long = "normalize-from", default_value = "X")]
    normalize_from: String,

    #[arg(long = "skip-normalized-layer", default_value_t = false)]
    skip_normalized_layer: bool,

    #[arg(long = "target-sum", default_value_t = 10_000.0)]
    target_sum: f32,

    #[arg(long = "no-log1p", default_value_t = false)]
    no_log1p: bool,

    #[arg(long = "skip-aggregation", default_value_t = false)]
    skip_aggregation: bool,

    #[arg(long = "aggregate-from", default_value = "normalized")]
    aggregate_from: String,

    #[arg(long = "nmf-components")]
    nmf_components: Option<usize>,

    #[arg(long = "nmf-max-iter", default_value_t = 100)]
    nmf_max_iter: usize,

    #[arg(long = "nmf-seed", default_value_t = 42)]
    nmf_seed: u64,

    #[arg(long = "overwrite-derived", default_value_t = false)]
    overwrite_derived: bool,

    #[arg(long = "viewer-json")]
    viewer_json: Option<PathBuf>,

    #[arg(long = "persist-analytics-in-h5ad", default_value_t = false)]
    persist_analytics_in_h5ad: bool,

    #[arg(long = "initial-color")]
    initial_color: Option<String>,

    #[arg(long = "viewer-genes")]
    viewer_genes: Option<String>,

    #[arg(long = "viewer-analytics-columns")]
    viewer_analytics_columns: Option<String>,

    #[arg(long = "skip-viewer-interaction-markers", default_value_t = false)]
    skip_viewer_interaction_markers: bool,

    #[arg(long = "viewer-gene-encoding", default_value = "auto")]
    viewer_gene_encoding: String,

    #[arg(long = "viewer-gene-zero-threshold", default_value_t = 0.8)]
    viewer_gene_zero_threshold: f32,

    #[arg(long = "viewer-pack-min-len", default_value_t = 1024)]
    viewer_pack_min_len: usize,

    #[arg(long = "no-viewer-pack-arrays", default_value_t = false)]
    no_viewer_pack_arrays: bool,

    #[arg(long = "viewer-gene-sparse-pack-min-nnz", default_value_t = 256)]
    viewer_gene_sparse_pack_min_nnz: usize,

    #[arg(long = "viewer-cluster-de-method", default_value = "t-test")]
    viewer_cluster_de_method: String,

    #[arg(long = "viewer-interaction-method", default_value = "t-test")]
    viewer_interaction_method: String,

    #[arg(long = "viewer-interaction-gene-limit", default_value_t = 512)]
    viewer_interaction_gene_limit: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Prepare(args) => {
            let graph_mode = match (args.radius, args.k) {
                (Some(radius), None) => GraphMode::Radius(radius),
                (None, Some(k)) => GraphMode::Knn(k),
                _ => anyhow::bail!("exactly one of --radius or --k must be provided"),
            };
            let viewer_precompute = if args.viewer_json.is_some() || args.persist_analytics_in_h5ad {
                Some(ViewerPrecomputeConfig {
                    output_json: args.viewer_json,
                    initial_color: args.initial_color,
                    genes: parse_csv_list(args.viewer_genes.as_deref()),
                    analytics_columns: parse_csv_list(args.viewer_analytics_columns.as_deref()),
                    include_interaction_markers: !args.skip_viewer_interaction_markers,
                    gene_encoding: GeneEncodingMode::parse(&args.viewer_gene_encoding)?,
                    gene_sparse_zero_threshold: args.viewer_gene_zero_threshold,
                    pack_arrays: !args.no_viewer_pack_arrays,
                    pack_arrays_min_len: args.viewer_pack_min_len,
                    gene_sparse_pack_min_nnz: args.viewer_gene_sparse_pack_min_nnz,
                    gene_correlation_top_n: 10,
                    cluster_means_n_genes: 500,
                    spatial_variable_genes_n: 200,
                    marker_genes_top_n: 30,
                    cluster_de_method: DeMethod::parse(&args.viewer_cluster_de_method)?,
                    cluster_de_top_n: 20,
                    cluster_de_min_cells: 20,
                    neighbor_stats_seed: 0,
                    interaction_markers_method: DeMethod::parse(&args.viewer_interaction_method)?,
                    interaction_markers_gene_limit: args.viewer_interaction_gene_limit,
                    interaction_markers_top_targets: 8,
                    interaction_markers_top_genes: 20,
                    interaction_markers_min_cells: 30,
                    interaction_markers_min_neighbors: 1,
                })
            } else {
                None
            };

            let config = PrepareConfig {
                input: args.input,
                output: args.output,
                graph_mode,
                groupby: args.groupby,
                composition_cell_type: args.composition_cell_type,
                normalize_from: MatrixSource::parse(&args.normalize_from)?,
                skip_normalized_layer: args.skip_normalized_layer,
                target_sum: args.target_sum,
                log1p: !args.no_log1p,
                skip_aggregation: args.skip_aggregation,
                aggregate_from: MatrixSource::parse(&args.aggregate_from)?,
                nmf_components: args.nmf_components,
                nmf_max_iter: args.nmf_max_iter,
                nmf_seed: args.nmf_seed,
                overwrite_derived: args.overwrite_derived,
                viewer_precompute,
            };

            let summary = prepare(&config)?;
            println!(
                "Wrote {} (cells={}, genes={}, undirected_edges={})",
                summary.output.display(),
                summary.n_cells,
                summary.n_genes,
                summary.n_edges
            );
            if let Some(viewer_json) = summary.viewer_json {
                println!("Viewer precompute JSON: {}", viewer_json.display());
            }
        }
    }

    Ok(())
}

fn parse_csv_list(raw: Option<&str>) -> Option<Vec<String>> {
    raw.map(|value| {
        value
            .split(',')
            .map(str::trim)
            .filter(|entry| !entry.is_empty())
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>()
    })
    .filter(|values| !values.is_empty())
}
