use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use anyhow::{bail, Context, Result};
use hdf5::types::VarLenUnicode;
use hdf5_metno as hdf5;
use ndarray::Array2;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixSource {
    X,
    Layer(String),
    Obsm(String),
    Normalized,
}

impl MatrixSource {
    pub fn parse(value: &str) -> Result<Self> {
        let raw = value.trim();
        if raw.eq_ignore_ascii_case("X") {
            return Ok(Self::X);
        }
        if raw.eq_ignore_ascii_case("normalized") {
            return Ok(Self::Normalized);
        }
        if let Some(name) = raw.strip_prefix("layer:") {
            if name.is_empty() {
                bail!("layer source cannot be empty");
            }
            return Ok(Self::Layer(name.to_string()));
        }
        if let Some(name) = raw.strip_prefix("obsm:") {
            if name.is_empty() {
                bail!("obsm source cannot be empty");
            }
            return Ok(Self::Obsm(name.to_string()));
        }
        bail!("unsupported matrix source '{raw}'");
    }

    pub fn as_cli_value(&self) -> String {
        match self {
            Self::X => "X".to_string(),
            Self::Layer(name) => format!("layer:{name}"),
            Self::Obsm(name) => format!("obsm:{name}"),
            Self::Normalized => "normalized".to_string(),
        }
    }
}

pub struct AnnDataRead {
    pub obs_names: Vec<String>,
    pub var_names: Vec<String>,
    pub coordinates: Array2<f64>,
    pub obs: HashMap<String, Vec<String>>,
    pub expression: Option<Array2<f32>>,
    pub embeddings: HashMap<String, Array2<f64>>,
}

#[derive(Debug, Clone)]
pub enum ObsColumnData {
    Categorical {
        values: Vec<String>,
        categories: Vec<String>,
    },
    Numeric {
        values: Vec<Option<f64>>,
    },
}

impl ObsColumnData {
    pub fn len(&self) -> usize {
        match self {
            Self::Categorical { values, .. } => values.len(),
            Self::Numeric { values } => values.len(),
        }
    }

    pub fn is_continuous(&self) -> bool {
        matches!(self, Self::Numeric { .. })
    }

    pub fn as_strings(&self) -> Vec<String> {
        match self {
            Self::Categorical { values, .. } => values.clone(),
            Self::Numeric { values } => values
                .iter()
                .map(|value| value.map_or_else(String::new, format_numeric_for_string))
                .collect(),
        }
    }

    pub fn categories(&self) -> Option<&[String]> {
        match self {
            Self::Categorical { categories, .. } => Some(categories.as_slice()),
            Self::Numeric { .. } => None,
        }
    }

    pub fn categorical_codes_f32(&self) -> Option<Vec<f32>> {
        match self {
            Self::Categorical { values, categories } => {
                let lookup: HashMap<&str, usize> = categories
                    .iter()
                    .enumerate()
                    .map(|(index, value)| (value.as_str(), index))
                    .collect();
                Some(
                    values
                        .iter()
                        .map(|value| {
                            if value.is_empty() {
                                f32::NAN
                            } else {
                                lookup
                                    .get(value.as_str())
                                    .map(|index| *index as f32)
                                    .unwrap_or(f32::NAN)
                            }
                        })
                        .collect(),
                )
            }
            Self::Numeric { .. } => None,
        }
    }

    pub fn numeric_values_f32(&self) -> Option<Vec<f32>> {
        match self {
            Self::Numeric { values } => Some(
                values
                    .iter()
                    .map(|value| value.map(|v| v as f32).unwrap_or(f32::NAN))
                    .collect(),
            ),
            Self::Categorical { .. } => None,
        }
    }
}

pub struct ViewerAnnDataRead {
    pub obs_names: Vec<String>,
    pub var_names: Vec<String>,
    pub coordinates: Array2<f64>,
    pub obs_order: Vec<String>,
    pub obs_columns: HashMap<String, ObsColumnData>,
    pub expression: Array2<f32>,
    pub expression_source: String,
    pub umap: Option<Array2<f64>>,
    pub graph: Option<CsrMatrix>,
    pub graph_key: Option<String>,
}

#[derive(Clone)]
pub struct DenseMatrixUpdate {
    pub group_path: String,
    pub name: String,
    pub data: Array2<f32>,
    pub overwrite: bool,
    pub require_flag: bool,
    pub encoding_type: String,
    pub encoding_version: String,
}

#[derive(Clone)]
pub struct CsrMatrix {
    pub data: Vec<f32>,
    pub indices: Vec<i32>,
    pub indptr: Vec<i32>,
    pub nrows: usize,
    pub ncols: usize,
}

#[derive(Clone)]
pub struct GraphMatrixUpdate {
    pub name: String,
    pub matrix: CsrMatrix,
    pub overwrite: bool,
    pub require_flag: bool,
}

#[derive(Default, Clone)]
pub struct MetadataScalars {
    pub strings: HashMap<String, String>,
    pub numbers: HashMap<String, f64>,
}

#[derive(Clone)]
pub struct StringArrayUpdate {
    pub path: String,
    pub values: Vec<String>,
    pub overwrite: bool,
    pub require_flag: bool,
}

pub struct AugmentedOutputs {
    pub dense_updates: Vec<DenseMatrixUpdate>,
    pub graph_updates: Vec<GraphMatrixUpdate>,
    pub metadata: MetadataScalars,
    pub string_arrays: Vec<StringArrayUpdate>,
}

pub fn read_h5ad(
    path: &Path,
    obs_cols: &[&str],
    obsm_keys: &[&str],
    load_expression: bool,
    expression_layer: Option<&str>,
) -> Result<AnnDataRead> {
    let file = hdf5::File::open(path).with_context(|| format!("cannot open {:?}", path))?;

    let obs_names = read_obs_names(&file)?;
    let var_names = read_var_names(&file)?;
    let coordinates = read_spatial(&file)?;
    ensure_row_count("spatial coordinates", coordinates.nrows(), obs_names.len())?;

    let mut obs = HashMap::new();
    for col in obs_cols {
        let values = read_obs_column(&file, col).with_context(|| format!("reading obs/{col}"))?;
        ensure_row_count(&format!("obs/{col}"), values.len(), obs_names.len())?;
        obs.insert((*col).to_string(), values);
    }

    let expression = if load_expression {
        let matrix = read_expression_matrix(&file, expression_layer)?;
        ensure_row_count("expression matrix", matrix.nrows(), obs_names.len())?;
        if matrix.ncols() != var_names.len() {
            bail!(
                "expression column count ({}) does not match var count ({})",
                matrix.ncols(),
                var_names.len()
            );
        }
        Some(matrix)
    } else {
        None
    };

    let mut embeddings = HashMap::new();
    for key in obsm_keys {
        let array =
            read_obsm_embedding(&file, key).with_context(|| format!("reading obsm/{key}"))?;
        ensure_row_count(&format!("obsm/{key}"), array.nrows(), obs_names.len())?;
        embeddings.insert((*key).to_string(), array);
    }

    Ok(AnnDataRead {
        obs_names,
        var_names,
        coordinates,
        obs,
        expression,
        embeddings,
    })
}

pub fn read_viewer_h5ad(path: &Path) -> Result<ViewerAnnDataRead> {
    let file = hdf5::File::open(path).with_context(|| format!("cannot open {:?}", path))?;

    let obs_names = read_obs_names(&file)?;
    let var_names = read_var_names(&file)?;
    let coordinates = read_spatial(&file)?;
    ensure_row_count("spatial coordinates", coordinates.nrows(), obs_names.len())?;

    let (obs_order, obs_columns) = read_all_obs_columns(&file)?;
    for (name, column) in &obs_columns {
        ensure_row_count(&format!("obs/{name}"), column.len(), obs_names.len())?;
    }

    let expression_layer = if layer_exists(&file, "normalized") {
        Some("normalized")
    } else {
        None
    };
    let expression = read_expression_matrix(&file, expression_layer)?;
    ensure_row_count("expression matrix", expression.nrows(), obs_names.len())?;
    if expression.ncols() != var_names.len() {
        bail!(
            "expression column count ({}) does not match var count ({})",
            expression.ncols(),
            var_names.len()
        );
    }

    let umap = match read_optional_obsm_embedding(&file, "X_umap")? {
        Some(array) => {
            ensure_row_count("obsm/X_umap", array.nrows(), obs_names.len())?;
            Some(slice_embedding_two_cols(array)?)
        }
        None => None,
    };

    let graph_candidates = [
        "obsp/spatial_connectivities",
        "obsp/connectivities",
        "obsp/neighbors",
        "obsp/neighbor_graph",
    ];
    let mut graph = None;
    let mut graph_key = None;
    for path in graph_candidates {
        if file.link_exists(path) {
            let matrix = read_csr_matrix_path(&file, path)?;
            if matrix.nrows == obs_names.len() && matrix.ncols == obs_names.len() {
                graph = Some(matrix);
                graph_key = path.rsplit('/').next().map(ToOwned::to_owned);
                break;
            }
        }
    }

    Ok(ViewerAnnDataRead {
        obs_names,
        var_names,
        coordinates,
        obs_order,
        obs_columns,
        expression,
        expression_source: expression_layer.unwrap_or("X").to_string(),
        umap,
        graph,
        graph_key,
    })
}

pub fn ensure_copy_of_input(input: &Path, output: &Path) -> Result<()> {
    fs::copy(input, output)
        .with_context(|| format!("copying input {:?} to {:?}", input, output))?;
    Ok(())
}

pub fn write_augmented_outputs(path: &Path, outputs: AugmentedOutputs) -> Result<()> {
    let file =
        hdf5::File::open_rw(path).with_context(|| format!("opening {:?} for update", path))?;

    ensure_dict_group(&file, "obsm")?;
    ensure_dict_group(&file, "obsp")?;
    ensure_dict_group(&file, "layers")?;
    ensure_dict_group(&file, "uns")?;
    ensure_dict_group(&file, "varm")?;
    ensure_dict_group(&file, "uns/karospace_companion")?;

    for update in outputs.dense_updates {
        let parent = ensure_dict_group(&file, &update.group_path)?;
        write_dense_matrix(
            &parent,
            &update.name,
            &update.data,
            update.overwrite,
            update.require_flag,
            &update.encoding_type,
            &update.encoding_version,
        )?;
    }

    let obsp = file.group("obsp").context("opening obsp group")?;
    for update in outputs.graph_updates {
        write_csr_matrix(
            &obsp,
            &update.name,
            &update.matrix,
            update.overwrite,
            update.require_flag,
        )?;
    }

    let metadata_group = ensure_dict_group(&file, "uns/karospace_companion")?;
    for (key, value) in outputs.metadata.strings {
        write_string_scalar(&metadata_group, &key, &value, true)?;
    }
    for (key, value) in outputs.metadata.numbers {
        write_numeric_scalar(&metadata_group, &key, value, true)?;
    }

    for update in outputs.string_arrays {
        write_string_array_path(
            &file,
            &update.path,
            &update.values,
            update.overwrite,
            update.require_flag,
        )?;
    }

    Ok(())
}

fn read_obs_names(file: &hdf5::File) -> Result<Vec<String>> {
    let obs = file.group("obs").context("no 'obs' group")?;
    read_index_column(&obs, "obs")
}

fn read_var_names(file: &hdf5::File) -> Result<Vec<String>> {
    let var = file.group("var").context("no 'var' group")?;
    read_index_column(&var, "var")
}

fn read_index_column(group: &hdf5::Group, label: &str) -> Result<Vec<String>> {
    if let Ok(ds) = group.dataset("_index") {
        return read_string_dataset(&ds).with_context(|| format!("reading {label}/_index"));
    }
    if let Ok(ds) = group.dataset("index") {
        return read_string_dataset(&ds).with_context(|| format!("reading {label}/index"));
    }
    if let Ok(attr) = group.attr("_index") {
        let name: VarLenUnicode = attr.read_scalar().context("reading _index attribute")?;
        let ds = group
            .dataset(name.as_str())
            .with_context(|| format!("opening {label}/{}", name.as_str()))?;
        return read_string_dataset(&ds);
    }
    bail!("could not locate {label} index")
}

fn read_spatial(file: &hdf5::File) -> Result<Array2<f64>> {
    let obsm = file.group("obsm").context("no 'obsm' group")?;
    for key in ["spatial", "X_spatial"] {
        if let Ok(ds) = obsm.dataset(key) {
            if let Ok(array) = ds.read_2d::<f64>() {
                return slice_two_cols(array);
            }
            if let Ok(array) = ds.read_2d::<f32>() {
                return slice_two_cols(array.mapv(|value| value as f64));
            }
            bail!("obsm/{key} is not a numeric matrix");
        }
    }
    bail!("no spatial coordinates found in obsm/spatial or obsm/X_spatial")
}

fn slice_two_cols(array: Array2<f64>) -> Result<Array2<f64>> {
    if array.ncols() < 2 {
        bail!("spatial coordinate matrix must have at least 2 columns");
    }
    Ok(array.slice(ndarray::s![.., 0..2]).to_owned())
}

fn read_expression_matrix(file: &hdf5::File, layer: Option<&str>) -> Result<Array2<f32>> {
    let group_key = match layer {
        Some(name) => format!("layers/{name}"),
        None => "X".to_string(),
    };

    if let Ok(group) = file.group(&group_key) {
        if let (Ok(data_ds), Ok(indices_ds), Ok(indptr_ds)) = (
            group.dataset("data"),
            group.dataset("indices"),
            group.dataset("indptr"),
        ) {
            let data = read_f32_vec(&data_ds)?;
            let indices = read_usize_vec(&indices_ds)?;
            let indptr = read_usize_vec(&indptr_ds)?;
            let (nrows, ncols) = read_csr_shape(&group)?;
            validate_csr_layout(nrows, ncols, data.len(), &indices, &indptr)?;
            let mut dense = Array2::<f32>::zeros((nrows, ncols));
            for row in 0..nrows {
                let start = indptr[row];
                let end = indptr[row + 1];
                for idx in start..end {
                    dense[[row, indices[idx]]] = data[idx];
                }
            }
            return Ok(dense);
        }
    }

    let ds = file
        .dataset(&group_key)
        .with_context(|| format!("opening dataset '{group_key}'"))?;
    if let Ok(array) = ds.read_2d::<f32>() {
        return Ok(array);
    }
    if let Ok(array) = ds.read_2d::<f64>() {
        return Ok(array.mapv(|value| value as f32));
    }
    bail!("expression matrix '{group_key}' is not readable as f32 or f64")
}

fn read_obsm_embedding(file: &hdf5::File, key: &str) -> Result<Array2<f64>> {
    let obsm = file.group("obsm").context("no 'obsm' group")?;
    let ds = obsm
        .dataset(key)
        .with_context(|| format!("obsm/{key} not found"))?;
    if let Ok(array) = ds.read_2d::<f64>() {
        return Ok(array);
    }
    if let Ok(array) = ds.read_2d::<f32>() {
        return Ok(array.mapv(|value| value as f64));
    }
    bail!("obsm/{key} is not a numeric matrix")
}

fn read_obs_column(file: &hdf5::File, col: &str) -> Result<Vec<String>> {
    let obs = file.group("obs").context("no 'obs' group")?;
    if let Ok(group) = obs.group(col) {
        if let (Ok(codes_ds), Ok(categories_ds)) =
            (group.dataset("codes"), group.dataset("categories"))
        {
            let categories = read_string_dataset(&categories_ds)?;
            let codes = read_category_codes(&codes_ds)?;
            return decode_categorical_codes(&categories, &codes, col);
        }
    }
    let ds = obs
        .dataset(col)
        .with_context(|| format!("obs/{col} not found"))?;
    if let Ok(values) = read_string_dataset(&ds) {
        return Ok(values);
    }
    if let Ok(values) = read_numeric_obs_values(&ds) {
        return Ok(values
            .into_iter()
            .map(|value| value.map_or_else(String::new, format_numeric_for_string))
            .collect());
    }
    if let Ok(values) = ds.read_1d::<bool>() {
        return Ok(values
            .iter()
            .map(|value| if *value { "true".to_string() } else { "false".to_string() })
            .collect());
    }
    bail!("obs/{col} is not readable as string, numeric, or bool")
}

fn read_all_obs_columns(file: &hdf5::File) -> Result<(Vec<String>, HashMap<String, ObsColumnData>)> {
    let obs = file.group("obs").context("no 'obs' group")?;
    let index_name = resolve_index_member_name(&obs)?;
    let mut order = Vec::new();
    let mut columns = HashMap::new();

    for name in obs.member_names().context("listing obs members")? {
        if name == "_index" || name == "index" || Some(name.as_str()) == index_name.as_deref() {
            continue;
        }
        let column = read_obs_column_data(&obs, &name)?;
        order.push(name.clone());
        columns.insert(name, column);
    }

    Ok((order, columns))
}

fn resolve_index_member_name(group: &hdf5::Group) -> Result<Option<String>> {
    if group.dataset("_index").is_ok() {
        return Ok(Some("_index".to_string()));
    }
    if group.dataset("index").is_ok() {
        return Ok(Some("index".to_string()));
    }
    if let Ok(attr) = group.attr("_index") {
        let name: VarLenUnicode = attr.read_scalar().context("reading _index attribute")?;
        return Ok(Some(name.to_string()));
    }
    Ok(None)
}

fn read_obs_column_data(obs: &hdf5::Group, col: &str) -> Result<ObsColumnData> {
    if let Ok(group) = obs.group(col) {
        if let (Ok(codes_ds), Ok(categories_ds)) =
            (group.dataset("codes"), group.dataset("categories"))
        {
            let categories = read_string_dataset(&categories_ds)?;
            let codes = read_category_codes(&codes_ds)?;
            let values = decode_categorical_codes(&categories, &codes, col)?;
            return Ok(ObsColumnData::Categorical { values, categories });
        }
    }

    let ds = obs
        .dataset(col)
        .with_context(|| format!("obs/{col} not found"))?;
    if let Ok(values) = read_string_dataset(&ds) {
        let categories = collect_categories_in_order(&values);
        return Ok(ObsColumnData::Categorical { values, categories });
    }
    if let Ok(values) = ds.read_1d::<bool>() {
        let strings: Vec<String> = values
            .iter()
            .map(|value| if *value { "true".to_string() } else { "false".to_string() })
            .collect();
        let categories = collect_categories_in_order(&strings);
        return Ok(ObsColumnData::Categorical {
            values: strings,
            categories,
        });
    }
    if let Ok(values) = read_numeric_obs_values(&ds) {
        return Ok(ObsColumnData::Numeric { values });
    }
    bail!("obs/{col} is not a supported categorical or numeric column")
}

fn collect_categories_in_order(values: &[String]) -> Vec<String> {
    let mut seen: HashMap<&str, ()> = HashMap::new();
    let mut categories = Vec::new();
    for value in values {
        if value.is_empty() || seen.contains_key(value.as_str()) {
            continue;
        }
        seen.insert(value.as_str(), ());
        categories.push(value.clone());
    }
    categories
}

fn read_string_dataset(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    ds.read_1d::<VarLenUnicode>()
        .map(|values| values.iter().map(|value| value.to_string()).collect())
        .context("reading string dataset")
}

fn read_category_codes(ds: &hdf5::Dataset) -> Result<Vec<i64>> {
    if let Ok(values) = ds.read_1d::<i8>() {
        return Ok(values.iter().map(|&value| value as i64).collect());
    }
    if let Ok(values) = ds.read_1d::<i16>() {
        return Ok(values.iter().map(|&value| value as i64).collect());
    }
    if let Ok(values) = ds.read_1d::<i32>() {
        return Ok(values.iter().map(|&value| value as i64).collect());
    }
    if let Ok(values) = ds.read_1d::<i64>() {
        return Ok(values.to_vec());
    }
    bail!("categorical codes are not stored as signed integers")
}

fn decode_categorical_codes(
    categories: &[String],
    codes: &[i64],
    col: &str,
) -> Result<Vec<String>> {
    codes
        .iter()
        .map(|&code| match code {
            -1 => Ok(String::new()),
            c if c < -1 => bail!("negative code {c} is invalid for obs/{col}"),
            c => categories
                .get(c as usize)
                .cloned()
                .with_context(|| format!("code {c} is out of range for obs/{col}")),
        })
        .collect()
}

fn read_f32_vec(ds: &hdf5::Dataset) -> Result<Vec<f32>> {
    if let Ok(values) = ds.read_1d::<f32>() {
        return Ok(values.to_vec());
    }
    if let Ok(values) = ds.read_1d::<f64>() {
        return Ok(values.iter().map(|&value| value as f32).collect());
    }
    bail!("floating-point dataset is not f32 or f64")
}

fn read_numeric_obs_values(ds: &hdf5::Dataset) -> Result<Vec<Option<f64>>> {
    if let Ok(values) = ds.read_1d::<f64>() {
        return Ok(values
            .iter()
            .map(|value| if value.is_finite() { Some(*value) } else { None })
            .collect());
    }
    if let Ok(values) = ds.read_1d::<f32>() {
        return Ok(values
            .iter()
            .map(|value| if value.is_finite() { Some(*value as f64) } else { None })
            .collect());
    }
    if let Ok(values) = ds.read_1d::<i64>() {
        return Ok(values.iter().map(|value| Some(*value as f64)).collect());
    }
    if let Ok(values) = ds.read_1d::<i32>() {
        return Ok(values.iter().map(|value| Some(*value as f64)).collect());
    }
    if let Ok(values) = ds.read_1d::<u64>() {
        return Ok(values.iter().map(|value| Some(*value as f64)).collect());
    }
    if let Ok(values) = ds.read_1d::<u32>() {
        return Ok(values.iter().map(|value| Some(*value as f64)).collect());
    }
    bail!("dataset is not a supported numeric vector")
}

fn format_numeric_for_string(value: f64) -> String {
    if (value - value.round()).abs() < 1e-9 {
        format!("{}", value.round() as i64)
    } else {
        value.to_string()
    }
}

fn read_usize_vec(ds: &hdf5::Dataset) -> Result<Vec<usize>> {
    if let Ok(values) = ds.read_1d::<i32>() {
        return values
            .iter()
            .map(|&value| checked_i64_to_usize(value as i64))
            .collect();
    }
    if let Ok(values) = ds.read_1d::<i64>() {
        return values
            .iter()
            .map(|&value| checked_i64_to_usize(value))
            .collect();
    }
    if let Ok(values) = ds.read_1d::<u32>() {
        return Ok(values.iter().map(|&value| value as usize).collect());
    }
    if let Ok(values) = ds.read_1d::<u64>() {
        return Ok(values.iter().map(|&value| value as usize).collect());
    }
    bail!("integer dataset is not a supported index type")
}

fn checked_i64_to_usize(value: i64) -> Result<usize> {
    if value < 0 {
        bail!("negative index value {value} is invalid");
    }
    Ok(value as usize)
}

fn read_csr_matrix_path(file: &hdf5::File, path: &str) -> Result<CsrMatrix> {
    let group = file
        .group(path)
        .with_context(|| format!("opening csr group '{path}'"))?;
    let data = read_f32_vec(&group.dataset("data").context("reading csr data dataset")?)?;
    let indices_u = read_usize_vec(&group.dataset("indices").context("reading csr indices dataset")?)?;
    let indptr_u = read_usize_vec(&group.dataset("indptr").context("reading csr indptr dataset")?)?;
    let (nrows, ncols) = read_csr_shape(&group)?;
    validate_csr_layout(nrows, ncols, data.len(), &indices_u, &indptr_u)?;
    Ok(CsrMatrix {
        data,
        indices: indices_u.into_iter().map(|value| value as i32).collect(),
        indptr: indptr_u.into_iter().map(|value| value as i32).collect(),
        nrows,
        ncols,
    })
}

fn layer_exists(file: &hdf5::File, name: &str) -> bool {
    match file.group("layers") {
        Ok(group) => group.link_exists(name),
        Err(_) => false,
    }
}

fn read_optional_obsm_embedding(file: &hdf5::File, key: &str) -> Result<Option<Array2<f64>>> {
    let obsm = match file.group("obsm") {
        Ok(group) => group,
        Err(_) => return Ok(None),
    };
    if !obsm.link_exists(key) {
        return Ok(None);
    }
    let ds = obsm
        .dataset(key)
        .with_context(|| format!("opening obsm/{key}"))?;
    if let Ok(array) = ds.read_2d::<f64>() {
        return Ok(Some(array));
    }
    if let Ok(array) = ds.read_2d::<f32>() {
        return Ok(Some(array.mapv(|value| value as f64)));
    }
    bail!("obsm/{key} is not a numeric matrix")
}

fn slice_embedding_two_cols(array: Array2<f64>) -> Result<Array2<f64>> {
    if array.ncols() < 2 {
        bail!("embedding matrix must have at least 2 columns");
    }
    Ok(array.slice(ndarray::s![.., 0..2]).to_owned())
}

fn read_csr_shape(group: &hdf5::Group) -> Result<(usize, usize)> {
    let attr = group
        .attr("shape")
        .context("csr group missing shape attribute")?;
    if let Ok(values) = attr.read_1d::<i64>() {
        let vec = values.to_vec();
        return parse_shape_pair(
            vec.iter()
                .map(|value| *value as i128)
                .collect::<Vec<_>>()
                .as_slice(),
        );
    }
    if let Ok(values) = attr.read_1d::<i32>() {
        let vec = values.to_vec();
        return parse_shape_pair(
            vec.iter()
                .map(|value| *value as i128)
                .collect::<Vec<_>>()
                .as_slice(),
        );
    }
    if let Ok(values) = attr.read_1d::<u64>() {
        let vec = values.to_vec();
        return parse_shape_pair(
            vec.iter()
                .map(|value| *value as i128)
                .collect::<Vec<_>>()
                .as_slice(),
        );
    }
    bail!("unsupported csr shape attribute type")
}

fn parse_shape_pair(values: &[i128]) -> Result<(usize, usize)> {
    if values.len() != 2 {
        bail!("csr shape attribute must have exactly 2 elements");
    }
    if values[0] < 0 || values[1] < 0 {
        bail!("csr shape contains negative values");
    }
    Ok((values[0] as usize, values[1] as usize))
}

fn validate_csr_layout(
    nrows: usize,
    ncols: usize,
    data_len: usize,
    indices: &[usize],
    indptr: &[usize],
) -> Result<()> {
    if indptr.len() != nrows + 1 {
        bail!(
            "csr indptr length ({}) must equal nrows + 1 ({})",
            indptr.len(),
            nrows + 1
        );
    }
    if indptr.first().copied() != Some(0) {
        bail!("csr indptr must start at 0");
    }
    if indptr[indptr.len() - 1] != data_len {
        bail!("csr indptr terminal value does not equal data length");
    }
    if indices.len() != data_len {
        bail!(
            "csr indices length ({}) does not match data length ({data_len})",
            indices.len()
        );
    }
    for window in indptr.windows(2) {
        if window[0] > window[1] {
            bail!("csr indptr must be monotonically non-decreasing");
        }
    }
    for &index in indices {
        if index >= ncols {
            bail!("csr column index {index} exceeds width {ncols}");
        }
    }
    Ok(())
}

fn ensure_row_count(label: &str, actual: usize, expected: usize) -> Result<()> {
    if actual != expected {
        bail!("{label} row count ({actual}) does not match expected count ({expected})");
    }
    Ok(())
}

fn ensure_dict_group(file: &hdf5::File, path: &str) -> Result<hdf5::Group> {
    if file.link_exists(path) {
        let group = file
            .group(path)
            .with_context(|| format!("opening group '{path}'"))?;
        set_string_attr(&*group, "encoding-type", "dict")?;
        set_string_attr(&*group, "encoding-version", "0.1.0")?;
        return Ok(group);
    }
    let group = file
        .create_group(path)
        .with_context(|| format!("creating group '{path}'"))?;
    set_string_attr(&*group, "encoding-type", "dict")?;
    set_string_attr(&*group, "encoding-version", "0.1.0")?;
    Ok(group)
}

fn write_dense_matrix(
    parent: &hdf5::Group,
    name: &str,
    data: &Array2<f32>,
    overwrite: bool,
    require_flag: bool,
    encoding_type: &str,
    encoding_version: &str,
) -> Result<()> {
    let contiguous = data.as_standard_layout().to_owned();
    if parent.link_exists(name) {
        if !overwrite && require_flag {
            bail!(
                "refusing to replace existing '{}/{}' without --overwrite-derived",
                parent.name(),
                name
            );
        }
        parent
            .unlink(name)
            .with_context(|| format!("removing '{}/{}'", parent.name(), name))?;
    }
    let dataset = parent
        .new_dataset_builder()
        .with_data(&contiguous)
        .create(name)
        .with_context(|| format!("creating dense matrix '{}/{}'", parent.name(), name))?;
    set_string_attr(&**dataset, "encoding-type", encoding_type)?;
    set_string_attr(&**dataset, "encoding-version", encoding_version)?;
    Ok(())
}

fn write_csr_matrix(
    parent: &hdf5::Group,
    name: &str,
    matrix: &CsrMatrix,
    overwrite: bool,
    require_flag: bool,
) -> Result<()> {
    if parent.link_exists(name) {
        if !overwrite && require_flag {
            bail!(
                "refusing to replace existing '{}/{}' without --overwrite-derived",
                parent.name(),
                name
            );
        }
        parent
            .unlink(name)
            .with_context(|| format!("removing '{}/{}'", parent.name(), name))?;
    }
    let group = parent
        .create_group(name)
        .with_context(|| format!("creating csr group '{}/{}'", parent.name(), name))?;
    set_string_attr(&*group, "encoding-type", "csr_matrix")?;
    set_string_attr(&*group, "encoding-version", "0.1.0")?;
    set_shape_attr(&group, matrix.nrows, matrix.ncols)?;
    group
        .new_dataset_builder()
        .with_data(&matrix.data)
        .create("data")
        .context("writing csr data")?;
    group
        .new_dataset_builder()
        .with_data(&matrix.indices)
        .create("indices")
        .context("writing csr indices")?;
    group
        .new_dataset_builder()
        .with_data(&matrix.indptr)
        .create("indptr")
        .context("writing csr indptr")?;
    Ok(())
}

fn set_shape_attr(group: &hdf5::Group, nrows: usize, ncols: usize) -> Result<()> {
    if group.attr("shape").is_ok() {
        group
            .delete_attr("shape")
            .context("deleting prior shape attribute")?;
    }
    group
        .new_attr_builder()
        .with_data(&[nrows as i64, ncols as i64])
        .create("shape")
        .context("writing shape attribute")?;
    Ok(())
}

fn write_string_scalar(
    parent: &hdf5::Group,
    name: &str,
    value: &str,
    overwrite: bool,
) -> Result<()> {
    if parent.link_exists(name) {
        if overwrite {
            parent
                .unlink(name)
                .with_context(|| format!("removing '{}/{}'", parent.name(), name))?;
        } else {
            bail!("dataset '{}/{}' already exists", parent.name(), name);
        }
    }
    let encoded = VarLenUnicode::from_str(value).context("encoding string scalar")?;
    let dataset = parent
        .new_dataset::<VarLenUnicode>()
        .shape(())
        .create(name)
        .with_context(|| format!("creating scalar string '{}/{}'", parent.name(), name))?;
    dataset
        .write_scalar(&encoded)
        .context("writing scalar string")?;
    set_string_attr(&**dataset, "encoding-type", "string")?;
    set_string_attr(&**dataset, "encoding-version", "0.2.0")?;
    Ok(())
}

fn write_numeric_scalar(
    parent: &hdf5::Group,
    name: &str,
    value: f64,
    overwrite: bool,
) -> Result<()> {
    if parent.link_exists(name) {
        if overwrite {
            parent
                .unlink(name)
                .with_context(|| format!("removing '{}/{}'", parent.name(), name))?;
        } else {
            bail!("dataset '{}/{}' already exists", parent.name(), name);
        }
    }
    let dataset = parent
        .new_dataset::<f64>()
        .shape(())
        .create(name)
        .with_context(|| format!("creating scalar numeric '{}/{}'", parent.name(), name))?;
    dataset
        .write_scalar(&value)
        .context("writing numeric scalar")?;
    set_string_attr(&**dataset, "encoding-type", "numeric-scalar")?;
    set_string_attr(&**dataset, "encoding-version", "0.2.0")?;
    Ok(())
}

fn write_string_array_path(
    file: &hdf5::File,
    path: &str,
    values: &[String],
    overwrite: bool,
    require_flag: bool,
) -> Result<()> {
    let (parent_path, name) = split_parent(path)?;
    let parent = ensure_dict_group(file, parent_path)?;
    if parent.link_exists(name) {
        if !overwrite && require_flag {
            bail!("refusing to replace existing '{path}' without --overwrite-derived");
        }
        parent
            .unlink(name)
            .with_context(|| format!("removing '{path}'"))?;
    }
    let encoded: Vec<VarLenUnicode> = values
        .iter()
        .map(|value| VarLenUnicode::from_str(value))
        .collect::<std::result::Result<Vec<_>, _>>()
        .context("encoding string array")?;
    let dataset = parent
        .new_dataset_builder()
        .with_data(&encoded)
        .create(name)
        .with_context(|| format!("creating string-array '{path}'"))?;
    set_string_attr(&**dataset, "encoding-type", "string-array")?;
    set_string_attr(&**dataset, "encoding-version", "0.2.0")?;
    Ok(())
}

fn split_parent(path: &str) -> Result<(&str, &str)> {
    let (parent, name) = path
        .rsplit_once('/')
        .with_context(|| format!("path '{path}' does not contain a parent group"))?;
    if parent.is_empty() || name.is_empty() {
        bail!("invalid group path '{path}'");
    }
    Ok((parent, name))
}

fn set_string_attr(location: &hdf5::Location, name: &str, value: &str) -> Result<()> {
    if location.attr(name).is_ok() {
        location
            .delete_attr(name)
            .with_context(|| format!("deleting attr '{name}'"))?;
    }
    let encoded = VarLenUnicode::from_str(value).context("encoding string attribute")?;
    let attr = location
        .new_attr::<VarLenUnicode>()
        .shape(())
        .create(name)
        .with_context(|| format!("creating attr '{name}'"))?;
    attr.write_scalar(&encoded)
        .with_context(|| format!("writing attr '{name}'"))?;
    Ok(())
}
