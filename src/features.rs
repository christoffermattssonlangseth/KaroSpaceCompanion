use std::collections::{BTreeSet, HashMap};

use ndarray::Array2;

use crate::anndata::ExpressionMatrix;

pub fn normalize_counts(matrix: &ExpressionMatrix, target_sum: f32, log1p: bool) -> Array2<f32> {
    match matrix {
        ExpressionMatrix::Dense(dense) => {
            let mut normalized = dense.clone();
            for mut row in normalized.rows_mut() {
                let total: f32 = row.iter().copied().sum();
                if total > 0.0 {
                    let scale = target_sum / total;
                    for value in &mut row {
                        *value *= scale;
                        if log1p {
                            *value = value.ln_1p();
                        }
                    }
                }
            }
            normalized
        }
        ExpressionMatrix::Sparse(sparse) => {
            let mut row_sums = vec![0.0f32; sparse.nrows];
            for (row, row_sum) in row_sums.iter_mut().enumerate() {
                sparse.for_each_nonzero_in_row(row, |_, value| {
                    *row_sum += value;
                });
            }

            let mut normalized = Array2::<f32>::zeros((sparse.nrows, sparse.ncols));
            for (row, total) in row_sums.into_iter().enumerate() {
                if total <= 0.0 {
                    continue;
                }
                let scale = target_sum / total;
                sparse.for_each_nonzero_in_row(row, |col, value| {
                    let mut scaled = value * scale;
                    if log1p {
                        scaled = scaled.ln_1p();
                    }
                    normalized[[row, col]] = scaled;
                });
            }
            normalized
        }
    }
}

pub fn aggregate_neighbors(neighbors: &[Vec<(usize, f32)>], source: &Array2<f32>) -> Array2<f32> {
    let mut output = Array2::<f32>::zeros((source.nrows(), source.ncols()));
    for (row_index, row_neighbors) in neighbors.iter().enumerate() {
        if row_neighbors.is_empty() {
            continue;
        }
        let weight = 1.0 / row_neighbors.len() as f32;
        for &(neighbor_index, _) in row_neighbors {
            for col in 0..source.ncols() {
                output[[row_index, col]] += source[[neighbor_index, col]] * weight;
            }
        }
    }
    output
}

pub fn compute_composition_matrix(
    neighbors: &[Vec<(usize, f32)>],
    cell_types: &[String],
) -> (Array2<f32>, Vec<String>) {
    let categories: Vec<String> = cell_types
        .iter()
        .filter(|value| !value.is_empty())
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let category_to_index: HashMap<&str, usize> = categories
        .iter()
        .enumerate()
        .map(|(index, value)| (value.as_str(), index))
        .collect();

    let mut output = Array2::<f32>::zeros((cell_types.len(), categories.len()));
    for (cell_index, row_neighbors) in neighbors.iter().enumerate() {
        if row_neighbors.is_empty() {
            continue;
        }
        let mut counts = vec![0usize; categories.len()];
        for &(neighbor_index, _) in row_neighbors {
            if let Some(&category_index) =
                category_to_index.get(cell_types[neighbor_index].as_str())
            {
                counts[category_index] += 1;
            }
        }
        let total: usize = counts.iter().sum();
        if total == 0 {
            continue;
        }
        for (category_index, count) in counts.into_iter().enumerate() {
            output[[cell_index, category_index]] = count as f32 / total as f32;
        }
    }

    (output, categories)
}

#[cfg(test)]
mod tests {
    use ndarray::arr2;

    use super::{aggregate_neighbors, compute_composition_matrix, normalize_counts};
    use crate::anndata::ExpressionMatrix;

    #[test]
    fn normalize_counts_log1p_scales_rows() {
        let matrix = arr2(&[[1.0f32, 1.0], [0.0, 4.0]]);
        let normalized = normalize_counts(&ExpressionMatrix::Dense(matrix), 10.0, false);
        assert_eq!(normalized[[0, 0]], 5.0);
        assert_eq!(normalized[[0, 1]], 5.0);
        assert_eq!(normalized[[1, 1]], 10.0);
    }

    #[test]
    fn aggregate_neighbors_returns_neighbor_mean() {
        let source = arr2(&[[1.0f32, 3.0], [5.0, 7.0]]);
        let neighbors = vec![vec![(1usize, 1.0)], vec![(0usize, 1.0)]];
        let aggregated = aggregate_neighbors(&neighbors, &source);
        assert_eq!(aggregated[[0, 0]], 5.0);
        assert_eq!(aggregated[[1, 1]], 3.0);
    }

    #[test]
    fn composition_matrix_uses_sorted_categories() {
        let neighbors = vec![vec![(1usize, 1.0), (2usize, 1.0)]];
        let cell_types = vec!["A".to_string(), "B".to_string(), "A".to_string()];
        let (matrix, categories) = compute_composition_matrix(&neighbors, &cell_types);
        assert_eq!(categories, vec!["A".to_string(), "B".to_string()]);
        assert_eq!(matrix[[0, 0]], 0.5);
        assert_eq!(matrix[[0, 1]], 0.5);
    }
}
