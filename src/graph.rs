use std::collections::HashMap;

use anyhow::{bail, Result};
use rstar::{PointDistance, RTree, RTreeObject, AABB};

use crate::anndata::CsrMatrix;

#[derive(Debug, Clone)]
pub enum GraphMode {
    Radius(f64),
    Knn(usize),
}

impl GraphMode {
    pub fn label(&self) -> &'static str {
        match self {
            Self::Radius(_) => "radius",
            Self::Knn(_) => "knn",
        }
    }
}

pub struct SpatialGraph {
    pub neighbors: Vec<Vec<(usize, f32)>>,
    pub connectivities: CsrMatrix,
    pub distances: CsrMatrix,
    pub n_undirected_edges: usize,
}

#[derive(Clone)]
struct IndexedPoint {
    coords: [f64; 2],
    local_index: usize,
}

impl RTreeObject for IndexedPoint {
    type Envelope = AABB<[f64; 2]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(self.coords)
    }
}

impl PointDistance for IndexedPoint {
    fn distance_2(&self, point: &[f64; 2]) -> f64 {
        let dx = self.coords[0] - point[0];
        let dy = self.coords[1] - point[1];
        dx * dx + dy * dy
    }
}

pub fn build_spatial_graph(
    coordinates: &ndarray::Array2<f64>,
    groups: &[String],
    mode: GraphMode,
) -> Result<SpatialGraph> {
    if coordinates.ncols() < 2 {
        bail!("coordinate matrix must have at least 2 columns");
    }
    if coordinates.nrows() != groups.len() {
        bail!(
            "coordinate row count ({}) does not match group count ({})",
            coordinates.nrows(),
            groups.len()
        );
    }
    match mode {
        GraphMode::Radius(radius) if !radius.is_finite() || radius <= 0.0 => {
            bail!("radius must be a finite value > 0")
        }
        GraphMode::Knn(k) if k == 0 => bail!("k must be greater than 0"),
        _ => {}
    }

    let coords: Vec<[f64; 2]> = coordinates
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1]])
        .collect();

    let mut grouped_indices: HashMap<&str, Vec<usize>> = HashMap::new();
    for (index, group) in groups.iter().enumerate() {
        grouped_indices
            .entry(group.as_str())
            .or_default()
            .push(index);
    }

    let mut undirected: HashMap<(usize, usize), f32> = HashMap::new();
    for indices in grouped_indices.values() {
        build_group_edges(indices, &coords, &mode, &mut undirected)?;
    }

    let mut neighbors = vec![Vec::<(usize, f32)>::new(); coords.len()];
    for ((i, j), distance) in &undirected {
        neighbors[*i].push((*j, *distance));
        neighbors[*j].push((*i, *distance));
    }
    for row in &mut neighbors {
        row.sort_by_key(|(index, _)| *index);
    }

    let connectivities = csr_from_neighbors(&neighbors, |_, _| 1.0);
    let distances = csr_from_neighbors(&neighbors, |_, distance| distance);

    Ok(SpatialGraph {
        neighbors,
        connectivities,
        distances,
        n_undirected_edges: undirected.len(),
    })
}

fn build_group_edges(
    indices: &[usize],
    coordinates: &[[f64; 2]],
    mode: &GraphMode,
    target: &mut HashMap<(usize, usize), f32>,
) -> Result<()> {
    if indices.len() < 2 {
        return Ok(());
    }

    let points: Vec<IndexedPoint> = indices
        .iter()
        .enumerate()
        .map(|(local_index, &global_index)| IndexedPoint {
            coords: coordinates[global_index],
            local_index,
        })
        .collect();
    let tree = RTree::bulk_load(points);

    match mode {
        GraphMode::Radius(radius) => {
            let radius_sq = radius * radius;
            for (local_i, &global_i) in indices.iter().enumerate() {
                let coord = coordinates[global_i];
                for point in tree.locate_within_distance(coord, radius_sq) {
                    if point.local_index <= local_i {
                        continue;
                    }
                    let global_j = indices[point.local_index];
                    let distance = point.distance_2(&coord).sqrt() as f32;
                    target.insert((global_i, global_j), distance);
                }
            }
        }
        GraphMode::Knn(k) => {
            for (local_i, &global_i) in indices.iter().enumerate() {
                let coord = coordinates[global_i];
                for (point, dist_sq) in tree
                    .nearest_neighbor_iter_with_distance_2(&coord)
                    .filter(|(point, _)| point.local_index != local_i)
                    .take(*k)
                {
                    let global_j = indices[point.local_index];
                    let pair = if global_i < global_j {
                        (global_i, global_j)
                    } else {
                        (global_j, global_i)
                    };
                    let distance = dist_sq.sqrt() as f32;
                    target
                        .entry(pair)
                        .and_modify(|existing| {
                            if distance < *existing {
                                *existing = distance;
                            }
                        })
                        .or_insert(distance);
                }
            }
        }
    }

    Ok(())
}

fn csr_from_neighbors(
    neighbors: &[Vec<(usize, f32)>],
    value_for_edge: impl Fn(usize, f32) -> f32,
) -> CsrMatrix {
    let n = neighbors.len();
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = Vec::with_capacity(n + 1);
    indptr.push(0i32);

    for row in neighbors {
        for &(col, distance) in row {
            data.push(value_for_edge(col, distance));
            indices.push(col as i32);
        }
        indptr.push(data.len() as i32);
    }

    CsrMatrix {
        data,
        indices,
        indptr,
        nrows: n,
        ncols: n,
    }
}

#[cfg(test)]
mod tests {
    use ndarray::arr2;

    use super::{build_spatial_graph, GraphMode};

    #[test]
    fn radius_graph_respects_group_boundaries() {
        let coords = arr2(&[[0.0, 0.0], [1.0, 0.0], [0.2, 0.0], [1.2, 0.0]]);
        let groups = vec![
            "A".to_string(),
            "A".to_string(),
            "B".to_string(),
            "B".to_string(),
        ];
        let graph = build_spatial_graph(&coords, &groups, GraphMode::Radius(1.1)).unwrap();
        assert_eq!(graph.n_undirected_edges, 2);
        assert_eq!(graph.neighbors[0].len(), 1);
        assert_eq!(graph.neighbors[2].len(), 1);
    }
}
