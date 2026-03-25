use std::collections::HashMap;

use anyhow::{bail, Result};
use rstar::{PointDistance, RTree, RTreeObject, AABB};
use spade::{DelaunayTriangulation, HasPosition, Point2, Triangulation};

use crate::anndata::CsrMatrix;

#[derive(Debug, Clone)]
pub enum GraphMode {
    Radius(f64),
    Knn(usize),
    Delaunay { remove_long_links_percentile: f64 },
}

impl GraphMode {
    pub fn label(&self) -> &'static str {
        match self {
            Self::Radius(_) => "radius",
            Self::Knn(_) => "knn",
            Self::Delaunay { .. } => "delaunay",
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

#[derive(Clone, Copy, Debug)]
struct IndexedVertex {
    position: Point2<f64>,
    global_index: usize,
}

impl HasPosition for IndexedVertex {
    type Scalar = f64;

    fn position(&self) -> Point2<Self::Scalar> {
        self.position
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
    validate_graph_mode(&mode)?;

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
    if let GraphMode::Delaunay {
        remove_long_links_percentile,
    } = &mode
    {
        trim_long_links(&mut undirected, *remove_long_links_percentile);
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

fn validate_graph_mode(mode: &GraphMode) -> Result<()> {
    match mode {
        GraphMode::Radius(radius) if !radius.is_finite() || *radius <= 0.0 => {
            bail!("radius must be a finite value > 0")
        }
        GraphMode::Knn(k) if *k == 0 => bail!("k must be greater than 0"),
        GraphMode::Delaunay {
            remove_long_links_percentile,
        } if !remove_long_links_percentile.is_finite()
            || !(0.0..=100.0).contains(remove_long_links_percentile) =>
        {
            bail!("remove-long-links-percentile must be a finite value between 0 and 100")
        }
        _ => Ok(()),
    }
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

    match mode {
        GraphMode::Radius(radius) => build_radius_edges(indices, coordinates, *radius, target),
        GraphMode::Knn(k) => build_knn_edges(indices, coordinates, *k, target),
        GraphMode::Delaunay { .. } => build_delaunay_edges(indices, coordinates, target),
    }
}

fn build_radius_edges(
    indices: &[usize],
    coordinates: &[[f64; 2]],
    radius: f64,
    target: &mut HashMap<(usize, usize), f32>,
) -> Result<()> {
    let points: Vec<IndexedPoint> = indices
        .iter()
        .enumerate()
        .map(|(local_index, &global_index)| IndexedPoint {
            coords: coordinates[global_index],
            local_index,
        })
        .collect();
    let tree = RTree::bulk_load(points);
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
    Ok(())
}

fn build_knn_edges(
    indices: &[usize],
    coordinates: &[[f64; 2]],
    k: usize,
    target: &mut HashMap<(usize, usize), f32>,
) -> Result<()> {
    let points: Vec<IndexedPoint> = indices
        .iter()
        .enumerate()
        .map(|(local_index, &global_index)| IndexedPoint {
            coords: coordinates[global_index],
            local_index,
        })
        .collect();
    let tree = RTree::bulk_load(points);
    for (local_i, &global_i) in indices.iter().enumerate() {
        let coord = coordinates[global_i];
        for (point, dist_sq) in tree
            .nearest_neighbor_iter_with_distance_2(&coord)
            .filter(|(point, _)| point.local_index != local_i)
            .take(k)
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
    Ok(())
}

fn build_delaunay_edges(
    indices: &[usize],
    coordinates: &[[f64; 2]],
    target: &mut HashMap<(usize, usize), f32>,
) -> Result<()> {
    let vertices: Vec<IndexedVertex> = indices
        .iter()
        .map(|&global_index| IndexedVertex {
            position: Point2::new(coordinates[global_index][0], coordinates[global_index][1]),
            global_index,
        })
        .collect();
    let triangulation = DelaunayTriangulation::<IndexedVertex>::bulk_load(vertices)?;
    for edge in triangulation.undirected_edges() {
        let [from_vertex, to_vertex] = edge.vertices();
        let from = from_vertex.data().global_index;
        let to = to_vertex.data().global_index;
        let pair = if from < to { (from, to) } else { (to, from) };
        let [from_pos, to_pos] = edge.positions();
        let distance = from_pos.distance_2(to_pos).sqrt() as f32;
        target
            .entry(pair)
            .and_modify(|existing| {
                if distance < *existing {
                    *existing = distance;
                }
            })
            .or_insert(distance);
    }
    Ok(())
}

fn trim_long_links(edges: &mut HashMap<(usize, usize), f32>, percentile: f64) {
    if edges.is_empty() {
        return;
    }

    let mut distances: Vec<f32> = edges
        .values()
        .copied()
        .filter(|value| *value > 0.0)
        .collect();
    if distances.is_empty() {
        return;
    }
    distances.sort_by(|left, right| left.total_cmp(right));
    let threshold = percentile_value(&distances, percentile);
    edges.retain(|_, distance| *distance <= threshold);
}

fn percentile_value(sorted: &[f32], percentile: f64) -> f32 {
    debug_assert!(!sorted.is_empty());
    if sorted.len() == 1 {
        return sorted[0];
    }

    let clamped = percentile.clamp(0.0, 100.0);
    let rank = (clamped / 100.0) * ((sorted.len() - 1) as f64);
    let lower = rank.floor() as usize;
    let upper = rank.ceil() as usize;
    if lower == upper {
        return sorted[lower];
    }
    let weight = (rank - lower as f64) as f32;
    sorted[lower] + (sorted[upper] - sorted[lower]) * weight
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

    #[test]
    fn delaunay_graph_respects_group_boundaries() {
        let coords = arr2(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [10.0, 0.0],
            [11.0, 0.0],
            [10.5, 1.0],
        ]);
        let groups = vec![
            "A".to_string(),
            "A".to_string(),
            "A".to_string(),
            "B".to_string(),
            "B".to_string(),
            "B".to_string(),
        ];
        let graph = build_spatial_graph(
            &coords,
            &groups,
            GraphMode::Delaunay {
                remove_long_links_percentile: 100.0,
            },
        )
        .unwrap();
        assert_eq!(graph.n_undirected_edges, 6);
        assert!(graph.neighbors[0].iter().all(|(idx, _)| *idx < 3));
        assert!(graph.neighbors[3].iter().all(|(idx, _)| *idx >= 3));
    }

    #[test]
    fn delaunay_long_link_trim_removes_percentile_tail() {
        let coords = arr2(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [1.5, 100.0]]);
        let groups = vec!["A".to_string(); coords.nrows()];

        let full = build_spatial_graph(
            &coords,
            &groups,
            GraphMode::Delaunay {
                remove_long_links_percentile: 100.0,
            },
        )
        .unwrap();
        let trimmed = build_spatial_graph(
            &coords,
            &groups,
            GraphMode::Delaunay {
                remove_long_links_percentile: 80.0,
            },
        )
        .unwrap();

        assert!(trimmed.n_undirected_edges < full.n_undirected_edges);
        let max_full = full
            .neighbors
            .iter()
            .flat_map(|row| row.iter().map(|(_, distance)| *distance))
            .fold(0.0f32, f32::max);
        let max_trimmed = trimmed
            .neighbors
            .iter()
            .flat_map(|row| row.iter().map(|(_, distance)| *distance))
            .fold(0.0f32, f32::max);
        assert!(max_trimmed < max_full);
    }
}
