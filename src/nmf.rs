use anyhow::{bail, Result};
use ndarray::{Array2, Zip};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

pub struct NmfConfig {
    pub n_components: usize,
    pub max_iter: usize,
    pub tol: f32,
    pub seed: u64,
    pub epsilon: f32,
}

pub struct NmfResult {
    pub w: Array2<f32>,
    pub h: Array2<f32>,
    pub n_iter: usize,
    pub final_error: f32,
}

pub fn run_nmf(x: &Array2<f32>, config: &NmfConfig) -> Result<NmfResult> {
    let (n_obs, n_var) = x.dim();
    if n_obs == 0 || n_var == 0 {
        bail!("expression matrix must be non-empty");
    }
    if config.n_components == 0 {
        bail!("n_components must be greater than 0");
    }
    if x.iter().any(|value| *value < 0.0) {
        bail!("NMF requires a non-negative matrix");
    }

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut w = Array2::from_shape_fn((n_obs, config.n_components), |_| {
        rng.random::<f32>().max(config.epsilon)
    });
    let mut h = Array2::from_shape_fn((config.n_components, n_var), |_| {
        rng.random::<f32>().max(config.epsilon)
    });

    let mut previous_error = f32::INFINITY;
    let mut final_error = f32::INFINITY;
    let mut n_iter = config.max_iter;

    for iteration in 0..config.max_iter {
        let wtx = w.t().dot(x);
        let wtwh = w.t().dot(&w).dot(&h);
        Zip::from(&mut h)
            .and(&wtx)
            .and(&wtwh)
            .for_each(|h_value, &numerator, &denominator| {
                *h_value *= numerator / (denominator + config.epsilon);
                if *h_value < config.epsilon {
                    *h_value = config.epsilon;
                }
            });

        let xht = x.dot(&h.t());
        let whht = w.dot(&h.dot(&h.t()));
        Zip::from(&mut w)
            .and(&xht)
            .and(&whht)
            .for_each(|w_value, &numerator, &denominator| {
                *w_value *= numerator / (denominator + config.epsilon);
                if *w_value < config.epsilon {
                    *w_value = config.epsilon;
                }
            });

        if (iteration + 1) % 10 == 0 || iteration + 1 == config.max_iter {
            final_error = frobenius_error(x, &w, &h);
            if (previous_error - final_error).abs() < config.tol {
                n_iter = iteration + 1;
                break;
            }
            previous_error = final_error;
        }
    }

    Ok(NmfResult {
        w,
        h,
        n_iter,
        final_error,
    })
}

fn frobenius_error(x: &Array2<f32>, w: &Array2<f32>, h: &Array2<f32>) -> f32 {
    let wh = w.dot(h);
    Zip::from(x)
        .and(&wh)
        .fold(0.0f32, |acc, &expected, &actual| {
            let diff = expected - actual;
            acc + diff * diff
        })
        .sqrt()
}

#[cfg(test)]
mod tests {
    use ndarray::arr2;

    use super::{run_nmf, NmfConfig};

    #[test]
    fn nmf_runs_on_simple_non_negative_input() {
        let matrix = arr2(&[[1.0f32, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        let config = NmfConfig {
            n_components: 2,
            max_iter: 20,
            tol: 1e-4,
            seed: 1,
            epsilon: 1e-12,
        };
        let result = run_nmf(&matrix, &config).unwrap();
        assert_eq!(result.w.nrows(), 3);
        assert_eq!(result.h.ncols(), 2);
    }
}
