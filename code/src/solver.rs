use crate::{
    config::Config,
    field::{C, Field, Grid, field_error, gaussian, intensity_error},
    propagation::Cavity,
};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct History {
    pub iteration: usize,
    pub delta: f64,
    pub intensity_difference: f64,
    pub gamma_real: f64,
    pub gamma_imag: f64,
    pub gamma_abs: f64,
    pub gamma_phase: f64,
    pub rho: f64,
    pub power_retention: f64,
    pub loss: f64,
    pub eigenvalue_difference: Option<f64>,
}
pub struct Solution {
    pub field: Field,
    pub gamma: C,
    pub iterations: usize,
    pub converged: bool,
    pub history: Vec<History>,
}
pub fn solve(mut f: Field, cavity: &Cavity, c: &Config) -> Solution {
    f.normalise();
    let mut history = Vec::new();
    let mut previous = None;
    let mut stable = 0usize;
    let mut gamma = C::new(0.0, 0.0);
    for iteration in 1..=c.fox_li.max_iterations {
        let raw = cavity.round_trip(&f);
        gamma = f.dot(&raw) / f.dot(&f);
        let rho = raw.norm() / f.norm();
        let mut next = raw;
        next.normalise();
        let delta = field_error(&f, &next);
        let di = intensity_error(&f, &next);
        let dg = previous.map(|old: C| (gamma - old).norm() / gamma.norm().max(f64::EPSILON));
        history.push(History {
            iteration,
            delta,
            intensity_difference: di,
            gamma_real: gamma.re,
            gamma_imag: gamma.im,
            gamma_abs: gamma.norm(),
            gamma_phase: gamma.arg(),
            rho,
            power_retention: gamma.norm_sqr(),
            loss: (1.0 - gamma.norm_sqr()).max(0.0),
            eigenvalue_difference: dg,
        });
        if delta < c.fox_li.convergence_threshold
            && di < c.fox_li.convergence_threshold
            && dg.is_some_and(|e| e < c.fox_li.eigenvalue_relative_threshold)
        {
            stable += 1;
        } else {
            stable = 0;
        }
        f = next;
        previous = Some(gamma);
        if stable >= c.fox_li.consecutive_converged {
            return Solution {
                field: f,
                gamma,
                iterations: iteration,
                converged: true,
                history,
            };
        }
    }
    Solution {
        field: f,
        gamma,
        iterations: c.fox_li.max_iterations,
        converged: false,
        history,
    }
}
pub fn initial(c: &Config, g: &Grid, two_d: bool, theta: f64) -> Field {
    let x = gaussian(c, g, 0.0, 0.0);
    let y = if two_d {
        let s = theta.sin();
        let (offset, slope) = if c.tilt.mirror == "M2" {
            (c.project.cavity_length_m * s, -s)
        } else {
            (0.0, s)
        };
        Some(gaussian(c, g, offset, slope))
    } else {
        None
    };
    Field { x, y }
}

pub struct DenseSolution {
    pub field: Vec<C>,
    pub gamma: C,
    pub iterations: usize,
    pub converged: bool,
    pub history: Vec<History>,
}
/// Arbitrary dense complex initial fields use the full 2-D propagation, never a rank-one approximation.
pub fn solve_dense(mut f: Vec<C>, cavity: &Cavity, c: &Config, two_d: bool) -> DenseSolution {
    crate::field::normalise(&mut f);
    let mut history = Vec::new();
    let mut previous = None;
    let mut stable = 0;
    let mut gamma = C::new(0.0, 0.0);
    for iteration in 1..=c.fox_li.max_iterations {
        let mut next = if two_d {
            cavity.full_round_trip(&f)
        } else {
            cavity
                .round_trip(&Field {
                    x: f.clone(),
                    y: None,
                })
                .x
        };
        gamma = crate::field::dot(&f, &next) / crate::field::dot(&f, &f);
        let rho = crate::field::norm(&next) / crate::field::norm(&f);
        crate::field::normalise(&mut next);
        let overlap = crate::field::dot(&f, &next);
        let phase = if overlap.norm() > f64::EPSILON {
            overlap / overlap.norm()
        } else {
            C::new(1.0, 0.0)
        };
        let delta = f
            .iter()
            .zip(&next)
            .map(|(a, b)| (*b - phase * *a).norm_sqr())
            .sum::<f64>()
            .sqrt()
            / crate::field::norm(&next);
        let di = f
            .iter()
            .zip(&next)
            .map(|(a, b)| (b.norm_sqr() - a.norm_sqr()).powi(2))
            .sum::<f64>()
            .sqrt()
            / next
                .iter()
                .map(|z| z.norm_sqr().powi(2))
                .sum::<f64>()
                .sqrt();
        let dg = previous.map(|old: C| (gamma - old).norm() / gamma.norm().max(f64::EPSILON));
        history.push(History {
            iteration,
            delta,
            intensity_difference: di,
            gamma_real: gamma.re,
            gamma_imag: gamma.im,
            gamma_abs: gamma.norm(),
            gamma_phase: gamma.arg(),
            rho,
            power_retention: gamma.norm_sqr(),
            loss: (1.0 - gamma.norm_sqr()).max(0.0),
            eigenvalue_difference: dg,
        });
        if delta < c.fox_li.convergence_threshold
            && di < c.fox_li.convergence_threshold
            && dg.is_some_and(|e| e < c.fox_li.eigenvalue_relative_threshold)
        {
            stable += 1
        } else {
            stable = 0
        }
        f = next;
        previous = Some(gamma);
        if stable >= c.fox_li.consecutive_converged {
            return DenseSolution {
                field: f,
                gamma,
                iterations: iteration,
                converged: true,
                history,
            };
        }
    }
    DenseSolution {
        field: f,
        gamma,
        iterations: c.fox_li.max_iterations,
        converged: false,
        history,
    }
}
pub fn random_field(n: usize, seed: u64) -> Vec<C> {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let r = (-2.0 * rng.gen_range(f64::MIN_POSITIVE..1.0).ln()).sqrt();
            let phi = std::f64::consts::TAU * rng.gen_range(0.0..1.0);
            C::from_polar(r, phi)
        })
        .collect()
}
