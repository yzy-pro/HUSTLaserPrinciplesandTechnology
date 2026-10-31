//! 陈恪瑾 / U202413925 / 光电信息科学与工程202409班 — Rust Fox–Li numerical core.
mod config;
mod experiment;
mod field;
mod propagation;
mod solver;

fn main() {
    if let Err(error) = experiment::run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{C, Field, Grid};
    #[test]
    fn propagation_matches_integral_and_dense_roundtrip() {
        let c = config::Config::load().unwrap();
        let mut spec = c.grid.two_d.clone();
        spec.n = c.validation.direct_grid_n;
        let g = Grid::new(&c, &spec);
        let cavity = propagation::Cavity::new(&g, &c, 0.001);
        let f = solver::initial(&c, &g, true, 0.001);
        let direct = propagation::direct_1d(&f.x, &g, &c);
        let fft = cavity.prop.axis(&f.x, true);
        let err = direct
            .iter()
            .zip(&fft)
            .map(|(a, b)| (*a - *b).norm_sqr())
            .sum::<f64>()
            .sqrt()
            / field::norm(&direct);
        assert!(err < 1e-8, "{err}");
        let full: Vec<C> = (0..g.n())
            .flat_map(|y| {
                (0..g.n()).map({
                    let f = &f;
                    move |x| f.value(y, x)
                })
            })
            .collect();
        let dense = cavity.full_round_trip(&full);
        let n = g.n();
        let factored = cavity.round_trip(&f);
        let err = (0..g.n())
            .flat_map(|y| {
                (0..g.n()).map({
                    let factored = &factored;
                    let dense = &dense;
                    move |x| (factored.value(y, x) - dense[y * n + x]).norm_sqr()
                })
            })
            .sum::<f64>()
            .sqrt()
            / field::norm(&dense);
        assert!(err < 1e-12, "{err}");
    }
    #[test]
    fn dense_random_solver_preserves_nonseparable_field() {
        let mut c = config::Config::load().unwrap();
        c.fox_li.max_iterations = 2;
        c.fox_li.consecutive_converged = 3;
        let mut spec = c.grid.two_d.clone();
        spec.n = 16;
        let g = Grid::new(&c, &spec);
        let cavity = propagation::Cavity::new(&g, &c, 0.0);
        let mut f = solver::random_field(g.n() * g.n(), c.fox_li.random_seed);
        field::normalise(&mut f);
        let mut expected = cavity.full_round_trip(&f);
        field::normalise(&mut expected);
        expected = cavity.full_round_trip(&expected);
        field::normalise(&mut expected);
        let result = solver::solve_dense(f, &cavity, &c, true);
        assert_eq!(result.iterations, 2);
        assert!(!result.converged);
        assert!(1.0 - field::overlap(&expected, &result.field) < 1e-12);
    }
    #[test]
    fn separable_metrics_match_dense_arrays_near_convergence() {
        let c = config::Config::load().unwrap();
        let g = Grid::new(&c, &c.grid.two_d);
        let a = solver::initial(&c, &g, true, 0.0);
        let b = solver::initial(&c, &g, true, 1e-10);
        let flat = |f: &Field| -> Vec<C> {
            (0..g.n())
                .flat_map(|y| (0..g.n()).map(move |x| f.value(y, x)))
                .collect()
        };
        let av = flat(&a);
        let bv = flat(&b);
        let ov = field::dot(&av, &bv);
        let phase = ov / ov.norm();
        let d = av
            .iter()
            .zip(&bv)
            .map(|(x, y)| (*y - phase * *x).norm_sqr())
            .sum::<f64>()
            .sqrt();
        assert!((field::field_error(&a, &b) - d).abs() < 1e-12);
    }
    #[test]
    fn rejects_invalid_configuration() {
        let mut c = config::Config::load().unwrap();
        c.project.mirror1_radius_m *= 2.0;
        assert!(c.validate().is_err());
        let mut c = config::Config::load().unwrap();
        c.tilt.angles_mrad.push(c.tilt.angles_mrad[0]);
        assert!(c.validate().is_err());
    }
}
