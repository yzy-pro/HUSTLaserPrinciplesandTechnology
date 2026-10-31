use crate::config::{Config, GridSpec};
pub type C = rustfft::num_complex::Complex64;
use std::f64::consts::PI;

#[derive(Clone)]
pub struct Grid {
    pub x: Vec<f64>,
    pub dx: f64,
    pub aperture: f64,
    pub half_width: f64,
    pub padding: usize,
}
impl Grid {
    pub fn new(c: &Config, s: &GridSpec) -> Self {
        let wm = (c.project.wavelength_m * c.project.cavity_length_m / PI).sqrt();
        let half_width = s.window_half_width_factor * wm;
        let dx = 2.0 * half_width / s.n as f64;
        Self {
            x: (0..s.n)
                .map(|j| (j as f64 - s.n as f64 / 2.0) * dx)
                .collect(),
            dx,
            aperture: s.aperture_half_width_factor * wm,
            half_width,
            padding: c.grid.propagation_padding_factor,
        }
    }
    pub fn n(&self) -> usize {
        self.x.len()
    }
}
pub fn dot(a: &[C], b: &[C]) -> C {
    a.iter().zip(b).map(|(x, y)| x.conj() * y).sum()
}
pub fn norm(a: &[C]) -> f64 {
    a.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt()
}
pub fn normalise(a: &mut [C]) {
    let n = norm(a);
    assert!(n.is_finite() && n > 0.0, "zero or non-finite field");
    for z in a {
        *z /= n;
    }
}
pub fn overlap(a: &[C], b: &[C]) -> f64 {
    (dot(a, b).norm_sqr() / (dot(a, a).re * dot(b, b).re)).clamp(0.0, 1.0)
}
pub fn gaussian(c: &Config, g: &Grid, offset: f64, slope: f64) -> Vec<C> {
    let wm2 = c.project.wavelength_m * c.project.cavity_length_m / PI;
    let k = 2.0 * PI / c.project.wavelength_m;
    let mut v: Vec<C> =
        g.x.iter()
            .map(|x| {
                let r = (x - offset).powi(2);
                C::from_polar(
                    (-r / wm2).exp(),
                    -k * r / (2.0 * c.project.cavity_length_m) + k * slope * x,
                )
            })
            .collect();
    normalise(&mut v);
    v
}
#[derive(Clone)]
pub struct Field {
    pub x: Vec<C>,
    pub y: Option<Vec<C>>,
}
impl Field {
    pub fn normalise(&mut self) {
        normalise(&mut self.x);
        if let Some(y) = &mut self.y {
            normalise(y)
        }
    }
    pub fn dot(&self, b: &Self) -> C {
        dot(&self.x, &b.x)
            * match (&self.y, &b.y) {
                (Some(a), Some(b)) => dot(a, b),
                (None, None) => C::new(1.0, 0.0),
                _ => panic!("dimension mismatch"),
            }
    }
    pub fn norm(&self) -> f64 {
        norm(&self.x) * self.y.as_ref().map_or(1.0, |y| norm(y))
    }
    pub fn overlap(&self, b: &Self) -> f64 {
        overlap(&self.x, &b.x)
            * match (&self.y, &b.y) {
                (Some(a), Some(b)) => overlap(a, b),
                (None, None) => 1.0,
                _ => panic!("dimension mismatch"),
            }
    }
    pub fn value(&self, row: usize, col: usize) -> C {
        self.x[col] * self.y.as_ref().map_or(C::new(1.0, 0.0), |y| y[row])
    }
}
fn aligned_error(a: &[C], b: &[C]) -> f64 {
    let d = dot(a, b);
    let phase = if d.norm() > f64::EPSILON {
        d / d.norm()
    } else {
        C::new(1.0, 0.0)
    };
    a.iter()
        .zip(b)
        .map(|(x, y)| (*y - phase * *x).norm_sqr())
        .sum::<f64>()
        .sqrt()
        / norm(b)
}
pub fn field_error(a: &Field, b: &Field) -> f64 {
    let dx = aligned_error(&a.x, &b.x);
    match (&a.y, &b.y) {
        (Some(a), Some(b)) => {
            let dy = aligned_error(a, b);
            (dx * dx + dy * dy - 0.5 * dx * dx * dy * dy)
                .max(0.0)
                .sqrt()
        }
        _ => dx,
    }
}
pub fn intensity_error(a: &Field, b: &Field) -> f64 {
    let ax: Vec<f64> = a.x.iter().map(|z| z.norm_sqr()).collect();
    let cx: Vec<f64> = b.x.iter().map(|z| z.norm_sqr()).collect();
    let (Some(ay), Some(cy)) = (&a.y, &b.y) else {
        return ax
            .iter()
            .zip(&cx)
            .map(|(x, y)| (y - x).powi(2))
            .sum::<f64>()
            .sqrt()
            / cx.iter().map(|x| x * x).sum::<f64>().sqrt();
    };
    // Orthogonal decomposition of the rank-one intensity difference avoids an N^2 temporary.
    let aa: f64 = ax.iter().map(|x| x * x).sum();
    let t = ax.iter().zip(&cx).map(|(x, y)| x * y).sum::<f64>() / aa;
    let orth: f64 = ax.iter().zip(&cx).map(|(x, y)| (y - t * x).powi(2)).sum();
    let dd: f64 = cy.iter().map(|z| z.norm_sqr().powi(2)).sum();
    let along: f64 = ay
        .iter()
        .zip(cy)
        .map(|(x, y)| (t * y.norm_sqr() - x.norm_sqr()).powi(2))
        .sum();
    let cc: f64 = cx.iter().map(|x| x * x).sum();
    ((orth * dd + aa * along) / (cc * dd)).sqrt()
}
pub fn moments(v: &[C], g: &Grid) -> (f64, f64) {
    let total: f64 = v.iter().map(|z| z.norm_sqr()).sum();
    let mean =
        g.x.iter()
            .zip(v)
            .map(|(x, z)| x * z.norm_sqr())
            .sum::<f64>()
            / total;
    let var =
        g.x.iter()
            .zip(v)
            .map(|(x, z)| (x - mean).powi(2) * z.norm_sqr())
            .sum::<f64>()
            / total;
    (mean, 2.0 * var.max(0.0).sqrt())
}
#[derive(Clone, Debug, serde::Serialize)]
pub struct Metrics {
    pub x_centroid_m: f64,
    pub y_centroid_m: f64,
    pub mode_radius_x_m: f64,
    pub mode_radius_y_m: f64,
    pub hg00_overlap: f64,
    pub zero_angle_overlap: Option<f64>,
    pub y_asymmetry: f64,
}
pub fn metrics(f: &Field, g: &Grid, reference: &Field, zero: Option<&Field>) -> Metrics {
    let (xc, wx) = moments(&f.x, g);
    let (yc, wy) = f.y.as_ref().map_or((0.0, 0.0), |y| moments(y, g));
    let asym = f.y.as_ref().map_or(0.0, |y| {
        let total: f64 = y.iter().map(|z| z.norm_sqr()).sum();
        (0..g.n())
            .map(|i| (y[i].norm_sqr() - y[(g.n() - i) % g.n()].norm_sqr()).abs())
            .sum::<f64>()
            / total
    });
    Metrics {
        x_centroid_m: xc,
        y_centroid_m: yc,
        mode_radius_x_m: wx,
        mode_radius_y_m: wy,
        hg00_overlap: f.overlap(reference),
        zero_angle_overlap: zero.map(|z| f.overlap(z)),
        y_asymmetry: asym,
    }
}
