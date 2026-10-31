//! Fresnel transfer function in native (unshifted) FFT ordering.
//! The longitudinal carrier occurs once per physical propagation, not per axis.
use crate::{
    config::Config,
    field::{C, Field, Grid},
};
use rustfft::{Fft, FftPlanner};
use std::{f64::consts::PI, sync::Arc};

pub struct Propagator {
    n: usize,
    padded: usize,
    start: usize,
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
    transverse: Vec<C>,
    carrier: C,
}
impl Propagator {
    pub fn new(g: &Grid, c: &Config) -> Self {
        let padded = g.n() * g.padding;
        let mut planner = FftPlanner::new();
        let transverse = (0..padded)
            .map(|q| {
                let index = if q < padded.div_ceil(2) {
                    q as isize
                } else {
                    q as isize - padded as isize
                };
                let frequency = index as f64 / (padded as f64 * g.dx);
                C::from_polar(
                    1.0,
                    -PI * c.project.wavelength_m
                        * c.project.cavity_length_m
                        * frequency
                        * frequency,
                )
            })
            .collect();
        Self {
            n: g.n(),
            padded,
            start: (padded - g.n()) / 2,
            forward: planner.plan_fft_forward(padded),
            inverse: planner.plan_fft_inverse(padded),
            transverse,
            carrier: C::from_polar(
                1.0,
                2.0 * PI * c.project.cavity_length_m / c.project.wavelength_m,
            ),
        }
    }
    pub fn axis(&self, input: &[C], carrier: bool) -> Vec<C> {
        assert_eq!(input.len(), self.n);
        let mut a = vec![C::new(0.0, 0.0); self.padded];
        a[self.start..self.start + self.n].copy_from_slice(input);
        self.forward.process(&mut a);
        for (z, h) in a.iter_mut().zip(&self.transverse) {
            *z *= h;
        }
        self.inverse.process(&mut a);
        let scale = if carrier {
            self.carrier
        } else {
            C::new(1.0, 0.0)
        } / self.padded as f64;
        a[self.start..self.start + self.n]
            .iter()
            .map(|z| z * scale)
            .collect()
    }
    pub fn propagate(&self, f: &Field) -> Field {
        Field {
            x: self.axis(&f.x, true),
            y: f.y.as_ref().map(|y| self.axis(y, false)),
        }
    }
    /// Dense row/column FFT, accepting arbitrary (including nonseparable) 2-D fields.
    pub fn full_2d(&self, input: &[C]) -> Vec<C> {
        assert_eq!(input.len(), self.n * self.n);
        let p = self.padded;
        let mut a = vec![C::new(0.0, 0.0); p * p];
        for y in 0..self.n {
            a[(self.start + y) * p + self.start..(self.start + y) * p + self.start + self.n]
                .copy_from_slice(&input[y * self.n..(y + 1) * self.n]);
        }
        for row in a.chunks_exact_mut(p) {
            self.forward.process(row);
        }
        let mut col = vec![C::new(0.0, 0.0); p];
        for x in 0..p {
            for y in 0..p {
                col[y] = a[y * p + x];
            }
            self.forward.process(&mut col);
            for y in 0..p {
                a[y * p + x] = col[y];
            }
        }
        for y in 0..p {
            for x in 0..p {
                a[y * p + x] *= self.carrier * self.transverse[y] * self.transverse[x];
            }
        }
        for row in a.chunks_exact_mut(p) {
            self.inverse.process(row);
        }
        for x in 0..p {
            for y in 0..p {
                col[y] = a[y * p + x];
            }
            self.inverse.process(&mut col);
            for y in 0..p {
                a[y * p + x] = col[y];
            }
        }
        let mut out = Vec::with_capacity(self.n * self.n);
        for y in 0..self.n {
            out.extend(
                a[(self.start + y) * p + self.start..(self.start + y) * p + self.start + self.n]
                    .iter()
                    .map(|z| z / (p * p) as f64),
            );
        }
        out
    }
}

pub struct Cavity<'a> {
    pub prop: Propagator,
    pub grid: &'a Grid,
    mirror1_x: Vec<C>,
    mirror1_y: Vec<C>,
    mirror2_x: Vec<C>,
    mirror2_y: Vec<C>,
}
impl<'a> Cavity<'a> {
    pub fn new(g: &'a Grid, c: &Config, theta: f64) -> Self {
        let k = 2.0 * PI / c.project.wavelength_m;
        let phase = |radius: f64, tilted: bool| -> Vec<C> {
            g.x.iter()
                .map(|x| {
                    if x.abs() <= g.aperture {
                        C::from_polar(
                            1.0,
                            -k * x * x / radius
                                + if tilted {
                                    2.0 * k * x * theta.sin()
                                } else {
                                    0.0
                                },
                        )
                    } else {
                        C::new(0.0, 0.0)
                    }
                })
                .collect()
        };
        Self {
            prop: Propagator::new(g, c),
            grid: g,
            mirror1_x: phase(c.project.mirror1_radius_m, false),
            mirror1_y: phase(c.project.mirror1_radius_m, c.tilt.mirror == "M1"),
            mirror2_x: phase(c.project.mirror2_radius_m, false),
            mirror2_y: phase(c.project.mirror2_radius_m, c.tilt.mirror == "M2"),
        }
    }
    pub fn mirror(&self, f: &mut Field, second: bool) {
        let (mx, my) = if second {
            (&self.mirror2_x, &self.mirror2_y)
        } else {
            (&self.mirror1_x, &self.mirror1_y)
        };
        for (z, m) in f.x.iter_mut().zip(mx) {
            *z *= m;
        }
        if let Some(y) = &mut f.y {
            for (z, m) in y.iter_mut().zip(my) {
                *z *= m;
            }
        }
    }
    pub fn m2(&self, f: &Field) -> Field {
        let mut v = self.prop.propagate(f);
        self.mirror(&mut v, true);
        v
    }
    pub fn round_trip(&self, f: &Field) -> Field {
        let m2 = self.m2(f);
        let mut v = self.prop.propagate(&m2);
        self.mirror(&mut v, false);
        v
    }
    pub fn full_mirror(&self, v: &mut [C], second: bool) {
        let n = self.grid.n();
        let (mx, my) = if second {
            (&self.mirror2_x, &self.mirror2_y)
        } else {
            (&self.mirror1_x, &self.mirror1_y)
        };
        for y in 0..n {
            for x in 0..n {
                v[y * n + x] *= mx[x] * my[y];
            }
        }
    }
    pub fn full_round_trip(&self, v: &[C]) -> Vec<C> {
        let mut a = self.prop.full_2d(v);
        self.full_mirror(&mut a, true);
        let mut b = self.prop.full_2d(&a);
        self.full_mirror(&mut b, false);
        b
    }
}

pub fn direct_1d(input: &[C], g: &Grid, c: &Config) -> Vec<C> {
    let lambda_l = c.project.wavelength_m * c.project.cavity_length_m;
    let carrier = C::from_polar(
        1.0,
        2.0 * PI * c.project.cavity_length_m / c.project.wavelength_m,
    );
    let prefactor = carrier / C::new(0.0, lambda_l).sqrt() * g.dx;
    g.x.iter()
        .map(|x| {
            prefactor
                * g.x
                    .iter()
                    .zip(input)
                    .map(|(source, u)| u * C::from_polar(1.0, PI * (x - source).powi(2) / lambda_l))
                    .sum::<C>()
        })
        .collect()
}
pub fn direct_2d(input: &[C], g: &Grid, c: &Config) -> Vec<C> {
    let n = g.n();
    let lambda_l = c.project.wavelength_m * c.project.cavity_length_m;
    let phase: Vec<C> =
        g.x.iter()
            .flat_map(|x| {
                g.x.iter()
                    .map(move |s| C::from_polar(1.0, PI * (x - s).powi(2) / lambda_l))
            })
            .collect();
    let mut temp = vec![C::new(0.0, 0.0); n * n];
    let mut out = temp.clone();
    for y in 0..n {
        for x in 0..n {
            temp[y * n + x] = (0..n).map(|s| input[y * n + s] * phase[x * n + s]).sum();
        }
    }
    let prefactor = C::from_polar(
        1.0,
        2.0 * PI * c.project.cavity_length_m / c.project.wavelength_m,
    ) / C::new(0.0, lambda_l)
        * g.dx
        * g.dx;
    for y in 0..n {
        for x in 0..n {
            out[y * n + x] = prefactor
                * (0..n)
                    .map(|s| phase[y * n + s] * temp[s * n + x])
                    .sum::<C>();
        }
    }
    out
}
