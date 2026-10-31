use serde::Deserialize;
use std::{
    error::Error,
    fs,
    path::{Path, PathBuf},
};

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Clone, Debug, Deserialize)]
pub struct Config {
    pub project: Project,
    pub grid: Grids,
    pub fox_li: FoxLi,
    pub tilt: Tilt,
    pub output: Output,
    pub validation: Validation,
}
#[derive(Clone, Debug, Deserialize)]
pub struct Project {
    pub name: String,
    pub wavelength_m: f64,
    pub cavity_length_m: f64,
    pub mirror1_radius_m: f64,
    pub mirror2_radius_m: f64,
}
#[derive(Clone, Debug, Deserialize)]
pub struct Grids {
    pub one_d: GridSpec,
    pub two_d: GridSpec,
    pub propagation_padding_factor: usize,
}
#[derive(Clone, Debug, Deserialize)]
pub struct GridSpec {
    pub n: usize,
    pub window_half_width_factor: f64,
    pub aperture_half_width_factor: f64,
}
#[derive(Clone, Debug, Deserialize)]
pub struct FoxLi {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub eigenvalue_relative_threshold: f64,
    pub consecutive_converged: usize,
    pub intensity_floor: f64,
    pub initializer: String,
    pub random_seed: u64,
}
#[derive(Clone, Debug, Deserialize)]
pub struct Tilt {
    pub mirror: String,
    pub angles_mrad: Vec<f64>,
    pub warm_start: bool,
    pub linear_fit_max_mrad: f64,
}
#[derive(Clone, Debug, Deserialize)]
pub struct Output {
    pub root: PathBuf,
}

#[derive(Clone, Debug, Deserialize)]
pub struct Validation {
    pub direct_grid_n: usize,
    pub grid_sizes: Vec<usize>,
    pub window_factors: Vec<f64>,
    pub tilt_angles_mrad: Vec<f64>,
    pub dense_check_iterations: usize,
    pub random_grid_n: usize,
    pub random_iterations: usize,
    pub perturbation_amplitude: f64,
    pub repeat_seed_offsets: Vec<u64>,
}

pub fn project_root() -> PathBuf {
    if let Ok(executable) = std::env::current_exe() {
        for candidate in executable.ancestors().skip(1) {
            if candidate.join("config/parameter.yaml").is_file() {
                return candidate.to_path_buf();
            }
        }
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}
pub fn parameter_path() -> PathBuf {
    project_root().join("config/parameter.yaml")
}
impl Config {
    pub fn load() -> Result<Self> {
        let path = parameter_path();
        let mut c: Self = serde_yaml::from_str(&fs::read_to_string(&path)?)?;
        c.tilt.mirror = c.tilt.mirror.to_uppercase();
        c.validate()?;
        Ok(c)
    }
    pub fn validate(&self) -> Result<()> {
        let p = &self.project;
        for value in [
            p.wavelength_m,
            p.cavity_length_m,
            p.mirror1_radius_m,
            p.mirror2_radius_m,
            self.fox_li.convergence_threshold,
            self.fox_li.eigenvalue_relative_threshold,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err("physical values and tolerances must be positive and finite".into());
            }
        }
        if (p.mirror1_radius_m - p.cavity_length_m).abs() > 1e-9 * p.cavity_length_m
            || (p.mirror2_radius_m - p.cavity_length_m).abs() > 1e-9 * p.cavity_length_m
        {
            return Err("confocal model requires R1=R2=L".into());
        }
        for g in [&self.grid.one_d, &self.grid.two_d] {
            if g.n < 8
                || !g.n.is_power_of_two()
                || !g.window_half_width_factor.is_finite()
                || !g.aperture_half_width_factor.is_finite()
                || g.aperture_half_width_factor <= 0.0
                || g.window_half_width_factor <= g.aperture_half_width_factor
            {
                return Err("grid must use power-of-two N>=8 and 0<aperture<window".into());
            }
            g.n.checked_mul(self.grid.propagation_padding_factor)
                .ok_or("grid size overflow")?;
        }
        if self.grid.propagation_padding_factor == 0
            || self.fox_li.max_iterations == 0
            || self.fox_li.consecutive_converged == 0
            || !self.fox_li.intensity_floor.is_finite()
            || !(0.0..1.0).contains(&self.fox_li.intensity_floor)
            || !["gaussian", "random"]
                .contains(&self.fox_li.initializer.to_ascii_lowercase().as_str())
        {
            return Err("invalid grid or Fox-Li settings".into());
        }
        if !["M1", "M2"].contains(&self.tilt.mirror.as_str())
            || self.tilt.angles_mrad.is_empty()
            || !self.tilt.linear_fit_max_mrad.is_finite()
            || self.tilt.linear_fit_max_mrad < 0.0
            || self.tilt.angles_mrad.iter().any(|a| !a.is_finite())
            || self
                .tilt
                .angles_mrad
                .iter()
                .enumerate()
                .any(|(i, a)| self.tilt.angles_mrad[..i].contains(a))
        {
            return Err("invalid tilt settings or duplicate/nonfinite angles".into());
        }
        if !self.grid.propagation_padding_factor.is_power_of_two() {
            return Err("padding factor must be a power of two".into());
        }
        if self.validation.direct_grid_n < 8
            || self.validation.random_grid_n < 8
            || !self.validation.direct_grid_n.is_power_of_two()
            || !self.validation.random_grid_n.is_power_of_two()
            || self.validation.dense_check_iterations == 0
            || self.validation.random_iterations == 0
            || self.validation.grid_sizes.is_empty()
            || self.validation.repeat_seed_offsets.is_empty()
            || !self.validation.perturbation_amplitude.is_finite()
            || self.validation.perturbation_amplitude < 0.0
            || self
                .validation
                .grid_sizes
                .iter()
                .any(|n| *n < 8 || !n.is_power_of_two())
            || self
                .validation
                .window_factors
                .iter()
                .any(|f| !f.is_finite() || *f <= self.grid.two_d.aperture_half_width_factor)
            || self
                .validation
                .tilt_angles_mrad
                .iter()
                .any(|a| !a.is_finite())
        {
            return Err("invalid explicit validation settings".into());
        }
        Ok(())
    }
}
