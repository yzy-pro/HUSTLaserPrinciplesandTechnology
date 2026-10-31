use crate::{
    config::{self, Config, Result},
    field::{self, C, Field, Grid, Metrics},
    propagation::{self, Cavity, Propagator},
    solver::{self, Solution},
};
use chrono::Utc;
use serde::{Serialize, de::DeserializeOwned};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};

fn write_json(path: &Path, data: &impl Serialize) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    serde_json::to_writer_pretty(BufWriter::new(File::create(path)?), data)?;
    Ok(())
}
fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T> {
    Ok(serde_json::from_reader(File::open(path)?)?)
}
fn write_csv<T: Serialize>(path: &Path, rows: &[T]) -> Result<()> {
    let mut w = csv::Writer::from_path(path)?;
    for row in rows {
        w.serialize(row)?;
    }
    w.flush()?;
    Ok(())
}
fn summary(s: &Solution, m1: &Metrics, m2: &Metrics, g: &Grid, c: &Config) -> Value {
    json!({"implementation":"rust","representation":if s.field.y.is_some(){"exact_separable_2d"}else{"1d"},"n":g.n(),"dx_m":g.dx,"dy_m":if s.field.y.is_some(){Some(g.dx)}else{None},"window_half_width_m":g.half_width,"aperture_half_width_m":g.aperture,"seed":c.fox_li.random_seed,"initialization":c.fox_li.initializer,"gamma_real":s.gamma.re,"gamma_imag":s.gamma.im,"gamma_abs":s.gamma.norm(),"gamma_phase_rad":s.gamma.arg(),"round_trip_loss":(1.0-s.gamma.norm_sqr()).max(0.0),"iterations":s.iterations,"converged":s.converged,"stop_reason":if s.converged{"all three criteria met for consecutive iterations"}else{"maximum_iterations_reached"},"final_errors":s.history.last(),"m1":m1,"m2":m2})
}
fn save_fields(path: &Path, g: &Grid, m1: &Field, m2: &Field) -> Result<()> {
    let mut w = csv::Writer::from_path(path)?;
    w.write_record([
        "coordinate_m",
        "m1_x_real",
        "m1_x_imag",
        "m2_x_real",
        "m2_x_imag",
        "m1_y_real",
        "m1_y_imag",
        "m2_y_real",
        "m2_y_imag",
    ])?;
    for i in 0..g.n() {
        let a = m1.x[i];
        let b = m2.x[i];
        let mut row = vec![
            g.x[i].to_string(),
            a.re.to_string(),
            a.im.to_string(),
            b.re.to_string(),
            b.im.to_string(),
        ];
        if let (Some(y1), Some(y2)) = (&m1.y, &m2.y) {
            row.extend([
                y1[i].re.to_string(),
                y1[i].im.to_string(),
                y2[i].re.to_string(),
                y2[i].im.to_string(),
            ]);
        } else {
            row.extend(["".into(), "".into(), "".into(), "".into()]);
        }
        w.write_record(row)?;
    }
    w.flush()?;
    Ok(())
}
fn save_dense(path: &Path, field: &[C]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for z in field {
        w.write_all(&z.re.to_le_bytes())?;
        w.write_all(&z.im.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}
fn materialise(f: &Field) -> Vec<C> {
    let ny = if f.y.is_some() { f.x.len() } else { 1 };
    (0..ny)
        .flat_map(|y| (0..f.x.len()).map(move |x| f.value(y, x)))
        .collect()
}

pub fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!(
            "confocal-fox-li [all|baseline|tilt|validate] [--mode 1d|2d|both]\nAll numerical parameters: config/parameter.yaml\nAuthor: 陈恪瑾 / U202413925 / 光电信息科学与工程202409班"
        );
        return Ok(());
    }
    let command = args.first().map(String::as_str).unwrap_or("all");
    if !["all", "baseline", "tilt", "validate"].contains(&command) {
        return Err("unknown command; use --help".into());
    }
    let mut mode = "both";
    if args.len() > 1 {
        if command != "baseline"
            || args.len() != 3
            || args[1] != "--mode"
            || !["1d", "2d", "both"].contains(&args[2].as_str())
        {
            return Err(
                "only baseline accepts --mode 1d|2d|both; physical CLI overrides are disabled"
                    .into(),
            );
        }
        mode = &args[2];
    }
    let c = Config::load()?;
    let root = config::project_root().join(&c.output.root);
    let stamp = Utc::now().format("%Y%m%dT%H%M%SZ").to_string();
    let mut run = root.join("runs").join(&stamp);
    let mut suffix = 0;
    while run.exists() {
        suffix += 1;
        run = root.join("runs").join(format!("{stamp}_{suffix:02}"));
    }
    fs::create_dir_all(&run)?;
    fs::copy(config::parameter_path(), run.join("parameter_used.yaml"))?;
    let timer = Instant::now();
    let mut results = BTreeMap::new();
    if command == "all" || command == "validate" {
        println!("Running numerical verification ...");
        let v = validate(&c, &run)?;
        write_json(&run.join("validation.json"), &v)?;
        results.insert(
            "rust_validation",
            json!({"files":{"summary":"validation.json"},"summary":v}),
        );
    }
    if command == "all" || command == "baseline" {
        for dim in [1, 2] {
            if mode != "both" && mode != format!("{dim}d") {
                continue;
            }
            let entry = baseline(&c, &run, dim)?;
            results.insert(
                if dim == 1 {
                    "rust_baseline_1d"
                } else {
                    "rust_baseline_2d"
                },
                entry,
            );
        }
    }
    if command == "all" || command == "tilt" {
        let entry = tilt(&c, &run)?;
        results.insert("rust_tilt", entry);
    }
    let parameter_value: Value =
        serde_yaml::from_str(&fs::read_to_string(config::parameter_path())?)?;
    let metadata = json!({"generated_at_utc":Utc::now().to_rfc3339(),"implementation":"rust","project_name":c.project.name,"authors":{"name":"陈恪瑾","student_id":"U202413925","class":"光电信息科学与工程202409班"},"command":args,"config":parameter_value,"elapsed_seconds":timer.elapsed().as_secs_f64(),"results":results});
    write_json(&run.join("manifest.json"), &metadata)?;
    let index_path = root.join("manifest.json");
    let mut index: Value = if index_path.exists() {
        read_json(&index_path)?
    } else {
        json!({"results":{}})
    };
    index["generated_at_utc"] = metadata["generated_at_utc"].clone();
    index["latest_run"] = json!(run.strip_prefix(&root)?.to_string_lossy());
    for (key, entry) in metadata["results"].as_object().unwrap() {
        let mut entry = entry.clone();
        entry["parameters"] = metadata["config"].clone();
        for path in entry["files"].as_object_mut().unwrap().values_mut() {
            *path = json!(
                run.strip_prefix(&root)?
                    .join(path.as_str().unwrap())
                    .to_string_lossy()
            );
        }
        index["results"][key] = entry;
    }
    write_json(&index_path, &index)?;
    println!(
        "Results written to {} ({:.3} s)",
        run.display(),
        timer.elapsed().as_secs_f64()
    );
    Ok(())
}
fn baseline(c: &Config, run: &Path, dim: usize) -> Result<Value> {
    let g = Grid::new(
        c,
        if dim == 1 {
            &c.grid.one_d
        } else {
            &c.grid.two_d
        },
    );
    let cavity = Cavity::new(&g, c, 0.0);
    let mut initial = solver::initial(c, &g, dim == 2, 0.0);
    if c.fox_li.initializer.eq_ignore_ascii_case("random") {
        return dense_baseline(c, run, dim, &g, &cavity);
    }
    if dim == 1 {
        initial.y = None;
    }
    let reference = initial.clone();
    let s = solver::solve(initial, &cavity, c);
    let m2 = cavity.m2(&s.field);
    let metric1 = field::metrics(&s.field, &g, &reference, None);
    let metric2 = field::metrics(&m2, &g, &reference, None);
    let sum = summary(&s, &metric1, &metric2, &g, c);
    let folder = format!("baseline/{dim}d");
    let out = run.join(&folder);
    fs::create_dir_all(&out)?;
    save_fields(&out.join("factors.csv"), &g, &s.field, &m2)?;
    write_csv(&out.join("history.csv"), &s.history)?;
    write_json(&out.join("summary.json"), &sum)?;
    println!(
        "{dim}D baseline: iterations={}, converged={}, |gamma|={:.15}",
        s.iterations,
        s.converged,
        s.gamma.norm()
    );
    Ok(
        json!({"summary":sum,"files":{"summary":format!("{folder}/summary.json"),"history":format!("{folder}/history.csv"),"factors":format!("{folder}/factors.csv")}}),
    )
}
#[derive(Serialize)]
#[allow(non_snake_case)]
struct Scan {
    theta_mrad: f64,
    theta_rad: f64,
    initialization: String,
    iterations: usize,
    converged: bool,
    gamma_abs_M1: f64,
    gamma_phase_M1_rad: f64,
    loss_M1: f64,
    loss_M2: f64,
    y_centroid_M1_mm: f64,
    y_centroid_M2_mm: f64,
    x_centroid_M1_mm: f64,
    x_centroid_M2_mm: f64,
    mode_radius_y_M1_mm: f64,
    mode_radius_y_M2_mm: f64,
    hg00_overlap_M1: f64,
    hg00_overlap_M2: f64,
    zero_angle_overlap_M1: Option<f64>,
    zero_angle_overlap_M2: Option<f64>,
    asymmetry_M1: f64,
    asymmetry_M2: f64,
    final_field_error: f64,
    final_intensity_error: f64,
    final_eigenvalue_error: Option<f64>,
}
fn tilt(c: &Config, run: &Path) -> Result<Value> {
    if !c.fox_li.initializer.eq_ignore_ascii_case("gaussian") {
        return Err("tilt scan uses the matched/geometric Gaussian seed; random initial fields are supported by baseline and the explicit validation diagnostics".into());
    }
    let g = Grid::new(c, &c.grid.two_d);
    let k = 2.0 * std::f64::consts::PI / c.project.wavelength_m;
    let phase_step = c
        .tilt
        .angles_mrad
        .iter()
        .map(|a| 2.0 * k * g.dx * (a * 1e-3).sin().abs())
        .fold(0.0, f64::max);
    if phase_step >= std::f64::consts::PI {
        return Err("tilt phase is undersampled: 2*k*dy*abs(sin(theta_max)) >= pi".into());
    }
    let out = run.join("tilt");
    fs::create_dir_all(&out)?;
    let reference = solver::initial(c, &g, true, 0.0);
    let mut fields = Vec::new();
    let mut previous = None;
    let mut rows = Vec::new();
    for (i, &angle) in c.tilt.angles_mrad.iter().enumerate() {
        let theta = angle * 1e-3;
        let cavity = Cavity::new(&g, c, theta);
        let need_geometry = previous.is_none()
            || (i > 0 && c.tilt.angles_mrad[i - 1].abs() < 1e-15 && angle.abs() >= 1e-15);
        let warm = c.tilt.warm_start && !need_geometry;
        let initial = if warm {
            previous.clone().unwrap()
        } else {
            solver::initial(c, &g, true, theta)
        };
        let s = solver::solve(initial, &cavity, c);
        let m2 = cavity.m2(&s.field);
        let a = field::metrics(&s.field, &g, &reference, None);
        let b = field::metrics(&m2, &g, &reference, None);
        let last = s.history.last().unwrap();
        rows.push(Scan {
            theta_mrad: angle,
            theta_rad: theta,
            initialization: if warm {
                "warm_start_previous_angle"
            } else if angle.abs() > 0.0 {
                "paraxial_geometric_fixed_point"
            } else {
                "gaussian"
            }
            .into(),
            iterations: s.iterations,
            converged: s.converged,
            gamma_abs_M1: s.gamma.norm(),
            gamma_phase_M1_rad: s.gamma.arg(),
            loss_M1: (1.0 - s.gamma.norm_sqr()).max(0.0),
            loss_M2: (1.0 - s.gamma.norm_sqr()).max(0.0),
            y_centroid_M1_mm: a.y_centroid_m * 1e3,
            y_centroid_M2_mm: b.y_centroid_m * 1e3,
            x_centroid_M1_mm: a.x_centroid_m * 1e3,
            x_centroid_M2_mm: b.x_centroid_m * 1e3,
            mode_radius_y_M1_mm: a.mode_radius_y_m * 1e3,
            mode_radius_y_M2_mm: b.mode_radius_y_m * 1e3,
            hg00_overlap_M1: a.hg00_overlap,
            hg00_overlap_M2: b.hg00_overlap,
            zero_angle_overlap_M1: None,
            zero_angle_overlap_M2: None,
            asymmetry_M1: a.y_asymmetry,
            asymmetry_M2: b.y_asymmetry,
            final_field_error: last.delta,
            final_intensity_error: last.intensity_difference,
            final_eigenvalue_error: last.eigenvalue_difference,
        });
        let name = format!("theta_{angle}mrad");
        save_fields(&out.join(format!("factors_{name}.csv")), &g, &s.field, &m2)?;
        write_csv(&out.join(format!("history_{name}.csv")), &s.history)?;
        let mut sum = summary(&s, &a, &b, &g, c);
        sum["initialization"] = json!(rows.last().unwrap().initialization);
        sum["theta_mrad"] = json!(angle);
        write_json(&out.join(format!("summary_{name}.json")), &sum)?;
        println!(
            "theta={angle:>4} mrad: iterations={:>4}, converged={}, delta={:.3e}",
            s.iterations, s.converged, last.delta
        );
        previous = Some(s.field.clone());
        fields.push((s.field, m2));
    }
    if let Some(i) = c.tilt.angles_mrad.iter().position(|a| a.abs() < 1e-15) {
        for (row, (a, b)) in rows.iter_mut().zip(&fields) {
            row.zero_angle_overlap_M1 = Some(a.overlap(&fields[i].0));
            row.zero_angle_overlap_M2 = Some(b.overlap(&fields[i].1));
        }
    }
    write_csv(&out.join("scan.csv"), &rows)?;
    let mut fits = BTreeMap::new();
    for plane in ["M1", "M2"] {
        let small: Vec<_> = rows
            .iter()
            .filter(|r| r.converged && r.theta_mrad.abs() <= c.tilt.linear_fit_max_mrad)
            .map(|r| {
                (
                    r.theta_mrad,
                    if plane == "M1" {
                        r.y_centroid_M1_mm
                    } else {
                        r.y_centroid_M2_mm
                    },
                )
            })
            .collect();
        if small.len() < 2 {
            continue;
        }
        let n = small.len() as f64;
        let meanx = small.iter().map(|p| p.0).sum::<f64>() / n;
        let meany = small.iter().map(|p| p.1).sum::<f64>() / n;
        let variance = small.iter().map(|p| (p.0 - meanx).powi(2)).sum::<f64>();
        if variance == 0.0 {
            continue;
        }
        let slope = small
            .iter()
            .map(|p| (p.0 - meanx) * (p.1 - meany))
            .sum::<f64>()
            / variance;
        let intercept = meany - slope * meanx;
        let ss = small.iter().map(|p| (p.1 - meany).powi(2)).sum::<f64>();
        let residual = small
            .iter()
            .map(|p| (p.1 - slope * p.0 - intercept).powi(2))
            .sum::<f64>();
        fits.insert(plane,json!({"slope_mm_per_mrad":slope,"intercept_mm":intercept,"r2":if ss>f64::EPSILON{Some(1.0-residual/ss)}else{None},"max_abs_theta_mrad":c.tilt.linear_fit_max_mrad,"n":small.len()}));
    }
    let sum = json!({"dimension":"2d","mirror":c.tilt.mirror,"phase_step_max_rad":phase_step,"linear_fit":fits,"rows":rows});
    write_json(&out.join("summary.json"), &sum)?;
    Ok(json!({"summary":sum,"files":{"scan":"tilt/scan.csv","summary":"tilt/summary.json"}}))
}

fn relative_error(a: &[C], b: &[C]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (*x - *y).norm_sqr())
        .sum::<f64>()
        .sqrt()
        / field::norm(b)
}
fn check(name: &str, value: f64, threshold: f64, checks: &mut Vec<Value>) -> Result<()> {
    let passed = value.is_finite() && value < threshold;
    checks.push(json!({"name":name,"value":value,"threshold":threshold,"passed":passed}));
    if !passed {
        return Err(format!("validation failed: {name} = {value:e} >= {threshold:e}").into());
    }
    Ok(())
}
fn validate(c: &Config, run: &Path) -> Result<Value> {
    let mut checks = Vec::new();
    let mut small = c.grid.two_d.clone();
    small.n = c.validation.direct_grid_n;
    let g = Grid::new(c, &small);
    let initial = solver::initial(c, &g, true, 0.0);
    let prop = Propagator::new(&g, c);
    let direct = propagation::direct_1d(&initial.x, &g, c);
    let fft = prop.axis(&initial.x, true);
    check(
        "direct_vs_fft_1d_relative",
        relative_error(&fft, &direct),
        1e-8,
        &mut checks,
    )?;
    let full = materialise(&initial);
    let direct = propagation::direct_2d(&full, &g, c);
    let dense = prop.full_2d(&full);
    check(
        "direct_vs_fft_2d_relative",
        relative_error(&dense, &direct),
        1e-8,
        &mut checks,
    )?;
    let factored = materialise(&prop.propagate(&initial));
    check(
        "dense_vs_separable_single_pass",
        relative_error(&factored, &dense),
        1e-12,
        &mut checks,
    )?;
    let cavity = Cavity::new(&g, c, 1e-3);
    let dense = cavity.full_round_trip(&full);
    let sep = materialise(&cavity.round_trip(&initial));
    check(
        "dense_vs_separable_round_trip",
        relative_error(&sep, &dense),
        1e-12,
        &mut checks,
    )?;
    let mut periodic = g.clone();
    periodic.padding = 1;
    let propagated = Propagator::new(&periodic, c).axis(&initial.x, true);
    check(
        "unapertured_periodic_power_error",
        (field::norm(&propagated).powi(2) - field::norm(&initial.x).powi(2)).abs(),
        1e-12,
        &mut checks,
    )?;
    let large = Grid::new(c, &c.grid.two_d);
    let cavity = Cavity::new(&large, c, 1e-3);
    let f = solver::initial(c, &large, true, 1e-3);
    let dense = cavity.full_round_trip(&materialise(&f));
    let sep = materialise(&cavity.round_trip(&f));
    check(
        "production_grid_dense_vs_separable",
        relative_error(&sep, &dense),
        1e-12,
        &mut checks,
    )?;
    let n = c.grid.two_d.n;
    let mut full = materialise(&f);
    let reference = &full;
    let prev = solver::initial(c, &large, true, 0.0);
    let dsep = field::field_error(&prev, &f);
    let d = field::dot(&materialise(&prev), reference);
    let phase = d / d.norm();
    let ddense = materialise(&prev)
        .iter()
        .zip(reference)
        .map(|(a, b)| (*b - phase * *a).norm_sqr())
        .sum::<f64>()
        .sqrt()
        / field::norm(reference);
    check(
        "rank_one_field_error_identity",
        (dsep - ddense).abs(),
        1e-12,
        &mut checks,
    )?;
    let idense = materialise(&prev)
        .iter()
        .zip(reference)
        .map(|(a, b)| (b.norm_sqr() - a.norm_sqr()).powi(2))
        .sum::<f64>()
        .sqrt()
        / reference
            .iter()
            .map(|z| z.norm_sqr().powi(2))
            .sum::<f64>()
            .sqrt();
    check(
        "rank_one_intensity_error_identity",
        (field::intensity_error(&prev, &f) - idense).abs(),
        1e-12,
        &mut checks,
    )?;
    // Verify the entire baseline iteration using dense 2-D arrays at the configured production grid.
    full = materialise(&solver::initial(c, &large, true, 0.0));
    let base_cavity = Cavity::new(&large, c, 0.0);
    let mut histories = Vec::new();
    for _ in 0..c.validation.dense_check_iterations {
        let mut next = base_cavity.full_round_trip(&full);
        let gamma = field::dot(&full, &next) / field::dot(&full, &full);
        field::normalise(&mut next);
        full = next;
        histories.push(gamma);
    }
    let mut step_config = c.clone();
    step_config.fox_li.max_iterations = c.validation.dense_check_iterations;
    step_config.fox_li.consecutive_converged = c.validation.dense_check_iterations + 1;
    let sep_result = solver::solve(
        solver::initial(c, &large, true, 0.0),
        &base_cavity,
        &step_config,
    );
    check(
        "dense_baseline_steps_overlap_error",
        1.0 - field::overlap(&full, &materialise(&sep_result.field)),
        1e-12,
        &mut checks,
    )?;
    check(
        "dense_baseline_gamma_error",
        (histories.last().unwrap() - sep_result.gamma).norm(),
        1e-12,
        &mut checks,
    )?;
    let mut convergence = Vec::new();
    for &grid_n in &c.validation.grid_sizes {
        let mut spec = c.grid.two_d.clone();
        spec.n = grid_n;
        let grid = Grid::new(c, &spec);
        let cavity = Cavity::new(&grid, c, 0.0);
        let reference = solver::initial(c, &grid, true, 0.0);
        let s = solver::solve(reference.clone(), &cavity, c);
        let m = field::metrics(&s.field, &grid, &reference, None);
        convergence.push(json!({"kind":"grid","n":spec.n,"window_factor":spec.window_half_width_factor,"dx_m":grid.dx,"width_m":m.mode_radius_y_m,"hg00_overlap":m.hg00_overlap,"iterations":s.iterations,"converged":s.converged}));
    }
    for &factor in &c.validation.window_factors {
        let mut spec = c.grid.two_d.clone();
        spec.window_half_width_factor = factor;
        let grid = Grid::new(c, &spec);
        let cavity = Cavity::new(&grid, c, 0.0);
        let reference = solver::initial(c, &grid, true, 0.0);
        let s = solver::solve(reference.clone(), &cavity, c);
        let m = field::metrics(&s.field, &grid, &reference, None);
        convergence.push(json!({"kind":"window","n":spec.n,"window_factor":factor,"dx_m":grid.dx,"width_m":m.mode_radius_y_m,"hg00_overlap":m.hg00_overlap,"iterations":s.iterations,"converged":s.converged}));
    }
    let theoretical =
        (c.project.wavelength_m * c.project.cavity_length_m / std::f64::consts::PI).sqrt();
    for row in &convergence {
        check(
            &format!(
                "baseline_radius_n{}_window{}",
                row["n"], row["window_factor"]
            ),
            (row["width_m"].as_f64().unwrap() / theoretical - 1.0).abs(),
            1e-6,
            &mut checks,
        )?;
    }
    let k = 2.0 * std::f64::consts::PI / c.project.wavelength_m;
    let phase_step = c
        .tilt
        .angles_mrad
        .iter()
        .map(|a| 2.0 * k * large.dx * (a * 1e-3).sin().abs())
        .fold(0.0, f64::max);
    check(
        "tilt_phase_sampling_rad",
        phase_step,
        std::f64::consts::PI,
        &mut checks,
    )?;
    write_json(
        &run.join("dense_baseline_gamma.json"),
        &histories
            .iter()
            .enumerate()
            .map(|(i, g)| json!({"iteration":i+1,"gamma_real":g.re,"gamma_imag":g.im}))
            .collect::<Vec<_>>(),
    )?;
    // Cross-resolution 1-D and 2-D profiles live on common coordinates at the default sampling ratio.
    let g1 = Grid::new(c, &c.grid.one_d);
    let ref1 = solver::initial(c, &g1, false, 0.0);
    let s1 = solver::solve(ref1, &Cavity::new(&g1, c, 0.0), c);
    if g1.n().is_multiple_of(large.n()) && (g1.half_width - large.half_width).abs() < 1e-15 {
        let stride = g1.n() / large.n();
        let a: Vec<C> = (0..large.n()).map(|i| s1.field.x[i * stride]).collect();
        check(
            "baseline_1d_2d_profile_overlap_error",
            1.0 - field::overlap(&a, &sep_result.field.x),
            1e-12,
            &mut checks,
        )?;
    }
    let mut tilt_diagnostics = Vec::new();
    for &angle in &c.validation.tilt_angles_mrad {
        for &grid_n in &c.validation.grid_sizes {
            let mut spec = c.grid.two_d.clone();
            spec.n = grid_n;
            let grid = Grid::new(c, &spec);
            let cavity = Cavity::new(&grid, c, angle * 1e-3);
            let init = solver::initial(c, &grid, true, angle * 1e-3);
            let result = solver::solve(init, &cavity, c);
            let m2 = cavity.m2(&result.field);
            let ref0 = solver::initial(c, &grid, true, 0.0);
            let a = field::metrics(&result.field, &grid, &ref0, None);
            let b = field::metrics(&m2, &grid, &ref0, None);
            tilt_diagnostics.push(json!({"theta_mrad":angle,"n":grid_n,"iterations":result.iterations,"converged":result.converged,"m1":a,"m2":b,"final_errors":result.history.last()}));
        }
    }
    let mut seed_diagnostics = Vec::new();
    let mut spec = c.grid.two_d.clone();
    spec.n = c.validation.random_grid_n;
    let grid = Grid::new(c, &spec);
    let cavity = Cavity::new(&grid, c, 0.0);
    let mut repeat_config = c.clone();
    repeat_config.fox_li.max_iterations = c.validation.random_iterations;
    for offset in &c.validation.repeat_seed_offsets {
        let seed = c.fox_li.random_seed + offset;
        let mut init = solver::initial(c, &grid, true, 0.0);
        let mut rng = solver::random_field(2 * grid.n(), seed);
        let mut y = rng.split_off(grid.n());
        for (a, r) in init.x.iter_mut().zip(rng) {
            *a += r * c.validation.perturbation_amplitude;
        }
        for (a, r) in init.y.as_mut().unwrap().iter_mut().zip(y.drain(..)) {
            *a += r * c.validation.perturbation_amplitude;
        }
        let result = solver::solve(init, &cavity, &repeat_config);
        let name = format!("validation_perturbed_seed_{seed}.csv");
        write_csv(&run.join(&name), &result.history)?;
        seed_diagnostics.push(json!({"seed":seed,"kind":"HG00_plus_complex_noise","noise_amplitude":c.validation.perturbation_amplitude,"n":grid.n(),"iterations":result.iterations,"converged":result.converged,"history":name,"final_errors":result.history.last(),"hg00_overlap":result.field.overlap(&solver::initial(c,&grid,true,0.0))}));
    }
    for &offset in &c.validation.repeat_seed_offsets {
        let seed = c.fox_li.random_seed + offset;
        let result = solver::solve_dense(
            solver::random_field(grid.n() * grid.n(), seed),
            &cavity,
            &repeat_config,
            true,
        );
        let name = format!("validation_dense_random_seed_{seed}.csv");
        write_csv(&run.join(&name), &result.history)?;
        seed_diagnostics.push(json!({"seed":seed,"kind":"dense_random_complex_gaussian","n":grid.n(),"iterations":result.iterations,"converged":result.converged,"history":name,"final_errors":result.history.last(),"hg00_overlap":field::overlap(&result.field,&materialise(&solver::initial(c,&grid,true,0.0)))}));
    }
    Ok(
        json!({"passed":true,"direct_grid_n":g.n(),"production_grid_n":n,"checks":checks,"grid_window_convergence":convergence,"tilt_grid_convergence":tilt_diagnostics,"initial_field_diagnostics":seed_diagnostics,"notes":"Diagnostic changes are explicitly configured in validation.*; un-converged initial-field and tilt diagnostics are findings, not successful eigenmode solves."}),
    )
}

fn dense_baseline(c: &Config, run: &Path, dim: usize, g: &Grid, cavity: &Cavity) -> Result<Value> {
    let count = if dim == 2 { g.n() * g.n() } else { g.n() };
    let initial = solver::random_field(count, c.fox_li.random_seed);
    let result = solver::solve_dense(initial, cavity, c, dim == 2);
    let m2 = if dim == 2 {
        let mut v = cavity.prop.full_2d(&result.field);
        cavity.full_mirror(&mut v, true);
        v
    } else {
        cavity
            .m2(&Field {
                x: result.field.clone(),
                y: None,
            })
            .x
    };
    let reference = materialise(&solver::initial(c, g, dim == 2, 0.0));
    let metrics = |f: &[C]| {
        let n = g.n();
        let total = f.iter().map(|z| z.norm_sqr()).sum::<f64>();
        let x = f
            .iter()
            .enumerate()
            .map(|(i, z)| g.x[i % n] * z.norm_sqr())
            .sum::<f64>()
            / total;
        let y = if dim == 2 {
            f.iter()
                .enumerate()
                .map(|(i, z)| g.x[i / n] * z.norm_sqr())
                .sum::<f64>()
                / total
        } else {
            0.0
        };
        let wx = 2.0
            * (f.iter()
                .enumerate()
                .map(|(i, z)| (g.x[i % n] - x).powi(2) * z.norm_sqr())
                .sum::<f64>()
                / total)
                .sqrt();
        let wy = if dim == 2 {
            2.0 * (f
                .iter()
                .enumerate()
                .map(|(i, z)| (g.x[i / n] - y).powi(2) * z.norm_sqr())
                .sum::<f64>()
                / total)
                .sqrt()
        } else {
            0.0
        };
        json!({"x_centroid_m":x,"y_centroid_m":y,"mode_radius_x_m":wx,"mode_radius_y_m":wy,"hg00_overlap":field::overlap(f,&reference)})
    };
    let sum = json!({"implementation":"rust","representation":"dense_complex128","dimension":format!("{dim}d"),"n":g.n(),"dx_m":g.dx,"window_half_width_m":g.half_width,"aperture_half_width_m":g.aperture,"seed":c.fox_li.random_seed,"rng":"rand::StdRng + Box-Muller; sequence differs from NumPy","initialization":"random","gamma_real":result.gamma.re,"gamma_imag":result.gamma.im,"gamma_abs":result.gamma.norm(),"round_trip_loss":(1.0-result.gamma.norm_sqr()).max(0.0),"iterations":result.iterations,"converged":result.converged,"m1":metrics(&result.field),"m2":metrics(&m2),"final_errors":result.history.last()});
    let folder = format!("baseline/{dim}d");
    let out = run.join(&folder);
    fs::create_dir_all(&out)?;
    save_dense(&out.join("field_m1.c128"), &result.field)?;
    save_dense(&out.join("field_m2.c128"), &m2)?;
    write_json(&out.join("summary.json"), &sum)?;
    write_csv(&out.join("history.csv"), &result.history)?;
    Ok(
        json!({"summary":sum,"files":{"summary":format!("{folder}/summary.json"),"history":format!("{folder}/history.csv"),"field_m1":format!("{folder}/field_m1.c128"),"field_m2":format!("{folder}/field_m2.c128")}}),
    )
}
