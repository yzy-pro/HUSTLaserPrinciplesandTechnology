"""Render saved Rust fields, generate LaTeX tables; no numerical solving."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/laser-matplotlib")
from .paths import project_root
ROOT = project_root()
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .plots import plot_1d_fields, plot_2d_fields, plot_convergence, plot_tilt_scan


def read_fields(folder, stem=""):
    suffix = f"_{stem}" if stem else ""
    summary = json.loads((folder / f"summary{suffix}.json").read_text())
    path = folder / f"factors{suffix}.csv"
    if path.exists():
        data = np.genfromtxt(path, delimiter=",", names=True)
        x = data["coordinate_m"]
        fields = []
        for plane in ("m1", "m2"):
            a = data[f"{plane}_x_real"] + 1j * data[f"{plane}_x_imag"]
            y = data[f"{plane}_y_real"] + 1j * data[f"{plane}_y_imag"]
            fields.append(a if np.isnan(y).all() else y[:, None] * a[None, :])
        return x, *fields, summary
    n = summary["n"]
    x = (np.arange(n) - n / 2) * summary["dx_m"]
    shape = (n, n) if summary["dimension"] == "2d" else (n,)
    fields = [np.fromfile(folder / f"field_{p}.c128", dtype="<c16").reshape(shape) for p in ("m1", "m2")]
    return x, *fields, summary


def history(path):
    return [{k: float(v) if v else np.inf for k, v in r.items()} for r in csv.DictReader(path.open())]


def tex_number(value):
    if value is None or not np.isfinite(value):
        return r"\textemdash"
    if value == 0:
        return r"\(0\)"
    if 0.001 <= abs(value) < 1000:
        return rf"\({value:.3g}\)"
    m, e = f"{value:.2e}".split("e")
    return rf"\({m}\times10^{{{int(e)}}}\)"


def table(path, columns, header, rows):
    path.write_text("\n".join([r"\begin{tabular}{" + columns + r"}\toprule", header + r"\\\midrule",
                              *[" & ".join(row) + r"\\" for row in rows], r"\bottomrule\end{tabular}"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    run = parser.parse_args().run_directory.resolve()
    manifest_path = run / "plots_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
    cfg = yaml.safe_load((run / "parameter_used.yaml").read_text())
    image_dir = run / "figures"
    image_dir.mkdir(exist_ok=True)
    floor = cfg["fox_li"]["intensity_floor"]
    project = cfg["project"]
    common = (f"lambda={project['wavelength_m']*1e9:g} nm, "
              f"L=R1=R2={project['cavity_length_m']:g} m, "
              f"a={cfg['grid']['two_d']['aperture_half_width_factor']:g} wm, "
              f"Xmax={cfg['grid']['two_d']['window_half_width_factor']:g} wm")
    records = []

    def record(name, sources):
        path = image_dir / name
        records.append(dict(file=path.relative_to(ROOT).as_posix(),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            source_data=[p.relative_to(ROOT).as_posix() for p in sources],
                            parameter_snapshot=(run / "parameter_used.yaml").relative_to(ROOT).as_posix()))

    baselines = {}
    for dim in (1, 2):
        folder = run / f"baseline/{dim}d"
        if not folder.exists():
            continue
        x, m1, m2, summary = read_fields(folder)
        baselines[dim] = (x, m1, m2, summary)
        sources = [folder / "summary.json", folder / "history.csv"]
        sources += [folder / "factors.csv"] if (folder / "factors.csv").exists() else [folder / "field_m1.c128", folder / "field_m2.c128"]
        title = (f"lambda={project['wavelength_m']*1e9:g} nm, "
                 f"L=R1=R2={project['cavity_length_m']:g} m, "
                 f"a={summary['aperture_half_width_m']*1e3:.4g} mm, "
                 f"Xmax={summary['window_half_width_m']*1e3:.4g} mm, "
                 f"N={len(x)}, theta=0 mrad | after reflection")
        name = f"baseline_{dim}d_fields.png"
        if dim == 1:
            plot_1d_fields(image_dir / name, x, m1, m2, title, floor)
        else:
            plot_2d_fields(image_dir / name, x, x, m1, m2, title, floor)
        record(name, sources)
        name = f"baseline_{dim}d_convergence.png"
        plot_convergence(image_dir / name, history(folder / "history.csv"), title)
        record(name, sources)
        np.savez_compressed(folder / "fields.npz", x=x, **({"y": x} if dim == 2 else {}), field_m1=m1, field_m2=m2)

    if len(baselines) == 2:
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        wm = np.sqrt(project["wavelength_m"] * project["cavity_length_m"] / np.pi)
        x, f, _, _ = baselines[1]
        ax.plot(x * 1e3, np.exp(-x*x/wm**2), "k-", label="HG00 theory")
        ax.plot(x * 1e3, abs(f)/abs(f).max(), "--", label="1-D numerical")
        xx, f2, _, _ = baselines[2]
        line = abs(f2[len(xx)//2])
        ax.plot(xx[::4]*1e3, line[::4]/line.max(), ".", label="2-D central slice")
        ax.set(xlabel="x (mm)", ylabel="amplitude / peak", xlim=(-.6,.6), title="M1 after reflection, theta=0")
        ax.grid(alpha=.25)
        ax.legend()
        fig.savefig(image_dir / "baseline_comparison.png", dpi=180)
        plt.close(fig)
        record("baseline_comparison.png", [run / "baseline/1d/fields.npz",run / "baseline/2d/fields.npz"])

    tilt_folder = run / "tilt"
    tilt = json.loads((tilt_folder / "summary.json").read_text()) if tilt_folder.exists() else None
    if tilt:
        rows = tilt["rows"]
        plot_tilt_scan(image_dir / "tilt_scan.png", rows, common + " | after reflection", tilt["linear_fit"])
        record("tilt_scan.png", [tilt_folder / "scan.csv",tilt_folder / "summary.json"])
        chosen = [rows[i]["theta_mrad"] for i in sorted({0,len(rows)//2,len(rows)-1})]
        for angle in chosen:
            stem = f"theta_{angle:g}mrad"
            x, a, b, summary = read_fields(tilt_folder, stem)
            title = common + f", N={len(x)}, theta={angle:g} mrad, converged={summary['converged']}"
            plot_2d_fields(image_dir / f"tilt_{angle:g}mrad_fields.png", x, x, a, b, title, floor)
            plot_convergence(image_dir / f"tilt_{angle:g}mrad_convergence.png", history(tilt_folder / f"history_{stem}.csv"), title)
            sources = [tilt_folder / f"factors_{stem}.csv",tilt_folder / f"history_{stem}.csv",tilt_folder / f"summary_{stem}.json"]
            for suffix in ("fields","convergence"):
                record(f"tilt_{angle:g}mrad_{suffix}.png", sources)
            np.savez_compressed(tilt_folder / f"fields_{stem}.npz", x=x,y=x,field_m1=a,field_m2=b)
        selected = [r for r in rows if r["converged"] and r["theta_mrad"] in (0,.5,1)]
        fig, axes = plt.subplots(1,2,figsize=(10,3.6),constrained_layout=True)
        for r in selected:
            x, a, b, _ = read_fields(tilt_folder,f"theta_{r['theta_mrad']:g}mrad")
            for ax, field in zip(axes,(a,b)):
                line = np.abs(field[:,len(x)//2])**2
                ax.plot(x*1e3,line/line.max(),label=f"{r['theta_mrad']:g} mrad")
        for ax, plane in zip(axes,("M1","M2")):
            ax.set(xlabel="y (mm)",ylabel="intensity / peak",title=plane+" after reflection",xlim=(-.5,.6))
            ax.grid(alpha=.25)
            ax.legend()
        fig.savefig(image_dir / "tilt_profiles.png",dpi=180)
        plt.close(fig)
        record("tilt_profiles.png",[tilt_folder / f"factors_theta_{r['theta_mrad']:g}mrad.csv" for r in selected])

    vpath = run / "validation.json"
    validation = json.loads(vpath.read_text()) if vpath.exists() else None
    if validation and validation["initial_field_diagnostics"]:
        fig, axes = plt.subplots(1,2,figsize=(10,3.6),constrained_layout=True)
        for item in validation["initial_field_diagnostics"]:
            ax = axes[0] if item["kind"].startswith("HG00") else axes[1]
            h = history(run / item["history"])
            ax.semilogy([r["iteration"] for r in h],[r["delta"] for r in h],label=f"seed={item['seed']}")
        for ax, title in zip(axes,("HG00 + complex perturbation","Dense random complex field")):
            ax.axhline(cfg["fox_li"]["convergence_threshold"],c="k",ls="--",label="threshold")
            ax.set(xlabel="iteration",ylabel="phase-aligned field error",title=title)
            ax.legend(fontsize=8)
            ax.grid(alpha=.25)
        fig.savefig(image_dir / "initial_field_diagnostics.png",dpi=180)
        plt.close(fig)
        record("initial_field_diagnostics.png",[vpath]+[run / r["history"] for r in validation["initial_field_diagnostics"]])

    gen = run / "report_generated"
    gen.mkdir(exist_ok=True)
    escaped_id = run.name.replace("_",r"\_")
    macros = ["% Generated from saved numerical data.",
              rf"\newcommand{{\RunID}}{{\texttt{{{escaped_id}}}}}",
              rf"\newcommand{{\RunPath}}{{\path{{results/runs/{run.name}}}}}"]
    for dim,(_,_,_,s) in baselines.items():
        label = "One" if dim == 1 else "Two"
        macros += [rf"\newcommand{{\{label}Iterations}}{{{s['iterations']}}}",
                   rf"\newcommand{{\{label}Gamma}}{{{s['gamma_real']:.9f}{s['gamma_imag']:+.9f}i}}",
                   rf"\newcommand{{\{label}Radius}}{{{s['m1']['mode_radius_x_m']*1e3:.9f}}}"]
    if tilt:
        slope=tilt["linear_fit"].get("M1",{}).get("slope_mm_per_mrad")
        if slope is not None:
            macros.append(rf"\newcommand{{\CentroidSlope}}{{{slope:.9f}}}")
        macros.append(rf"\newcommand{{\MaxPhaseStep}}{{{tilt['phase_step_max_rad']:.6f}}}")
        body=[]
        for r in rows:
            values=([f"{r['theta_mrad']:g}",f"{r['y_centroid_M1_mm']:.5f}",f"{r['y_centroid_M2_mm']:.2e}",f"{r['mode_radius_y_M1_mm']:.6f}",f"{r['zero_angle_overlap_M1']:.6f}",str(r["iterations"]),"是"] if r["converged"] else [f"{r['theta_mrad']:g}"]+[r"\textemdash"]*4+[str(r["iterations"]),"否"])
            body.append(values)
        table(gen / "scan_table.tex","rrrrrrl",r"$\theta$/mrad & $\bar y_1$/mm & $\bar y_2$/mm & $w_{y,1}$/mm & $\eta_0$ & 轮数 & 收敛",body)
        body=[[f"{r['theta_mrad']:g}",tex_number(r["final_field_error"]),tex_number(r["final_intensity_error"]),tex_number(r["final_eigenvalue_error"]),tex_number(r["loss_M1"])] for r in rows if not r["converged"]]
        table(gen / "unconverged_table.tex","rrrrr",r"$\theta$/mrad & $\delta_U$ & $\delta_I$ & $\delta_\gamma$ & 末态投影损耗",body)
    if validation:
        selected={"direct_vs_fft_1d_relative":"一维直接积分/FFT 复场相对差","direct_vs_fft_2d_relative":"二维直接积分/FFT 复场相对差","production_grid_dense_vs_separable":"生产网格完整二维/分离往返相对差","unapertured_periodic_power_error":"无孔径周期传播功率误差","dense_baseline_gamma_error":"完整二维/分离基模本征值差","baseline_1d_2d_profile_overlap_error":"1D/2D 中心截面重叠误差"}
        body=[[selected[r["name"]],tex_number(r["value"]),tex_number(r["threshold"])] for r in validation["checks"] if r["name"] in selected]
        table(gen / "verification_table.tex","lrr",r"检查项 & 实际误差 & 阈值",body)
        body=[[str(r["n"]),f"{r['window_factor']:g}",f"{r['dx_m']*1e6:.4f}",f"{r['width_m']*1e3:.9f}",f"{r['hg00_overlap']:.12f}",str(r["iterations"])] for r in validation["grid_window_convergence"]]
        table(gen / "grid_table.tex","rrrrrr",r"$N$ & 窗口半宽/$w_m$ & $\Delta y/\mu$m & $w_y$/mm & HG00 重叠 & 轮数",body)
        body=[[f"{r['theta_mrad']:g}",str(r["n"]),f"{r['m1']['y_centroid_m']*1e3:.6f}",f"{r['m1']['mode_radius_y_m']*1e3:.6f}",tex_number(r["final_errors"]["delta"]),"是" if r["converged"] else "否"] for r in validation["tilt_grid_convergence"]]
        table(gen / "tilt_grid_table.tex","rrrrrl",r"$\theta$/mrad & $N$ & $\bar y_1$/mm & $w_{y,1}$/mm & $\delta_U$ & 收敛",body)
    (gen / "numbers.tex").write_text("\n".join(macros)+"\n")
    result=dict(source_run=run.relative_to(ROOT).as_posix(),images=records,software=dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__),generated_tables=[p.relative_to(ROOT).as_posix() for p in sorted(gen.glob("*.tex"))])
    manifest_path.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(f"Rendered {len(records)} figures; tables: {gen}")


if __name__=="__main__":
    main()
