"""Cross-check saved Rust production fields against the original dense Python solver.

Uses the saved YAML snapshot; writes a new verification file (never replaces it).
The production physics and iteration limits are not overridden.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/laser-matplotlib")
from .paths import project_root
ROOT = project_root()
import numpy as np
import yaml
from Python.cavity import Cavity1D, Cavity2D
from Python.fox_li import run_fox_li
from Python.grid import make_grid_1d, make_grid_2d
from Python.metrics import mode_overlap


def factors(path: Path, dim: int):
    table = np.genfromtxt(path, delimiter=",", names=True)
    fields = []
    for plane in ("m1", "m2"):
        x = table[f"{plane}_x_real"] + 1j * table[f"{plane}_x_imag"]
        if dim == 2:
            y = table[f"{plane}_y_real"] + 1j * table[f"{plane}_y_imag"]
            fields.append(y[:, None] * x[None, :])
        else:
            fields.append(x)
    return table["coordinate_m"], *fields


def initial(parameters, grid, dim: int, angle=0.0):
    wavelength = parameters["project"]["wavelength_m"]
    length = parameters["project"]["cavity_length_m"]
    k = 2 * np.pi / wavelength
    wm2 = wavelength * length / np.pi
    if dim == 1:
        field = np.exp(-grid.x**2 / wm2 - 1j * k * grid.x**2 / (2 * length))
    else:
        xx, yy = np.meshgrid(grid.x, grid.y, indexing="xy")
        if parameters["tilt"]["mirror"] == "M2":
            offset, slope = length * np.sin(angle), -np.sin(angle)
        else:
            offset, slope = 0.0, np.sin(angle)
        r2 = xx**2 + (yy - offset)**2
        field = np.exp(-r2 / wm2 - 1j * k * r2 / (2 * length) + 1j * k * slope * yy)
    return field / np.linalg.norm(field.ravel())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    run = args.run_directory.resolve()
    target = run / "python_verification.json"
    if target.exists():
        raise FileExistsError(target)
    cfg = yaml.safe_load((run / "parameter_used.yaml").read_text())
    checks = []
    artifacts = []

    def check(name, error, limit):
        passed = bool(np.isfinite(error) and error < limit)
        checks.append(dict(name=name, value=float(error), threshold=limit, passed=passed))
        print(name, f"{error:.3e}", "PASS" if passed else "FAIL", flush=True)

    for dim in (1, 2):
        grid = make_grid_1d(cfg) if dim == 1 else make_grid_2d(cfg)
        cavity = Cavity1D(cfg, grid) if dim == 1 else Cavity2D(cfg, grid)
        py = run_fox_li(cavity.round_trip, initial(cfg, grid, dim), cfg, area_element=grid.dx**dim)
        _, rust, _ = factors(run / f"baseline/{dim}d/factors.csv", dim)
        summary = json.loads((run / f"baseline/{dim}d/summary.json").read_text())
        gamma = complex(summary["gamma_real"], summary["gamma_imag"])
        check(f"baseline_{dim}d_overlap_error", max(0.0, 1 - mode_overlap(py.field, rust)), 1e-12)
        check(f"baseline_{dim}d_gamma_error", abs(py.gamma - gamma), 1e-10)
        check(f"baseline_{dim}d_iterations", abs(py.iterations - summary["iterations"]), 0.5)
        artifacts.append(dict(dimension=dim, iterations=py.iterations, converged=py.converged,
                              gamma_real=py.gamma.real, gamma_imag=py.gamma.imag))
        np.savez_compressed(run / f"baseline/{dim}d/python_verification_fields.npz",
                            x=grid.x, field_m1=py.field)

    # Full Python solve at every configured converged scan point; for un-converged
    # points verify a single dense round trip at the saved terminal field.
    grid = make_grid_2d(cfg)
    for angle in cfg["tilt"]["angles_mrad"]:
        theta = angle * 1e-3
        cavity = Cavity2D(cfg, grid, theta=theta)
        name = f"theta_{angle:g}mrad"
        _, rust, rust_m2 = factors(run / f"tilt/factors_{name}.csv", 2)
        summary = json.loads((run / f"tilt/summary_{name}.json").read_text())
        raw = cavity.round_trip(rust)
        gamma = np.vdot(rust, raw) / np.vdot(rust, rust)
        check(f"{name}_terminal_gamma_error", abs(gamma - complex(summary["gamma_real"], summary["gamma_imag"])),
              1e-6 if not summary["converged"] else 1e-10)
        from Python.propagation import fresnel_propagate_2d
        py_m2 = cavity.mirror2(fresnel_propagate_2d(rust, grid, cfg["project"]["wavelength_m"], cfg["project"]["cavity_length_m"]))
        check(f"{name}_m2_relative_error", np.linalg.norm(py_m2-rust_m2)/np.linalg.norm(py_m2), 1e-9)
        if summary["converged"]:
            py = run_fox_li(cavity.round_trip, initial(cfg, grid, 2, theta), cfg, area_element=grid.dx**2)
            check(f"{name}_solve_overlap_error", max(0.0, 1-mode_overlap(py.field, rust)), 1e-12)
            check(f"{name}_solve_iterations", abs(py.iterations-summary["iterations"]), 0.5)

    result = dict(passed=all(c["passed"] for c in checks), checks=checks, baseline=artifacts,
                  numpy_version=np.__version__, notes="Un-converged tilt terminal fields were checked by one dense round trip, not independently iterated for 3000 dense steps.")
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2))
    if not result["passed"]:
        raise SystemExit("Python/Rust verification failed")


if __name__ == "__main__":
    main()
