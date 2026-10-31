"""Check report assets and portable executables, then save submission metadata."""
from pathlib import Path
import hashlib
import json
import platform
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((ROOT / "report/manifest.json").read_text())
    run = ROOT / manifest["source_run"]
    validation = json.loads((run / "validation.json").read_text())
    if not validation["passed"] or not all(item["passed"] for item in validation["checks"]):
        raise ValueError("Rust validation failed")
    for image in manifest["images"]:
        if sha256(ROOT / image["file"]) != image["sha256"]:
            raise ValueError(f"Figure hash mismatch: {image['file']}")
        for source in [*image["source_data"], image["parameter_snapshot"], image["source_image"]]:
            if not (ROOT / source).is_file():
                raise FileNotFoundError(source)
    for table in manifest["generated_tables"]:
        if not (ROOT / table).is_file():
            raise FileNotFoundError(table)
    text = (ROOT / "report/report.tex").read_text()
    forbidden = "待核查|待核验|待人工|人工审阅|尚须|草稿|重新实现|重实现|后续改进|仍需|需要继续"
    pdf_text = subprocess.check_output(["pdftotext", str(ROOT / "report/report.pdf"), "-"], text=True)
    if re.search(forbidden, text + pdf_text):
        raise ValueError("Report contains draft or pending-review wording")
    figures = re.findall(r"\\fig\{([^}]+)\}", text)
    if not all((ROOT / "assets" / name).is_file() for name in figures):
        raise ValueError("Missing report figure")
    git_checked = (ROOT / ".git").exists()
    if git_checked:
        tracked = subprocess.check_output(["git", "ls-files", "bin", "results", "submission_checks.json"], cwd=ROOT, text=True)
        if tracked:
            raise ValueError("Generated executables or results are tracked by Git")
    with tempfile.TemporaryDirectory(prefix="laser-relocated-") as temp:
        moved = Path(temp) / "submission"
        shutil.copytree(ROOT / "bin", moved / "bin")
        shutil.copytree(ROOT / "config", moved / "config")
        for name in ("confocal-fox-li", "laser-tools"):
            subprocess.run([str(moved / "bin" / name), "--help"], cwd=temp, check=True, capture_output=True)
        subprocess.run([str(moved / "bin/confocal-fox-li"), "all"], cwd=temp, check=True)
        moved_run = next((moved / "results/runs").iterdir())
        subprocess.run([str(moved / "bin/laser-tools"), "plot", str(moved_run)], cwd=temp, check=True)
        if len(list((moved_run / "figures").glob("*.png"))) != 14:
            raise ValueError("Relocated plotting output is incomplete")
        subprocess.run([str(moved / "bin/laser-tools"), "prepare", str(moved_run)], cwd=temp, check=True)
        if not (moved / "assets/manifest.json").is_file():
            raise ValueError("Relocated report preparation failed")
    info = subprocess.check_output(["pdfinfo", str(ROOT / "report/report.pdf")], text=True)
    pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    record = {
        "platform": platform.platform(), "architecture": platform.machine(), "libc": platform.libc_ver(),
        "source_run": manifest["source_run"], "numerical_implementation": "Rust",
        "python_role": "plotting, saved-field reconstruction and report asset preparation",
        "numerical_validation": {"passed": True, "count": len(validation["checks"])},
        "relocated_executables": {"passed": True, "tasks": ["all", "plot", "prepare"], "plotted_figures": 14},
        "report": {"pages": pages, "figure_references": len(figures), "asset_manifest_images": len(manifest["images"]),
                   "source_hash_checks_passed": True, "draft_or_pending_review_wording_present": False},
        "git_tracking_checked": git_checked,
        "generated_files_tracked": False if git_checked else None,
        "executables": [{"path": f"bin/{p.name}", "size_bytes": p.stat().st_size, "sha256": sha256(p)}
                        for p in sorted((ROOT / "bin").iterdir()) if p.is_file()],
    }
    (ROOT / "submission_checks.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    print(f"PASS: {len(validation['checks'])} Rust checks, {len(figures)} report figures, portable computation and plotting")


if __name__ == "__main__":
    main()
