"""Publish saved figures to assets/ and saved tables to report/generated/."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from .paths import project_root
ROOT = project_root()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    run = parser.parse_args().run_directory.resolve()
    source = json.loads((run / "plots_manifest.json").read_text())
    validation_path = run / "validation.json"
    if validation_path.exists():
        validation = json.loads(validation_path.read_text())
        if not validation["passed"]:
            raise ValueError("Rust numerical validation has not passed")
    # Check every source before updating final report assets.
    for entry in source["images"]:
        original = ROOT / entry["file"]
        if hashlib.sha256(original.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"source hash mismatch: {original}")
    assets = ROOT / "assets"
    assets.mkdir(exist_ok=True)
    generated = ROOT / "report/generated"
    generated.mkdir(parents=True, exist_ok=True)
    frozen = []
    for entry in source["images"]:
        original = ROOT / entry["file"]
        target = assets / original.name
        shutil.copy2(original, target)
        frozen.append(dict(entry, file=target.relative_to(ROOT).as_posix(),
                           source_image=entry["file"]))
    for original in (run / "report_generated").glob("*.tex"):
        shutil.copy2(original, generated / original.name)
    (generated / "paths.tex").write_text("\\graphicspath{{../assets/}}\n")
    manifest = dict(report="report/report.tex", source_run=run.relative_to(ROOT).as_posix(),
                    images=frozen,
                    validation=validation_path.relative_to(ROOT).as_posix() if validation_path.exists() else None,
                    generated_tables=[p.relative_to(ROOT).as_posix() for p in sorted(generated.glob("*.tex"))])
    for path in (ROOT / "report/manifest.json", assets / "manifest.json"):
        path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(f"Report assets: {len(frozen)} figures in assets/")


if __name__ == "__main__":
    main()
