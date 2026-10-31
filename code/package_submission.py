"""Create the submission archive without development caches."""
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
NAME = "陈恪瑾_U202413925_光学谐振腔仿真_linux-x86_64"
SKIP_DIRS = {"target", "__pycache__", ".venv", ".git", "build", "dist"}
SKIP_SUFFIXES = {".pyc", ".aux", ".log", ".out", ".toc", ".fls", ".fdb_latexmk", ".spec"}


def main():
    output = ROOT / f"{NAME}.zip"
    files = []
    for name in ("README.md", "LICENSE", "pyproject.toml", "uv.lock", "submission_checks.json",
                 ".github", "bin", "code", "config", "assets", "report", "results"):
        path = ROOT / name
        if not path.exists():
            raise FileNotFoundError(path)
        candidates = path.rglob("*") if path.is_dir() else [path]
        files.extend(p for p in candidates if p.is_file()
                     and not (set(p.relative_to(ROOT).parts) & SKIP_DIRS)
                     and p.suffix not in SKIP_SUFFIXES)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(files):
            archive.write(path, Path(NAME) / path.relative_to(ROOT))
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise ValueError("Archive CRC check failed")
    print(f"{output.name}: {len(files)} files, {output.stat().st_size / 1024**2:.1f} MiB")


if __name__ == "__main__":
    main()
