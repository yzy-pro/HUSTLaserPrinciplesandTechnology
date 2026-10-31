#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
uv sync --locked --group build
cargo build --release --locked --manifest-path code/Cargo.toml
mkdir -p bin
install -m 755 code/target/release/confocal-fox-li bin/confocal-fox-li
uv run --locked --group build pyinstaller --noconfirm --clean --onefile \
  --name laser-tools --paths code/python --distpath bin \
  --workpath build/pyinstaller --specpath build \
  --collect-submodules matplotlib.backends.backend_agg \
  code/python/launcher.py
