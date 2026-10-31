"""Find the submission directory from source or a relocated executable."""
from pathlib import Path
import sys


def project_root() -> Path:
    start = Path(sys.executable) if getattr(sys, "frozen", False) else Path(__file__)
    for candidate in start.resolve().parents:
        if (candidate / "config/parameter.yaml").is_file():
            return candidate
    raise FileNotFoundError("Keep config/parameter.yaml beside the bin directory")
