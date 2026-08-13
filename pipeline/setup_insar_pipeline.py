"""Create project virtual environments and install all workflow packages."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECTS = (
    "InSAR_KML_Cropper",
    "InSAR_Filtering",
    "InSAR_Unwrapping",
    "InSAR_Detrending",
    "InSAR_Timeseries",
)


def main() -> None:
    for index, name in enumerate(PROJECTS, 1):
        project = ROOT / name
        environment = project / ".venv"
        python = environment / "Scripts" / "python.exe"
        print(f"\n[{index}/{len(PROJECTS)}] Setting up {name}", flush=True)
        if not python.exists():
            subprocess.run([sys.executable, "-m", "venv", str(environment)], check=True)
        subprocess.run([str(python), "-m", "pip", "install", "-e", str(project)], check=True)
    print("\nAll InSAR environments are ready.")


if __name__ == "__main__":
    main()
