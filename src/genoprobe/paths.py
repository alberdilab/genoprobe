"""Output path management for genoprobe."""

from __future__ import annotations

from pathlib import Path


def resolve_output_root(output: str | Path) -> Path:
    return Path(output).expanduser().resolve()


def targets_dir(output: Path) -> Path:
    return output / "targets"


def probes_dir(output: Path) -> Path:
    return output / "probes"


def screen_dir(output: Path) -> Path:
    return output / "screen"


def panels_dir(output: Path) -> Path:
    return output / "panels"


def index_dir(output: Path) -> Path:
    return output / "index"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
