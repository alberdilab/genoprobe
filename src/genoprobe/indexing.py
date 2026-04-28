"""Optional Fulgor-based off-target index management."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

DEFAULT_KMER_LENGTH: int = 31
DEFAULT_MINIMIZER_LENGTH: int = 20
DEFAULT_FULGOR_THREADS: int = 4
DEFAULT_INDEX_PENALTY_WEIGHT: float = 1.0


def find_fulgor() -> Path | None:
    exe = shutil.which("fulgor")
    return Path(exe) if exe else None


def find_ggcat() -> Path | None:
    exe = shutil.which("ggcat")
    return Path(exe) if exe else None


def require_fulgor() -> Path:
    path = find_fulgor()
    if path is None:
        raise RuntimeError(
            "Fulgor executable not found. Install Fulgor and ensure it is on PATH, "
            "or omit --index to skip secondary off-target screening."
        )
    return path


def build_fulgor_index(
    fasta_paths: list[Path],
    output_dir: Path,
    *,
    kmer_length: int = DEFAULT_KMER_LENGTH,
    minimizer_length: int = DEFAULT_MINIMIZER_LENGTH,
    threads: int = DEFAULT_FULGOR_THREADS,
) -> Path:
    """Build a Fulgor index from a list of FASTA files. Returns the index path."""
    fulgor = require_fulgor()
    output_dir.mkdir(parents=True, exist_ok=True)
    filelist = output_dir / "genome_list.txt"
    filelist.write_text("\n".join(str(p.resolve()) for p in fasta_paths) + "\n")
    index_path = output_dir / "index.fur"
    cmd = [
        str(fulgor), "build",
        "-i", str(filelist),
        "-o", str(index_path),
        "-k", str(kmer_length),
        "-m", str(minimizer_length),
        "-t", str(threads),
    ]
    subprocess.run(cmd, check=True)
    return index_path


def query_fulgor_index(
    index_path: Path,
    probe_sequences: list[str],
    *,
    threads: int = DEFAULT_FULGOR_THREADS,
) -> dict[str, list[int]]:
    """Query a Fulgor index with probe sequences. Returns probe → list of genome indices."""
    fulgor = require_fulgor()
    # Write probes to a temporary FASTA
    tmp_fasta = index_path.parent / "_query_probes.fa"
    with tmp_fasta.open("w") as fh:
        for i, seq in enumerate(probe_sequences):
            fh.write(f">probe_{i}\n{seq}\n")
    tmp_out = index_path.parent / "_query_results.txt"
    cmd = [
        str(fulgor), "query",
        "-i", str(index_path),
        "-q", str(tmp_fasta),
        "-o", str(tmp_out),
        "-t", str(threads),
    ]
    subprocess.run(cmd, check=True)
    results: dict[str, list[int]] = {}
    if tmp_out.exists():
        for line in tmp_out.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                name = parts[0].lstrip(">")
                hits = [int(x) for x in parts[1:] if x.strip().isdigit()]
                results[name] = hits
    tmp_fasta.unlink(missing_ok=True)
    tmp_out.unlink(missing_ok=True)
    return results
