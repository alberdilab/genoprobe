"""FASTA genome handling via pyfaidx."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from pyfaidx import Fasta, FastaIndexingError


@dataclass(slots=True)
class GenomeRecord:
    name: str
    sequence: str
    source_file: Path

    @property
    def length(self) -> int:
        return len(self.sequence)


def load_fasta(path: str | Path) -> Fasta:
    p = Path(path)
    try:
        return Fasta(str(p), as_raw=True, sequence_always_upper=True)
    except FastaIndexingError as exc:
        raise ValueError(f"Cannot index FASTA at {p}: {exc}") from exc


def iter_sequences(fasta: Fasta, source: Path) -> Iterator[GenomeRecord]:
    for name in fasta.keys():
        yield GenomeRecord(
            name=name,
            sequence=str(fasta[name]),
            source_file=source,
        )


def load_genomes(paths: list[str | Path]) -> list[GenomeRecord]:
    """Load all sequences from one or more FASTA files."""
    records: list[GenomeRecord] = []
    for p in paths:
        p = Path(p)
        fasta = load_fasta(p)
        records.extend(iter_sequences(fasta, p))
    return records


def extract_region(
    fasta: Fasta,
    seqid: str,
    start: int,
    end: int,
) -> str:
    """Extract sequence for a genomic region (1-based, inclusive)."""
    return str(fasta[seqid][start - 1 : end]).upper()


def count_kmers(sequence: str, k: int = 18) -> dict[str, int]:
    """Count all k-mers in a sequence."""
    counts: dict[str, int] = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        if "N" not in kmer:
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def build_kmer_index(records: list[GenomeRecord], k: int = 18) -> dict[str, int]:
    """Build a genome-wide k-mer frequency index across all records."""
    index: dict[str, int] = {}
    for rec in records:
        for kmer, count in count_kmers(rec.sequence, k).items():
            index[kmer] = index.get(kmer, 0) + count
    return index


def query_kmer_frequency(index: dict[str, int], sequence: str, k: int = 18) -> int:
    """Return the maximum k-mer frequency for any k-mer in the sequence."""
    max_freq = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        max_freq = max(max_freq, index.get(kmer, 0))
    return max_freq
