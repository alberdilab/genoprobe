"""GFF3 and GTF annotation parser with auto-format detection."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


@dataclass(slots=True)
class AnnotationRecord:
    seqid: str
    source: str
    feature: str
    start: int        # 1-based, inclusive
    end: int          # 1-based, inclusive
    score: str
    strand: str
    phase: str
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def get_attr(self, key: str, default: str = "") -> str:
        return self.attributes.get(key, default)


def detect_format(path: Path) -> str:
    """Return 'gff3' or 'gtf' based on file extension and content."""
    suffix = path.suffix.lower()
    if suffix in {".gff", ".gff3"}:
        return "gff3"
    if suffix == ".gtf":
        return "gtf"
    # Fall back to inspecting the first non-comment line
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # GTF uses gene_id "value"; GFF3 uses ID=value
            if 'gene_id "' in line or 'transcript_id "' in line:
                return "gtf"
            return "gff3"
    return "gff3"


def _parse_gff3_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, _, value = part.partition("=")
            attrs[key.strip()] = value.strip()
    return attrs


def _parse_gtf_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        tokens = part.split(None, 1)
        if len(tokens) == 2:
            key = tokens[0]
            value = tokens[1].strip().strip('"')
            attrs[key] = value
    return attrs


def _iter_records(path: Path, fmt: str) -> Iterator[AnnotationRecord]:
    parse_attrs = _parse_gff3_attributes if fmt == "gff3" else _parse_gtf_attributes
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            yield AnnotationRecord(
                seqid=cols[0],
                source=cols[1],
                feature=cols[2],
                start=start,
                end=end,
                score=cols[5],
                strand=cols[6],
                phase=cols[7],
                attributes=parse_attrs(cols[8]),
            )


def load_annotation(
    path: str | Path,
    *,
    features: list[str] | None = None,
) -> list[AnnotationRecord]:
    """Load all annotation records, optionally filtered to specific feature types."""
    p = Path(path)
    fmt = detect_format(p)
    records = [r for r in _iter_records(p, fmt)]
    if features:
        feature_set = {f.lower() for f in features}
        records = [r for r in records if r.feature.lower() in feature_set]
    return records


def summarize_annotation(records: list[AnnotationRecord]) -> dict[str, int]:
    """Return feature type → count mapping for a list of records."""
    counts: dict[str, int] = {}
    for r in records:
        counts[r.feature] = counts.get(r.feature, 0) + 1
    return dict(sorted(counts.items()))
