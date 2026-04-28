"""Pure-Python mismatch-tolerant sequence matching for off-target screening."""

from __future__ import annotations

from dataclasses import dataclass, field

DEFAULT_MAX_MISMATCHES: int = 2
DEFAULT_MISMATCH_PENALTY_PER_MISMATCH: float = 0.5


@dataclass(slots=True)
class MatchHit:
    subject_name: str
    position: int       # 0-based start on subject
    mismatches: int
    strand: str         # '+' or '-'


@dataclass
class MatchSummary:
    query: str
    on_target_hits: list[MatchHit] = field(default_factory=list)
    off_target_hits: list[MatchHit] = field(default_factory=list)

    @property
    def on_target_score(self) -> float:
        if not self.on_target_hits:
            return 0.0
        best = min(h.mismatches for h in self.on_target_hits)
        return mismatch_support(best, DEFAULT_MISMATCH_PENALTY_PER_MISMATCH)

    @property
    def off_target_score(self) -> float:
        return sum(
            mismatch_support(h.mismatches, DEFAULT_MISMATCH_PENALTY_PER_MISMATCH)
            for h in self.off_target_hits
        )


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def count_mismatches(query: str, subject_window: str) -> int:
    return sum(a != b for a, b in zip(query, subject_window))


def find_hits(
    query: str,
    subject: str,
    subject_name: str,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> list[MatchHit]:
    """Find all approximate matches of query in subject (both strands)."""
    hits: list[MatchHit] = []
    qlen = len(query)
    rc_query = reverse_complement(query)
    for i in range(len(subject) - qlen + 1):
        window = subject[i : i + qlen]
        mm = count_mismatches(query, window)
        if mm <= max_mismatches:
            hits.append(MatchHit(subject_name, i, mm, "+"))
        mm_rc = count_mismatches(rc_query, window)
        if mm_rc <= max_mismatches:
            hits.append(MatchHit(subject_name, i, mm_rc, "-"))
    return hits


def mismatch_support(mismatches: int, penalty_per: float) -> float:
    return max(0.0, 1.0 - penalty_per * mismatches)
