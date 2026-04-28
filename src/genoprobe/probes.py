"""Probe candidate generation — blockParse logic ported to Python 3."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator

from genoprobe.thermo import (
    DEFAULT_MONOVALENT_SALT_MM,
    DEFAULT_PROBE_CONC_NM,
    DEFAULT_THERMO_BACKEND,
    DEFAULT_THERMO_PENALTY_WEIGHT,
    calc_tm,
    gc_content,
    hairpin_tm,
    homodimer_tm,
    self_complementarity,
    sequence_entropy,
)

DEFAULT_MIN_PROBE_LENGTH: int = 36
DEFAULT_MAX_PROBE_LENGTH: int = 41
DEFAULT_MIN_TM: float = 42.0
DEFAULT_MAX_TM: float = 47.0
DEFAULT_MIN_GC: float = 20.0
DEFAULT_MAX_GC: float = 75.0
DEFAULT_MAX_HOMOPOLYMER_RUN: int = 4
DEFAULT_MAX_DINUCLEOTIDE_RUN: int | None = 4
DEFAULT_MIN_SEQUENCE_ENTROPY: float | None = 0.6
DEFAULT_MAX_SELF_COMPLEMENTARITY: int | None = 14
DEFAULT_MAX_CONTIGUOUS_SELF_COMPLEMENTARITY: int | None = 6
DEFAULT_MAX_HAIRPIN_TM: float | None = None
DEFAULT_MAX_HOMODIMER_TM: float | None = None
DEFAULT_PROBE_SPACING: int = 0
DEFAULT_PROHIBITED_SEQUENCES: tuple[str, ...] = ("AAAAA", "TTTTT", "CCCCC", "GGGGG")


@dataclass(slots=True)
class ProbeCandidate:
    sequence: str
    target_name: str
    start: int          # 0-based start on the target sequence
    end: int            # 0-based exclusive end
    tm: float
    gc: float
    entropy: float
    self_comp_total: int
    self_comp_run: int
    hairpin_tm: float
    homodimer_tm: float
    score: float = 0.0
    passes_filters: bool = True
    fail_reason: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class ProbeDesignConfig:
    min_probe_length: int = DEFAULT_MIN_PROBE_LENGTH
    max_probe_length: int = DEFAULT_MAX_PROBE_LENGTH
    min_tm: float = DEFAULT_MIN_TM
    max_tm: float = DEFAULT_MAX_TM
    min_gc: float = DEFAULT_MIN_GC
    max_gc: float = DEFAULT_MAX_GC
    max_homopolymer_run: int = DEFAULT_MAX_HOMOPOLYMER_RUN
    max_dinucleotide_run: int | None = DEFAULT_MAX_DINUCLEOTIDE_RUN
    min_sequence_entropy: float | None = DEFAULT_MIN_SEQUENCE_ENTROPY
    max_self_complementarity: int | None = DEFAULT_MAX_SELF_COMPLEMENTARITY
    max_contiguous_self_complementarity: int | None = DEFAULT_MAX_CONTIGUOUS_SELF_COMPLEMENTARITY
    max_hairpin_tm: float | None = DEFAULT_MAX_HAIRPIN_TM
    max_homodimer_tm: float | None = DEFAULT_MAX_HOMODIMER_TM
    probe_spacing: int = DEFAULT_PROBE_SPACING
    prohibited_sequences: tuple[str, ...] = DEFAULT_PROHIBITED_SEQUENCES
    thermo_penalty_weight: float = DEFAULT_THERMO_PENALTY_WEIGHT
    thermo_backend: str = DEFAULT_THERMO_BACKEND
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM


def _max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    max_run = 1
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 1
    return max_run


def _max_dinucleotide_run(seq: str) -> int:
    if len(seq) < 4:
        return 0
    max_run = 0
    i = 0
    while i < len(seq) - 1:
        dinuc = seq[i : i + 2]
        count = 1
        j = i + 2
        while j + 1 < len(seq) and seq[j : j + 2] == dinuc:
            count += 1
            j += 2
        max_run = max(max_run, count)
        i += 1
    return max_run


def _score_probe(candidate: ProbeCandidate, config: ProbeDesignConfig) -> float:
    """Compute a composite score in [0, 1]; higher is better."""
    # Tm proximity to midpoint of [min_tm, max_tm]
    tm_mid = (config.min_tm + config.max_tm) / 2.0
    tm_range = max(1.0, (config.max_tm - config.min_tm) / 2.0)
    tm_score = max(0.0, 1.0 - abs(candidate.tm - tm_mid) / tm_range)

    # GC proximity to 50 %
    gc_score = max(0.0, 1.0 - abs(candidate.gc - 50.0) / 30.0)

    base_score = 0.5 * tm_score + 0.3 * gc_score + 0.2 * candidate.entropy

    if config.thermo_penalty_weight == 0.0:
        return base_score

    # Thermodynamic penalties
    penalty = 0.0
    if config.max_self_complementarity is not None and candidate.self_comp_total > 0:
        penalty += min(1.0, candidate.self_comp_total / max(1, config.max_self_complementarity))
    if config.max_hairpin_tm is not None and candidate.hairpin_tm > 0:
        penalty += min(1.0, candidate.hairpin_tm / max(1.0, config.max_hairpin_tm))
    if config.max_homodimer_tm is not None and candidate.homodimer_tm > 0:
        penalty += min(1.0, candidate.homodimer_tm / max(1.0, config.max_homodimer_tm))

    return base_score * (1.0 - config.thermo_penalty_weight * min(1.0, penalty / 3.0))


def generate_candidates(
    sequence: str,
    target_name: str,
    config: ProbeDesignConfig,
) -> list[ProbeCandidate]:
    """Sliding-window probe candidate generation with thermodynamic filtering."""
    sequence = sequence.upper()
    n = len(sequence)
    candidates: list[ProbeCandidate] = []
    last_accepted_end: int = -1

    for start in range(n):
        for probe_len in range(config.min_probe_length, config.max_probe_length + 1):
            end = start + probe_len
            if end > n:
                break

            seq = sequence[start:end]

            # Hard filters — cheap checks first
            if "N" in seq:
                break
            if any(p in seq for p in config.prohibited_sequences):
                continue

            gc = gc_content(seq) * 100.0
            if gc < config.min_gc or gc > config.max_gc:
                continue

            if _max_homopolymer(seq) > config.max_homopolymer_run:
                continue

            if config.max_dinucleotide_run is not None:
                if _max_dinucleotide_run(seq) > config.max_dinucleotide_run:
                    continue

            # Thermodynamic calculations
            tm = calc_tm(seq, config.thermo_backend, config.monovalent_mm, config.probe_conc_nm)
            if tm < config.min_tm or tm > config.max_tm:
                continue

            entropy = sequence_entropy(seq)
            if config.min_sequence_entropy is not None and entropy < config.min_sequence_entropy:
                continue

            sc_total, sc_run = self_complementarity(seq)
            if config.max_self_complementarity is not None and sc_total > config.max_self_complementarity:
                continue
            if (
                config.max_contiguous_self_complementarity is not None
                and sc_run > config.max_contiguous_self_complementarity
            ):
                continue

            h_tm = hairpin_tm(seq, config.monovalent_mm, config.probe_conc_nm)
            if config.max_hairpin_tm is not None and h_tm > config.max_hairpin_tm:
                continue

            hd_tm = homodimer_tm(seq, config.monovalent_mm, config.probe_conc_nm)
            if config.max_homodimer_tm is not None and hd_tm > config.max_homodimer_tm:
                continue

            # Spacing filter
            if config.probe_spacing > 0 and last_accepted_end > 0:
                if start < last_accepted_end + config.probe_spacing:
                    continue

            candidate = ProbeCandidate(
                sequence=seq,
                target_name=target_name,
                start=start,
                end=end,
                tm=tm,
                gc=gc,
                entropy=entropy,
                self_comp_total=sc_total,
                self_comp_run=sc_run,
                hairpin_tm=h_tm,
                homodimer_tm=hd_tm,
            )
            candidate.score = _score_probe(candidate, config)
            candidates.append(candidate)
            last_accepted_end = end
            # Accept the first passing length for this start position
            break

    return sorted(candidates, key=lambda c: c.score, reverse=True)
