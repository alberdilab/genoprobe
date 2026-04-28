"""Thermodynamic heuristics for probe scoring."""

from __future__ import annotations

import math

DEFAULT_MONOVALENT_SALT_MM: float = 390.0
DEFAULT_DIVALENT_SALT_MM: float = 0.0
DEFAULT_DNTP_MM: float = 0.0
DEFAULT_PROBE_CONC_NM: float = 50.0    # total strand concentration (probe + target), Ct/4 = 12.5 nM
DEFAULT_FORMAMIDE_PCT: float = 50.0    # % formamide; correction = -0.72 °C / %
DEFAULT_HYB_TEMP_C: float = 42.0
DEFAULT_THERMO_BACKEND: str = "basic"
DEFAULT_THERMO_PENALTY_WEIGHT: float = 0.15
DEFAULT_PANEL_COMPLEMENTARITY_PENALTY_WEIGHT: float = 0.15
DEFAULT_PANEL_HETERODIMER_PENALTY_WEIGHT: float = 0.0

# Nearest-neighbour parameters (SantaLucia 1998, unified)
_NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_INIT_AT: tuple[float, float] = (2.3, 4.1)
_INIT_GC: tuple[float, float] = (0.1, -2.8)
_R: float = 1.987  # cal / (mol·K)


def gc_content(seq: str) -> float:
    seq = seq.upper()
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq) if seq else 0.0


def basic_tm(
    seq: str,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
) -> float:
    """Salt-adjusted Tm using the Wallace rule for short oligos."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    at = n - gc
    tm = 81.5 + 16.6 * math.log10(monovalent_mm / 1000.0) + 41.0 * (gc / n) - (500.0 / n)
    return tm


def nn_tm(
    seq: str,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM,
) -> float:
    """Nearest-neighbour Tm (SantaLucia 1998) with method-5 salt correction.

    probe_conc_nm is the *total* strand concentration (probe + target).
    Salt correction follows SantaLucia 1998 / Biopython method 5:
    ΔS_salt = 0.368 × (n−1) × ln([Na+] in M), added directly to ΔS.
    """
    seq = seq.upper()
    n = len(seq)
    if n < 2:
        return 0.0

    dh: float = 0.0
    ds: float = 0.0

    # Initiation parameters (SantaLucia 1998, Table 2)
    if seq[0] in "AT":
        dh += _INIT_AT[0]; ds += _INIT_AT[1]
    else:
        dh += _INIT_GC[0]; ds += _INIT_GC[1]
    if seq[-1] in "AT":
        dh += _INIT_AT[0]; ds += _INIT_AT[1]
    else:
        dh += _INIT_GC[0]; ds += _INIT_GC[1]

    for i in range(n - 1):
        dinuc = seq[i : i + 2]
        params = _NN_PARAMS.get(dinuc)
        if params:
            dh += params[0]; ds += params[1]

    # Method-5 salt correction: additive to ΔS (cal/mol/K)
    salt_corr = 0.368 * (n - 1) * math.log(monovalent_mm / 1000.0)
    ds_total = ds + salt_corr

    # Ct/4 for non-self-complementary duplex
    ct = probe_conc_nm * 1e-9
    tm_k = (dh * 1000.0) / (ds_total + _R * math.log(ct / 4.0)) - 273.15
    return tm_k


def sequence_entropy(seq: str) -> float:
    """Shannon entropy of nucleotide composition, normalised 0–1."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    counts = {b: seq.count(b) for b in "ACGT"}
    entropy = 0.0
    for c in counts.values():
        if c > 0:
            p = c / n
            entropy -= p * math.log2(p)
    return entropy / 2.0  # max entropy for 4 bases = 2 bits


def self_complementarity(seq: str) -> tuple[int, int]:
    """Return (total_complementary_bases, max_contiguous_run) vs reverse complement."""
    rc = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    total = sum(a == b for a, b in zip(seq, rc))
    max_run = 0
    run = 0
    for a, b in zip(seq, rc):
        if a == b:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return total, max_run


def hairpin_tm(
    seq: str,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM,
) -> float:
    """Estimate hairpin melting temperature from the most stable stem."""
    seq = seq.upper()
    n = len(seq)
    best_stem: list[str] = []
    for stem_len in range(n // 2, 2, -1):
        arm5 = seq[:stem_len]
        arm3 = seq[n - stem_len :]
        rc3 = arm3.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if arm5 == rc3:
            best_stem = list(arm5)
            break
    if not best_stem:
        return 0.0
    stem_seq = "".join(best_stem)
    return nn_tm(stem_seq, monovalent_mm, probe_conc_nm)


def homodimer_tm(
    seq: str,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM,
) -> float:
    """Estimate homodimer Tm from self-complementarity."""
    total, _run = self_complementarity(seq)
    if total < 4:
        return 0.0
    # Use the complementary portion as a proxy duplex
    rc = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    paired = "".join(a for a, b in zip(seq, rc) if a == b)
    if len(paired) < 2:
        return 0.0
    return nn_tm(paired, monovalent_mm, probe_conc_nm)


def heterodimer_tm(
    seq_a: str,
    seq_b: str,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM,
) -> float:
    """Estimate worst-case heterodimer Tm between two probe sequences."""
    rc_b = seq_b.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    best_run = 0
    run_seq: list[str] = []
    current: list[str] = []
    # Sliding alignment of seq_a against rc_b
    for offset in range(-(len(seq_b) - 1), len(seq_a)):
        run = 0
        for i in range(len(seq_a)):
            j = i - offset
            if 0 <= j < len(rc_b) and seq_a[i] == rc_b[j]:
                run += 1
            else:
                run = 0
            if run > best_run:
                best_run = run
    if best_run < 4:
        return 0.0
    proxy = "GC" * (best_run // 2) + ("G" if best_run % 2 else "")
    return nn_tm(proxy, monovalent_mm, probe_conc_nm)


def calc_tm(
    seq: str,
    backend: str = DEFAULT_THERMO_BACKEND,
    monovalent_mm: float = DEFAULT_MONOVALENT_SALT_MM,
    probe_conc_nm: float = DEFAULT_PROBE_CONC_NM,
    formamide_pct: float = DEFAULT_FORMAMIDE_PCT,
) -> float:
    """Calculate Tm with salt and formamide correction.

    The formamide correction (-0.72 °C per % formamide) is applied after the
    NN Tm, matching the blockParse / OligoMiner convention.
    """
    if backend == "primer3":
        try:
            import primer3
            tm = float(primer3.calc_tm(
                seq,
                mv_conc=monovalent_mm,
                dv_conc=DEFAULT_DIVALENT_SALT_MM,
                dntp_conc=DEFAULT_DNTP_MM,
                dna_conc=probe_conc_nm,
            ))
        except ImportError:
            tm = nn_tm(seq, monovalent_mm, probe_conc_nm)
    else:
        tm = nn_tm(seq, monovalent_mm, probe_conc_nm)
    return tm - 0.72 * formamide_pct
