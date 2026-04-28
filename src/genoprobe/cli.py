"""Command-line interface for genoprobe."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import csv
import json
from pathlib import Path
import sys
from typing import Any

from genoprobe import __version__
from genoprobe.annotation import load_annotation, summarize_annotation
from genoprobe.genome import (
    build_kmer_index,
    extract_region,
    load_fasta,
    load_genomes,
    query_kmer_frequency,
)
from genoprobe.indexing import (
    DEFAULT_FULGOR_THREADS,
    DEFAULT_KMER_LENGTH,
    DEFAULT_MINIMIZER_LENGTH,
    build_fulgor_index,
    find_fulgor,
)
from genoprobe.matching import (
    DEFAULT_MAX_MISMATCHES,
    DEFAULT_MISMATCH_PENALTY_PER_MISMATCH,
    MatchSummary,
    find_hits,
)
from genoprobe.panels import (
    DEFAULT_MAX_PANEL_CONTIGUOUS_COMPLEMENTARITY,
    DEFAULT_MAX_PANEL_HETERODIMER_TM,
    DEFAULT_MAX_PROBES_PER_TARGET,
    DEFAULT_PANEL_COMPLEMENTARITY_PENALTY_WEIGHT,
    DEFAULT_PANEL_HETERODIMER_PENALTY_WEIGHT,
    PanelConfig,
    assemble_panel,
)
from genoprobe.parallel import DEFAULT_AUTO_WORKERS, resolve_worker_count
from genoprobe.paths import (
    ensure_dir,
    index_dir,
    panels_dir,
    probes_dir,
    resolve_output_root,
    screen_dir,
    targets_dir,
)
from genoprobe.probes import (
    DEFAULT_MAX_GC,
    DEFAULT_MAX_HOMOPOLYMER_RUN,
    DEFAULT_MAX_PROBE_LENGTH,
    DEFAULT_MAX_TM,
    DEFAULT_MIN_GC,
    DEFAULT_MIN_PROBE_LENGTH,
    DEFAULT_MIN_TM,
    DEFAULT_PROBE_SPACING,
    ProbeDesignConfig,
    generate_candidates,
)
from genoprobe.profile_defaults import get_stage_profile_defaults, get_stage_profile_names
from genoprobe.reports import (
    write_panels_report,
    write_probes_report,
    write_screen_report,
    write_targets_report,
)
from genoprobe.thermo import (
    DEFAULT_FORMAMIDE_PCT,
    DEFAULT_HYB_TEMP_C,
    DEFAULT_MONOVALENT_SALT_MM,
    DEFAULT_PROBE_CONC_NM,
    DEFAULT_THERMO_BACKEND,
    DEFAULT_THERMO_PENALTY_WEIGHT,
)

try:
    from rich.console import Console
    from rich.table import Table
    _RICH = True
except ImportError:
    _RICH = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print(msg: str) -> None:
    print(msg, file=sys.stderr)


def _die(msg: str, code: int = 1) -> None:
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(code)


def _resolve_profile(stage: str, profile: str | None, default: str = "balanced") -> dict[str, Any]:
    name = (profile or default).strip().lower()
    try:
        return get_stage_profile_defaults(stage, name)
    except ValueError as exc:
        _die(str(exc))


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# targets command
# ---------------------------------------------------------------------------

def cmd_targets(args: argparse.Namespace) -> int:
    """Extract or define target regions and write FASTA + BED."""
    output = resolve_output_root(args.output)
    out_dir = ensure_dir(targets_dir(output))

    genomes = args.genomes
    mode = args.mode.lower()

    if mode not in {"genome", "gene"}:
        _die(f"--mode must be 'genome' or 'gene', got '{mode}'.")

    if mode == "gene" and not args.annotation:
        _die("--annotation is required in gene mode.")

    _print(f"[genoprobe targets] mode={mode}  genomes={len(genomes)}")

    target_records: list[dict[str, Any]] = []
    fasta_path = out_dir / "targets.fa"
    bed_path = out_dir / "targets.bed"

    with fasta_path.open("w") as fa_out, bed_path.open("w") as bed_out:

        if mode == "genome":
            for genome_path in genomes:
                fasta = load_fasta(genome_path)
                for seq_name in fasta.keys():
                    seq = str(fasta[seq_name]).upper()
                    if args.region:
                        # Parse region string: seqid:start-end (1-based)
                        try:
                            rid, coords = args.region.split(":")
                            rstart, rend = (int(x) for x in coords.split("-"))
                            if seq_name != rid:
                                continue
                            seq = seq[rstart - 1 : rend]
                            label = f"{seq_name}:{rstart}-{rend}"
                            bed_out.write(f"{seq_name}\t{rstart - 1}\t{rend}\t{label}\n")
                        except (ValueError, KeyError):
                            _die(f"Cannot parse --region '{args.region}'. Use seqid:start-end.")
                    else:
                        label = seq_name
                        bed_out.write(f"{seq_name}\t0\t{len(seq)}\t{label}\n")
                    fa_out.write(f">{label}\n{seq}\n")
                    target_records.append({"target": label, "length": len(seq)})

        else:  # gene mode
            records = load_annotation(
                args.annotation,
                features=args.feature if args.feature else None,
            )
            _print(f"  Loaded {len(records)} annotation records.")
            summary = summarize_annotation(records)
            _print(f"  Feature counts: {summary}")

            for genome_path in genomes:
                fasta = load_fasta(genome_path)
                for rec in records:
                    if rec.seqid not in fasta.keys():
                        continue
                    seq = extract_region(fasta, rec.seqid, rec.start, rec.end)
                    attr_id = rec.get_attr("ID") or rec.get_attr("gene_id") or rec.get_attr("Name")
                    label = attr_id or f"{rec.seqid}:{rec.start}-{rec.end}"
                    fa_out.write(f">{label}\n{seq}\n")
                    bed_out.write(
                        f"{rec.seqid}\t{rec.start - 1}\t{rec.end}\t{label}\t"
                        f"{rec.strand}\t{rec.feature}\n"
                    )
                    target_records.append({
                        "target": label,
                        "seqid": rec.seqid,
                        "start": rec.start,
                        "end": rec.end,
                        "strand": rec.strand,
                        "feature": rec.feature,
                        "length": rec.length,
                    })

    summary_path = out_dir / "targets_summary.json"
    summary_path.write_text(
        json.dumps({"mode": mode, "target_count": len(target_records)}, indent=2),
        encoding="utf-8",
    )
    write_targets_report(out_dir, {r["target"]: r["length"] for r in target_records})

    _print(f"  Wrote {len(target_records)} targets → {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# probes command
# ---------------------------------------------------------------------------

def cmd_probes(args: argparse.Namespace) -> int:
    """Generate probe candidates from target sequences."""
    output = resolve_output_root(args.output)
    targets_path = targets_dir(output) / "targets.fa"
    if not targets_path.exists():
        _die(f"targets.fa not found at {targets_path}. Run 'genoprobe targets' first.")

    out_dir = ensure_dir(probes_dir(output))
    profile = _resolve_profile("probes", args.profile)

    config = ProbeDesignConfig(
        min_probe_length=args.min_length or profile["min_probe_length"],
        max_probe_length=args.max_length or profile["max_probe_length"],
        min_tm=args.min_tm or profile["min_tm"],
        max_tm=args.max_tm or profile["max_tm"],
        min_gc=args.min_gc or profile["min_gc"],
        max_gc=args.max_gc or profile["max_gc"],
        max_homopolymer_run=profile["max_homopolymer_run"],
        max_dinucleotide_run=profile["max_dinucleotide_run"],
        min_sequence_entropy=profile["min_sequence_entropy"],
        max_self_complementarity=profile["max_self_complementarity"],
        max_contiguous_self_complementarity=profile["max_contiguous_self_complementarity"],
        max_hairpin_tm=profile["max_hairpin_tm"],
        max_homodimer_tm=profile["max_homodimer_tm"],
        probe_spacing=args.probe_spacing if args.probe_spacing is not None else profile["probe_spacing"],
        thermo_penalty_weight=profile["thermo_penalty_weight"],
        thermo_backend=args.thermo_backend or DEFAULT_THERMO_BACKEND,
        monovalent_mm=args.monovalent_mm or DEFAULT_MONOVALENT_SALT_MM,
        probe_conc_nm=args.probe_conc or DEFAULT_PROBE_CONC_NM,
        formamide_pct=args.formamide if args.formamide is not None else DEFAULT_FORMAMIDE_PCT,
    )

    workers = resolve_worker_count(args.workers)
    _print(f"[genoprobe probes] profile={args.profile or 'balanced'}  workers={workers}")

    from pyfaidx import Fasta
    fasta = Fasta(str(targets_path), as_raw=True, sequence_always_upper=True)

    all_rows: list[dict[str, Any]] = []
    target_counts: dict[str, int] = {}

    for target_name in fasta.keys():
        seq = str(fasta[target_name])
        candidates = generate_candidates(seq, target_name, config)
        target_counts[target_name] = len(candidates)
        for c in candidates:
            all_rows.append({
                "target": c.target_name,
                "start": c.start,
                "end": c.end,
                "sequence": c.sequence,
                "length": c.length,
                "tm": round(c.tm, 2),
                "gc": round(c.gc, 2),
                "entropy": round(c.entropy, 4),
                "self_comp_total": c.self_comp_total,
                "self_comp_run": c.self_comp_run,
                "hairpin_tm": round(c.hairpin_tm, 2),
                "homodimer_tm": round(c.homodimer_tm, 2),
                "score": round(c.score, 4),
            })

    candidates_path = out_dir / "candidates.tsv"
    _write_tsv(candidates_path, all_rows)
    summary_path = out_dir / "probes_summary.json"
    summary_path.write_text(
        json.dumps({"total_candidates": len(all_rows), "per_target": target_counts}, indent=2),
        encoding="utf-8",
    )
    write_probes_report(out_dir, target_counts)

    _print(f"  Generated {len(all_rows)} candidates across {len(target_counts)} targets → {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# screen command
# ---------------------------------------------------------------------------

def cmd_screen(args: argparse.Namespace) -> int:
    """Off-target screening via mismatch matching (+ optional Fulgor index)."""
    output = resolve_output_root(args.output)
    candidates_path = probes_dir(output) / "candidates.tsv"
    targets_bed = targets_dir(output) / "targets.bed"

    if not candidates_path.exists():
        _die(f"candidates.tsv not found at {candidates_path}. Run 'genoprobe probes' first.")

    out_dir = ensure_dir(screen_dir(output))
    profile = _resolve_profile("screen", args.profile)

    max_mm = args.max_mismatches if args.max_mismatches is not None else profile["max_mismatches"]
    max_kmer = args.max_kmer_frequency if args.max_kmer_frequency is not None else profile["max_kmer_frequency"]
    offtarget_pw = profile["offtarget_penalty_weight"]
    min_ontarget = profile["min_ontarget_fraction"]

    _print(f"[genoprobe screen] profile={args.profile or 'balanced'}  max_mismatches={max_mm}")

    # Load candidate probes
    with candidates_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        candidates = list(reader)

    # Load genome sequences for matching
    genome_records = load_genomes(args.genomes)
    _print(f"  Loaded {len(genome_records)} genome sequences for off-target screening.")

    # Build k-mer index if filtering enabled
    kmer_index: dict[str, int] = {}
    if max_kmer is not None:
        _print(f"  Building k-mer index (k=18)...")
        kmer_index = build_kmer_index(genome_records, k=18)

    # Parse target BED to know on-target intervals
    target_intervals: dict[str, list[tuple[int, int]]] = {}
    if targets_bed.exists():
        for line in targets_bed.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                seqid, start, end = parts[0], int(parts[1]), int(parts[2])
                target_intervals.setdefault(seqid, []).append((start, end))

    screened_rows: list[dict[str, Any]] = []
    target_passing: dict[str, int] = {}

    for row in candidates:
        seq = row["sequence"]
        target_name = row["target"]

        # k-mer frequency filter
        if max_kmer is not None and kmer_index:
            max_freq = query_kmer_frequency(kmer_index, seq)
            if max_freq > max_kmer:
                continue

        # Mismatch matching across all genome sequences
        on_hits = []
        off_hits = []
        for grec in genome_records:
            hits = find_hits(seq, grec.sequence, grec.name, max_mismatches=max_mm)
            for hit in hits:
                intervals = target_intervals.get(grec.name, [])
                is_on = any(s <= hit.position <= e for s, e in intervals)
                (on_hits if is_on else off_hits).append(hit)

        summary = MatchSummary(query=seq, on_target_hits=on_hits, off_target_hits=off_hits)

        if summary.on_target_score < min_ontarget:
            continue

        off_score = summary.off_target_score * offtarget_pw
        final_score = max(0.0, float(row["score"]) - off_score / max(1, len(off_hits) + 1))

        screened_rows.append({
            **row,
            "on_target_hits": len(on_hits),
            "off_target_hits": len(off_hits),
            "on_target_score": round(summary.on_target_score, 4),
            "off_target_score": round(summary.off_target_score, 4),
            "final_score": round(final_score, 4),
        })
        target_passing[target_name] = target_passing.get(target_name, 0) + 1

    screened_path = out_dir / "screened.tsv"
    _write_tsv(screened_path, screened_rows)
    summary_path = out_dir / "screen_summary.json"
    summary_path.write_text(
        json.dumps({"total_passing": len(screened_rows), "per_target": target_passing}, indent=2),
        encoding="utf-8",
    )
    write_screen_report(out_dir, target_passing)

    _print(f"  {len(screened_rows)} probes passed screening → {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# panels command
# ---------------------------------------------------------------------------

def cmd_panels(args: argparse.Namespace) -> int:
    """Assemble final probe panels from screened candidates."""
    output = resolve_output_root(args.output)

    # Prefer screened candidates if available, fall back to raw candidates
    screened_path = screen_dir(output) / "screened.tsv"
    candidates_path = probes_dir(output) / "candidates.tsv"
    source = screened_path if screened_path.exists() else candidates_path
    if not source.exists():
        _die("No screened.tsv or candidates.tsv found. Run 'genoprobe probes' first.")

    out_dir = ensure_dir(panels_dir(output))
    profile = _resolve_profile("panels", args.profile)

    config = PanelConfig(
        max_probes_per_target=args.max_probes or profile["max_probes_per_target"],
        max_panel_contiguous_complementarity=profile["max_panel_contiguous_complementarity"],
        max_panel_heterodimer_tm=profile["max_panel_heterodimer_tm"],
        panel_complementarity_penalty_weight=profile["panel_complementarity_penalty_weight"],
        panel_heterodimer_penalty_weight=profile["panel_heterodimer_penalty_weight"],
    )

    _print(f"[genoprobe panels] profile={args.profile or 'balanced'}  source={source.name}")

    with source.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        all_candidates = list(reader)

    # Group by target
    by_target: dict[str, list[dict[str, Any]]] = {}
    for row in all_candidates:
        by_target.setdefault(row["target"], []).append(row)

    from genoprobe.probes import ProbeCandidate
    final_rows: list[dict[str, Any]] = []
    panel_counts: dict[str, int] = {}

    for target_name, rows in by_target.items():
        # Reconstruct ProbeCandidate objects for panel assembly
        probe_objs: list[ProbeCandidate] = []
        for r in rows:
            try:
                probe_objs.append(ProbeCandidate(
                    sequence=r["sequence"],
                    target_name=target_name,
                    start=int(r["start"]),
                    end=int(r["end"]),
                    tm=float(r["tm"]),
                    gc=float(r["gc"]),
                    entropy=float(r["entropy"]),
                    self_comp_total=int(r["self_comp_total"]),
                    self_comp_run=int(r["self_comp_run"]),
                    hairpin_tm=float(r["hairpin_tm"]),
                    homodimer_tm=float(r["homodimer_tm"]),
                    score=float(r.get("final_score") or r.get("score") or 0),
                ))
            except (KeyError, ValueError):
                continue

        # Sort by score descending before assembly
        probe_objs.sort(key=lambda p: p.score, reverse=True)
        result = assemble_panel(probe_objs, target_name, config)
        panel_counts[target_name] = result.probe_count

        for p in result.selected_probes:
            final_rows.append({
                "target": p.target_name,
                "start": p.start,
                "end": p.end,
                "sequence": p.sequence,
                "length": p.length,
                "tm": round(p.tm, 2),
                "gc": round(p.gc, 2),
                "score": round(p.score, 4),
            })

    final_path = out_dir / "final_probes.tsv"
    _write_tsv(final_path, final_rows)
    summary_path = out_dir / "panels_summary.json"
    summary_path.write_text(
        json.dumps({"total_probes": len(final_rows), "per_target": panel_counts}, indent=2),
        encoding="utf-8",
    )
    write_panels_report(out_dir, panel_counts)

    _print(f"  Assembled {len(final_rows)} panel probes across {len(panel_counts)} targets → {out_dir}")
    return 0


# ---------------------------------------------------------------------------
# index command
# ---------------------------------------------------------------------------

def cmd_index(args: argparse.Namespace) -> int:
    """Build an optional Fulgor index for secondary off-target screening."""
    if find_fulgor() is None:
        _die("Fulgor not found on PATH. Install Fulgor before running 'genoprobe index'.")
    output = resolve_output_root(args.output)
    out_dir = ensure_dir(index_dir(output))
    fasta_paths = [Path(p) for p in args.genomes]
    _print(f"[genoprobe index] Building Fulgor index from {len(fasta_paths)} FASTA file(s)...")
    index_path = build_fulgor_index(
        fasta_paths,
        out_dir,
        kmer_length=args.kmer_length or DEFAULT_KMER_LENGTH,
        minimizer_length=args.minimizer_length or DEFAULT_MINIMIZER_LENGTH,
        threads=args.threads or DEFAULT_FULGOR_THREADS,
    )
    _print(f"  Index written → {index_path}")
    return 0


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genoprobe",
        description="Genome- and annotation-driven probe design toolkit.",
    )
    parser.add_argument("--version", action="version", version=f"genoprobe {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # ---- shared arguments --------------------------------------------------
    def _add_genomes(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--genomes", "-g", nargs="+", required=True, metavar="FASTA",
            help="One or more genome FASTA files.",
        )

    def _add_output(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--output", "-o", required=True, metavar="DIR",
            help="Output directory.",
        )

    def _add_profile(p: argparse.ArgumentParser, stage: str) -> None:
        names = get_stage_profile_names(stage)
        p.add_argument(
            "--profile", metavar="PROFILE", default=None,
            help=f"Parameter profile: {', '.join(names)}. Default: balanced.",
        )

    def _add_workers(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--workers", "-w", type=int, default=None, metavar="N",
            help=f"Worker processes (default: auto, max {DEFAULT_AUTO_WORKERS}).",
        )

    # ---- targets -----------------------------------------------------------
    p_targets = sub.add_parser("targets", help="Define probe-design target regions.")
    _add_genomes(p_targets)
    _add_output(p_targets)
    p_targets.add_argument(
        "--mode", "-m", choices=["genome", "gene"], required=True,
        help="'genome': tile full genome(s); 'gene': use annotated features.",
    )
    p_targets.add_argument(
        "--annotation", "-a", metavar="GFF_GTF",
        help="Annotation file (GFF3 or GTF). Required for --mode gene.",
    )
    p_targets.add_argument(
        "--region", metavar="SEQID:START-END",
        help="Restrict genome mode to a specific region (1-based, inclusive).",
    )
    p_targets.add_argument(
        "--feature", nargs="+", metavar="TYPE",
        help="Filter gene mode to specific feature types (e.g. gene CDS exon).",
    )
    p_targets.set_defaults(func=cmd_targets)

    # ---- probes ------------------------------------------------------------
    p_probes = sub.add_parser("probes", help="Generate probe candidates from targets.")
    _add_output(p_probes)
    _add_profile(p_probes, "probes")
    _add_workers(p_probes)
    p_probes.add_argument("--min-length", type=int, default=None, metavar="N",
        help=f"Minimum probe length (default: {DEFAULT_MIN_PROBE_LENGTH}).")
    p_probes.add_argument("--max-length", type=int, default=None, metavar="N",
        help=f"Maximum probe length (default: {DEFAULT_MAX_PROBE_LENGTH}).")
    p_probes.add_argument("--min-tm", type=float, default=None, metavar="T",
        help=f"Minimum melting temperature °C (default: {DEFAULT_MIN_TM}).")
    p_probes.add_argument("--max-tm", type=float, default=None, metavar="T",
        help=f"Maximum melting temperature °C (default: {DEFAULT_MAX_TM}).")
    p_probes.add_argument("--min-gc", type=float, default=None, metavar="PCT",
        help=f"Minimum GC %% (default: {DEFAULT_MIN_GC}).")
    p_probes.add_argument("--max-gc", type=float, default=None, metavar="PCT",
        help=f"Maximum GC %% (default: {DEFAULT_MAX_GC}).")
    p_probes.add_argument("--probe-spacing", type=int, default=None, metavar="N",
        help=f"Minimum bp spacing between adjacent probes (default: {DEFAULT_PROBE_SPACING}).")
    p_probes.add_argument("--thermo-backend", choices=["basic", "primer3"], default=None,
        help="Thermodynamic backend (default: basic).")
    p_probes.add_argument("--monovalent-mm", type=float, default=None, metavar="MM",
        help=f"Monovalent salt concentration mM (default: {DEFAULT_MONOVALENT_SALT_MM}).")
    p_probes.add_argument("--probe-conc", type=float, default=None, metavar="NM",
        help=f"Total strand concentration nM (default: {DEFAULT_PROBE_CONC_NM}).")
    p_probes.add_argument("--formamide", type=float, default=None, metavar="PCT",
        help=f"Formamide %% for Tm correction (default: {DEFAULT_FORMAMIDE_PCT}).")
    p_probes.set_defaults(func=cmd_probes)

    # ---- screen ------------------------------------------------------------
    p_screen = sub.add_parser("screen", help="Off-target screen probe candidates.")
    _add_genomes(p_screen)
    _add_output(p_screen)
    _add_profile(p_screen, "screen")
    _add_workers(p_screen)
    p_screen.add_argument("--max-mismatches", type=int, default=None, metavar="N",
        help=f"Max mismatches for off-target hit (default: {DEFAULT_MAX_MISMATCHES}).")
    p_screen.add_argument("--max-kmer-frequency", type=int, default=None, metavar="N",
        help="Discard probes where any 18-mer exceeds this genome frequency.")
    p_screen.add_argument("--index", metavar="DIR",
        help="Fulgor index directory for secondary off-target screening.")
    p_screen.set_defaults(func=cmd_screen)

    # ---- panels ------------------------------------------------------------
    p_panels = sub.add_parser("panels", help="Assemble final probe panels.")
    _add_output(p_panels)
    _add_profile(p_panels, "panels")
    p_panels.add_argument("--max-probes", type=int, default=None, metavar="N",
        help=f"Max probes per target (default: {DEFAULT_MAX_PROBES_PER_TARGET}).")
    p_panels.set_defaults(func=cmd_panels)

    # ---- index -------------------------------------------------------------
    p_index = sub.add_parser("index", help="Build optional Fulgor index for secondary off-target screening.")
    _add_genomes(p_index)
    _add_output(p_index)
    p_index.add_argument("--kmer-length", type=int, default=None, metavar="K",
        help=f"k-mer length (default: {DEFAULT_KMER_LENGTH}).")
    p_index.add_argument("--minimizer-length", type=int, default=None, metavar="M",
        help=f"Minimizer length (default: {DEFAULT_MINIMIZER_LENGTH}).")
    p_index.add_argument("--threads", type=int, default=None, metavar="N",
        help=f"Threads for Fulgor (default: {DEFAULT_FULGOR_THREADS}).")
    p_index.set_defaults(func=cmd_index)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
