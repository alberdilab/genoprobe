"""Command-line interface for genoprobe."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import csv
import json
import lzma
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


def _write_tsv_xz(path: Path, rows: list[dict[str, Any]]) -> None:
    with lzma.open(path, "wt", encoding="utf-8") as fh:
        if not rows:
            return
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _open_tsv(path: Path):
    """Open a plain or xz-compressed TSV for reading."""
    if path.suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


# ---------------------------------------------------------------------------
# targets command
# ---------------------------------------------------------------------------

def _parse_targets_file(path: Path) -> tuple[list[dict[str, str]], bool]:
    """Parse a TSV/CSV batch file into a list of genome/annotation/output/group dicts.

    Accepts files with or without a header row.  With headers the recognised
    column names are ``genomes`` / ``genome``, ``annotations`` / ``annotation``,
    ``output``, and ``group`` / ``groups``.  Without headers the columns are
    positional: genome (1), annotation (2, optional), output name (3, optional),
    group (4, optional).

    Genomes sharing the same non-empty ``group`` value are combined into a single
    joint target in genome mode; the group name becomes the output subdirectory.

    Returns a tuple of (pairs, has_annotation_col).  ``has_annotation_col`` is
    True when a named annotation column is present in the header, or when at
    least one row contains a non-empty annotation value in a headerless file.
    """
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        _die(f"Targets file is empty: {path}")

    lines = text.splitlines()
    delimiter = "\t" if "\t" in lines[0] else ","

    reader = csv.reader(lines, delimiter=delimiter)
    rows = [r for r in reader if any(c.strip() for c in r)]

    if not rows:
        _die(f"No data rows found in targets file: {path}")

    first_lower = [c.strip().lower() for c in rows[0]]
    has_header = first_lower[0] in {"genomes", "genome"}

    if has_header:
        header = first_lower
        data_rows = rows[1:]
        genome_idx = next(i for i, h in enumerate(header) if h in {"genomes", "genome"})
        annotation_idx = next(
            (i for i, h in enumerate(header) if h in {"annotations", "annotation"}), None
        )
        output_idx = next((i for i, h in enumerate(header) if h in {"outputs", "output"}), None)
        group_idx = next((i for i, h in enumerate(header) if h in {"group", "groups"}), None)
        has_annotation_col = annotation_idx is not None
    else:
        data_rows = rows
        genome_idx = 0
        annotation_idx = 1
        output_idx = 2
        group_idx = 3
        has_annotation_col = False  # determined after parsing

    pairs: list[dict[str, str]] = []
    for row in data_rows:
        genome = row[genome_idx].strip() if genome_idx < len(row) else ""
        if not genome:
            continue
        annotation = (
            row[annotation_idx].strip()
            if annotation_idx is not None and annotation_idx < len(row)
            else ""
        )
        output_name = (
            row[output_idx].strip()
            if output_idx is not None and output_idx < len(row)
            else ""
        )
        group = (
            row[group_idx].strip()
            if group_idx is not None and group_idx < len(row)
            else ""
        )
        pairs.append({"genome": genome, "annotation": annotation, "output": output_name, "group": group})

    if not pairs:
        _die(f"No valid genome entries found in targets file: {path}")

    if not has_header:
        has_annotation_col = any(p["annotation"] for p in pairs)

    return pairs, has_annotation_col


def _run_targets(
    genome: str,
    annotation: str | None,
    out_dir: Path,
    name: str,
    mode: str,
    region: str | None,
    features: list[str] | None,
) -> list[dict[str, Any]]:
    """Process a single genome-annotation pair and write named output files.

    out_dir should already be <output>/targets/<name>/.
    Writes <name>_targets.fa, <name>_targets.bed, <name>_targets_summary.json.
    Returns the list of target record dicts.
    """
    target_records: list[dict[str, Any]] = []
    fasta_path = out_dir / f"{name}_targets.fa"
    bed_path = out_dir / f"{name}_targets.bed"

    with fasta_path.open("w") as fa_out, bed_path.open("w") as bed_out:

        if mode == "genome":
            fasta = load_fasta(genome)
            for seq_name in fasta.keys():
                seq = str(fasta[seq_name]).upper()
                if region:
                    try:
                        rid, coords = region.split(":")
                        rstart, rend = (int(x) for x in coords.split("-"))
                        if seq_name != rid:
                            continue
                        seq = seq[rstart - 1 : rend]
                        label = f"{seq_name}:{rstart}-{rend}"
                        bed_out.write(f"{seq_name}\t{rstart - 1}\t{rend}\t{label}\n")
                    except (ValueError, KeyError):
                        _die(f"Cannot parse --region '{region}'. Use seqid:start-end.")
                else:
                    label = seq_name
                    bed_out.write(f"{seq_name}\t0\t{len(seq)}\t{label}\n")
                fa_out.write(f">{label}\n{seq}\n")
                target_records.append({"target": label, "length": len(seq)})

        else:  # gene mode
            all_records = load_annotation(annotation)
            full_summary = summarize_annotation(all_records)
            _print(f"  Loaded {len(all_records)} annotation records. Feature types: {full_summary}")
            if features:
                feature_set = {f.lower() for f in features}
                records = [r for r in all_records if r.feature.lower() in feature_set]
                _print(f"  After --feature {features} filter: {len(records)} records.")
            else:
                records = all_records

            fasta = load_fasta(genome)
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

    summary_path = out_dir / f"{name}_targets_summary.json"
    summary_path.write_text(
        json.dumps({"name": name, "mode": mode, "target_count": len(target_records)}, indent=2),
        encoding="utf-8",
    )

    return target_records


def _run_joint_targets(
    members: list[tuple[str, str | None]],
    out_dir: Path,
    name: str,
) -> list[dict[str, Any]]:
    """Concatenate multiple genomes into one joint target (genome mode only).

    Sequence headers are tagged as ``<genome_stem>|<seqid>`` so that sequences
    from different constituent genomes with identical contig names remain
    distinguishable throughout downstream stages.
    """
    target_records: list[dict[str, Any]] = []
    fasta_path = out_dir / f"{name}_targets.fa"
    bed_path = out_dir / f"{name}_targets.bed"

    with fasta_path.open("w") as fa_out, bed_path.open("w") as bed_out:
        for genome_path, _ in members:
            genome_stem = Path(genome_path).stem
            fasta = load_fasta(genome_path)
            for seq_name in fasta.keys():
                seq = str(fasta[seq_name]).upper()
                tagged = f"{genome_stem}|{seq_name}"
                fa_out.write(f">{tagged}\n{seq}\n")
                bed_out.write(f"{tagged}\t0\t{len(seq)}\t{tagged}\n")
                target_records.append({"target": tagged, "length": len(seq)})

    (out_dir / f"{name}_targets_summary.json").write_text(
        json.dumps({
            "name": name,
            "mode": "genome",
            "joint": True,
            "members": [Path(g).stem for g, _ in members],
            "target_count": len(target_records),
        }, indent=2),
        encoding="utf-8",
    )
    return target_records


def cmd_targets(args: argparse.Namespace) -> int:
    """Extract or define target regions and write FASTA + BED."""
    base_output = resolve_output_root(args.output)
    explicit_mode = args.mode  # None when not provided on the command line

    targets_file = getattr(args, "file", None)

    # Build a unified list of (genome_path, annotation_path, name) entries.
    if targets_file:
        if args.genomes:
            _die("--file and --genomes are mutually exclusive.")
        pairs, has_annotation_col = _parse_targets_file(Path(targets_file))
        entries = [
            (p["genome"], p["annotation"] or None, p["output"] or Path(p["genome"]).stem, p.get("group", ""))
            for p in pairs
        ]
        if explicit_mode is None:
            mode = "gene" if has_annotation_col else "genome"
        else:
            mode = explicit_mode.lower()
        _print(f"[genoprobe targets] mode={mode}  pairs={len(entries)}  (batch)")
    else:
        if not args.genomes:
            _die("Either --genomes or --file is required.")
        if explicit_mode is None:
            mode = "gene" if args.annotation else "genome"
        else:
            mode = explicit_mode.lower()
        if mode == "gene" and not args.annotation:
            _die("--annotation is required in gene mode.")
        entries = [
            (g, getattr(args, "annotation", None), Path(g).stem, "")
            for g in args.genomes
        ]
        _print(f"[genoprobe targets] mode={mode}  genomes={len(entries)}")

    base_targets = targets_dir(base_output)

    # Separate joint-group entries (same group name → shared output) from solo entries.
    joint_groups: dict[str, list[tuple[str, str | None]]] = {}
    solo_entries: list[tuple[str, str | None, str]] = []
    for genome_path, annotation_path, name, group in entries:
        if group:
            if group not in joint_groups:
                joint_groups[group] = []
            joint_groups[group].append((genome_path, annotation_path))
        else:
            solo_entries.append((genome_path, annotation_path, name))

    for group_name, members in joint_groups.items():
        out_dir = ensure_dir(base_targets / group_name)
        sentinel = out_dir / f"{group_name}_targets.fa"
        if not args.overwrite and sentinel.exists():
            _print(f"  [{group_name}] Skipping (outputs exist; use --overwrite to redo).")
            continue
        _print(f"  [{group_name}] joint target  members={[Path(g).stem for g, _ in members]}")
        records = _run_joint_targets(members=members, out_dir=out_dir, name=group_name)
        _print(f"    Wrote {len(records)} targets → {out_dir}")

    for genome_path, annotation_path, name in solo_entries:
        if mode == "gene" and not annotation_path:
            _die(
                f"--mode gene requires an annotation for genome '{genome_path}'. "
                "Add an 'annotations' column to your targets file."
            )
        out_dir = ensure_dir(base_targets / name)
        sentinel = out_dir / f"{name}_targets.fa"
        if not args.overwrite and sentinel.exists():
            _print(f"  [{name}] Skipping (outputs exist; use --overwrite to redo).")
            continue
        _print(f"  [{name}] genome={genome_path}  annotation={annotation_path or '—'}")

        records = _run_targets(
            genome=genome_path,
            annotation=annotation_path,
            out_dir=out_dir,
            name=name,
            mode=mode,
            region=getattr(args, "region", None),
            features=getattr(args, "feature", None),
        )
        _print(f"    Wrote {len(records)} targets → {out_dir}")

    return 0


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _read_project_mode(base_targets: Path, names: list[str]) -> str:
    """Read analysis mode from the first available targets summary JSON.

    Falls back to 'genome' when the file is absent or lacks a mode field.
    """
    for name in names:
        summary_path = base_targets / name / f"{name}_targets_summary.json"
        if summary_path.exists():
            try:
                data = json.loads(summary_path.read_text(encoding="utf-8"))
                return data.get("mode", "genome")
            except (json.JSONDecodeError, OSError):
                pass
    return "genome"


def _discover_names(stage_dir: Path, stage_label: str) -> list[str]:
    """Return sorted genome names from per-name subdirs within a stage directory."""
    if not stage_dir.exists():
        _die(f"{stage_label}/ directory not found: {stage_dir}. Run the preceding stage first.")
    names = sorted(d.name for d in stage_dir.iterdir() if d.is_dir())
    if not names:
        _die(f"No genome subfolders found in {stage_dir}. Run the preceding stage first.")
    return names


# ---------------------------------------------------------------------------
# probes command
# ---------------------------------------------------------------------------

def _process_genome_probes(
    name: str,
    base_targets: Path,
    base_probes: Path,
    config: ProbeDesignConfig,
) -> tuple[str, list[str], int]:
    """Process one genome's probe candidates. Runs in a worker process.

    Returns (name, log_lines, total_candidates) so the caller can flush
    output atomically and avoid interleaving across parallel workers.
    """
    from pyfaidx import Fasta

    log: list[str] = []
    targets_path = base_targets / name / f"{name}_targets.fa"

    if not targets_path.exists():
        log.append(f"  [{name}] WARNING: {targets_path.name} not found, skipping.")
        return name, log, 0

    out_dir = ensure_dir(base_probes / name)
    fasta = Fasta(str(targets_path), as_raw=True, sequence_always_upper=True)

    all_rows: list[dict[str, Any]] = []
    target_counts: dict[str, int] = {}

    for target_name in fasta.keys():
        seq = str(fasta[target_name])
        candidates, stats = generate_candidates(seq, target_name, config)
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

    _write_tsv_xz(out_dir / f"{name}_candidates.tsv.xz", all_rows)
    (out_dir / f"{name}_probes_summary.json").write_text(
        json.dumps({"name": name, "total_candidates": len(all_rows), "per_target": target_counts}, indent=2),
        encoding="utf-8",
    )
    write_probes_report(out_dir, target_counts)
    log.append(f"  [{name}] Generated {len(all_rows)} candidates across {len(target_counts)} targets → {out_dir}")

    return name, log, len(all_rows)


def cmd_probes(args: argparse.Namespace) -> int:
    """Generate probe candidates from target sequences."""
    output = resolve_output_root(args.output)
    base_targets = targets_dir(output)
    base_probes = probes_dir(output)

    names = _discover_names(base_targets, "targets")
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
    _print(f"[genoprobe probes] profile={args.profile or 'balanced'}  workers={workers}  genomes={len(names)}")

    names_to_run: list[str] = []
    for name in names:
        sentinel = base_probes / name / f"{name}_candidates.tsv.xz"
        if not args.overwrite and sentinel.exists():
            _print(f"  [{name}] Skipping (outputs exist; use --overwrite to redo).")
        else:
            names_to_run.append(name)

    if workers <= 1 or len(names_to_run) <= 1:
        for name in names_to_run:
            _, log_lines, _ = _process_genome_probes(name, base_targets, base_probes, config)
            for line in log_lines:
                _print(line)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_process_genome_probes, name, base_targets, base_probes, config): name
                for name in names_to_run
            }
            for future in as_completed(futures):
                _, log_lines, _ = future.result()
                for line in log_lines:
                    _print(line)

    return 0


# ---------------------------------------------------------------------------
# screen command
# ---------------------------------------------------------------------------

# Module-level globals loaded once per worker process via _init_screen_worker.
_g_kmer_index: dict[str, int] = {}
_g_project_records: dict[str, list] = {}
_g_external_records: list = []


def _init_screen_worker(
    project_paths: dict[str, str],
    external_paths: list[str],
    kmer_paths: list[str],
    use_kmer: bool,
) -> None:
    """ProcessPoolExecutor initializer — loads shared read-only data from disk once per worker.

    Accepts file paths rather than pre-loaded data to avoid pickling large dicts through
    IPC, which causes BrokenProcessPool via OOM when there are many genomes or workers.
    """
    global _g_kmer_index, _g_project_records, _g_external_records
    for pname, path_str in project_paths.items():
        p = Path(path_str)
        if p.exists():
            _g_project_records[pname] = load_genomes([p])
    if external_paths:
        _g_external_records = load_genomes([Path(p) for p in external_paths])
    if use_kmer:
        project_path_set = set(project_paths.values())
        extra_paths = [Path(p) for p in kmer_paths if p not in project_path_set]
        all_recs: list = [r for recs in _g_project_records.values() for r in recs]
        if extra_paths:
            all_recs.extend(load_genomes(extra_paths))
        all_recs.extend(_g_external_records)
        _g_kmer_index = build_kmer_index(all_recs, k=18)


def _screen_genome(
    name: str,
    base_targets: Path,
    base_probes: Path,
    base_screen: Path,
    mode: str,
    max_mm: int,
    max_kmer: int | None,
    offtarget_pw: float,
    min_ontarget: float,
) -> tuple[str, list[str], int]:
    """Screen one focal genome's probe candidates. Runs in a worker process.

    Returns (name, log_lines, total_passing) so the caller can flush output
    atomically and avoid interleaving across parallel workers.
    """
    kmer_index = _g_kmer_index
    project_records = _g_project_records
    external_records = _g_external_records
    log: list[str] = []

    candidates_path = base_probes / name / f"{name}_candidates.tsv.xz"
    if not candidates_path.exists():
        log.append(f"  [{name}] WARNING: {candidates_path.name} not found, skipping.")
        return name, log, 0

    targets_bed = base_targets / name / f"{name}_targets.bed"
    out_dir = ensure_dir(base_screen / name)

    with _open_tsv(candidates_path) as fh:
        candidates = list(csv.DictReader(fh, delimiter="\t"))

    target_intervals: dict[str, list[tuple[int, int]]] = {}
    if mode == "genome" and targets_bed.exists():
        for line in targets_bed.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                seqid, start, end = parts[0], int(parts[1]), int(parts[2])
                target_intervals.setdefault(seqid, []).append((start, end))

    focal_fa = base_targets / name / f"{name}_targets.fa"
    focal_recs = project_records.get(name) or (load_genomes([focal_fa]) if focal_fa.exists() else [])

    other_recs: list = []
    for pname, recs in project_records.items():
        if pname != name:
            other_recs.extend(recs)

    genome_records = focal_recs + other_recs + external_records

    screened_rows: list[dict[str, Any]] = []
    target_passing: dict[str, int] = {}

    for row in candidates:
        seq = row["sequence"]
        target_name = row["target"]

        if max_kmer is not None and kmer_index:
            if query_kmer_frequency(kmer_index, seq) > max_kmer:
                continue

        on_hits = []
        off_hits = []
        for grec in genome_records:
            hits = find_hits(seq, grec.sequence, grec.name, max_mismatches=max_mm)
            for hit in hits:
                if mode == "gene":
                    is_on = (grec.name == target_name)
                else:
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

    _write_tsv(out_dir / f"{name}_screened.tsv", screened_rows)
    (out_dir / f"{name}_screen_summary.json").write_text(
        json.dumps({"name": name, "total_passing": len(screened_rows), "per_target": target_passing}, indent=2),
        encoding="utf-8",
    )
    write_screen_report(out_dir, target_passing)
    log.append(f"  [{name}] {len(screened_rows)} probes passed screening → {out_dir}")

    return name, log, len(screened_rows)


def cmd_screen(args: argparse.Namespace) -> int:
    """Off-target screening via mismatch matching (+ optional Fulgor index)."""
    output = resolve_output_root(args.output)
    base_probes = probes_dir(output)
    base_targets = targets_dir(output)
    base_screen = screen_dir(output)

    names = _discover_names(base_probes, "probes")
    profile = _resolve_profile("screen", args.profile)

    max_mm = args.max_mismatches if args.max_mismatches is not None else profile["max_mismatches"]
    max_kmer = args.max_kmer_frequency if args.max_kmer_frequency is not None else profile["max_kmer_frequency"]
    offtarget_pw = profile["offtarget_penalty_weight"]
    min_ontarget = profile["min_ontarget_fraction"]

    # Mode is recorded in the targets summary written by genoprobe targets.
    mode = _read_project_mode(base_targets, names)

    # Determine which project genomes participate in off-target screening.
    all_project_names = _discover_names(base_targets, "targets")
    include_set: set[str] | None = None
    exclude_set: set[str] = set()
    if getattr(args, "include", None):
        include_set = {n.strip() for n in args.include.split(",")}
    if getattr(args, "exclude", None):
        exclude_set = {n.strip() for n in args.exclude.split(",")}
    if include_set is not None:
        screened_project_names = set(include_set) & set(all_project_names)
    else:
        screened_project_names = set(all_project_names) - exclude_set

    # Collect file paths for project genomes — sequences are loaded on-demand
    # (in workers or in the main process for single-worker runs) to avoid
    # holding large data in the parent before workers are spawned.
    project_paths: dict[str, str] = {}
    for pname in screened_project_names:
        fa = base_targets / pname / f"{pname}_targets.fa"
        if fa.exists():
            project_paths[pname] = str(fa)

    # Load external genomes now so we can report the count; paths are also
    # forwarded to workers so they can reload from the OS page-cached files.
    external_records: list = []
    external_paths: list[str] = []
    if getattr(args, "external", None):
        external_paths = [str(p) for p in args.external]
        external_records = load_genomes(args.external)
        _print(f"  Loaded {len(external_records)} external genome sequences.")

    workers = resolve_worker_count(args.workers)
    _print(
        f"[genoprobe screen] mode={mode}  profile={args.profile or 'balanced'}  "
        f"max_mismatches={max_mm}  project_genomes={len(project_paths)}  "
        f"external_sequences={len(external_records)}  workers={workers}"
    )

    # Collect all target FASTA paths used for k-mer index building.
    all_kmer_paths: list[str] = list(project_paths.values())
    for n in names:
        if n not in screened_project_names:
            fa = base_targets / n / f"{n}_targets.fa"
            if fa.exists():
                all_kmer_paths.append(str(fa))

    names_to_run: list[str] = []
    for name in names:
        sentinel = base_screen / name / f"{name}_screened.tsv"
        if not args.overwrite and sentinel.exists():
            _print(f"  [{name}] Skipping (outputs exist; use --overwrite to redo).")
        else:
            names_to_run.append(name)

    global _g_kmer_index, _g_project_records, _g_external_records

    if workers <= 1 or len(names_to_run) <= 1:
        # Single-worker path: build all data in the main process.
        project_records: dict[str, list] = {
            pname: load_genomes([Path(p)]) for pname, p in project_paths.items()
        }
        kmer_index: dict[str, int] = {}
        if max_kmer is not None:
            _print("  Building k-mer index (k=18)...")
            kmer_index = build_kmer_index(
                load_genomes([Path(p) for p in all_kmer_paths]), k=18
            )
            if external_records:
                ext_index = build_kmer_index(external_records, k=18)
                for kmer, count in ext_index.items():
                    kmer_index[kmer] = kmer_index.get(kmer, 0) + count
        _g_kmer_index = kmer_index
        _g_project_records = project_records
        _g_external_records = external_records
        for name in names_to_run:
            _, log_lines, _ = _screen_genome(
                name, base_targets, base_probes, base_screen,
                mode, max_mm, max_kmer, offtarget_pw, min_ontarget,
            )
            for line in log_lines:
                _print(line)
    else:
        # Multi-worker path: build all shared data in the main process once,
        # then fork workers so they inherit a single COW copy.  This avoids
        # both large-object IPC pickling and per-worker rebuilds, which caused
        # OOM → BrokenProcessPool when many genomes and workers were combined.
        project_records_mp: dict[str, list] = {
            pname: load_genomes([Path(p)])
            for pname, p in project_paths.items()
            if Path(p).exists()
        }
        kmer_index_mp: dict[str, int] = {}
        if max_kmer is not None:
            _print("  Building k-mer index (k=18)...")
            kmer_index_mp = build_kmer_index(
                load_genomes([Path(p) for p in all_kmer_paths]), k=18
            )
            if external_records:
                ext_index = build_kmer_index(external_records, k=18)
                for kmer, count in ext_index.items():
                    kmer_index_mp[kmer] = kmer_index_mp.get(kmer, 0) + count
        _g_kmer_index = kmer_index_mp
        _g_project_records = project_records_mp
        _g_external_records = external_records

        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Fork on Unix: workers inherit pre-built globals via COW (one physical
        # copy shared across all workers).  Windows lacks fork; fall back to
        # spawn with the initializer-based loading path.
        pool_kwargs: dict[str, Any] = {"max_workers": workers}
        if sys.platform != "win32":
            pool_kwargs["mp_context"] = mp.get_context("fork")
        else:
            pool_kwargs["initializer"] = _init_screen_worker
            pool_kwargs["initargs"] = (project_paths, external_paths, all_kmer_paths, max_kmer is not None)

        with ProcessPoolExecutor(**pool_kwargs) as pool:
            futures = {
                pool.submit(
                    _screen_genome,
                    name, base_targets, base_probes, base_screen,
                    mode, max_mm, max_kmer, offtarget_pw, min_ontarget,
                ): name
                for name in names_to_run
            }
            for future in as_completed(futures):
                _, log_lines, _ = future.result()
                for line in log_lines:
                    _print(line)

    return 0


# ---------------------------------------------------------------------------
# panels command
# ---------------------------------------------------------------------------

def cmd_panels(args: argparse.Namespace) -> int:
    """Assemble final probe panels from screened candidates."""
    output = resolve_output_root(args.output)
    base_probes = probes_dir(output)

    names = _discover_names(base_probes, "probes")
    profile = _resolve_profile("panels", args.profile)

    config = PanelConfig(
        max_probes_per_target=args.max_probes or profile["max_probes_per_target"],
        max_panel_contiguous_complementarity=profile["max_panel_contiguous_complementarity"],
        max_panel_heterodimer_tm=profile["max_panel_heterodimer_tm"],
        panel_complementarity_penalty_weight=profile["panel_complementarity_penalty_weight"],
        panel_heterodimer_penalty_weight=profile["panel_heterodimer_penalty_weight"],
    )

    _print(f"[genoprobe panels] profile={args.profile or 'balanced'}  genomes={len(names)}")

    from genoprobe.probes import ProbeCandidate

    for name in names:
        # Prefer screened candidates; fall back to raw candidates
        screened_path = screen_dir(output) / name / f"{name}_screened.tsv"
        candidates_path = base_probes / name / f"{name}_candidates.tsv.xz"

        if screened_path.exists():
            source_path = screened_path
            source_label = "screened"
        elif candidates_path.exists():
            source_path = candidates_path
            source_label = "candidates"
        else:
            _print(f"  [{name}] WARNING: no candidates found, skipping.")
            continue

        out_dir = ensure_dir(panels_dir(output) / name)
        sentinel = out_dir / f"{name}_final_probes.tsv"
        if not args.overwrite and sentinel.exists():
            _print(f"  [{name}] Skipping (outputs exist; use --overwrite to redo).")
            continue
        _print(f"  [{name}] source={source_label}")

        with _open_tsv(source_path) as fh:
            all_candidates = list(csv.DictReader(fh, delimiter="\t"))

        by_target: dict[str, list[dict[str, Any]]] = {}
        for row in all_candidates:
            by_target.setdefault(row["target"], []).append(row)

        final_rows: list[dict[str, Any]] = []
        panel_counts: dict[str, int] = {}

        for target_name, rows in by_target.items():
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

        _write_tsv(out_dir / f"{name}_final_probes.tsv", final_rows)
        (out_dir / f"{name}_panels_summary.json").write_text(
            json.dumps({"name": name, "total_probes": len(final_rows), "per_target": panel_counts}, indent=2),
            encoding="utf-8",
        )
        write_panels_report(out_dir, panel_counts)
        _print(f"  [{name}] Assembled {len(final_rows)} panel probes → {out_dir}")

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
    sentinel = out_dir / "index.fur"
    if not args.overwrite and sentinel.exists():
        _print(f"[genoprobe index] Skipping (index exists at {sentinel}; use --overwrite to redo).")
        return 0
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

    def _add_overwrite(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--overwrite", action="store_true", default=False,
            help="Overwrite existing output files (default: skip if outputs already exist).",
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
    _add_output(p_targets)
    p_targets.add_argument(
        "--genomes", "-g", nargs="+", default=None, metavar="FASTA",
        help="One or more genome FASTA files. Mutually exclusive with --file.",
    )
    p_targets.add_argument(
        "--file", "-f", metavar="TSV_CSV",
        help=(
            "TSV/CSV file with one genome per row. "
            "Recognised column headers: 'genomes', 'annotations', 'output' (subfolder name), "
            "'group' (joint-target group name). "
            "Genomes sharing the same group name are combined into a single joint target "
            "(genome mode only); the group name becomes the output subdirectory. "
            "Headerless files use columns in that order. "
            "Mutually exclusive with --genomes."
        ),
    )
    p_targets.add_argument(
        "--mode", "-m", choices=["genome", "gene"], default=None,
        help=(
            "'genome': tile full genome(s); 'gene': use annotated features. "
            "Auto-detected from --annotation or file contents when omitted "
            "(gene if annotations are present, genome otherwise)."
        ),
    )
    p_targets.add_argument(
        "--annotation", "-a", metavar="GFF_GTF",
        help="Annotation file (GFF3 or GTF). Required for --mode gene (single-genome mode).",
    )
    p_targets.add_argument(
        "--region", metavar="SEQID:START-END",
        help="Restrict genome mode to a specific region (1-based, inclusive).",
    )
    p_targets.add_argument(
        "--feature", nargs="+", metavar="TYPE",
        help="Filter gene mode to specific feature types (e.g. gene CDS exon).",
    )
    _add_overwrite(p_targets)
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
    _add_overwrite(p_probes)
    p_probes.set_defaults(func=cmd_probes)

    # ---- screen ------------------------------------------------------------
    p_screen = sub.add_parser("screen", help="Off-target screen probe candidates.")
    _add_output(p_screen)
    _add_profile(p_screen, "screen")
    _add_workers(p_screen)
    p_screen.add_argument(
        "--external", "-x", nargs="+", default=None, metavar="FASTA",
        help=(
            "External genome FASTA files to screen against in addition to project genomes. "
            "Use for genomes outside the current project."
        ),
    )
    p_screen.add_argument(
        "--include", "-i", metavar="NAMES",
        help=(
            "Comma-separated project genome names to include in off-target screening. "
            "When set, only these project genomes are used (always from within the project)."
        ),
    )
    p_screen.add_argument(
        "--exclude", "-e", metavar="NAMES",
        help=(
            "Comma-separated project genome names to skip during off-target screening. "
            "Mutually exclusive with --include."
        ),
    )
    p_screen.add_argument("--max-mismatches", type=int, default=None, metavar="N",
        help=f"Max mismatches for off-target hit (default: {DEFAULT_MAX_MISMATCHES}).")
    p_screen.add_argument("--max-kmer-frequency", type=int, default=None, metavar="N",
        help="Discard probes where any 18-mer exceeds this genome frequency.")
    p_screen.add_argument("--index", metavar="DIR",
        help="Fulgor index directory for secondary off-target screening.")
    _add_overwrite(p_screen)
    p_screen.set_defaults(func=cmd_screen)

    # ---- panels ------------------------------------------------------------
    p_panels = sub.add_parser("panels", help="Assemble final probe panels.")
    _add_output(p_panels)
    _add_profile(p_panels, "panels")
    p_panels.add_argument("--max-probes", type=int, default=None, metavar="N",
        help=f"Max probes per target (default: {DEFAULT_MAX_PROBES_PER_TARGET}).")
    _add_overwrite(p_panels)
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
    _add_overwrite(p_index)
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
