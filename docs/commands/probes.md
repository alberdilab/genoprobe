# `probes` — Generate probe candidates

The `probes` command tiles target sequences with a sliding window and applies a
series of thermodynamic and sequence-complexity filters to produce a list of
candidate probes.

## Usage

```
genoprobe probes --output DIR [OPTIONS]
```

**Requires:** `targets` has already been run (reads `<output>/targets/<name>/`).

## Required arguments

| Argument | Description |
|----------|-------------|
| `--output / -o` | Output directory used in the preceding `targets` run. |

## Optional arguments

### Profile

| Argument | Default | Description |
|----------|---------|-------------|
| `--profile` | `balanced` | Parameter preset: `fast`, `balanced`, or `strict`. |

See [Profiles](../profiles.md) for the full parameter tables.

### Probe length

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-length N` | 36 | Minimum probe length in nucleotides. |
| `--max-length N` | 41 | Maximum probe length in nucleotides. |

### Melting temperature

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-tm T` | 42.0 | Minimum Tm in °C (after salt and formamide correction). |
| `--max-tm T` | 47.0 | Maximum Tm in °C. |

### GC content

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-gc PCT` | 20.0 | Minimum GC percentage. |
| `--max-gc PCT` | 75.0 | Maximum GC percentage. |

### Probe spacing

| Argument | Default | Description |
|----------|---------|-------------|
| `--probe-spacing N` | 0 | Minimum gap in bp between the end of one accepted probe and the start of the next on the same target. Set to `0` to disable. |

### Thermodynamic backend

| Argument | Default | Description |
|----------|---------|-------------|
| `--thermo-backend` | `basic` | `basic` — built-in nearest-neighbour model; `primer3` — [primer3-py](https://github.com/libnano/primer3-py) (requires optional install). |
| `--monovalent-mm MM` | 390.0 | Monovalent salt concentration in mM. |
| `--probe-conc NM` | 50.0 | Total strand concentration in nM (probe + target). |
| `--formamide PCT` | 50.0 | Formamide percentage. Tm is corrected by −0.72 °C per % formamide. |

### Parallelism

| Argument | Default | Description |
|----------|---------|-------------|
| `--workers / -w N` | auto | Number of worker processes. Defaults to `min(CPU_count, 4)`. Each genome is processed in a separate worker, so speedup scales with genome count up to `--workers`. |

## Output files

All outputs are written to `<output>/probes/<name>/`, where `<name>` matches the genome subfolder created by `targets`.

| File | Description |
|------|-------------|
| `<name>_candidates.tsv` | Tab-separated table of all candidate probes passing all filters. |
| `<name>_probes_summary.json` | JSON summary: total candidate count and per-target counts. |

### `candidates.tsv` columns

| Column | Description |
|--------|-------------|
| `target` | Target name (from `targets.fa`). |
| `start` | 0-based start position on the target sequence. |
| `end` | 0-based exclusive end position. |
| `sequence` | Probe sequence (5′→3′). |
| `length` | Probe length in nt. |
| `tm` | Calculated Tm in °C (after formamide correction). |
| `gc` | GC content in %. |
| `entropy` | Normalised Shannon entropy (0–1). |
| `self_comp_total` | Total self-complementary base count vs reverse complement. |
| `self_comp_run` | Longest contiguous self-complementary run. |
| `hairpin_tm` | Estimated hairpin Tm in °C (0 if no stable stem found). |
| `homodimer_tm` | Estimated homodimer Tm in °C (0 if self-complementarity < 4). |
| `score` | Composite quality score in \[0, 1\]; higher is better. |

## Scoring

Each candidate receives a composite score before thermodynamic penalties:

```
base_score = 0.5 × tm_score + 0.3 × gc_score + 0.2 × entropy
```

where `tm_score` and `gc_score` measure proximity to the midpoint of the allowed
Tm window and 50% GC respectively. Thermodynamic penalties (self-complementarity,
hairpin, homodimer) are applied with weight `thermo_penalty_weight` from the
profile. The final score is in [0, 1]; higher values indicate better candidates.
