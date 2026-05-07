# `screen` — Off-target screening

The `screen` command filters probe candidates by searching for off-target matches
across one or more genome sequences using a mismatch-tolerant string search. It
also supports k-mer frequency pre-filtering to quickly reject probes that match
repetitive regions.

## Usage

```
genoprobe screen --genomes FASTA [FASTA ...] --output DIR [OPTIONS]
```

**Requires:** `probes` has already been run (reads `<output>/probes/candidates.tsv`).

## Required arguments

| Argument | Description |
|----------|-------------|
| `--genomes / -g` | One or more genome FASTA files to screen against. Typically the same genome(s) used in `targets`. |
| `--output / -o` | Output directory used in the preceding pipeline stages. |

## Optional arguments

### Profile

| Argument | Default | Description |
|----------|---------|-------------|
| `--profile` | `balanced` | Parameter preset: `fast`, `balanced`, or `strict`. |

See [Profiles](../profiles.md) for profile-specific parameter values.

### Mismatch matching

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-mismatches N` | 2 | Maximum number of mismatches allowed when calling an off-target hit. |

### k-mer frequency filter

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-kmer-frequency N` | 100 | Discard any probe whose most frequent 18-mer appears more than N times in the genome. Set to omit k-mer filtering. |

### Fulgor secondary screening

| Argument | Description |
|----------|-------------|
| `--index DIR` | Path to a pre-built Fulgor index directory (see [`index`](index_cmd.md)). When provided, Fulgor is used for a secondary k-mer presence query after the primary mismatch search. |

### Parallelism

| Argument | Default | Description |
|----------|---------|-------------|
| `--workers / -w N` | auto | Number of worker processes. |

## How screening works

1. **k-mer frequency filter** (optional): A 18-mer index is built from all genome
   sequences. Any candidate whose most frequent 18-mer exceeds `--max-kmer-frequency`
   is discarded immediately.

2. **Mismatch matching**: Each surviving candidate is aligned against every genome
   sequence. Hits within the probe's own target BED interval are counted as
   *on-target*; all other hits are counted as *off-target*.

3. **On-target fraction filter**: Candidates where the on-target alignment score
   falls below `min_ontarget_fraction` (profile parameter) are discarded.

4. **Off-target score adjustment**: The candidate's `score` from the `probes` stage
   is penalised proportionally to the number and strength of off-target hits, using
   the `offtarget_penalty_weight` profile parameter.

## Output files

All outputs are written to `<output>/screen/`.

| File | Description |
|------|-------------|
| `screened.tsv` | All candidates passing screening, with added off-target columns. |
| `screen_summary.json` | JSON summary: total passing count and per-target counts. |

### `screened.tsv` additional columns

The `screened.tsv` file contains all columns from `candidates.tsv` plus:

| Column | Description |
|--------|-------------|
| `on_target_hits` | Number of hits within the target BED interval. |
| `off_target_hits` | Number of hits outside target intervals. |
| `on_target_score` | Fraction of on-target alignment strength. |
| `off_target_score` | Aggregate off-target penalty score. |
| `final_score` | Adjusted composite score after off-target penalty. |

## Notes

- Screening is performed in-memory. For very large genomes (> 5 Gbp), memory usage
  can be substantial. Consider splitting the genome into chromosomes and screening
  in batches.
- If `screened.tsv` is present when `panels` is run, it is preferred over the raw
  `candidates.tsv`.
