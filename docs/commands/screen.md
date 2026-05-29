# `screen` — Off-target screening

The `screen` command filters probe candidates by searching for off-target matches
across genome sequences using a mismatch-tolerant string search. It automatically
screens every probe against all other project genomes (those processed by `targets`)
and optionally against external genomes supplied with `--external`. It also supports
k-mer frequency pre-filtering to quickly reject probes that match repetitive regions.

The analysis mode (`gene` or `genome`) is read automatically from the
`<output>/targets/<name>/<name>_targets_summary.json` files written by
`genoprobe targets` — no manual mode flag is needed at the screening stage.

## Usage

```
genoprobe screen --output DIR [OPTIONS]
```

**Requires:** `probes` has already been run (reads `<output>/probes/<name>/`).

## Required arguments

| Argument | Description |
|----------|-------------|
| `--output / -o` | Output directory used in the preceding pipeline stages. |

## Optional arguments

### Project genome selection

By default all project genomes found in `<output>/targets/` are used for
cross-genome off-target screening. Two mutually exclusive flags narrow this set:

| Argument | Description |
|----------|-------------|
| `--include / -i NAMES` | Comma-separated genome names to include. Only these project genomes are screened against (the focal genome is always included regardless). |
| `--exclude / -e NAMES` | Comma-separated genome names to skip. All other project genomes are still used. |

### External genomes

| Argument | Description |
|----------|-------------|
| `--external / -x FASTA [FASTA ...]` | One or more FASTA files from outside the current project. All hits against these sequences are counted as off-target. |

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
| `--max-kmer-frequency N` | 100 | Discard any probe whose most frequent 18-mer appears more than N times across all screening sequences. |

### Fulgor secondary screening

| Argument | Description |
|----------|-------------|
| `--index DIR` | Path to a pre-built Fulgor index directory (see [`index`](index_cmd.md)). When provided, Fulgor is used for a secondary k-mer presence query after the primary mismatch search. |

### Parallelism

| Argument | Default | Description |
|----------|---------|-------------|
| `--workers / -w N` | auto | Number of worker processes. |

## How screening works

### Genome sources

For each focal genome being screened, the sequences searched are:

1. **Focal genome** — always included. The `<name>_targets.fa` written by `genoprobe targets` is used directly, so no separate genome FASTA is needed.
2. **Other project genomes** — all genomes in `<output>/targets/` except the focal one, filtered by `--include` / `--exclude`. Their `targets.fa` files are loaded as screening sequences.
3. **External genomes** — any FASTA files supplied via `--external`. All hits here count as off-target.

In **gene mode**, `targets.fa` holds the extracted gene sequences, so cross-genome
screening is gene-vs-gene. In **genome mode**, `targets.fa` holds the full genome
sequences (all contigs/chromosomes), so cross-genome screening is genome-vs-genome.

### Per-probe classification

1. **k-mer frequency filter** (optional): An 18-mer index is built from all screening
   sequences. Any candidate whose most frequent 18-mer exceeds `--max-kmer-frequency`
   is discarded immediately.

2. **Mismatch matching**: Each surviving candidate is aligned against every screening
   sequence. Hit classification depends on mode:
   - **Gene mode**: a hit is *on-target* only when the sequence name exactly matches
     the probe's target gene. Hits in any other gene sequence — including other genes
     in the same genome — count as off-target.
   - **Genome mode**: a hit is *on-target* when its position falls within the target
     BED intervals of the focal genome. All other hits are off-target.

3. **On-target fraction filter**: Candidates where the on-target alignment score
   falls below `min_ontarget_fraction` (profile parameter) are discarded.

4. **Off-target score adjustment**: The candidate's `score` from the `probes` stage
   is penalised proportionally to the number and strength of off-target hits, using
   the `offtarget_penalty_weight` profile parameter.

## Output files

All outputs are written to `<output>/screen/<name>/`.

| File | Description |
|------|-------------|
| `<name>_screened.tsv` | All candidates passing screening, with added off-target columns. |
| `<name>_screen_summary.json` | JSON summary: total passing count and per-target counts. |

### `<name>_screened.tsv` additional columns

The file contains all columns from `<name>_candidates.tsv` plus:

| Column | Description |
|--------|-------------|
| `on_target_hits` | Number of hits classified as on-target. |
| `off_target_hits` | Number of hits classified as off-target. |
| `on_target_score` | Fraction of on-target alignment strength. |
| `off_target_score` | Aggregate off-target penalty score. |
| `final_score` | Adjusted composite score after off-target penalty. |

## Notes

- No FASTA files need to be supplied at the `screen` stage unless you want to include
  external genomes via `--external`. All project genomes are loaded automatically from
  the `targets/` directory.
- Screening is performed in-memory. For very large genomes (> 5 Gbp total across all
  project genomes), memory usage can be substantial.
- If `<name>_screened.tsv` is present when `panels` is run, it is preferred over the
  raw `<name>_candidates.tsv`.
