# Quick start

This guide walks through a complete probe-design run from raw genome sequences to a
final panel, using the four-stage pipeline.

## Overview

```
genome.fa  annotation.gff
     │
     ▼
┌─────────┐   targets.fa / targets.bed
│ targets │──────────────────────────────┐
└─────────┘                              │
                                         ▼
                                  ┌──────────┐   candidates.tsv
                                  │  probes  │──────────────────┐
                                  └──────────┘                  │
                                                                 ▼
                                                         ┌────────────┐   screened.tsv
                                                         │   screen   │──────────────┐
                                                         └────────────┘              │
                                                                                     ▼
                                                                             ┌────────────┐
                                                                             │   panels   │
                                                                             └────────────┘
                                                                              final_probes.tsv
```

## Example: gene-mode probe design

### 1. Define target regions

Extract gene sequences from an annotated genome:

```bash
genoprobe targets \
    --genomes genome.fa \
    --annotation genes.gff \
    --mode gene \
    --feature gene \
    --output my_run/
```

Output files written to `my_run/targets/`:

| File | Description |
|------|-------------|
| `targets.fa` | Target sequences in FASTA format |
| `targets.bed` | Genomic coordinates of each target |
| `targets_summary.json` | Target count summary |

### 2. Generate probe candidates

Tile each target sequence with candidate probes:

```bash
genoprobe probes \
    --output my_run/ \
    --profile balanced
```

Output files written to `my_run/probes/`:

| File | Description |
|------|-------------|
| `candidates.tsv` | All probes passing thermodynamic filters |
| `probes_summary.json` | Candidate counts per target |

### 3. Screen for off-target binding

Filter candidates against the full genome:

```bash
genoprobe screen \
    --genomes genome.fa \
    --output my_run/ \
    --profile balanced
```

Output files written to `my_run/screen/`:

| File | Description |
|------|-------------|
| `screened.tsv` | Candidates with on/off-target hit counts and scores |
| `screen_summary.json` | Passing counts per target |

### 4. Assemble final panels

Select the best non-cross-reactive probes per target:

```bash
genoprobe panels \
    --output my_run/ \
    --profile balanced
```

Output files written to `my_run/panels/`:

| File | Description |
|------|-------------|
| `final_probes.tsv` | Final probe sequences, ready for synthesis |
| `panels_summary.json` | Probe counts per target |

---

## Example: genome-tiling mode

To tile an entire genome (e.g. for capture-seq), use `--mode genome`:

```bash
genoprobe targets \
    --genomes genome.fa \
    --mode genome \
    --output tiling_run/

genoprobe probes --output tiling_run/
genoprobe screen --genomes genome.fa --output tiling_run/
genoprobe panels --output tiling_run/
```

To restrict tiling to a single region:

```bash
genoprobe targets \
    --genomes genome.fa \
    --mode genome \
    --region chr1:1000000-2000000 \
    --output region_run/
```

---

## Choosing a profile

All commands accept a `--profile` flag with three presets:

| Profile | When to use |
|---------|-------------|
| `fast` | Quick exploration — relaxed filters, no thermodynamic penalty |
| `balanced` | **Default** — good selectivity for most applications |
| `strict` | Demanding assays — tightest Tm window and cross-reactivity filters |

See [Profiles](profiles.md) for a full parameter breakdown.
