# `targets` — Define probe-design target regions

The `targets` command extracts or defines genomic regions to be used as probe
design templates. It writes a FASTA file of target sequences and a BED file of
their genomic coordinates.

## Usage

**Single genome:**
```
genoprobe targets --genomes FASTA [FASTA ...] --output DIR [--mode {genome,gene}] [OPTIONS]
```

**Batch (multiple genome-annotation pairs):**
```
genoprobe targets --file FILE --output DIR [--mode {genome,gene}] [OPTIONS]
```

## Required arguments

| Argument | Description |
|----------|-------------|
| `--genomes / -g` | One or more genome FASTA files. Mutually exclusive with `--file`. |
| `--file / -f` | TSV/CSV batch file of genome-annotation pairs. Mutually exclusive with `--genomes`. |
| `--output / -o` | Output directory. Sub-directories are created automatically. |

## Mode selection (`--mode / -m`)

`--mode` is **optional**. When omitted, it is inferred automatically:

| Situation | Auto-detected mode |
|-----------|-------------------|
| `--annotation` is provided (single-genome) | `gene` |
| Batch file header contains an `annotation`/`annotations` column with at least one value | `gene` |
| Neither of the above | `genome` |

An explicit `--mode genome` or `--mode gene` always overrides auto-detection.

| Value | Behaviour |
|-------|-----------|
| `genome` | Tile the full genome sequence(s). |
| `gene` | Extract sequences for annotated features from a GFF3/GTF file. |

## Mode: `genome`

Tiles every sequence in the input FASTA file(s). Use `--region` to restrict to a
single genomic interval.

```bash
# mode auto-detected as genome (no --annotation given)
genoprobe targets \
    --genomes genome.fa \
    --output run/
```

### `--region SEQID:START-END`

Restrict genome-mode tiling to one region. Coordinates are **1-based, inclusive**.

```bash
genoprobe targets \
    --genomes genome.fa \
    --region chr1:500000-1500000 \
    --output run/
```

## Mode: `gene`

Extracts sequences for annotated features from a GFF3 or GTF file.
`--annotation` is required in single-genome mode; in batch mode each row supplies
its own annotation path.

```bash
# mode auto-detected as gene because --annotation is provided
genoprobe targets \
    --genomes genome.fa \
    --annotation genes.gff \
    --output run/
```

### `--annotation / -a GFF_GTF`

Path to the annotation file (single-genome mode). Both GFF3 and GTF formats are
supported.

### `--feature TYPE [TYPE ...]`

Filter to specific feature types. If omitted, all features in the annotation are
used.

```bash
# Use only gene-level features
genoprobe targets --mode gene --annotation genes.gff --feature gene ...

# Use CDS and exon features
genoprobe targets --mode gene --annotation genes.gff --feature CDS exon ...
```

## Batch mode: `--file`

Use `--file / -f` to process multiple genome-annotation pairs in one command.
Each pair is processed independently and written to its own subfolder under
`<output>/`.

### File format

The file can be TSV (tab-separated) or CSV (comma-separated); the delimiter is
detected automatically.

**With column headers** — recognised header names:

| Column | Header names accepted | Required |
|--------|----------------------|----------|
| Genome FASTA | `genomes`, `genome` | Yes |
| Annotation file | `annotations`, `annotation` | No |
| Output subfolder | `outputs`, `output` | No |

```
genomes	annotations	output
/data/org1.fa	/data/org1.gff	org1
/data/org2.fa	/data/org2.gff	org2
/data/org3.fa		org3
```

**Without headers** — columns are positional:

1. Genome FASTA (required)
2. Annotation file (optional)
3. Output subfolder name (optional)

```
/data/org1.fa	/data/org1.gff	org1
/data/org2.fa	/data/org2.gff
/data/org3.fa
```

If the `output` column is absent or empty for a row, the subfolder name defaults
to the **genome filename stem** (e.g. `org1.fa` → `org1`).

### Example

```bash
# mode auto-detected as gene because the file has an annotations column
genoprobe targets \
    --file organisms.tsv \
    --output results/
```

Produces:
```
results/
  targets/
    org1/
      org1_targets.fa
      org1_targets.bed
      org1_targets_summary.json
    org2/
      org2_targets.fa
      org2_targets.bed
      org2_targets_summary.json
```

## Output files

All outputs are written to `<output>/targets/<name>/`, where `<name>` is the
genome filename stem (or the value from the `output` column in a batch file).
This layout keeps the pipeline stage (`targets/`) as the parent and organises
all genomes within it.

| File | Description |
|------|-------------|
| `<name>_targets.fa` | Target sequences in FASTA format. One entry per target region. |
| `<name>_targets.bed` | Genomic coordinates (BED format). Strand and feature type are included for gene mode. |
| `<name>_targets_summary.json` | JSON summary: `name`, `mode`, and `target_count`. |
| `report.html` | HTML report listing each target and its length. |

### Example layout

```
results/
  targets/
    org1/
      org1_targets.fa
      org1_targets.bed
      org1_targets_summary.json
      report.html
    org2/
      org2_targets.fa
      org2_targets.bed
      org2_targets_summary.json
      report.html
```

## Notes

- Target FASTA sequence names become the identifiers used throughout all
  subsequent pipeline stages.
- In gene mode, the target identifier is taken from the annotation `ID`, `gene_id`,
  or `Name` attribute (first found). If none is present, it defaults to
  `seqid:start-end`.
- Sequences containing `N` are not filtered at this stage; individual probe candidates
  that span `N` positions are dropped during the `probes` step.
