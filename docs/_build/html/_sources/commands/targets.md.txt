# `targets` — Define probe-design target regions

The `targets` command extracts or defines genomic regions to be used as probe
design templates. It writes a FASTA file of target sequences and a BED file of
their genomic coordinates.

## Usage

```
genoprobe targets --genomes FASTA [FASTA ...] --output DIR --mode {genome,gene} [OPTIONS]
```

## Required arguments

| Argument | Description |
|----------|-------------|
| `--genomes / -g` | One or more genome FASTA files. |
| `--output / -o` | Output directory. Sub-directories are created automatically. |
| `--mode / -m` | `genome` — tile the full genome(s); `gene` — use annotated features. |

## Mode: `genome`

Tiles every sequence in the input FASTA file(s). Use `--region` to restrict to a
single genomic interval.

```bash
genoprobe targets \
    --genomes genome.fa \
    --mode genome \
    --output run/
```

### `--region SEQID:START-END`

Restrict genome-mode tiling to one region. Coordinates are **1-based, inclusive**.

```bash
genoprobe targets \
    --genomes genome.fa \
    --mode genome \
    --region chr1:500000-1500000 \
    --output run/
```

## Mode: `gene`

Extracts sequences for annotated features from a GFF3 or GTF file.
`--annotation` is required in this mode.

```bash
genoprobe targets \
    --genomes genome.fa \
    --annotation genes.gff \
    --mode gene \
    --output run/
```

### `--annotation / -a GFF_GTF`

Path to the annotation file. Both GFF3 and GTF formats are supported.

### `--feature TYPE [TYPE ...]`

Filter to specific feature types. If omitted, all features in the annotation are
used.

```bash
# Use only gene-level features
genoprobe targets --mode gene --annotation genes.gff --feature gene ...

# Use CDS and exon features
genoprobe targets --mode gene --annotation genes.gff --feature CDS exon ...
```

## Output files

All outputs are written to `<output>/targets/`.

| File | Description |
|------|-------------|
| `targets.fa` | Target sequences in FASTA format. One entry per target region. |
| `targets.bed` | Genomic coordinates (BED format). Strand and feature type are included for gene mode. |
| `targets_summary.json` | JSON summary: `mode` and `target_count`. |

## Notes

- Target FASTA sequence names become the identifiers used throughout all
  subsequent pipeline stages.
- In gene mode, the target identifier is taken from the annotation `ID`, `gene_id`,
  or `Name` attribute (first found). If none is present, it defaults to
  `seqid:start-end`.
- Sequences containing `N` are not filtered at this stage; individual probe candidates
  that span `N` positions are dropped during the `probes` step.
