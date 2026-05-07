# `index` — Build a Fulgor index

The `index` command builds an optional [Fulgor](https://github.com/jermp/fulgor)
coloured compacted de Bruijn graph index from one or more genome FASTA files. The
index can then be passed to the `screen` command via `--index` for fast, exact
k-mer presence queries across large collections of genomes.

## Usage

```
genoprobe index --genomes FASTA [FASTA ...] --output DIR [OPTIONS]
```

**Requires:** Fulgor must be installed and on `PATH`.

## Required arguments

| Argument | Description |
|----------|-------------|
| `--genomes / -g` | One or more genome FASTA files to include in the index. |
| `--output / -o` | Output directory. The index is written to `<output>/index/`. |

## Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--kmer-length K` | 31 | k-mer length used to build the de Bruijn graph. |
| `--minimizer-length M` | 19 | Minimizer length for the coloured graph. |
| `--threads N` | 4 | Number of threads passed to Fulgor. |

## Installing Fulgor

Fulgor is a separate C++ tool and is not bundled with genoprobe. Installation
instructions are available at the [Fulgor GitHub repository](https://github.com/jermp/fulgor).

A typical install from source:

```bash
git clone --recursive https://github.com/jermp/fulgor.git
cd fulgor
cmake -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build
# add the build/ directory to PATH or copy the binary to /usr/local/bin
```

## Example

```bash
# Build an index from two genome assemblies
genoprobe index \
    --genomes genome_A.fa genome_B.fa \
    --output my_run/ \
    --kmer-length 31 \
    --threads 8
```

The index is written to `my_run/index/`. Pass this directory to `screen`:

```bash
genoprobe screen \
    --genomes genome_A.fa genome_B.fa \
    --output my_run/ \
    --index my_run/index/
```

## Notes

- Building a Fulgor index is optional. The `screen` command performs its own
  in-memory k-mer frequency calculation without Fulgor; the index adds a secondary,
  exact presence-or-absence layer across multiple genomes simultaneously.
- For single-genome workflows the k-mer frequency filter (`--max-kmer-frequency`)
  in `screen` is usually sufficient.
