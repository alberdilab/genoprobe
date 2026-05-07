# Installation

## Requirements

- Python ≥ 3.10
- [pyfaidx](https://github.com/mdshw5/pyfaidx) ≥ 0.8
- [rich](https://github.com/Textualize/rich) ≥ 13.7
- [typer](https://typer.tiangolo.com/) ≥ 0.12

## Installing from PyPI

```bash
pip install genoprobe
```

## Installing from source

```bash
git clone https://github.com/anttonalberdi/genoprobe.git
cd genoprobe
pip install -e .
```

## Optional dependencies

### Thermodynamic backend (`primer3`)

The default Tm calculation uses the built-in nearest-neighbour model (no extra
dependencies). For higher-accuracy calculations, install the
[primer3-py](https://github.com/libnano/primer3-py) backend:

```bash
pip install "genoprobe[thermo]"
```

Then pass `--thermo-backend primer3` to the `probes` command.

### Conda environment

A conda `environment.yml` is provided for a fully reproducible environment:

```bash
conda env create -f environment.yml
conda activate genoprobe
```

## Optional: Fulgor index

The `screen` command can use a [Fulgor](https://github.com/jermp/fulgor) coloured
compacted de Bruijn graph index for fast, genome-wide k-mer lookups. Fulgor must
be installed separately and available on `PATH`. The `index` command builds the
index; see [Building a Fulgor index](commands/index_cmd.md).

## Verifying the installation

```bash
genoprobe --version
```

You should see output such as `genoprobe 0.1.0`.
