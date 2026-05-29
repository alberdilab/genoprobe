# Changelog

## [Unreleased]

### Added

- No unreleased changes yet.

## [0.1.3] - 2026-05-29

### Changed

- `targets` command: `--targets-file / -t` renamed to `--file / -f` for brevity.
- `targets` command: `--mode / -m` is now optional. When omitted, mode is auto-detected — `gene` if `--annotation` is provided (single-genome mode) or the batch file header contains an `annotations` / `annotation` column with at least one value; `genome` otherwise. An explicit `--mode` always takes precedence.
- `targets` batch file: `outputs` is now accepted as an alias for the `output` column header, making all column names consistently available in both singular and plural forms (`genome`/`genomes`, `annotation`/`annotations`, `output`/`outputs`).

## [0.1.2] - 2026-05-29

### Added

- `targets` command: new `--targets-file / -t` argument accepting a TSV or CSV file of genome-annotation pairs for batch processing. Each row is run independently. Delimiter and header row are auto-detected.
- Consistent per-genome output structure across all pipeline stages. Each stage directory (`targets/`, `probes/`, `screen/`, `panels/`) now contains a named subfolder per genome, with all files prefixed by that name: e.g. `targets/<name>/<name>_targets.fa`, `probes/<name>/<name>_candidates.tsv`, `screen/<name>/<name>_screened.tsv`, `panels/<name>/<name>_final_probes.tsv`. The name defaults to the genome filename stem and can be overridden via the `output` column in a batch file.
- `probes`, `screen`, and `panels` commands now auto-discover genome names by scanning the preceding stage's subdirectories, so no per-genome arguments are needed after `targets`.
- `panels` per-genome source selection: prefers `<name>_screened.tsv` when available, falls back to `<name>_candidates.tsv` on a per-genome basis.

## [0.1.1] - 2026-05-07

### Added

- Sphinx documentation site compatible with ReadTheDocs (`docs/`), covering installation, quick-start guide, per-command reference (targets, probes, screen, panels, index), profile parameter tables, thermodynamic model background, and auto-generated Python API reference via autodoc.
- ReadTheDocs build configuration (`.readthedocs.yaml`).
- GitHub Actions release workflow (`.github/workflows/release.yml`) and release helper script (`scripts/release.py`).

## [0.1.0] - 2026-05-07

### Added

- Initial skeleton of the GenoProbe genome- and annotation-driven probe design toolkit.
