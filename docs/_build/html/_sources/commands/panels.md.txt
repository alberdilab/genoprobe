# `panels` — Assemble final probe panels

The `panels` command selects a final set of non-cross-reactive probes for each
target from the screened (or raw) candidate list. It enforces inter-probe
complementarity and heterodimer constraints to minimise probe–probe interactions
within the same hybridisation panel.

## Usage

```
genoprobe panels --output DIR [OPTIONS]
```

**Requires:** `probes` has been run. `screen` is optional but recommended — if
`<output>/screen/screened.tsv` exists it is used; otherwise `probes/candidates.tsv`
is used as input.

## Required arguments

| Argument | Description |
|----------|-------------|
| `--output / -o` | Output directory used in the preceding pipeline stages. |

## Optional arguments

### Profile

| Argument | Default | Description |
|----------|---------|-------------|
| `--profile` | `balanced` | Parameter preset: `fast`, `balanced`, or `strict`. |

See [Profiles](../profiles.md) for profile-specific parameter values.

### Panel size

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-probes N` | 50 | Maximum number of probes selected per target. |

## How panel assembly works

Probes are sorted by descending `final_score` (or `score` if screening was
skipped). A greedy selection algorithm then adds probes one by one, rejecting any
candidate that:

- would bring the selected set above `--max-probes` for that target; or
- has a contiguous self-complementary overlap with an already-selected probe longer
  than `max_panel_contiguous_complementarity` (profile parameter); or
- forms a heterodimer with an already-selected probe whose estimated Tm exceeds
  `max_panel_heterodimer_tm` (profile parameter, if set).

The complementarity and heterodimer penalties are also used to adjust scores during
selection when the corresponding profile penalty weights are non-zero.

## Output files

All outputs are written to `<output>/panels/`.

| File | Description |
|------|-------------|
| `final_probes.tsv` | Final probe sequences, one row per probe. Ready for synthesis ordering. |
| `panels_summary.json` | JSON summary: total probe count and per-target counts. |

### `final_probes.tsv` columns

| Column | Description |
|--------|-------------|
| `target` | Target name. |
| `start` | 0-based start on the target sequence. |
| `end` | 0-based exclusive end. |
| `sequence` | Probe sequence (5′→3′). |
| `length` | Length in nt. |
| `tm` | Calculated Tm in °C. |
| `gc` | GC content in %. |
| `score` | Final quality score (higher is better). |

## Notes

- To obtain more probes per target, increase `--max-probes` or switch to the `fast`
  profile which has a higher default cap.
- To reduce cross-reactivity between probes, switch to the `strict` profile or
  lower `max_panel_contiguous_complementarity` via a custom profile override.
