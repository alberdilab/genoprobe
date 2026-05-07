# Profiles

Every pipeline stage (`probes`, `screen`, `panels`) accepts a `--profile` flag
that selects a named preset of parameter values. Three built-in profiles are
provided:

| Profile | Use case |
|---------|----------|
| `fast` | Rapid exploration — relaxed filters, no thermodynamic penalty weight. |
| `balanced` | **Default** — good specificity for most FISH / capture-seq experiments. |
| `strict` | Demanding assays — tight Tm window, strong self-complementarity and cross-reactivity filters. |

Profile values can be overridden on the command line (e.g. `--min-tm 44`). The
command-line value always takes precedence over the profile value.

---

## `probes` profiles

| Parameter | `fast` | `balanced` | `strict` |
|-----------|--------|-----------|---------|
| `min_probe_length` | 36 | 36 | 38 |
| `max_probe_length` | 41 | 41 | 40 |
| `min_tm` (°C) | 40.0 | 42.0 | 44.0 |
| `max_tm` (°C) | 50.0 | 47.0 | 47.0 |
| `min_gc` (%) | 20.0 | 20.0 | 20.0 |
| `max_gc` (%) | 80.0 | 75.0 | 70.0 |
| `max_homopolymer_run` | 5 | 4 | 4 |
| `max_dinucleotide_run` | — | 4 | 3 |
| `min_sequence_entropy` | — | 0.6 | 0.7 |
| `max_self_complementarity` | — | 14 | 12 |
| `max_contiguous_self_complementarity` | — | 6 | 5 |
| `max_hairpin_tm` (°C) | — | — | 40.0 |
| `max_homodimer_tm` (°C) | — | — | 40.0 |
| `probe_spacing` | 0 | 0 | 0 |
| `thermo_penalty_weight` | 0.0 | 0.15 | 0.30 |

`—` means the filter is disabled (no upper bound applied).

---

## `screen` profiles

| Parameter | `fast` | `balanced` | `strict` |
|-----------|--------|-----------|---------|
| `max_mismatches` | 3 | 2 | 1 |
| `max_kmer_frequency` | — | 100 | 50 |
| `offtarget_penalty_weight` | 0.5 | 1.0 | 2.0 |
| `min_ontarget_fraction` | 0.50 | 0.70 | 0.85 |

---

## `panels` profiles

| Parameter | `fast` | `balanced` | `strict` |
|-----------|--------|-----------|---------|
| `max_probes_per_target` | 100 | 50 | 20 |
| `max_panel_contiguous_complementarity` | — | 7 | 6 |
| `max_panel_heterodimer_tm` (°C) | — | — | 40.0 |
| `panel_complementarity_penalty_weight` | 0.0 | 0.15 | 0.25 |
| `panel_heterodimer_penalty_weight` | 0.0 | 0.0 | 0.30 |

---

## Parameter glossary

**`thermo_penalty_weight`**
: Weight applied to the thermodynamic component of the probe score. At `0.0` the
  score is based purely on Tm proximity and GC content. Higher values penalise
  probes with high self-complementarity, hairpin, or homodimer Tm.

**`offtarget_penalty_weight`**
: Multiplier for the off-target hit penalty applied during `screen`. Higher values
  more aggressively reduce the score of probes with off-target matches.

**`min_ontarget_fraction`**
: Minimum fraction of alignment strength that must come from on-target hits. Probes
  below this threshold are discarded entirely during screening.

**`max_panel_contiguous_complementarity`**
: During panel assembly, any probe pair whose contiguous complementary run exceeds
  this value is considered incompatible. The second probe is not added to the panel
  if the first is already selected.

**`panel_complementarity_penalty_weight` / `panel_heterodimer_penalty_weight`**
: Weights used to soft-penalise the candidate score during greedy panel assembly
  when complementarity or heterodimer criteria are partially violated.
