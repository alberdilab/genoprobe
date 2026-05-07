# Thermodynamics

genoprobe uses a nearest-neighbour (NN) thermodynamic model to calculate probe
melting temperatures (Tm) and to estimate the stability of secondary structures
(hairpins, homodimers, heterodimers).

## Melting temperature (Tm)

### Nearest-neighbour model

The default Tm calculation implements the unified nearest-neighbour parameters of
[SantaLucia (1998)](https://doi.org/10.1073/pnas.95.4.1460):

$$
T_m = \frac{\Delta H \times 1000}{\Delta S + R \ln(C_T / 4)} - 273.15
$$

where:
- $\Delta H$ and $\Delta S$ are the sum of NN enthalpy and entropy parameters over
  all dinucleotide steps, plus initiation corrections for terminal AT/GC pairs.
- $C_T$ is the total strand concentration (probe + target, in mol/L). The factor
  of 4 assumes a non-self-complementary duplex ($C_T / 4 \approx 12.5$ nM at the
  default 50 nM total concentration).
- $R = 1.987$ cal/(mol·K).

### Salt correction

The salt correction follows SantaLucia (1998) method 5 — an additive correction
to $\Delta S$:

$$
\Delta S_{salt} = 0.368 \times (n-1) \times \ln([\text{Na}^+] / \text{M})
$$

where $n$ is the probe length and $[\text{Na}^+]$ is the monovalent salt
concentration (default: 390 mM, representative of 2× SSC).

### Formamide correction

Formamide destabilises duplex formation. After NN Tm calculation, a linear
correction is applied:

$$
T_m^{corrected} = T_m^{NN} - 0.72 \times [\text{formamide}]\%
$$

The default formamide concentration is 50%, matching the OligoMiner / blockParse
convention and a typical FISH hybridisation buffer.

### primer3 backend

When `--thermo-backend primer3` is used, Tm is calculated using the
[primer3-py](https://github.com/libnano/primer3-py) library. The formamide
correction is still applied afterwards. Install the optional dependency with:

```bash
pip install "genoprobe[thermo]"
```

## Secondary structure estimates

### Self-complementarity

`self_complementarity(seq)` counts positions where the probe sequence matches its
own reverse complement at the same index (i.e. where `seq[i] == rc_seq[i]`). It
returns both the total matching count and the longest contiguous matching run.
These are used as filters and penalty terms during candidate scoring.

### Hairpin Tm

The hairpin algorithm scans for the longest stem in which the 5′ arm of the probe
is the reverse complement of the 3′ arm. If a stable stem is found, the stem
sequence is passed to `nn_tm` to estimate its melting temperature. A value of 0
is returned when no stem is found.

### Homodimer Tm

The homodimer estimate uses the self-complementary paired bases as a proxy duplex.
If fewer than 4 bases are self-complementary, the homodimer Tm is reported as 0.

### Heterodimer Tm

Used during panel assembly. A sliding-alignment search finds the longest
complementary run between probe A and the reverse complement of probe B. A
poly-GC proxy of that length is then passed to `nn_tm` as a conservative
upper-bound estimate.

## Default thermodynamic parameters

| Parameter | Default | Flag |
|-----------|---------|------|
| Monovalent salt (mM) | 390 | `--monovalent-mm` |
| Total strand concentration (nM) | 50 | `--probe-conc` |
| Formamide (%) | 50 | `--formamide` |
| Thermodynamic backend | `basic` | `--thermo-backend` |

## References

- SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide
  DNA nearest-neighbor thermodynamics. *PNAS*, 95(4), 1460–1465.
  <https://doi.org/10.1073/pnas.95.4.1460>
- Beliveau, B.J. et al. (2018). OligoMiner provides a rapid, flexible environment
  for the design of genome-scale oligonucleotide in situ hybridization probes.
  *PNAS*, 115(10), E2183–E2192. <https://doi.org/10.1073/pnas.1714530115>
