# genoprobe

**genoprobe** is a genome- and annotation-driven probe design toolkit for fluorescence in-situ hybridisation (FISH) and related capture-hybridisation assays.

It implements the core logic of the [OligoMiner / blockParse](https://github.com/beliveau-lab/OligoMiner) pipeline in a modern Python package with a clean command-line interface and a four-stage workflow:

| Stage | Command | Purpose |
|-------|---------|---------|
| 1 | `targets` | Extract target regions from genome(s) |
| 2 | `probes`  | Generate and filter probe candidates |
| 3 | `screen`  | Screen candidates for off-target binding |
| 4 | `panels`  | Assemble final, non-cross-reactive probe panels |

---

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
quickstart
```

```{toctree}
:maxdepth: 2
:caption: Command reference

commands/targets
commands/probes
commands/screen
commands/panels
commands/index_cmd
```

```{toctree}
:maxdepth: 2
:caption: Background

profiles
thermodynamics
```

```{toctree}
:maxdepth: 1
:caption: Developer reference

api/index
```
