"""Final probe selection and panel assembly."""

from __future__ import annotations

from dataclasses import dataclass, field

from genoprobe.probes import ProbeCandidate
from genoprobe.thermo import (
    DEFAULT_PANEL_COMPLEMENTARITY_PENALTY_WEIGHT,
    DEFAULT_PANEL_HETERODIMER_PENALTY_WEIGHT,
    heterodimer_tm,
    self_complementarity,
)

DEFAULT_MAX_PROBES_PER_TARGET: int = 50
DEFAULT_MAX_PANEL_CONTIGUOUS_COMPLEMENTARITY: int | None = 7
DEFAULT_MAX_PANEL_HETERODIMER_TM: float | None = None


@dataclass
class PanelConfig:
    max_probes_per_target: int = DEFAULT_MAX_PROBES_PER_TARGET
    max_panel_contiguous_complementarity: int | None = DEFAULT_MAX_PANEL_CONTIGUOUS_COMPLEMENTARITY
    max_panel_heterodimer_tm: float | None = DEFAULT_MAX_PANEL_HETERODIMER_TM
    panel_complementarity_penalty_weight: float = DEFAULT_PANEL_COMPLEMENTARITY_PENALTY_WEIGHT
    panel_heterodimer_penalty_weight: float = DEFAULT_PANEL_HETERODIMER_PENALTY_WEIGHT


@dataclass
class PanelResult:
    target_name: str
    selected_probes: list[ProbeCandidate]
    config: PanelConfig
    screened: bool = False

    @property
    def probe_count(self) -> int:
        return len(self.selected_probes)


def _panel_penalty(probe: ProbeCandidate, panel: list[ProbeCandidate], config: PanelConfig) -> float:
    if not panel:
        return 0.0
    penalty = 0.0
    for existing in panel:
        _total, run = self_complementarity(probe.sequence)
        if (
            config.max_panel_contiguous_complementarity is not None
            and run > config.max_panel_contiguous_complementarity
        ):
            penalty += config.panel_complementarity_penalty_weight

        if config.panel_heterodimer_penalty_weight > 0 or config.max_panel_heterodimer_tm is not None:
            ht = heterodimer_tm(probe.sequence, existing.sequence)
            if config.max_panel_heterodimer_tm is not None and ht > config.max_panel_heterodimer_tm:
                return 1.0  # hard reject
            if config.panel_heterodimer_penalty_weight > 0 and ht > 0:
                penalty += config.panel_heterodimer_penalty_weight * min(1.0, ht / 60.0)

    return min(1.0, penalty)


def assemble_panel(
    candidates: list[ProbeCandidate],
    target_name: str,
    config: PanelConfig,
) -> PanelResult:
    """Greedy panel assembly: iteratively add highest-scoring non-penalised probes."""
    panel: list[ProbeCandidate] = []
    for probe in candidates:
        if len(panel) >= config.max_probes_per_target:
            break
        penalty = _panel_penalty(probe, panel, config)
        if penalty >= 1.0:
            continue
        panel.append(probe)
    return PanelResult(
        target_name=target_name,
        selected_probes=panel,
        config=config,
    )
