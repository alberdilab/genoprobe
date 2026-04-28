"""Load bundled stage-profile defaults."""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files
import json
from typing import Any


@lru_cache(maxsize=1)
def load_profile_catalog() -> dict[str, Any]:
    resource = files("genoprobe").joinpath("data/profile_defaults.json")
    with resource.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Profile catalog must be a JSON object.")
    return payload


def get_stage_default_profile(stage: str) -> str:
    section = load_profile_catalog().get(stage)
    if not isinstance(section, dict):
        raise ValueError(f"Profile catalog is missing stage '{stage}'.")
    profile = section.get("default_profile")
    if not isinstance(profile, str) or not profile.strip():
        raise ValueError(f"Profile catalog stage '{stage}' is missing a valid default_profile.")
    return profile.strip().lower()


def get_stage_profile_names(stage: str) -> tuple[str, ...]:
    section = load_profile_catalog().get(stage)
    if not isinstance(section, dict):
        raise ValueError(f"Profile catalog is missing stage '{stage}'.")
    profiles = section.get("profiles")
    if not isinstance(profiles, dict) or not profiles:
        raise ValueError(f"Profile catalog stage '{stage}' is missing profile definitions.")
    return tuple(sorted(str(name).strip().lower() for name in profiles))


def get_stage_profile_defaults(stage: str, profile: str) -> dict[str, Any]:
    section = load_profile_catalog().get(stage)
    if not isinstance(section, dict):
        raise ValueError(f"Profile catalog is missing stage '{stage}'.")
    profiles = section.get("profiles")
    if not isinstance(profiles, dict) or not profiles:
        raise ValueError(f"Profile catalog stage '{stage}' is missing profile definitions.")
    normalized = str(profile).strip().lower()
    if normalized not in profiles:
        valid = ", ".join(sorted(str(name) for name in profiles))
        raise ValueError(
            f"Profile catalog stage '{stage}' does not define profile '{profile}'. "
            f"Choose one of: {valid}."
        )
    values = profiles[normalized]
    if not isinstance(values, dict):
        raise ValueError(
            f"Profile catalog stage '{stage}' profile '{profile}' must map to a JSON object."
        )
    return dict(values)
