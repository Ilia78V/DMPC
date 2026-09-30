"""Strict configuration resolution. No executable expressions are accepted."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import yaml

DEFAULTS: dict[str, Any] = {
    "schema_version": 1,
    "experiment": {"name": "tanks", "seed": 42},
    "plant": {
        "model": "coupled_water_tanks",
        "agents": 4,
        "initial_level_m": 0.5,
        "reference_level_m": 2.0,
        "coupling_model": "legacy_polynomial_4",
        "area_m2": 0.1,
        "coupling_gain": 1.0,
        "endpoint_gain": 6.0,
        "state_weight": 1.0,
        "endpoint_state_weight": 0.1,
        "input_weight": 0.1,
        "terminal_weight": 1.0,
    },
    "controller": {
        "type": "matlab_dmpc",
        "horizon_steps": 20,
        "sample_time_s": 0.1,
        "optimizer": "ipopt",
        "admm": {
            "max_iterations": 120,
            "convergence_tolerance": 0.001,
            "initial_penalty": 20.0,
            "adaptive_penalty": False,
            "approximation": "none",
        },
    },
    "scenario": {"duration_s": 3.0, "integrator": "forward_euler"},
    "constraints": {
        "level_m": {"min": 0.0, "max": 3.0},
        "input": {"min": 0.0, "max": 1000.0},
        "tolerance": 1e-6,
    },
    "requirements": [],
    "execution": {"timeout_s": 600.0},
}


def merge(defaults: dict, supplied: dict, path: str = "") -> dict:
    if not isinstance(supplied, dict):
        raise ValueError(f"{path or 'configuration'} must be a mapping")
    unknown = supplied.keys() - defaults.keys()
    if unknown:
        raise ValueError(f"Unknown keys at {path or 'root'}: {sorted(unknown)}")
    result = copy.deepcopy(defaults)
    for key, value in supplied.items():
        result[key] = (
            merge(defaults[key], value, f"{path}.{key}")
            if isinstance(defaults[key], dict)
            else value
        )
    return result


def number(value: Any, name: str, minimum: float | None = None, integer: bool = False) -> None:
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite numeric data")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    if integer and not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")


def resolve(raw: dict) -> dict:
    c = merge(DEFAULTS, raw)
    if c["schema_version"] != 1 or isinstance(c["schema_version"], bool):
        raise ValueError("Unsupported schema_version")
    if not isinstance(c["experiment"]["name"], str) or not c["experiment"]["name"].strip():
        raise ValueError("experiment.name must be nonempty")
    number(c["experiment"]["seed"], "seed", 0, True)
    if c["experiment"]["seed"] > 2**32 - 1:
        raise ValueError("seed must fit a uint32")
    p, ctl, s = c["plant"], c["controller"], c["scenario"]
    if p["model"] not in ("coupled_water_tanks", "linear_chain"):
        raise ValueError("Supported plants: coupled_water_tanks, linear_chain")
    if p["coupling_model"] not in ("legacy_polynomial_4", "linear"):
        raise ValueError("Unsupported coupling_model")
    if (p["model"] == "linear_chain") != (p["coupling_model"] == "linear"):
        raise ValueError("linear_chain requires linear coupling; tanks require legacy_polynomial_4")
    number(p["agents"], "agents", 2, True)
    for k in ("initial_level_m", "reference_level_m"):
        number(p[k], k)
    for k in (
        "coupling_gain",
        "endpoint_gain",
        "state_weight",
        "endpoint_state_weight",
        "input_weight",
        "terminal_weight",
    ):
        number(p[k], k, 0)
    number(p["area_m2"], "area_m2", 1e-12)
    if ctl["type"] not in ("matlab_dmpc", "matlab_mpc") or ctl["optimizer"] != "ipopt":
        raise ValueError("Supported controllers: matlab_dmpc, matlab_mpc with ipopt")
    number(ctl["horizon_steps"], "horizon_steps", 2, True)
    number(ctl["sample_time_s"], "sample_time_s", 1e-12)
    number(s["duration_s"], "duration_s", 1e-12)
    ratio = s["duration_s"] / ctl["sample_time_s"]
    if not math.isclose(ratio, round(ratio), abs_tol=1e-9, rel_tol=0) or ratio < 1:
        raise ValueError("duration_s must be an integral number of sample intervals")
    if s["integrator"] != "forward_euler":
        raise ValueError("Only forward_euler is supported")
    a = ctl["admm"]
    number(a["max_iterations"], "max_iterations", 1, True)
    number(a["convergence_tolerance"], "convergence_tolerance", 1e-12)
    number(a["initial_penalty"], "initial_penalty", 1e-12)
    if not isinstance(a["adaptive_penalty"], bool):
        raise ValueError("adaptive_penalty must be boolean")
    if a["approximation"] not in ("none", "all"):
        raise ValueError("approximation must be none or all")
    number(c["execution"]["timeout_s"], "timeout_s", 1)
    number(c["constraints"]["tolerance"], "tolerance", 0)
    for key in ("level_m", "input"):
        bounds = c["constraints"][key]
        number(bounds["min"], f"{key}.min")
        number(bounds["max"], f"{key}.max")
        if bounds["min"] > bounds["max"]:
            raise ValueError(f"Reversed {key} bounds")
    from .analysis import metric_names

    ids: set[str] = set()
    if not isinstance(c["requirements"], list):
        raise ValueError("requirements must be a list")
    for r in c["requirements"]:
        if not isinstance(r, dict) or set(r) != {"id", "metric", "operator", "limit"}:
            raise ValueError("Requirement needs exactly id, metric, operator, limit")
        if not isinstance(r["id"], str) or not r["id"] or r["id"] in ids:
            raise ValueError("Requirement IDs must be unique nonempty strings")
        ids.add(r["id"])
        if r["metric"] not in metric_names(p["agents"]):
            raise ValueError(f"Unknown metric: {r['metric']}")
        if r["operator"] not in ("le", "ge", "eq"):
            raise ValueError("Unsupported requirement operator")
        number(r["limit"], "requirement.limit")
    return c


class ConfigLoader(yaml.SafeLoader):
    """Safe loader with scientific notation and duplicate-key rejection."""


def unique_mapping(loader: ConfigLoader, node: yaml.MappingNode) -> dict:
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node)
        if not isinstance(key, str) or key in mapping:
            raise ValueError("Configuration keys must be unique strings")
        mapping[key] = loader.construct_object(value_node)
    return mapping


ConfigLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique_mapping)
ConfigLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(r"^[-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)[eE][-+]?[0-9]+$"),
    list("-+0123456789."),
)


def document(path: Path):
    return yaml.load(path.read_text(encoding="utf-8"), Loader=ConfigLoader)


def load(path: Path) -> dict:
    return resolve(document(path))


def digest(config: dict) -> str:
    return hashlib.sha256(json.dumps(config, sort_keys=True, allow_nan=False).encode()).hexdigest()
