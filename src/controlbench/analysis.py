"""Versioned analysis registry; metrics never decide engineering acceptance."""

from __future__ import annotations

from importlib.metadata import entry_points
from typing import Protocol

import numpy as np

from .data import ExperimentResult


class Analyzer(Protocol):
    name: str
    version: str

    def evaluate(self, result: ExperimentResult) -> dict[str, dict]: ...


def metric(value: float, unit: str) -> dict:
    if not np.isfinite(value):
        raise ValueError("Nonfinite metric")
    return {"value": float(value), "unit": unit}


class Tracking:
    name, version = "tracking", "1"

    def evaluate(self, r: ExperimentResult) -> dict[str, dict]:
        error = r.states - r.references
        result = {
            "tracking_rmse": metric(np.sqrt(np.mean(error**2)), "m"),
            "max_tracking_error": metric(np.max(np.abs(error)), "m"),
        }
        for a in range(error.shape[1]):
            result[f"agent_{a + 1}_tracking_rmse"] = metric(np.sqrt(np.mean(error[:, a] ** 2)), "m")
        return result


class Constraints:
    name, version = "constraints", "1"

    def evaluate(self, r: ExperimentResult) -> dict[str, dict]:
        result = {}
        for name, data, key, unit in (
            ("state", r.states, "level_m", "m"),
            ("input", r.inputs, "input", "m3/s"),
        ):
            b = r.config["constraints"][key]
            v = np.maximum(0, np.maximum(b["min"] - data, data - b["max"]))
            mask = v > r.config["constraints"]["tolerance"]
            result[f"{name}_violation_count"] = metric(float(np.count_nonzero(mask)), "count")
            result[f"{name}_violation_steps"] = metric(
                float(np.count_nonzero(mask.any(axis=1))), "count"
            )
            result[f"max_{name}_violation"] = metric(np.max(v), unit)
        return result


class Solver:
    name, version = "solver", "1"

    def evaluate(self, r: ExperimentResult) -> dict[str, dict]:
        t = 1000 * r.solve_times
        return {
            "mean_solve_time_ms": metric(np.mean(t), "ms"),
            "p95_solve_time_ms": metric(np.quantile(t, 0.95, method="linear"), "ms"),
            "max_solve_time_ms": metric(np.max(t), "ms"),
            "mean_iterations": metric(np.mean(r.iterations), "count"),
            "failed_local_solves": metric(
                sum(int(s["problem"]) != 0 for s in r.local_solves), "count"
            ),
            "converged_step_fraction": metric(np.mean(r.converged), "fraction"),
        }


class Effort:
    name, version = "effort", "1"

    def evaluate(self, r: ExperimentResult) -> dict[str, dict]:
        dt = r.config["controller"]["sample_time_s"]
        p = r.config["plant"]
        w = np.full(r.states.shape[1], p["state_weight"], dtype=float)
        w[[0, -1]] = p["endpoint_state_weight"]
        err = (r.states - r.references) ** 2
        effort = dt * np.sum(r.inputs**2)
        physical = (
            dt * np.sum(err[:-1] * w)
            + p["input_weight"] * effort
            + p["terminal_weight"] * np.sum(err[-1] * w)
        )
        return {
            "control_effort": metric(effort, "m6/s"),
            "physical_cost": metric(physical, "weighted_cost"),
        }


class Residuals:
    name, version = "residuals", "1"

    def evaluate(self, r: ExperimentResult) -> dict[str, dict]:
        if not r.residuals:
            return {}
        final = [
            float(v["normalized_primal"])
            for v in r.residuals
            if int(v["iteration"]) == r.iterations[int(v["step"])]
        ]
        return {
            "peak_normalized_primal_residual": metric(
                max(float(v["normalized_primal"]) for v in r.residuals), "mixed_norm"
            ),
            "max_final_normalized_primal_residual": metric(max(final), "mixed_norm"),
            "peak_dual_residual": metric(max(float(v["dual"]) for v in r.residuals), "mixed_norm"),
        }


REGISTRY: dict[str, Analyzer] = {
    a.name: a for a in (Tracking(), Constraints(), Solver(), Effort(), Residuals())
}
EXTRA_METRICS: dict[str, str] = {}
_discovered = False


def register(analyzer: Analyzer) -> None:
    if analyzer.name in REGISTRY:
        raise ValueError(f"Duplicate analyzer: {analyzer.name}")
    advertised = getattr(analyzer, "metric_units", {})
    if not isinstance(advertised, dict) or not advertised:
        raise ValueError("Custom analyzers must declare metric_units")
    if any(not isinstance(k, str) or not isinstance(v, str) for k, v in advertised.items()):
        raise ValueError("Plugin metric names and units must be strings")
    if set(advertised) & metric_names(0, discover=False):
        raise ValueError("Duplicate plugin metric identifier")
    if any(k.startswith("agent_") for k in advertised):
        raise ValueError("agent_ is a reserved metric prefix")
    EXTRA_METRICS.update(advertised)
    REGISTRY[analyzer.name] = analyzer


def discover_plugins() -> None:
    global _discovered
    if not _discovered:
        original_registry, original_metrics = dict(REGISTRY), dict(EXTRA_METRICS)
        try:
            for entry in entry_points(group="controlbench.analyzers"):
                register(entry.load()())
        except Exception:
            REGISTRY.clear()
            REGISTRY.update(original_registry)
            EXTRA_METRICS.clear()
            EXTRA_METRICS.update(original_metrics)
            raise
        _discovered = True


def analyze(result: ExperimentResult) -> dict:
    discover_plugins()
    metrics: dict = {}
    for analyzer in REGISTRY.values():
        added = analyzer.evaluate(result)
        if metrics.keys() & added.keys():
            raise ValueError("Duplicate metric identifier")
        advertised = getattr(analyzer, "metric_units", None)
        if advertised is not None and added.keys() != advertised.keys():
            raise ValueError("Plugin did not return its declared metrics")
        for key, value in added.items():
            if not np.isfinite(value["value"]) or not isinstance(value["unit"], str):
                raise ValueError("Invalid analyzer metric")
            if advertised is not None and value["unit"] != advertised[key]:
                raise ValueError("Plugin metric unit mismatch")
        metrics.update(
            {
                k: {**v, "analyzer": analyzer.name, "version": analyzer.version}
                for k, v in added.items()
            }
        )
    return metrics


def metric_names(agents: int, discover: bool = True) -> set[str]:
    if discover:
        discover_plugins()
    return {
        "tracking_rmse",
        "max_tracking_error",
        "control_effort",
        "physical_cost",
        "peak_normalized_primal_residual",
        "max_final_normalized_primal_residual",
        "peak_dual_residual",
        *EXTRA_METRICS,
        "mean_solve_time_ms",
        "p95_solve_time_ms",
        "max_solve_time_ms",
        "mean_iterations",
        "failed_local_solves",
        "converged_step_fraction",
        *[f"agent_{a}_tracking_rmse" for a in range(1, agents + 1)],
        *[
            f"{kind}_{m}"
            for kind in ("state", "input")
            for m in ("violation_count", "violation_steps")
        ],
        "max_state_violation",
        "max_input_violation",
    }
