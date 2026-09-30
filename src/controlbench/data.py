"""Portable result contract with strict sample and diagnostic validation."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .config import resolve


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    required = {
        "states.csv": "step time_s agent_id value reference lower upper state_name unit",
        "inputs.csv": "step time_s agent_id commanded applied lower upper input_name unit",
        "steps.csv": (
            "step solve_time_s plant_time_s iterations converged termination stopping_criterion"
        ),
        "local_solves.csv": "step iteration agent_id problem duration_s message",
        "residuals.csv": "step iteration agent_id primal dual normalized_primal",
        "events.csv": "step type message",
    }
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if path.name in required and set(reader.fieldnames or []) != set(
            required[path.name].split()
        ):
            raise ValueError(f"Invalid columns in {path.name}")
        values = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in values):
            raise ValueError(f"Incomplete row in {path.name}")
        return values


def finite(value: str) -> float:
    number = float(value)
    if not np.isfinite(number):
        raise ValueError("Nonfinite result value")
    return number


@dataclass(frozen=True)
class ExperimentResult:
    directory: Path
    config: dict
    states: NDArray[np.float64]
    references: NDArray[np.float64]
    inputs: NDArray[np.float64]
    solve_times: NDArray[np.float64]
    iterations: NDArray[np.float64]
    converged: NDArray[np.float64]
    local_solves: list[dict[str, str]]
    residuals: list[dict[str, str]]

    @classmethod
    def load(cls, directory: Path) -> ExperimentResult:
        c = resolve(read_json(directory / "config.resolved.json"))
        n = c["plant"]["agents"]
        dt = c["controller"]["sample_time_s"]
        k = round(c["scenario"]["duration_s"] / dt)

        def tensor(
            filename: str, count: int, field: str, component: str, expected_name: str, unit: str
        ) -> NDArray[np.float64]:
            table = rows(directory / filename)
            if len(table) != count * n:
                raise ValueError(f"Wrong sample count in {filename}")
            out = np.full((count, n), np.nan)
            for row in table:
                step, agent = int(row["step"]), int(row["agent_id"]) - 1
                if not 0 <= step < count or not 0 <= agent < n:
                    raise ValueError(f"Invalid sample key in {filename}")
                if np.isfinite(out[step, agent]):
                    raise ValueError(f"Duplicate sample in {filename}")
                if row[component] != expected_name or row["unit"] != unit:
                    raise ValueError(f"Invalid component/unit in {filename}")
                if not np.isclose(finite(row["time_s"]), step * dt, rtol=0, atol=1e-9):
                    raise ValueError(f"Time grid mismatch in {filename}")
                bound = c["constraints"]["level_m" if filename == "states.csv" else "input"]
                if finite(row["lower"]) != bound["min"] or finite(row["upper"]) != bound["max"]:
                    raise ValueError("Exported bounds differ from configuration")
                out[step, agent] = finite(row[field])
            return out

        x = tensor("states.csv", k + 1, "value", "state_name", "level", "m")
        ref = tensor("states.csv", k + 1, "reference", "state_name", "level", "m")
        if not np.all(ref == c["plant"]["reference_level_m"]):
            raise ValueError("Reference differs from configuration")
        if not np.allclose(x[0], c["plant"]["initial_level_m"], rtol=0, atol=1e-12):
            raise ValueError("Initial state differs from configuration")
        u = tensor("inputs.csv", k, "applied", "input_name", "flow", "m3/s")
        commanded = tensor("inputs.csv", k, "commanded", "input_name", "flow", "m3/s")
        if not np.array_equal(u, commanded):
            raise ValueError("v1 has no actuator transform: applied must equal commanded")
        table = rows(directory / "steps.csv")
        if [int(r["step"]) for r in table] != list(range(k)):
            raise ValueError("Incomplete or unordered step diagnostics")
        times = np.array([finite(r["solve_time_s"]) for r in table])
        its = np.array([finite(r["iterations"]) for r in table])
        conv = np.array([finite(r["converged"]) for r in table])
        if np.any(times < 0) or np.any(its < 1) or np.any(its != np.floor(its)):
            raise ValueError("Invalid solve times/iterations")
        if not np.isin(conv, [0, 1]).all():
            raise ValueError("Convergence flags must be 0 or 1")
        for i, row in enumerate(table):
            expected_criterion = (
                "normalized_primal_v1"
                if c["controller"]["type"] == "matlab_dmpc"
                else "local_optimizer_success"
            )
            if row["stopping_criterion"] != expected_criterion:
                raise ValueError("Unknown stopping criterion")
            maximum = (
                c["controller"]["admm"]["max_iterations"]
                if c["controller"]["type"] == "matlab_dmpc"
                else 1
            )
            if its[i] > maximum or (not conv[i] and its[i] != maximum):
                raise ValueError("Invalid iteration-limit termination")
            if finite(row["plant_time_s"]) < 0:
                raise ValueError("Negative plant time")
            if row["termination"] != ("converged" if conv[i] else "iteration_limit"):
                raise ValueError("Inconsistent termination diagnostics")
        local = rows(directory / "local_solves.csv")
        expected = (
            {
                (s, q, a)
                for s in range(k)
                for q in range(1, int(its[s]) + 1)
                for a in range(1, n + 1)
            }
            if c["controller"]["type"] == "matlab_dmpc"
            else {(s, 1, 0) for s in range(k)}
        )
        keys = [(int(r["step"]), int(r["iteration"]), int(r["agent_id"])) for r in local]
        if len(set(keys)) != len(keys) or set(keys) != expected:
            raise ValueError("Incomplete local solve diagnostics")
        for r in local:
            if finite(r["duration_s"]) < 0 or int(r["problem"]) != 0:
                raise ValueError("Failed or invalid local solve")
        residual = rows(directory / "residuals.csv")
        if c["controller"]["type"] == "matlab_dmpc":
            keys = [(int(r["step"]), int(r["iteration"]), int(r["agent_id"])) for r in residual]
            if len(set(keys)) != len(keys) or set(keys) != expected:
                raise ValueError("Incomplete residual diagnostics")
            for r in residual:
                for f in ("primal", "dual", "normalized_primal"):
                    if finite(r[f]) < 0:
                        raise ValueError("Negative residual")
            tol = c["controller"]["admm"]["convergence_tolerance"]
            for s in range(k):
                final = [
                    finite(r["normalized_primal"])
                    for r in residual
                    if int(r["step"]) == s and int(r["iteration"]) == its[s]
                ]
                if bool(conv[s]) != all(v <= tol for v in final):
                    raise ValueError("Convergence disagrees with recorded stopping criterion")
        if rows(directory / "events.csv"):
            raise ValueError("Exceptional events present; inspect execution logs")
        return cls(directory, c, x, ref, u, times, its, conv, local, residual)
