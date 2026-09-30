import csv

import pytest

from controlbench.config import resolve
from controlbench.data import write_json


@pytest.fixture
def exported(tmp_path):
    c = resolve({"plant": {"agents": 2}, "scenario": {"duration_s": 0.2}})
    write_json(tmp_path / "config.resolved.json", c)

    def table(name, fields, data):
        with (tmp_path / name).open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(fields)
            w.writerows(data)

    table(
        "states.csv",
        [
            "step",
            "time_s",
            "agent_id",
            "value",
            "reference",
            "lower",
            "upper",
            "state_name",
            "unit",
        ],
        [[s, s * 0.1, a, 0.5 + s * 0.5, 2, 0, 3, "level", "m"] for s in range(3) for a in (1, 2)],
    )
    table(
        "inputs.csv",
        [
            "step",
            "time_s",
            "agent_id",
            "commanded",
            "applied",
            "lower",
            "upper",
            "input_name",
            "unit",
        ],
        [[s, s * 0.1, a, 2, 2, 0, 1000, "flow", "m3/s"] for s in range(2) for a in (1, 2)],
    )
    table(
        "steps.csv",
        [
            "step",
            "solve_time_s",
            "plant_time_s",
            "iterations",
            "converged",
            "termination",
            "stopping_criterion",
        ],
        [[s, 0.01 * (s + 1), 0.001, 1, 1, "converged", "normalized_primal_v1"] for s in range(2)],
    )
    table(
        "local_solves.csv",
        ["step", "iteration", "agent_id", "problem", "duration_s", "message"],
        [[s, 1, a, 0, 0.001, "ok"] for s in range(2) for a in (1, 2)],
    )
    table(
        "residuals.csv",
        ["step", "iteration", "agent_id", "primal", "dual", "normalized_primal"],
        [[s, 1, a, 0, 0, 0] for s in range(2) for a in (1, 2)],
    )
    table("events.csv", ["step", "type", "message"], [])
    return tmp_path
