import math

import numpy as np
import pytest

from controlbench.analysis import analyze
from controlbench.campaigns import expand
from controlbench.config import resolve
from controlbench.data import ExperimentResult
from controlbench.validation import evaluate, regression


@pytest.mark.parametrize(
    "bad",
    [
        {"typo": 1},
        {"controller": {"horizon_steps": 1}},
        {"scenario": {"duration_s": 0.15}},
        {"plant": {"agents": True}},
        {"constraints": {"input": {"min": 2, "max": 1}}},
        {"execution": {"timeout_s": float("nan")}},
        {"requirements": [{"id": "a", "metric": "unknown", "operator": "le", "limit": 0}]},
        {"controller": {"admm": {"adaptive_penalty": "true"}}},
        {"plant": {"model": "vehicle_platoon"}},
    ],
)
def test_invalid_configuration(bad):
    with pytest.raises(ValueError):
        resolve(bad)


def test_resolved_defaults_are_independent():
    a = resolve({})
    a["plant"]["agents"] = 8
    assert resolve({})["plant"]["agents"] == 4


def test_independent_metric_values(exported):
    result = ExperimentResult.load(exported)
    metrics = analyze(result)
    assert result.states.shape == (3, 2)
    assert result.inputs.shape == (2, 2)
    assert metrics["tracking_rmse"]["value"] == pytest.approx(math.sqrt(3.5 / 3))
    assert metrics["control_effort"]["value"] == pytest.approx(1.6)
    assert metrics["mean_solve_time_ms"]["value"] == pytest.approx(15)
    assert metrics["p95_solve_time_ms"]["value"] == pytest.approx(19.5)
    assert metrics["state_violation_count"]["value"] == 0


@pytest.mark.parametrize(
    "filename,old,new",
    [
        ("states.csv", "1.0,2,0,3", "nan,2,0,3"),
        ("states.csv", "level,m", "level,cm"),
        ("states.csv", "0.2,1", "0.3,1"),
        ("local_solves.csv", "1,1,0,0.001", "1,1,1,0.001"),
        ("steps.csv", "converged", "iteration_limit"),
    ],
)
def test_corrupt_exports_fail(exported, filename, old, new):
    p = exported / filename
    text = p.read_text()
    assert old in text
    p.write_text(text.replace(old, new))
    with pytest.raises(ValueError):
        ExperimentResult.load(exported)


def test_missing_samples_fail(exported):
    p = exported / "states.csv"
    p.write_text("\n".join(p.read_text().splitlines()[:-1]))
    with pytest.raises(ValueError, match="sample count"):
        ExperimentResult.load(exported)


def test_constraint_violation_visible(exported):
    p = exported / "states.csv"
    p.write_text(p.read_text().replace("1.5,2,0,3", "3.5,2,0,3"))
    m = analyze(ExperimentResult.load(exported))
    assert m["state_violation_count"]["value"] == 2
    assert m["max_state_violation"]["value"] == 0.5


def test_requirements_fail_closed():
    rule = {"id": "x", "metric": "m", "operator": "le", "limit": 1}
    assert evaluate({}, [rule])["overall"] == "ERROR"
    assert evaluate({"m": {"value": np.nan}}, [rule])["overall"] == "ERROR"
    assert evaluate({"m": {"value": 2}}, [rule])["overall"] == "FAIL"
    assert evaluate({"m": {"value": 1}}, [rule])["overall"] == "PASS"
    assert evaluate({}, [])["overall"] == "NOT_EVALUATED"


@pytest.mark.parametrize(
    "b,c,direction,expected",
    [
        (0, 1, "lower", "FAIL"),
        (10, 11, "lower", "PASS"),
        (10, 12, "lower", "FAIL"),
        (10, 8, "higher", "FAIL"),
    ],
)
def test_regression_direction_and_zero_baseline(b, c, direction, expected):
    def m(v):
        return {"x": {"value": v, "unit": "m", "version": "1"}}

    rules = [{"metric": "x", "direction": direction, "absolute": 0, "relative": 0.1}]
    assert regression(m(b), m(c), rules)["overall"] == expected
    assert regression({}, m(c), rules)["overall"] == "ERROR"


def test_cartesian_sweep():
    configs = expand(
        {"base": {}, "sweep": {"plant.agents": [2, 4], "controller.horizon_steps": [3, 5, 10]}}
    )
    assert len(configs) == 6
    assert configs[0]["plant"]["agents"] == 2
    assert configs[-1]["controller"]["horizon_steps"] == 10


@pytest.mark.parametrize(
    "raw",
    [
        {"base": {}, "sweep": {"plant.agents": [2, 2]}},
        {"base": {}, "sweep": {"typo": [1]}},
        {"base": {}, "sweep": {"plant.agents": []}},
        {"base": {}, "sweep": {"plant.agents": [2, 3]}, "max_runs": 1},
    ],
)
def test_invalid_campaign(raw):
    with pytest.raises(ValueError):
        expand(raw)
