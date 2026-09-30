import shutil

import pytest

from controlbench import execution
from controlbench.data import read_json, write_json
from scripts.ci.stage_simulation_results import stage_run


def test_stage_run_accepts_only_verified_pass(exported, monkeypatch, tmp_path):
    def fake_matlab(root, out, entry, timeout, config_path=None):
        for path in exported.glob("*.csv"):
            shutil.copy2(path, out / path.name)
        write_json(out / "matlab.json", {"matlab": "test", "yalmip": "test"})
        (out / "stdout.log").write_text("fixture")
        (out / "stderr.log").write_text("")

    monkeypatch.setattr(execution, "matlab_process", fake_matlab)
    config = read_json(exported / "config.resolved.json")
    config["requirements"] = [
        {"id": "valid", "metric": "state_violation_count", "operator": "le", "limit": 0}
    ]
    results = tmp_path / "results"
    run, code = execution.execute(config, results)
    assert code == 0

    staging = tmp_path / "staging"
    staged = stage_run(results, staging)

    metadata = read_json(staging / "ci-metadata.json")
    assert staged.name == run.name == metadata["run_id"]
    assert metadata["validation"] == "PASS"
    assert (staged / "report.html").is_file()


def test_stage_run_rejects_failed_validation(exported, monkeypatch, tmp_path):
    def fake_matlab(root, out, entry, timeout, config_path=None):
        for path in exported.glob("*.csv"):
            shutil.copy2(path, out / path.name)
        write_json(out / "matlab.json", {"matlab": "test", "yalmip": "test"})
        (out / "stdout.log").write_text("fixture")
        (out / "stderr.log").write_text("")

    monkeypatch.setattr(execution, "matlab_process", fake_matlab)
    config = read_json(exported / "config.resolved.json")
    config["requirements"] = [
        {"id": "impossible", "metric": "tracking_rmse", "operator": "le", "limit": 0}
    ]
    results = tmp_path / "results"
    _, code = execution.execute(config, results)
    assert code == 1

    with pytest.raises(ValueError, match="Refusing to publish"):
        stage_run(results, tmp_path / "staging")
