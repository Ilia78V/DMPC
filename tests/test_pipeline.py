import shutil
import subprocess
from pathlib import Path

import pytest

from controlbench import benchmark as bench
from controlbench import campaigns, execution
from controlbench.cli import main
from controlbench.config import resolve
from controlbench.data import read_json, write_json
from controlbench.reporting import report


@pytest.fixture
def fake_matlab(exported, monkeypatch):
    def run(root, out, entry, timeout, config_path=None):
        for p in exported.glob("*.csv"):
            shutil.copy2(p, out / p.name)
        write_json(out / "matlab.json", {"matlab": "test", "yalmip": "test", "toolboxes": []})
        (out / "stdout.log").write_text("synthetic test fixture, not a MATLAB benchmark")
        (out / "stderr.log").write_text("")

    monkeypatch.setattr(execution, "matlab_process", run)
    return read_json(exported / "config.resolved.json")


def test_complete_lifecycle(fake_matlab, tmp_path):
    directory, code = execution.execute(fake_matlab, tmp_path)
    assert code == 0
    assert execution.verify_artifacts(directory)["status"] == "complete"
    assert (directory / "report.html").exists()
    assert (directory / "states.png").exists()
    assert (directory / "source" / "ADMM_Solver.m").exists()
    repeated, code = execution.reproduce(directory, tmp_path)
    assert code == 0
    assert read_json(repeated / "manifest.json")["parent"] == directory.name
    assert main(["analyze", str(directory)]) == 0
    assert execution.verify_artifacts(directory)["run_id"] == directory.name
    (directory / "states.csv").write_text("corruption")
    with pytest.raises(ValueError, match="checksum"):
        execution.verify_artifacts(directory)


def test_failure_retains_logs(tmp_path, monkeypatch):
    def fail(*args):
        raise RuntimeError("intentional solver failure")

    monkeypatch.setattr(execution, "matlab_process", fail)
    directory, code = execution.execute(resolve({}), tmp_path)
    assert code == 2
    assert read_json(directory / "manifest.json")["status"] == "error"
    assert "intentional solver failure" in (directory / "report.html").read_text()
    with pytest.raises(ValueError, match="did not complete"):
        execution.verify_artifacts(directory)


def test_requirement_failure_is_distinct(fake_matlab, tmp_path):
    fake_matlab["requirements"] = [
        {"id": "strict", "metric": "tracking_rmse", "operator": "le", "limit": 0}
    ]
    directory, code = execution.execute(fake_matlab, tmp_path)
    assert code == 1
    assert read_json(directory / "validation.json")["overall"] == "FAIL"
    assert execution.verify_artifacts(directory)["status"] == "complete"


def test_provenance_mismatch_rejected(fake_matlab, tmp_path):
    directory, _ = execution.execute(fake_matlab, tmp_path)
    m = read_json(directory / "manifest.json")
    m["provenance"]["python"] = "different"
    write_json(directory / "manifest.json", m)
    with pytest.raises(ValueError, match="Provenance"):
        execution.reproduce(directory, tmp_path)
    _, code = execution.reproduce(directory, tmp_path, allow_mismatch=True)
    assert code == 0


def test_campaign_resume_and_comparison(fake_matlab, tmp_path):
    raw = {"base": fake_matlab, "sweep": {"controller.horizon_steps": [2, 3]}}
    directory, code = campaigns.campaign(raw, tmp_path)
    assert code == 0
    original = read_json(directory / "campaign.json")
    assert len(original["runs"]) == 2
    assert (directory / "comparison.csv").exists()
    assert campaigns.campaign(raw, tmp_path, directory) == (directory, 0)
    assert read_json(directory / "campaign.json") == original
    run = Path(next(iter(original["runs"].values()))["path"])
    rules = [{"metric": "tracking_rmse", "direction": "lower", "absolute": 0, "relative": 0}]
    assert campaigns.compare(run, run, rules)["overall"] == "PASS"
    with pytest.raises(ValueError, match="timing"):
        campaigns.compare(run, run, [{**rules[0], "metric": "mean_solve_time_ms"}])
    with pytest.raises(ValueError, match="changed campaign"):
        campaigns.campaign({**raw, "max_runs": 3}, tmp_path, directory)


def test_campaign_error_propagates(fake_matlab, tmp_path, monkeypatch):
    original = campaigns.execute

    def fail_second(config, output):
        out, code = original(config, output)
        return out, 2 if config["controller"]["horizon_steps"] == 3 else code

    monkeypatch.setattr(campaigns, "execute", fail_second)
    _, code = campaigns.campaign(
        {"base": fake_matlab, "sweep": {"controller.horizon_steps": [2, 3]}}, tmp_path
    )
    assert code == 2


def test_report_escapes_labels(tmp_path):
    write_json(tmp_path / "config.resolved.json", {"name": "<script>bad()</script>"})
    report(tmp_path, {}, {"overall": "ERROR", "checks": []}, error="<script>bad()</script>")
    content = (tmp_path / "report.html").read_text()
    assert "<script>" not in content
    assert "&lt;script&gt;" in content


def test_cli_run_and_error(fake_matlab, tmp_path, capsys):
    config = tmp_path / "config.yaml"
    write_json(config, fake_matlab)
    assert main(["run", str(config), "--output", str(tmp_path / "runs")]) == 0
    assert "Report:" in capsys.readouterr().out
    config.write_text("unsupported: true")
    assert main(["run", str(config)]) == 2
    assert "Unknown keys" in capsys.readouterr().err


def test_missing_matlab(tmp_path, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="not found"):
        execution.matlab_process(tmp_path, tmp_path, "doctor", 1)


@pytest.mark.parametrize("result", [0, 1, "timeout"])
def test_process_boundary(tmp_path, monkeypatch, result):
    observed = {}

    class Child:
        pid = 12345
        calls = 0

        def wait(self, timeout=None):
            self.calls += 1
            if result == "timeout" and self.calls == 1:
                raise subprocess.TimeoutExpired("matlab", 1)
            return result if isinstance(result, int) else -1

    def popen(args, **kwargs):
        observed.update(kwargs)
        observed["args"] = args
        return Child()

    monkeypatch.setattr(shutil, "which", lambda _: "matlab")
    monkeypatch.setattr(subprocess, "Popen", popen)
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: None)
    # POSIX branch also uses a fake signal target, never a real process.
    monkeypatch.setattr(execution.os, "killpg", lambda *a: None, raising=False)
    if result == 0:
        execution.matlab_process(tmp_path, tmp_path, "doctor", 1)
        assert observed["args"][-1].endswith("cb_doctor()")
        assert "shell" not in observed
    else:
        with pytest.raises((RuntimeError, subprocess.TimeoutExpired)):
            execution.matlab_process(tmp_path, tmp_path, "doctor", 1)


def test_benchmark_repetitions(fake_matlab, tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "report", lambda *args, **kwargs: None)
    directory, code = bench.benchmark(fake_matlab, tmp_path, 5)
    assert code == 0
    summary = read_json(directory / "benchmark.json")
    assert len(summary["runs"]) == 6
    assert summary["metrics"]["median_mean_solve_time_ms"]["value"] == pytest.approx(15)
    assert bench.compare_benchmarks(directory, directory)["overall"] == "PASS"
    with pytest.raises(ValueError, match="five"):
        bench.benchmark(fake_matlab, tmp_path, 4)
    summary["profile"]["hostname"] = "other-host"
    other = tmp_path / "other"
    other.mkdir()
    write_json(other / "benchmark.json", summary)
    with pytest.raises(ValueError, match="profile mismatch"):
        bench.compare_benchmarks(directory, other)


def test_config_checksum(fake_matlab, tmp_path):
    directory, _ = execution.execute(fake_matlab, tmp_path)
    write_json(directory / "config.resolved.json", resolve({}))
    with pytest.raises(ValueError, match="Configuration checksum"):
        execution.verify_artifacts(directory)


def test_manifest_path_escape(fake_matlab, tmp_path):
    directory, _ = execution.execute(fake_matlab, tmp_path)
    m = read_json(directory / "manifest.json")
    m["checksums"]["../config.resolved.json"] = "invalid"
    write_json(directory / "manifest.json", m)
    with pytest.raises(ValueError, match="checksum"):
        execution.verify_artifacts(directory)
