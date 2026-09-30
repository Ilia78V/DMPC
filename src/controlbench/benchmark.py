"""Repeated timing experiments on one controlled host; never single-run timing gates."""

from __future__ import annotations

import uuid
from pathlib import Path

import numpy as np

from .config import digest, resolve
from .data import read_json, write_json
from .execution import execute, project_root, provenance, verify_artifacts
from .reporting import report
from .validation import regression


def benchmark(config: dict, output_root: Path, repetitions: int = 5) -> tuple[Path, int]:
    if repetitions < 5:
        raise ValueError("Timing benchmarks require at least five measured repetitions")
    c = resolve(config)
    directory = output_root.resolve() / ("benchmark-" + str(uuid.uuid4()))
    directory.mkdir(parents=True)
    profile = provenance(project_root())
    summary: dict = {
        "config": c,
        "profile": profile,
        "policy": "one discarded complete warm-up run; fresh MATLAB process per run",
        "runs": [],
        "status": "running",
    }
    values = []
    codes = []
    for index in range(repetitions + 1):
        path, code = execute(c, directory / "runs")
        codes.append(code)
        summary["runs"].append({"path": str(path), "warmup": index == 0, "exit_code": code})
        write_json(directory / "benchmark.json", summary)
        if code == 2:
            summary["status"] = "error"
            write_json(directory / "benchmark.json", summary)
            return directory, 2
        if index > 0:
            values.append(read_json(path / "metrics.json")["mean_solve_time_ms"]["value"])
    summary["status"] = "complete"
    summary["metrics"] = {
        "median_mean_solve_time_ms": {
            "value": float(np.median(values)),
            "unit": "ms",
            "version": "1",
        }
    }
    summary["spread_ms"] = {
        "min": min(values),
        "max": max(values),
        "p25": float(np.quantile(values, 0.25)),
        "p75": float(np.quantile(values, 0.75)),
    }
    write_json(directory / "benchmark.json", summary)
    write_json(directory / "config.resolved.json", c)
    report(
        directory,
        summary["metrics"],
        {
            "overall": "FAIL" if max(codes) else "PASS" if c["requirements"] else "NOT_EVALUATED",
            "checks": [],
        },
        error=f"Five or more measured repetitions. Spread (ms): {summary['spread_ms']}",
    )
    return directory, max(codes)


def compare_benchmarks(baseline: Path, candidate: Path, relative: float = 0.15) -> dict:
    b, c = (read_json(p / "benchmark.json") for p in (baseline, candidate))
    for data in (b, c):
        measured = [r for r in data["runs"] if not r["warmup"]]
        if data["status"] != "complete" or len(measured) < 5:
            raise ValueError("Incomplete timing benchmark")
        if len({run["path"] for run in data["runs"]}) != len(data["runs"]):
            raise ValueError("Repeated run IDs are not independent timing repetitions")
        values = []
        runtime = None
        for run in data["runs"]:
            manifest = verify_artifacts(Path(run["path"]))
            if manifest["config_hash"] != digest(data["config"]):
                raise ValueError("Timing child configuration mismatch")
            for field in (
                "platform",
                "processor",
                "machine",
                "hostname",
                "dependencies",
                "threads",
            ):
                if manifest["provenance"][field] != data["profile"][field]:
                    raise ValueError(f"Timing profile mismatch: {field}")
            actual_runtime = read_json(Path(run["path"]) / "matlab.json")
            identity = {
                k: actual_runtime.get(k) for k in ("matlab", "yalmip", "toolboxes", "ipopt_sha256")
            }
            if runtime is not None and identity != runtime:
                raise ValueError("Runtime changed within timing repetitions")
            runtime = identity
            if not run["warmup"]:
                values.append(
                    read_json(Path(run["path"]) / "metrics.json")["mean_solve_time_ms"]["value"]
                )
        if data["metrics"]["median_mean_solve_time_ms"]["value"] != float(np.median(values)):
            raise ValueError("Timing summary differs from measured artifacts")
    for field in ("platform", "processor", "machine", "hostname", "dependencies", "threads"):
        if b["profile"][field] != c["profile"][field]:
            raise ValueError(f"Timing profile mismatch: {field}")
    if b["config"] != c["config"] or b["policy"] != c["policy"]:
        raise ValueError("Timing configuration/policy mismatch")
    for field in ("matlab", "yalmip", "toolboxes", "ipopt_sha256"):
        bm = read_json(Path(b["runs"][0]["path"]) / "matlab.json")
        cm = read_json(Path(c["runs"][0]["path"]) / "matlab.json")
        if bm.get(field) != cm.get(field):
            raise ValueError(f"Timing MATLAB profile mismatch: {field}")
    return regression(
        b["metrics"],
        c["metrics"],
        [
            {
                "metric": "median_mean_solve_time_ms",
                "direction": "lower",
                "absolute": 0,
                "relative": relative,
            }
        ],
    )
