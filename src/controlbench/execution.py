"""MATLAB process boundary, immutable artifacts, and reproducibility metadata."""

from __future__ import annotations

import importlib.metadata
import os
import platform
import shutil
import signal
import subprocess
import sys
import uuid
from datetime import UTC, datetime
from pathlib import Path

from .analysis import analyze
from .config import digest, resolve
from .data import ExperimentResult, checksum, read_json, write_json
from .reporting import report
from .validation import evaluate


def project_root() -> Path:
    root = Path(os.environ.get("CONTROLBENCH_ROOT", Path(__file__).resolve().parents[2]))
    if not (root / "ADMM_Solver.m").is_file():
        raise ValueError("Set CONTROLBENCH_ROOT to the DMPC checkout")
    return root.resolve()


def source_files(root: Path) -> list[Path]:
    return sorted(
        [
            *root.glob("*.m"),
            *root.glob("matlab/controlbench/*.m"),
            *root.glob("src/controlbench/*.py"),
            root / "pyproject.toml",
        ]
    )


def provenance(root: Path) -> dict:
    def git(*args: str) -> str:
        p = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=False)
        return p.stdout.strip() if p.returncode == 0 else "unavailable"

    return {
        "git_sha": git("rev-parse", "HEAD"),
        "git_status": git("status", "--short"),
        "sources": {p.relative_to(root).as_posix(): checksum(p) for p in source_files(root)},
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "hostname": platform.node(),
        "dependencies": {
            name: importlib.metadata.version(name) for name in ("numpy", "PyYAML", "matplotlib")
        },
        "threads": "MATLAB -singleCompThread; sequential agent solves",
    }


def matlab_process(
    root: Path, out: Path, entry: str, timeout: float, config_path: Path | None = None
) -> None:
    executable = shutil.which(os.environ.get("CONTROLBENCH_MATLAB", "matlab"))
    if executable is None:
        raise RuntimeError("MATLAB executable not found; set CONTROLBENCH_MATLAB")
    env = os.environ.copy()
    # Execute the saved MATLAB sources, so edits in the checkout cannot change
    # an in-flight experiment after provenance was captured.
    saved_root = out / "source"
    run_root = saved_root if (saved_root / "ADMM_Solver.m").is_file() else root
    env.setdefault("CONTROLBENCH_YALMIP", str(root / "YALMIP-master"))
    env.setdefault("CONTROLBENCH_OPTI", str(root / "OPTI-master"))
    env["CONTROLBENCH_MATLAB_DIR"] = str(run_root / "matlab" / "controlbench")
    env["CONTROLBENCH_OUTPUT"] = str(out.resolve())
    env["CONTROLBENCH_CONFIG"] = str(config_path.resolve()) if config_path else ""
    expressions = {
        "doctor": "cb_doctor()",
        "run": "cb_run(getenv('CONTROLBENCH_CONFIG'),getenv('CONTROLBENCH_OUTPUT'))",
    }
    command = "addpath(getenv('CONTROLBENCH_MATLAB_DIR'));" + expressions[entry]
    with (
        (out / "stdout.log").open("w", encoding="utf-8") as stdout,
        (out / "stderr.log").open("w", encoding="utf-8") as stderr,
    ):
        kwargs = (
            {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}  # type: ignore[attr-defined]
            if os.name == "nt"
            else {"start_new_session": True}
        )
        process = subprocess.Popen(
            [executable, "-singleCompThread", "-batch", command],
            cwd=run_root,
            env=env,
            stdout=stdout,
            stderr=stderr,
            **kwargs,
        )
        try:
            code = process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            if os.name == "nt":
                subprocess.run(
                    ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                    capture_output=True,
                    check=False,
                )
            else:
                # These POSIX-only symbols are absent from Windows type stubs.
                os.killpg(process.pid, signal.SIGKILL)  # type: ignore[attr-defined]
            process.wait()
            raise
    if code:
        raise RuntimeError(f"MATLAB exited with code {code}; inspect {out / 'stderr.log'}")


def verify_artifacts(directory: Path) -> dict:
    m = read_json(directory / "manifest.json")
    if m["status"] != "complete":
        raise ValueError("Run did not complete")
    if m["config_hash"] != digest(read_json(directory / "config.resolved.json")):
        raise ValueError("Configuration checksum mismatch")
    required = {
        "config.resolved.json",
        "states.csv",
        "inputs.csv",
        "steps.csv",
        "local_solves.csv",
        "residuals.csv",
        "events.csv",
        "matlab.json",
        "metrics.json",
        "validation.json",
    }
    if not required <= m["checksums"].keys():
        raise ValueError("Manifest missing required artifact checksums")
    for name, expected in m["checksums"].items():
        p = (directory / name).resolve()
        if not p.is_relative_to(directory.resolve()) or checksum(p) != expected:
            raise ValueError(f"Artifact checksum mismatch: {name}")
    return m


def evaluate_run(directory: Path) -> dict:
    result = ExperimentResult.load(directory)
    metrics = analyze(result)
    validation = evaluate(metrics, result.config["requirements"])
    write_json(directory / "metrics.json", metrics)
    write_json(directory / "validation.json", validation)
    report(directory, metrics, validation, result)
    return validation


def execute(
    config: dict, output_root: Path, root: Path | None = None, parent: str | None = None
) -> tuple[Path, int]:
    c = resolve(config)
    root = root or project_root()
    identifier = str(uuid.uuid4())
    out = output_root.resolve() / identifier
    out.mkdir(parents=True, exist_ok=False)
    write_json(out / "config.resolved.json", c)
    meta: dict = {
        "schema_version": 1,
        "run_id": identifier,
        "status": "running",
        "started_at": datetime.now(UTC).isoformat(),
        "config_hash": digest(c),
        "parent": parent,
        "provenance": provenance(root),
    }
    snapshot = out / "source"
    for p in source_files(root):
        target = snapshot / p.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
    write_json(out / "manifest.json", meta)
    try:
        matlab_process(root, out, "run", c["execution"]["timeout_s"], out / "config.resolved.json")
        validation = evaluate_run(out)
        meta["status"] = "complete"
        meta["validation"] = validation["overall"]
        code = (
            1 if validation["overall"] == "FAIL" else 2 if validation["overall"] == "ERROR" else 0
        )
    except (Exception, KeyboardInterrupt) as exc:
        meta["status"] = "error"
        meta["error"] = f"{type(exc).__name__}: {exc}"
        validation = {"overall": "ERROR", "checks": []}
        write_json(out / "validation.json", validation)
        report(out, {}, validation, error=meta["error"])
        code = 2
    meta["finished_at"] = datetime.now(UTC).isoformat()
    meta["checksums"] = {
        p.relative_to(out).as_posix(): checksum(p)
        for p in out.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }
    write_json(out / "manifest.json", meta)
    return out, code


def reproduce(directory: Path, output_root: Path, allow_mismatch: bool = False) -> tuple[Path, int]:
    previous = verify_artifacts(directory)
    root = project_root()
    current = provenance(root)
    fields = ("sources", "dependencies", "python", "platform")
    mismatch = [f for f in fields if previous["provenance"][f] != current[f]]
    if mismatch and not allow_mismatch:
        raise ValueError(f"Provenance differs: {mismatch}; use --allow-mismatch for a linked run")
    preflight = output_root.resolve() / ("preflight-" + str(uuid.uuid4()))
    preflight.mkdir(parents=True)
    matlab_process(root, preflight, "doctor", 120)
    old = read_json(directory / "matlab.json")
    fresh = read_json(preflight / "matlab.json")
    runtime_fields = ("matlab", "yalmip", "toolboxes", "ipopt_sha256")
    changed = [k for k in runtime_fields if old.get(k) != fresh.get(k)]
    if changed and not allow_mismatch:
        raise ValueError(f"MATLAB provenance differs: {changed}; use --allow-mismatch")
    out, code = execute(
        read_json(directory / "config.resolved.json"), output_root, root, previous["run_id"]
    )
    if code != 2:
        old = read_json(directory / "matlab.json")
        new = read_json(out / "matlab.json")
        changed = [k for k in runtime_fields if old.get(k) != new.get(k)]
        if changed and not allow_mismatch:
            m = read_json(out / "manifest.json")
            m["status"] = "error"
            m["error"] = f"MATLAB provenance mismatch: {changed}"
            v = {"overall": "ERROR", "checks": []}
            write_json(out / "validation.json", v)
            report(out, {}, v, error=m["error"])
            m["checksums"] = {
                p.relative_to(out).as_posix(): checksum(p)
                for p in out.rglob("*")
                if p.is_file() and p.name != "manifest.json"
            }
            write_json(out / "manifest.json", m)
            return out, 2
    return out, code
