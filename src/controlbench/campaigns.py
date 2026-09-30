"""Deterministic Cartesian campaigns and compatibility-checked comparisons."""

from __future__ import annotations

import copy
import csv
import html
import itertools
import math
import uuid
from pathlib import Path

from .config import digest, resolve
from .data import read_json, write_json
from .execution import execute, verify_artifacts
from .validation import regression


def expand(raw: dict) -> list[dict]:
    if set(raw) - {"base", "sweep", "max_runs"} or not {"base", "sweep"} <= set(raw):
        raise ValueError("Campaign needs base, sweep, optional max_runs")
    base = resolve(raw["base"])
    axes = raw["sweep"]
    if not isinstance(axes, dict) or not axes:
        raise ValueError("sweep must be a nonempty mapping")
    limit = raw.get("max_runs", 100)
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 1:
        raise ValueError("max_runs must be a positive integer")
    for values in axes.values():
        if not isinstance(values, list) or not values:
            raise ValueError("Sweep axes must be nonempty lists")
    if math.prod(len(v) for v in axes.values()) > limit:
        raise ValueError("Campaign exceeds max_runs")
    configs = []
    for values in itertools.product(*axes.values()):
        c = copy.deepcopy(base)
        for path, value in zip(axes, values, strict=True):
            keys = path.split(".")
            target = c
            for key in keys[:-1]:
                if key not in target or not isinstance(target[key], dict):
                    raise ValueError(f"Unknown sweep path: {path}")
                target = target[key]
            if keys[-1] not in target:
                raise ValueError(f"Unknown sweep path: {path}")
            target[keys[-1]] = value
        configs.append(resolve(c))
    if len({digest(c) for c in configs}) != len(configs):
        raise ValueError("Duplicate sweep configurations")
    return configs


def campaign(raw: dict, output_root: Path, resume: Path | None = None) -> tuple[Path, int]:
    configs = expand(raw)
    directory = resume or output_root.resolve() / ("campaign-" + str(uuid.uuid4()))
    if resume:
        saved = read_json(directory / "campaign.json")
        if saved["definition"] != raw:
            raise ValueError("Cannot resume a changed campaign")
    else:
        directory.mkdir(parents=True)
        saved = {"definition": raw, "runs": {}}
    for config in configs:
        key = digest(config)
        old = saved["runs"].get(key)
        if old:
            try:
                manifest = verify_artifacts(Path(old["path"]))
                if manifest["config_hash"] == key:
                    continue
            except (OSError, ValueError, KeyError):
                pass
        path, code = execute(config, directory / "runs")
        saved["runs"][key] = {"path": str(path), "exit_code": code}
        write_json(directory / "campaign.json", saved)
    comparison = []
    for key, item in saved["runs"].items():
        path = Path(item["path"])
        row = {"config_hash": key, "run": path.name, "exit_code": item["exit_code"]}
        if (path / "metrics.json").exists():
            row.update(
                {name: val["value"] for name, val in read_json(path / "metrics.json").items()}
            )
        comparison.append(row)
    columns = list(dict.fromkeys(k for row in comparison for k in row))
    with (directory / "comparison.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(comparison)
    table = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(k, '')))}</td>" for k in columns) + "</tr>"
        for row in comparison
    )
    (directory / "report.html").write_text(
        '<!doctype html><meta charset="utf-8">'
        '<title>ControlBench campaign</title><h1>Campaign comparison</h1><table border="1"><tr>'
        + "".join(f"<th>{html.escape(k)}</th>" for k in columns)
        + "</tr>"
        + table
        + "</table>",
        encoding="utf-8",
    )
    return directory, max(item["exit_code"] for item in saved["runs"].values())


def compare(baseline: Path, candidate: Path, rules: list[dict]) -> dict:
    verify_artifacts(baseline)
    verify_artifacts(candidate)
    b = read_json(baseline / "config.resolved.json")
    c = read_json(candidate / "config.resolved.json")
    for key in ("plant", "scenario", "constraints"):
        if b[key] != c[key]:
            raise ValueError(f"Incompatible {key}")
    for key in ("horizon_steps", "sample_time_s"):
        if b["controller"][key] != c["controller"][key]:
            raise ValueError(f"Incompatible {key}")
    if b["experiment"]["seed"] != c["experiment"]["seed"]:
        raise ValueError("Incompatible seeds")
    if any("time_ms" in rule["metric"] for rule in rules):
        raise ValueError("Single-run timing gates are unsupported; use repeated benchmark data")
    return regression(
        read_json(baseline / "metrics.json"), read_json(candidate / "metrics.json"), rules
    )
