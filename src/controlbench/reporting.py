"""Offline, escaped HTML reports and explicit physical-signal plots."""

from __future__ import annotations

import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .data import ExperimentResult


def report(
    directory: Path,
    metrics: dict,
    validation: dict,
    result: ExperimentResult | None = None,
    error: str = "",
) -> None:
    pictures = []
    if result is not None:
        r = result
        t = np.arange(len(r.states)) * r.config["controller"]["sample_time_s"]
        for name, data, times, ylabel, bounds in (
            ("states", r.states, t, "Level (m)", r.config["constraints"]["level_m"]),
            ("inputs", r.inputs, t[:-1], "Applied flow (m³/s)", r.config["constraints"]["input"]),
            ("solver", r.solve_times * 1000, t[:-1], "Controller solve time (ms)", None),
        ):
            fig, ax = plt.subplots(figsize=(9, 3.6), layout="constrained")
            if data.ndim == 2:
                for a in range(data.shape[1]):
                    ax.plot(times, data[:, a], label=f"Agent {a + 1}")
            else:
                ax.plot(times, data)
            if name == "states":
                ax.plot(t, r.references[:, 0], "k--", label="Reference")
            if bounds:
                ax.axhline(bounds["min"], color="red", linestyle=":", label="Bounds")
                ax.axhline(bounds["max"], color="red", linestyle=":")
            if name == "inputs" and bounds is not None:
                low, high = min(float(data.min()), bounds["min"]), float(data.max())
                padding = max(high - low, 0.01) * 0.15
                ax.set_ylim(low - padding, high + padding)
                if bounds["max"] > high + padding:
                    ax.text(
                        0.99,
                        0.95,
                        f"Upper bound: {bounds['max']:g} m³/s (outside view)",
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=8,
                    )
            if data.ndim == 2:
                ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=6, fontsize=8)
            ax.set(xlabel="Time (s)", ylabel=ylabel)
            ax.grid(alpha=0.2)
            path = directory / f"{name}.png"
            fig.savefig(path, dpi=130)
            plt.close(fig)
            pictures.append(f'<img alt="{html.escape(ylabel)}" src="{path.name}">')
        b = r.config["constraints"]["level_m"]
        margin = np.minimum(r.states - b["min"], b["max"] - r.states)
        fig, ax = plt.subplots(figsize=(9, 3.6), layout="constrained")
        for a in range(margin.shape[1]):
            ax.plot(t, margin[:, a], label=f"Agent {a + 1}")
        ax.axhline(0, color="red", linestyle=":")
        ax.set(xlabel="Time (s)", ylabel="Nearest state-bound margin (m)")
        ax.legend(ncol=4)
        fig.savefig(directory / "margins.png", dpi=130)
        plt.close(fig)
        pictures.append('<img alt="State constraint margins" src="margins.png">')
        if r.residuals:
            fig, ax = plt.subplots(figsize=(9, 3.6), layout="constrained")
            for a in range(1, r.config["plant"]["agents"] + 1):
                values = [
                    float(v["normalized_primal"]) for v in r.residuals if int(v["agent_id"]) == a
                ]
                ax.semilogy(np.maximum(values, 1e-16), label=f"Agent {a}")
            ax.set(xlabel="Accumulated ADMM iterations", ylabel="Normalized primal residual")
            ax.legend()
            fig.savefig(directory / "residuals.png", dpi=130)
            plt.close(fig)
            pictures.append('<img alt="Primal residuals" src="residuals.png">')
    esc = html.escape
    metric_rows = "".join(
        f"<tr><td>{esc(k)}</td><td>{v['value']:.6g}</td><td>{esc(v['unit'])}</td></tr>"
        for k, v in metrics.items()
    )
    checks = "".join(
        f"<tr><td>{esc(c.get('id', c.get('metric', '')))}</td>"
        f"<td>{esc(c['status'])}</td><td>{esc(json.dumps(c))}</td></tr>"
        for c in validation.get("checks", [])
    )
    config = (
        (directory / "config.resolved.json").read_text()
        if (directory / "config.resolved.json").exists()
        else "{}"
    )
    document = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>ControlBench experiment report</title><style>
body{{font:16px system-ui;max-width:1050px;margin:40px auto;padding:20px;color:#182332}}
h1{{color:#155e75}}table{{border-collapse:collapse;width:100%;margin:20px 0}}
td,th{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}pre{{white-space:pre-wrap}}
img{{max-width:100%}}.status{{font-size:24px;font-weight:bold}}
</style><h1>ControlBench experiment report</h1>
<p class="status">{esc(validation["overall"])}</p><pre>{esc(error)}</pre>
<p>Physical simulation uses simultaneous, unclipped Euler propagation.
ADMM convergence refers to the recorded normalized primal criterion, not stability.</p>
<h2>Metrics</h2><table><tr><th>Metric</th><th>Actual</th><th>Unit</th></tr>{metric_rows}</table>
<h2>Requirements</h2><table>{checks}</table>{"".join(pictures)}
<h2>Resolved configuration</h2><pre>{esc(config)}</pre>
<p><a href="manifest.json">Provenance</a> · <a href="stdout.log">MATLAB output</a> ·
<a href="stderr.log">Error log</a> · <a href="matlab.json">MATLAB/model metadata</a></p></html>"""
    (directory / "report.html").write_text(document, encoding="utf-8")
