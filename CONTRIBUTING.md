# Contributing

Use a `codex/` feature branch for new work. Keep each change reviewable and use the
issue/PR templates to explain its problem, resulting behavior, tests and limitations.
Do not invent previous releases or reviews to create a development history.

Install `.[dev]` in a virtual environment. Run Ruff, its formatting check, mypy and
pytest before a Python change is ready. Numerical MATLAB changes also require the
MATLAB unit/integration suite and a saved experiment comparison. Keep changes to
algorithm behavior separate from changes to diagnostics where practical.

Never accept a new baseline solely because current tests fail. Investigate and
document the engineering reason, preserve the old artifacts, and make the baseline
choice explicit. Never count skipped or unavailable MATLAB tests as release evidence.
Run strict performance gates only on a controlled host with repeated measurements.

AI-assisted code must pass the same review and verification as manually written
code. Distinguish synthetic pipeline fixtures from real solver results. Do not
claim stability, global nonlinear optimality, real-time performance or parallel
speedup without the corresponding evidence.

Preserve uncommitted user changes and third-party code. This project does not assign
a new license to pre-existing dependencies. Do not commit licensed binaries,
machine-specific dependency paths, environment secrets or generated result folders.
