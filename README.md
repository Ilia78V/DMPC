# ControlBench for MATLAB DMPC

Run repeatable control experiments, check engineering requirements, and compare
results using the existing MATLAB `ADMM_Solver` framework. Python coordinates runs
and evaluates exported signals; MATLAB/YALMIP/IPOPT perform optimization.

## Quick start on Windows

From this repository in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\controlbench.exe doctor
.\.venv\Scripts\controlbench.exe run experiments/tanks/smoke.yaml
```

Each run prints its result directory and HTML report path. Open `report.html` to
see metrics, plots, requirements, and provenance. The smoke case is intentionally
short; `experiments/tanks/nominal.yaml` runs the complete 30-step tank scenario.

Python 3.12 or newer is required. MATLAB must be licensed and on PATH. Tested locally
with MATLAB R2026a Update 4 and the existing IPOPT/YALMIP installation. Dependency
installation is confined to `.venv`; MATLAB settings are not saved globally.

If dependencies are elsewhere, set these environment variables before running:

```powershell
$env:CONTROLBENCH_MATLAB = 'D:\Program Files\MATLAB\R2026a\bin\matlab.exe'
$env:CONTROLBENCH_YALMIP = 'D:\tools\YALMIP'
$env:CONTROLBENCH_OPTI = 'D:\tools\OPTI'
```

`CONTROLBENCH_ROOT` identifies the DMPC checkout for an installed package outside
the repository. By default, the runner uses this checkout's `YALMIP-master` and
`OPTI-master`. Those ignored third-party directories are not included in a clone;
install compatible dependencies separately. `doctor` solves an actual IPOPT problem.

## Available workflows

```powershell
controlbench run experiments/tanks/nominal.yaml
controlbench run experiments/tanks/centralized.yaml
controlbench sweep experiments/tanks/horizon-sweep.yaml
controlbench sweep experiments/tanks/horizon-sweep.yaml --resume results/<campaign-id>
controlbench analyze results/<run-id>
controlbench compare results/<baseline-id> results/<candidate-id> --rules experiments/regression-rules.yaml
controlbench reproduce results/<run-id>
controlbench benchmark experiments/tanks/smoke.yaml --repetitions 5
controlbench compare-benchmarks results/<baseline-benchmark> results/<candidate-benchmark>
```

Activate `.venv` or use the executable's full path as in Quick start. Replace angle
bracket placeholders with actual run IDs. `analyze` creates a new analysis directory;
it does not overwrite the original. Reproduction refuses changed Python/source
provenance unless `--allow-mismatch` is supplied. MATLAB versions are also checked.
Reproduction reruns a configuration; it does not install historical dependencies.

Exit codes: **0** means execution completed without failed requirements; **1** means
an engineering requirement failed; **2** means an execution/configuration/data error.
No requirements produces **NOT_EVALUATED**, never PASS. A failed local solver halts
the run and retains its diagnostics. ADMM iteration exhaustion remains a recorded
nonconverged step; the nominal configuration explicitly requires convergence.

## What is implemented

- Strict YAML/JSON configuration with resolved defaults and deterministic seeds.
- Existing DMPC controller and a centralized MPC baseline using the same Euler
  dynamics, costs, bounds, horizon, and sampling grid.
- Scalar coupled water tanks, generalized N-tank chains, and a convex linear-chain
  verification case. Exact double-integrator propagation is an independent MATLAB
  unit-test model, not a separate CLI controller family.
- Simultaneous, unclipped plant propagation independent of optimizer predictions.
- CSV signals, per-local-solve diagnostics, per-iteration residuals, checksums,
  environment metadata, and actual source snapshots including uncommitted code.
- Registered tracking, constraint, effort, solver, and residual analyzers; extension
  through installed Python entry points.
- Requirement checks, immutable baseline comparisons, sweeps with resume, repeated
  timing benchmarks, and offline HTML/plots.
- Python tests/lint/type checks, MATLAB tests, CI definitions and contribution templates.
- An Azure DevOps pipeline that runs the solver-backed smoke simulation, publishes
  only verified result bundles, and uploads them to private Azure Blob Storage.

## Engineering interpretation

The tank model deliberately retains the example's cubic coupling polynomial and
sixfold gain on the final pipe. Endpoint agents receive inflow; intermediate agents
have zero direct actuation. `area_m2` maps input flow in m³/s to level rate in m/s.
These are the platform's explicit model conventions; the original script's comments
are not treated as a separately validated physical installation specification.
State/input weights are resolved and saved. Unlike the original example, every
agent's unused input is penalized to avoid unconstrained dummy inputs in comparisons.

Time arrays contain K+1 states and K applied inputs. State bounds are checked on the
physical, unprojected trajectory. The original example still uses its legacy shift
method; ControlBench does not silently convert its historical output into validated
results. No stability or safety certification is inferred from a passing simulation.

ADMM convergence uses the existing normalized-primal stopping criterion. Dual
residuals are exported separately. Nonlinear IPOPT solutions need not be global
optima. Centralized/DMPC equality is tested on a small convex linear-chain case;
nonlinear water-tank results are empirical comparisons.

Solve times are measured in **milliseconds**, excluding MATLAB startup and plant
integration. This implementation solves agents sequentially on one host, so it makes
no parallel-speedup or real-time claim. Timing gates require at least five measured
repetitions after a discarded warm-up run and compatible host/runtime metadata.
Do not compare strict wall-clock gates across shared CI hosts.

## Verification

```powershell
.\.venv\Scripts\ruff.exe check src tests
.\.venv\Scripts\ruff.exe format --check src tests
.\.venv\Scripts\mypy.exe src/controlbench
.\.venv\Scripts\python.exe -m pytest --cov=controlbench --cov-fail-under=85
matlab -singleCompThread -batch "r=runtests('matlab/tests'); assertSuccess(r);"
```

Python fixtures are clearly synthetic and test the pipeline without requiring
MATLAB; they are not controller performance evidence. MATLAB integration tests use
real IPOPT solves, including infeasibility, iteration limits, and convex controller
agreement. See [verification notes](docs/verification.md) for local evidence.

The Azure DevOps MATLAB stage uses an owner-configured `ControlBench` self-hosted
runner with the licensed solver stack already installed. See
[Azure DevOps CI and Blob publication](docs/azure-devops.md) for the trust boundary,
least-privilege service connection, and artifact layout. The GitHub MATLAB workflow
remains a manual alternative.
The Docker image supports analysis of mounted result folders only; it does not
contain MATLAB or grant solver/dependency licenses.

## Project documents

- [Specification and acceptance criteria](docs/controlbench-v1-specification.md)
- [Architecture, data contract and plugins](docs/architecture.md)
- [Contribution workflow](CONTRIBUTING.md)

Cloud execution, communication faults, sensor/actuator fault injection, vehicle
platoons, and live playback remain deferred as specified. No new license is assigned
to the pre-existing or vendored code by this implementation.
