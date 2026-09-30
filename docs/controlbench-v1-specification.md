# ControlBench v1.0 — software specification

Status: proposed implementation contract, grounded in repository inspection on 2026-09-29.
This document specifies future work; it does not claim that the platform or its
acceptance tests already exist. No MATLAB execution was performed for this specification.

## 1. Purpose and release boundary

Build an experiment, validation, and benchmarking platform around the existing
MATLAB distributed model predictive control framework. MATLAB remains responsible
for controller optimization and plant simulation; Python orchestrates runs and
processes their exported results. Do not port or replace the DMPC algorithm.

The first release answers: how do horizon, ADMM settings, coupling approximation,
and system size affect tracking, constraint satisfaction, convergence, and runtime
for the coupled water-tank system?

v1.0 includes configuration-driven MATLAB runs, a documented result contract,
analysis plugins, engineering requirements, parameter sweeps, regression comparison,
reproducibility metadata, HTML reports, automated tests, and CI. A centralized MPC
baseline uses the same plant, horizon, costs, and bounds. A double-integrator case
provides a small numerical verification problem. Water tanks are the principal
application; a vehicle platoon is not required to make this project substantial.

## 2. Existing implementation and integration findings

| Component | Observed behavior | Integration requirement |
|---|---|---|
| `Agent.m`, `Neighbor.m` | Agent models, costs, constraints, directed neighbor registration | Reuse through model factories |
| `Agent_data.m`, `Neighbor_data.m` | YALMIP variables and mutable optimizer state | Create fresh instances for each experiment |
| `ADMM_Solver.m` | Local optimization, exchanges, residuals, penalty adaptation | Add structured per-step diagnostics without silently changing the algorithm |
| `Solution.m` | Stores time, state, first applied input, and optimizer cost | Export normalized arrays; compute physical cost separately |
| `coupled_water_tank_chain.m` | Four-tank script, IPOPT, hard-coded settings and plots | Preserve as a legacy example; extract an independently callable runner |
| `tank_chain/GRAMPCD/` | Existing comparison data | Import only with documented provenance and signal definitions |

Specific findings that influence the design:

1. `solve` prints convergence and warns at the iteration limit but returns no
   structured termination result. Capture termination reason and iteration count.
2. Local YALMIP diagnostics are printed. Preserve numeric status and diagnostic
   text for every local solve; a failed local solve cannot become a successful run.
3. `check_convergence` currently checks normalized primal residuals. Export that
   exact criterion under a versioned name and also record dual residuals. Do not
   label it a primal-and-dual stopping test or a proof of closed-loop stability.
4. `shift` propagates with forward Euler and projects out-of-bound states. Export
   pre-projection values and projection events for legacy runs. Benchmark mode must
   integrate the plant separately and must not conceal violations by clipping state.
5. `shift` updates agents inside a loop. Verify simultaneous plant propagation
   from an immutable snapshot before using this path as a benchmark reference.
6. The example loops over `0:dt_sample:T_sim` and advances at each iteration.
   The new runner must execute exactly K steps, producing K+1 state samples.
7. The existing timing combines `solve` and `shift`. New measurements must separate
   optimization, plant integration, initialization, and total process time.
8. `Solution.cost` contains the recorded optimization objective, which may include
   ADMM terms. It must not be reported as comparable physical control cost.
9. The current example contains a local modification. Preserve it during migration.

These are static observations and verification targets, not experimentally confirmed
performance findings. Existing vendored directories do not prove solver availability.

## 3. Supported use cases and command contract

Commands below are planned interfaces, not currently available commands.

| Command | Outcome |
|---|---|
| `controlbench doctor` | Checks Python, MATLAB, YALMIP, IPOPT and a small solver problem |
| `controlbench run experiments/tanks/nominal.yaml` | Executes, exports, analyzes, validates, and reports one experiment |
| `controlbench analyze results/<run-id>` | Recalculates metrics from recorded data without MATLAB |
| `controlbench sweep experiments/tanks/horizon-sweep.yaml` | Runs a deterministic Cartesian campaign and produces a comparison |
| `controlbench compare results/<baseline-id> results/<candidate-id>` | Evaluates explicit numerical/performance regression rules |
| `controlbench reproduce results/<run-id>` | Checks provenance and reruns the saved resolved configuration |

Exit codes: 0 = completed and all required checks passed; 1 = completed but a
requirement/regression failed; 2 = configuration, dependency, execution, or data
error. A campaign returns 2 if any execution errors occurred, otherwise 1 if any
required check failed. A partial campaign is never reported as a complete pass.

## 4. Configuration contract

Use YAML with safe parsing, typed validation, forbidden unknown keys, and an
explicit schema version. Resolve defaults before hashing or starting MATLAB.
No configuration may contain executable MATLAB or Python expressions.

```yaml
schema_version: 1
experiment:
  name: tanks_nominal
  seed: 42
plant:
  model: coupled_water_tanks
  agents: 4
  initial_level_m: 0.5
  reference_level_m: 2.0
  coupling_model: legacy_polynomial_4
controller:
  type: matlab_dmpc
  horizon_steps: 20
  sample_time_s: 0.1
  optimizer: ipopt
  admm:
    max_iterations: 120
    convergence_tolerance: 0.001
    initial_penalty: 20.0
    adaptive_penalty: false
    approximation: none
scenario:
  duration_s: 3.0
  integrator: forward_euler
constraints:
  level_m: {min: 0.0, max: 3.0}
  input: {min: 0.0, max: 1000.0}
requirements:
  - id: tank_bounds
    metric: state_violation_count
    operator: le
    limit: 0
  - id: successful_local_solves
    metric: failed_local_solves
    operator: le
    limit: 0
  - id: admm_termination
    metric: converged_step_fraction
    operator: ge
    limit: 1.0
execution:
  timeout_s: 600
```

This is a proposed schema, not a certified nominal case. The state/input bounds
and settings come from the example; a passing result is not assumed. Input units,
actuator placement, tank areas, coupling coefficients, per-agent cost weights,
and directed edges must appear explicitly in the resolved model definition.
The factory must not silently reproduce contradictory comments from the script.

Horizon steps H means H control intervals and H+1 state nodes; map to existing
`N=H+1` and `T=H*sample_time_s`. Initially require prediction and sampling intervals
to match. Duration must be a positive integral multiple of sampling time. Reject
invalid dimensions, nonfinite values, reversed bounds, nonexistent agents, invalid
edges, unsupported approximation modes, and requirements for unknown metrics.

Model factories support the exact legacy four-tank arrangement and an explicitly
defined N-tank chain. Generalizing the graph must preserve documented endpoint
actuation, symmetric pipe physics, and cost allocation; it is not merely changing
an agent-count field. Archive the complete resolved graph in every run.

## 5. Architecture and interfaces

```text
YAML -> validated/resolved ExperimentConfig -> LocalMatlabRunner
                                                 |
                                      MATLAB batch entry point
                                                 |
                                  model + controller + simulator
                                                 |
                               manifest.json + normalized CSV files
                                                 |
                        Python ExperimentResult -> analysis registry
                                                 |
                              metrics -> requirement/regression engine
                                                 |
                                    JSON + CSV + HTML/plots
```

Python modules: `config`, `execution`, `data`, `analysis`, `validation`,
`campaigns`, `reporting`, and `cli`. Keep these responsibilities separate without
creating empty modules for speculative features.

MATLAB additions live under `matlab/controlbench/`: batch entry point, tank and
double-integrator factories, controller adapters, simulator, and exporter.
Existing classes remain in their current locations.

Interfaces:

- Runner: `run(config, output_directory) -> ArtifactPaths`.
- Controller adapter: initialize; compute control from state/reference and time;
  return input plus diagnostics; update warm-start state. DMPC wraps existing classes.
- Plant: compute next state from all current states and applied inputs, independent
  of controller prediction variables. Snapshot all agents before advancing any agent.
- Analyzer: stable identifier and version, required signals, and
  `evaluate(result) -> mapping of metric names to values/units/status`.
- Requirement evaluator: consumes metrics and declarative rules only; never launches
  simulations or changes results to meet a threshold.

Use `matlab -batch` through a subprocess argument list with a fixed entry point;
pass paths through environment variables or safely encoded arguments, never inject
configuration strings into executable code. Capture stdout/stderr, enforce timeout,
and terminate the child process tree on cancellation. MATLAB Engine is deferred.
Dependency paths come from local settings/environment, not committed absolute paths.

## 6. Data contract and artifact lifecycle

Use JSON metadata and UTF-8 CSV for portable v1.0 interchange. MATLAB `.mat` files
may be supplemental raw diagnostics, not the sole cross-language contract.

| Artifact | Required contents |
|---|---|
| `manifest.json` | Schema, UUID, status, timestamps, resolved configuration hash, artifact checksums |
| `config.resolved.json` | All defaults, units, graph, model coefficients, costs, requirements |
| `states.csv` | step, time_s, agent_id, state_name, value, reference, lower/upper bound |
| `inputs.csv` | step, time_s, agent_id, input_name, commanded, applied, lower/upper bound |
| `steps.csv` | step, solve_time_s, plant_time_s, ADMM iterations, termination, convergence |
| `local_solves.csv` | step, ADMM iteration, agent_id, diagnostic code, duration and message |
| `residuals.csv` | step, iteration, agent_id, raw/normalized primal and dual residuals |
| `events.csv` | projection, solver failure, cancellation or other exceptional events |
| `metrics.json`, `validation.json` | Versioned metrics and individual requirement outcomes |
| `report.html`, `plots/` | Offline-readable report and figures |
| `stdout.log`, `stderr.log` | Execution diagnostics |

Canonical state tensors are `(K+1, agents, state_components)` and input tensors
`(K, agents, input_components)` for homogeneous models. References match states.
Tables retain named components and explicit agent IDs. Missing values are errors for
required signals; never coerce NaN, absent diagnostics, or missing intervals to zero.
Validate monotonic times, exact sample counts, unique keys, units, finite values,
and consistency with the resolved configuration before computing metrics.

Create a UUID output directory atomically, mark it running, and finalize the manifest
only after exports validate. Failed runs retain diagnostics and partial artifacts,
with failure status. Never overwrite a previous run. Hash inputs and outputs; record
Git SHA, dirty status and source snapshot/hash, Python/MATLAB/toolbox/solver versions,
OS, machine/CPU, random seed, thread settings, and controller/model revisions.
Dirty code must be identifiable; a commit SHA alone is insufficient provenance.

Reproduce checks source/configuration/dependency provenance before running. Mismatches
are reported; an explicit allow-mismatch option creates a linked new run. Numerical
equivalence uses tolerances; bitwise equivalence across machines is not promised.

## 7. Analysis and validation semantics

Built-in analyzers register through an explicit registry. Optional installed plugins
may register through a documented Python entry-point group; never load arbitrary
files named in a YAML configuration. Reject duplicate metric identifiers.

| Metric family | Definition |
|---|---|
| Tracking | Per-agent/component RMSE and max absolute error against saved references; include initial and terminal state |
| Constraints | Count scalar samples outside bounds beyond a declared tolerance; report maximum violation magnitude and affected steps |
| Effort | Sum of `dt * applied_input^2`, per actuator; keep units explicit |
| Physical cost | Common configured stage/terminal cost, independently recomputed; exclude ADMM augmentation |
| Solver | Mean, p95 and maximum measured whole-controller solve time; failed local solves; iterations; converged-step fraction |
| Residuals | Final and peak residuals per agent/step, with normalization and stopping-rule identifier |

Do not aggregate heterogeneous units into one RMSE. For homogeneous tanks, publish
both per-agent values and a clearly defined pooled RMSE. Report all measured steps;
warm-up exclusions belong to an explicit separate benchmark policy. Quantile method
is fixed and versioned. Numerical tolerances are independent of engineering limits.

Requirements have unique IDs, metric keys, comparison operators (`le`, `ge`, `eq`),
limits, units, and optional selectors. Outcomes are PASS, FAIL, or ERROR. Missing
metrics, invalid units, analyzer exceptions, and nonfinite results are ERROR and block
an overall pass. An experiment with no requirements is NOT_EVALUATED, not PASS.

Legacy projected trajectories cannot establish satisfaction of physical state bounds.
For them, required pre-projection signals must exist or constraint validation is ERROR.
No finite simulation, bounded trace, or converged optimizer is a stability proof.
Communication latency is not measured by in-process neighbor method calls.

## 8. Campaigns, baselines, and regression

Sweep only validated dotted configuration paths over explicit finite lists. Expand
the Cartesian product deterministically; enforce a configurable maximum run count
before launch. Store campaign configuration, run IDs, seeds, and per-run status.
Default execution is sequential to avoid MATLAB license and CPU oversubscription.
Resume only runs whose completed artifacts and configuration hashes validate.

Initial campaigns: horizon 5/10/20; ADMM penalty and tolerance; approximation modes;
and 2/4/8 tanks after graph generalization passes tests. Larger agent counts require
measured resource budgets. Never call sequential local optimization parallel speedup.

Centralized MPC and DMPC comparisons require identical dynamics, integration,
constraints, objective, initial conditions, reference, and sample/horizon grids.
Label differences explicitly and reject incompatible automated comparisons.

Baselines are immutable result artifacts chosen explicitly, never automatically
replaced following a failure. Each metric declares degradation direction, absolute
allowance and relative allowance. For a lower-is-better metric, failure occurs when
`candidate > baseline + absolute_allowance + relative_allowance*abs(baseline)`;
reverse the direction for higher-is-better metrics. A zero baseline uses absolute
allowance; never divide by zero. Missing metrics invalidate the comparison.

Timing regressions require the same recorded machine/runtime profile, warm-up policy,
and at least five repetitions. Report median and spread across repetitions. Ordinary
shared CI runners check numerical behavior, not strict wall-clock performance gates.

## 9. Reports

Produce standalone HTML with local/static plots: state/reference/bounds, applied
inputs/bounds, solve times, residual convergence, violations, and campaign comparisons.
Include configuration, provenance, requirement limits and actual values, execution
status, and links to logs. Escape all user-supplied text. Failures must remain visible
even if a plot cannot be produced. Display measured values only; the illustrative
performance tables in the initial brief are not project results.

## 10. Verification and release acceptance

| ID | Acceptance test |
|---|---|
| A01 | Unknown configuration keys, invalid horizons/durations, invalid bounds and missing dependencies fail before a long run |
| A02 | A 3-second, 0.1-second run yields exactly 30 inputs and 31 state samples per tank |
| A03 | Round-trip MATLAB export/Python import preserves values, IDs, dimensions, units and references |
| A04 | Independent hand-calculated fixtures verify RMSE, effort, boundary tolerance, violation count, and timing quantiles |
| A05 | Deliberately infeasible local optimization produces an execution failure with diagnostics; no PASS report |
| A06 | ADMM iteration exhaustion is distinct from local solver failure and is reflected in convergence requirements |
| A07 | An injected out-of-bound state is detected before any legacy projection; absent raw values cause ERROR |
| A08 | Simultaneous plant propagation is invariant to agent iteration order within numerical tolerance |
| A09 | A double-integrator constant-input trajectory agrees with its specified discretization; equilibrium remains stationary |
| A10 | A nominal four-tank DMPC run finishes on a licensed MATLAB host and produces all contract artifacts |
| A11 | Centralized and distributed formulations agree on a small convex case within declared tolerances after convergence |
| A12 | A 2-by-3 sweep produces six unique configs/runs; a failed child remains visible and affects campaign exit status |
| A13 | Regression fixtures detect both increased error and runtime degradation, including zero-baseline and missing-metric cases |
| A14 | Reproduction preserves resolved configuration and seeds; incompatible provenance is rejected unless explicitly allowed |
| A15 | Reports open offline; malformed exports, interrupted processes and paths containing spaces are covered |
| A16 | Python checks and licensed MATLAB integration checks both pass for the release commit; skipped MATLAB tests cannot satisfy release readiness |

Use pytest and coverage for Python, Ruff for lint/format checks, mypy for typed domain
interfaces, and MATLAB unit tests for factories/export/controller integration.
Use analytically justified or independently computed fixtures, not only snapshots
generated by the implementation under test. Freeze numerical tolerances before
accepting a baseline. Coverage target: at least 85% of Python domain/validation code,
with explicit tests for every failure path above; coverage alone is not correctness.

CI tiers: PRs run Python checks and portable fixtures; trusted licensed runners run
short MATLAB integration tests. Main-branch/manual runs execute broader numerical
campaigns; controlled-host manual runs assess timing. Do not execute untrusted fork
code on a privileged self-hosted MATLAB runner. Archive test/report artifacts.

Python dependencies: NumPy, pandas, PyYAML, Pydantic, Matplotlib, Jinja2, and a CLI
library; development dependencies as above. Select and lock tested versions during
implementation. MATLAB uses existing YALMIP/IPOPT integration. Document supported
MATLAB release and solver versions after an actual environment check. A Python-only
container may analyze saved results; it does not provide a licensed MATLAB solver.

## 11. Delivery sequence

1. Repository/environment audit and baseline capture: doctor command design,
   dependency inventory, legacy-run evidence, confirmed model equations/units.
2. MATLAB batch runner and diagnostics: configuration factory, simultaneous plant
   simulation, exact time grid, explicit termination, portable export; A02/A03/A05–A09.
3. Python ingestion and analysis: typed contract, registry, metrics and requirements;
   A01/A04/A07 with synthetic and MATLAB-produced fixtures.
4. CLI and reporting: single-run lifecycle, timeout/error handling, provenance,
   offline HTML; A10/A14/A15.
5. Campaigns and regression: sweep/resume, comparison compatibility, immutable
   baselines and controlled timing; A12/A13.
6. Centralized baseline and release verification: convex reference case, fair tank
   comparisons, CI, contribution guide, real reports; A11/A16.

Each increment should be a reviewable change with tests and limitations. Add issue
and PR templates locally, then use real issues/branches/PRs as development occurs.
Do not fabricate historical releases, reviews, or benchmark results. Do not choose
a license for existing or vendored code without checking ownership and permissions.

## 12. Deferred scope

Later releases may add sensor noise/bias, disturbances and actuator faults, with
separate true/measured/commanded/applied signals and deterministic seeds. Communication
delay/packet loss requires an explicit exchange transport abstraction and defined
stale-data policy; it cannot be implemented honestly by adding a delay field alone.

Vehicle platoons, pendulums, LQR comparisons, live playback, distributed worker
execution, Azure deployment, and a fully containerized MATLAB runtime are deferred.
Cloud infrastructure and paid resources are not prerequisites for v1.0. Stability
certification and safety certification are outside scope.

## 13. Decisions to resolve through implementation evidence

- Establish actual MATLAB/YALMIP/IPOPT availability and a tested version matrix.
- Verify tank input units, intended endpoint actuation, polynomial validity range,
  and whether the prediction model and physical plant should intentionally differ.
- Capture a legacy baseline before altering solver behavior; isolate instrumentation
  changes from numerical corrections and quantify resulting differences.
- Calibrate tracking and timing requirements using measured baselines. Never select
  a limit simply to make a report green; document the engineering rationale.
- Establish a licensed CI host before declaring complete release verification.

These checks can proceed without inventing user preferences or rewriting the solver.
