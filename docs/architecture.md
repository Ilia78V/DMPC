# Architecture and data contract

## Execution boundary

`config.py` resolves defaults and checks types, dimensions, limits, unique requirement
IDs, known metrics, and supported methods. `execution.py` launches a fixed MATLAB
entry point in a child process. Configuration and output paths pass through
environment variables; YAML cannot contain executable model expressions.

`cb_build_dmpc` constructs the existing `Agent`, `Agent_data`, `Neighbor`, and
`Neighbor_data` objects. The added optional `ControlBenchDiagnostics` handle observes
each local solve and residual update in `ADMM_Solver`. It also converts local solver
failure into an explicit batch error. With no observer, the legacy error-printing
behavior remains. The approximation external-influence sum now starts at numerical
zero; previously a free YALMIP variable allowed an unmodeled contribution.

The batch runner bypasses legacy `shift`: it applies the first controls to a separate
plant using a single snapshot of all current states, then shifts optimization state.
No state projection is performed. Initialization and warm starts remain fresh for
each experiment; only successive steps within a run share optimizer state.
The fixed initial-state bound admits the configured numerical validation tolerance;
future prediction bounds do not. IPOPT bound relaxation is disabled in ControlBench
runs. This avoids treating boundary roundoff as a new infeasible physical state.
MATLAB executes the saved source snapshot, not mutable files in the live checkout.

## Result files (schema 1)

- `states.csv`: step, time_s, agent_id, value, reference, lower, upper, state_name, unit.
- `inputs.csv`: step, time_s, agent_id, commanded, applied, lower, upper, input_name, unit.
- `steps.csv`: step, solve_time_s, plant_time_s, iterations, converged, termination,
  stopping_criterion.
- `local_solves.csv`: step, iteration, agent_id, problem, duration_s, message.
- `residuals.csv`: step, iteration, agent_id, primal, dual, normalized_primal.
- `events.csv`: step, type, message. Exceptional events prevent successful analysis.
- `matlab.json`: runtime/toolboxes, solver location, explicit graph, coefficients,
  actuation and cost weights; independent simulation convention.
- `manifest.json`: run identity, lifecycle status, configuration digest, checksums,
  source/dependency/machine provenance, parent run and timestamps.
- `source/`: the files hashed at launch, including dirty working-copy files.

Only homogeneous scalar agents are supported in this release. Accordingly Python
arrays are `(K+1, agents)` for states/references and `(K, agents)` for inputs; the
singleton component axis from the proposed general tensor contract is omitted.
Named CSV components preserve a future migration path. Arbitrary multidimensional
agents are not silently flattened. Time and agent keys must be complete and unique.

Inputs use zero-order hold; physical cost sums rectangle-rule stage cost and one
terminal cost, independently of augmented ADMM objectives. Constraint counts are
scalar samples violating bounds beyond the configured absolute tolerance; both
affected steps and largest raw violation are also reported. Tracking includes the
initial sample and terminal sample. Quantiles use NumPy's linear interpolation.

The current scalar models share meters and m³/s, including the linearized chain.
Normalized primal/dual residuals mix consensus coordinates; their metric unit is
`mixed_norm`, not a claim of a homogeneous physical unit. MPC has no ADMM residual
metrics; a requirement asking for one is an error, not an invented zero.

## Analysis extension

Built-ins use the same `evaluate(ExperimentResult)` interface as external analyzers.
An installed Python package can register a no-argument factory:

```toml
[project.entry-points."controlbench.analyzers"]
terminal_error = "my_analysis:TerminalError"
```

```python
class TerminalError:
    name = "terminal_error"
    version = "1"
    metric_units = {"terminal_max_error": "m"}

    def evaluate(self, result):
        error = abs(result.states[-1] - result.references[-1]).max()
        return {"terminal_max_error": {"value": float(error), "unit": "m"}}
```

Installed plugins are trusted executable software, discovered through package
metadata; configurations cannot load arbitrary files. Declared metric names become
valid requirement keys. Duplicate analyzers/metrics, incorrect units, nonfinite
metrics, and mismatches between declared and returned metrics fail closed. The
`agent_` prefix belongs to built-in per-agent metrics.

## Reproducibility and benchmarking

Runs never overwrite each other. Reports and analysis outputs are derived artifacts;
reanalysis makes a new copy with a new ID. Existing artifact hashes are checked first.
`reproduce` links a new simulation to the original and checks provenance, but does
not automatically check out old code or install dependencies. Keep original artifacts.

Sweeps enumerate finite Cartesian products in declared order. Resume skips only
complete, checksum-valid runs with matching resolved configurations. Execution is
sequential, with a finite timeout and process-tree termination on timeout/interruption.

Numerical regression rules explicitly declare degradation direction and both absolute
and relative allowances. Single-run wall-clock regression checks are rejected.
The benchmark command discards one full warm-up run, measures at least five fresh
process runs, and aggregates their mean solve times with median and interquartile
spread. Timings include only controller work, not process startup. This policy
does not eliminate all machine-load variability.
