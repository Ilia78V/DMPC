# Verification record — 2026-09-29

This distinguishes completed measurements from release checks still awaiting a
working MATLAB license connection. Generated artifacts stay in the ignored `results/`
directory; they are real local files, not example benchmark numbers.

## Completed MATLAB checks

MATLAB R2026a Update 4 successfully solved the dependency smoke problem through
YALMIP/IPOPT. The first full test suite completed with **10 passed, 0 failed,
0 incomplete**. It covered exact double-integrator propagation, equilibrium,
invalid intervals, conservative coupling, agent-order independence, no state
projection, sample counts, infeasibility, iteration exhaustion, and convex
centralized/DMPC control agreement.

Extending the same convex comparison to `approximation: all` exposed an existing
external-influence defect: its summation began with a free optimization variable
instead of zero. Before correction it exhausted 400 iterations and differed from
the centralized first input by about 0.04965 (approximately 5.3%). After correction,
the unchanged comparison passed for both modes: 118 iterations without approximation
and 38 with approximation, with absolute control tolerance 0.005. These are
specific test results, not general convergence guarantees.

## Completed end-to-end experiment

Run `c0be877f-b963-4ccd-b59c-3e52dccdaee3` ran the 4-tank, 20-interval horizon,
30-control-step nominal scenario and produced 31 states and 30 inputs per agent.
Its four configured requirements passed: no state violations, no input violations,
no failed local solves, and convergence at all control steps.

| Metric | Observed value |
|---|---:|
| Pooled tracking RMSE, including initial sample | 0.690498 m |
| State violation count | 0 |
| Input violation count | 0 |
| Converged-step fraction | 1.0 |
| Mean controller solve time | 2848.207 ms |
| Maximum controller solve time | 23660.203 ms |
| Independent physical objective | 5.404527 |

These timings are a single development-host run, with other verification work
occurring during development. They are not a controlled performance baseline or
a real-time claim. This saved source snapshot predates the later boundary-tolerance
and solver-bound-relaxation changes; rerun before attributing these measurements
to the final working-tree version. The original artifacts and their hashes remain
unchanged.

## Failures retained as evidence

The full centralized tank case initially stopped with a local IPOPT infeasibility
result. Partial export exposed a fixed initial level of `3.0000000292857` against
an upper bound of 3.0. The implemented correction keeps that actual measured state,
admits the configured `1e-6` tolerance only at the fixed initial-state constraint,
and disables IPOPT bound relaxation for future decision variables. Both controller
adapters use this convention. A dedicated boundary-roundoff test was added.

Subsequent MATLAB launches failed before executing project code with MathWorks
license-server error **-15.2**. Therefore the final boundary correction, source-snapshot
execution path, solver-binary provenance, full centralized baseline, and final
MATLAB suite need a rerun once the university license connection is available.
The MATLAB connector also could not attach; earlier successful tests used the
installed MATLAB batch interface.

## Python and software checks

The automated Python suite covers independent metric values, corrupt/missing data,
configuration validation, installed analysis plugins, requirement and regression
failures, process errors/timeouts, campaign expansion/resume, provenance checks,
report escaping, reanalysis, reproduction, and repeated timing aggregation. Tests
that replace MATLAB with a synthetic exporter are explicitly pipeline tests.

The latest completed run recorded **59 passing Python tests and 89.83% coverage**.
An installable Python wheel also built successfully. Run
`pytest --cov=controlbench --cov-fail-under=85` to refresh the count/coverage.
Ruff lint/format and mypy must also pass. Type checking is exercised with Windows
and Linux platform settings. CI configuration is provided but has not been run on
GitHub, and no licensed self-hosted runner has been provisioned here.

## Finish release verification

After restoring MATLAB's license-server connection:

1. Run `controlbench doctor`.
2. Run `matlab -singleCompThread -batch "r=runtests('matlab/tests'); assertSuccess(r);"`.
3. Run the nominal, centralized and approximation YAML experiments.
4. Run the linear comparison campaign and horizon sweep; retain their reports.
5. Reproduce one final run and compare its numerical metrics.
6. Run timing benchmarks only on an otherwise idle controlled host if runtime
   regression evidence is needed.

Until these MATLAB checks finish, this is a working implementation with substantial
local verification, not a fully verified v1.0 release. No GitHub PR, historical
release, review or cloud deployment has been fabricated.
