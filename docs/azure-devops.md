# Azure DevOps CI and Blob publication

`azure-pipelines.yml` implements three gates:

1. Portable Python quality, contract, and artifact-integrity tests.
2. Solver-backed MATLAB tests plus the real `experiments/tanks/smoke.yaml` DMPC run.
3. Publication of the verified bundle to private Azure Blob Storage.

The MATLAB job intentionally targets the self-hosted `ControlBench` pool. Its agent
must advertise `controlbench.dependencies=true` and provide MATLAB, YALMIP, IPOPT,
and Python 3.12. This avoids pretending that a hosted image contains the project's
solver stack or university MATLAB license. Pull-request code reaches that trusted
machine only after the portable stage passes.

The publication stage uses the Azure Resource Manager service connection
`controlbench-blob-wif`. It should use workload identity federation and receive only
`Storage Blob Data Contributor` on the result storage account. Do not grant it
subscription-wide contributor access. The pipeline variable `AZURE_STORAGE_ACCOUNT`
contains the account name, not a key or connection string.

Results are uploaded beneath:

```text
simulation-results/<git-commit>/<azure-build-id>/
```

The container remains private. Every upload contains `ci-metadata.json`, the HTML
report, numerical CSV files, diagnostics, environment metadata, source snapshot,
and the ControlBench manifest/checksums. `stage_simulation_results.py` refuses to
stage incomplete runs or runs whose engineering requirements do not pass.

`infra/main.bicep` defines the storage account and private container. Deploy it to
a dedicated resource group, then grant the service connection identity the
`Storage Blob Data Contributor` role at the storage-account scope.

For local parity:

```powershell
python -m pip install -e ".[dev]"
ruff check src tests scripts
ruff format --check src tests scripts
mypy src/controlbench scripts/ci
pytest --cov=controlbench --cov-fail-under=85
matlab -singleCompThread -batch "buildtool ci"
controlbench run experiments/tanks/smoke.yaml --output results/ci
python scripts/ci/stage_simulation_results.py --results results/ci --staging artifacts/simulation-results
```
