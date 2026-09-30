import pytest
import yaml

from controlbench import analysis
from controlbench.config import document, resolve
from controlbench.data import ExperimentResult


def test_azure_pipeline_keeps_validation_before_blob_publication():
    pipeline = yaml.safe_load(open("azure-pipelines.yml", encoding="utf-8"))
    stages = [stage["stage"] for stage in pipeline["stages"]]
    assert stages == ["PortableValidation", "DMPCSimulation", "PublishBlob"]
    text = open("azure-pipelines.yml", encoding="utf-8").read()
    assert "stage_simulation_results.py" in text
    assert "--auth-mode login" in text
    assert "connection-string" not in text.lower()
    assert "account-key" not in text.lower()


@pytest.fixture
def registry(monkeypatch):
    monkeypatch.setattr(analysis, "REGISTRY", dict(analysis.REGISTRY))
    monkeypatch.setattr(analysis, "EXTRA_METRICS", {})


class Terminal:
    name, version = "terminal", "1"
    metric_units = {"terminal_max": "m"}

    def evaluate(self, result):
        return {
            "terminal_max": {
                "value": float(abs(result.states[-1] - result.references[-1]).max()),
                "unit": "m",
            }
        }


def test_plugin_registration_and_requirements(registry, exported):
    analysis.register(Terminal())
    config = resolve(
        {
            "requirements": [
                {"id": "terminal", "metric": "terminal_max", "operator": "le", "limit": 1}
            ]
        }
    )
    assert config["requirements"][0]["metric"] == "terminal_max"
    assert analysis.analyze(ExperimentResult.load(exported))["terminal_max"]["value"] == 0.5
    with pytest.raises(ValueError, match="Duplicate analyzer"):
        analysis.register(Terminal())


@pytest.mark.parametrize("units", [{}, {"tracking_rmse": "m"}, {"agent_1_custom": "m"}, {1: "m"}])
def test_invalid_plugin_names(registry, units):
    plugin = Terminal()
    plugin.metric_units = units
    with pytest.raises(ValueError):
        analysis.register(plugin)


def test_plugin_invalid_units(registry, exported):
    plugin = Terminal()
    plugin.metric_units = {"terminal_max": "cm"}
    analysis.register(plugin)
    with pytest.raises(ValueError, match="unit mismatch"):
        analysis.analyze(ExperimentResult.load(exported))


def test_plugin_discovery(registry, monkeypatch):
    class Entry:
        def load(self):
            return Terminal

    monkeypatch.setattr(analysis, "_discovered", False)
    monkeypatch.setattr(analysis, "entry_points", lambda **kwargs: [Entry()])
    analysis.discover_plugins()
    assert "terminal_max" in analysis.metric_names(2)
    analysis.discover_plugins()  # Discovery happens once.


def test_failed_discovery_can_never_skip_broken_plugin(registry, monkeypatch):
    class Entry:
        def load(self):
            raise RuntimeError("broken installation")

    monkeypatch.setattr(analysis, "_discovered", False)
    monkeypatch.setattr(analysis, "entry_points", lambda **kwargs: [Entry()])
    with pytest.raises(RuntimeError, match="broken installation"):
        analysis.discover_plugins()
    assert not analysis._discovered
    with pytest.raises(RuntimeError, match="broken installation"):
        analysis.discover_plugins()


def test_scientific_notation_and_duplicate_keys(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("constraints:\n  tolerance: 1e-6\n")
    assert resolve(document(path))["constraints"]["tolerance"] == 1e-6
    path.write_text("plant:\n  agents: 2\n  agents: 4\n")
    with pytest.raises(ValueError, match="unique"):
        document(path)


@pytest.mark.parametrize(
    "file,old,new",
    [
        ("states.csv", "0.1,1", "0.1,2"),
        ("states.csv", "1.0,2,0,3", "1.0,2,0,4"),
        ("states.csv", "0.5,2,0,3", "0.6,2,0,3"),
        ("steps.csv", "normalized_primal_v1", "dual_only"),
        ("inputs.csv", ",2,2,0,1000", ",2,3,0,1000"),
        ("residuals.csv", ",0,0,0", ",0,0,1"),
    ],
)
def test_data_invariants(exported, file, old, new):
    path = exported / file
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new))
    with pytest.raises(ValueError):
        ExperimentResult.load(exported)


def test_missing_diagnostics(exported):
    path = exported / "local_solves.csv"
    path.write_text("\n".join(path.read_text().splitlines()[:-1]))
    with pytest.raises(ValueError, match="diagnostics"):
        ExperimentResult.load(exported)
