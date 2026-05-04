from __future__ import annotations

import subprocess

from lpmo_pipeline.utils.manifest import ToolVersionFetcher


def test_tool_version_fetcher_uses_privateer_sif_helper(monkeypatch) -> None:
    from lpmo_pipeline.utils import manifest as mod

    seen_commands: list[list[str]] = []
    real_run = subprocess.run

    class Proc:
        def __init__(self, stdout: str = "") -> None:
            self.stdout = stdout

    def _fake_run(command, *args, **kwargs):
        seen_commands.append(command)
        if command[:2] == ["placer", "--version"]:
            return Proc("placer 1.2.3")
        return real_run(command, *args, **kwargs)

    monkeypatch.setattr(mod, "get_privateer_version", lambda: "MKV")
    monkeypatch.setattr(mod.subprocess, "run", _fake_run)

    versions = ToolVersionFetcher.get_versions()

    assert versions["privateer"] == "MKV"
    assert ["privateer", "-V"] not in seen_commands


def test_tool_probe_fetcher_runs_all_analysis_probes(monkeypatch) -> None:
    from lpmo_pipeline.utils import manifest as mod

    monkeypatch.setattr(mod, "_probe_gemmi", lambda: "ok")
    monkeypatch.setattr(mod, "_probe_mdanalysis", lambda: "ok")
    monkeypatch.setattr(mod, "_probe_prolif", lambda: "ok")
    monkeypatch.setattr(mod, "_probe_hdbscan", lambda: "ok")

    statuses = ToolVersionFetcher.get_analysis_dependency_statuses()

    assert statuses == {
        "gemmi": "ok",
        "mdanalysis": "ok",
        "prolif": "ok",
        "hdbscan": "ok",
    }


def test_tool_probe_fetcher_keeps_probe_failure_details(monkeypatch) -> None:
    from lpmo_pipeline.utils import manifest as mod

    monkeypatch.setattr(mod, "_probe_gemmi", lambda: "missing_required_api")
    monkeypatch.setattr(mod, "_probe_mdanalysis", lambda: "ValueError: broken wheel")
    monkeypatch.setattr(mod, "_probe_prolif", lambda: "ok")
    monkeypatch.setattr(mod, "_probe_hdbscan", lambda: "RuntimeError: probe failed")

    statuses = ToolVersionFetcher.get_analysis_dependency_statuses()

    assert statuses["gemmi"] == "missing_required_api"
    assert statuses["mdanalysis"] == "ValueError: broken wheel"
    assert statuses["prolif"] == "ok"
    assert statuses["hdbscan"] == "RuntimeError: probe failed"
