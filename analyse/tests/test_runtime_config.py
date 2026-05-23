from __future__ import annotations

import importlib
from pathlib import Path

from lpmo_pipeline.config import clear_runtime_paths_cache, load_runtime_paths_config


def _write_runtime_paths_config(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    project_root = tmp_path / "analyse"
    configs_dir = project_root / "configs"
    schemas_dir = project_root / "schemas"
    configs_dir.mkdir(parents=True, exist_ok=True)
    schemas_dir.mkdir(parents=True, exist_ok=True)

    for file_name in [
        "thresholds.yaml",
        "prolif_features.yaml",
        "defaults.yaml",
        "geometry_rules.yaml",
        "residue_rules.yaml",
        "cv_hierarchy.yaml",
    ]:
        (configs_dir / file_name).write_text("{}\n")
    (schemas_dir / "qc_report_schema.json").write_text("{}\n")

    privateer_sif = tmp_path / "tools" / "privateer.sif"
    posebusters_sif = tmp_path / "tools" / "posebusters.sif"
    privateer_sif.parent.mkdir(parents=True, exist_ok=True)
    privateer_sif.write_text("fake\n")
    posebusters_sif.write_text("fake\n")

    runtime_config = configs_dir / "runtime_paths.yaml"
    runtime_config.write_text(
        """
runtime_paths:
  project_root: ".."
  external_tools:
    apptainer_executable: "/usr/bin/apptainer"
    privateer_sif_candidates:
      - "../tools/privateer.sif"
    posebusters_sif_candidates:
      - "../tools/posebusters.sif"
  runtime_settings:
    privateer_default_mode: "custom-mode"
  pipeline_assets:
    thresholds_config: "configs/thresholds.yaml"
    prolif_features_config: "configs/prolif_features.yaml"
    defaults_config: "configs/defaults.yaml"
    geometry_rules_config: "configs/geometry_rules.yaml"
    residue_rules_config: "configs/residue_rules.yaml"
    cv_hierarchy_config: "configs/cv_hierarchy.yaml"
    qc_report_schema: "schemas/qc_report_schema.json"
""".lstrip()
    )
    return runtime_config, project_root, privateer_sif, posebusters_sif, schemas_dir


def test_load_runtime_paths_config_resolves_project_relative_assets(tmp_path: Path) -> None:
    runtime_config, project_root, privateer_sif, posebusters_sif, schemas_dir = _write_runtime_paths_config(tmp_path)

    clear_runtime_paths_cache()
    config = load_runtime_paths_config(runtime_config)

    assert config.project_root == project_root.resolve()
    assert config.external_tools.apptainer_executable == "/usr/bin/apptainer"
    assert config.external_tools.privateer_sif_candidates == (privateer_sif.resolve(),)
    assert config.external_tools.posebusters_sif_candidates == (posebusters_sif.resolve(),)
    assert config.runtime_settings.privateer_default_mode == "custom-mode"
    assert config.pipeline_assets.thresholds_config == (project_root / "configs" / "thresholds.yaml").resolve()
    assert config.pipeline_assets.qc_report_schema == (schemas_dir / "qc_report_schema.json").resolve()


def test_runtime_paths_config_controls_tool_wrappers(monkeypatch, tmp_path: Path) -> None:
    runtime_config, _, privateer_sif, posebusters_sif, _ = _write_runtime_paths_config(tmp_path)

    import lpmo_pipeline.config as config_mod
    import lpmo_pipeline.qc.posebusters_runner as posebusters_mod
    import lpmo_pipeline.qc.privateer_runner as privateer_mod

    monkeypatch.setenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", str(runtime_config))

    try:
        config_mod.clear_runtime_paths_cache()
        config_mod = importlib.reload(config_mod)
        privateer_mod = importlib.reload(privateer_mod)
        posebusters_mod = importlib.reload(posebusters_mod)

        assert privateer_mod.APPTAINER_EXECUTABLE == "/usr/bin/apptainer"
        assert privateer_mod.PRIVATEER_SIF_CANDIDATES == (privateer_sif.resolve(),)
        assert privateer_mod.PRIVATEER_DEFAULT_MODE == "custom-mode"
        assert posebusters_mod.APPTAINER_EXECUTABLE == "/usr/bin/apptainer"
        assert posebusters_mod.POSEBUSTERS_SIF_CANDIDATES == (posebusters_sif.resolve(),)
    finally:
        monkeypatch.delenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", raising=False)
        config_mod.clear_runtime_paths_cache()
        importlib.reload(config_mod)
        importlib.reload(privateer_mod)
        importlib.reload(posebusters_mod)


def test_runtime_paths_config_controls_default_asset_paths(monkeypatch, tmp_path: Path) -> None:
    runtime_config, project_root, _, _, _ = _write_runtime_paths_config(tmp_path)

    (project_root / "configs" / "thresholds.yaml").write_text(
        """
hdbscan:
  min_cluster_size: 7
  metric: dice
  cluster_selection_method: leaf
  locked: true
""".lstrip()
    )
    (project_root / "configs" / "prolif_features.yaml").write_text(
        """
active_interaction_types:
  - ImplicitHBAcceptor
  - ImplicitHBDonor
  - VdWContact
""".lstrip()
    )

    import lpmo_pipeline.analysis.clustering_hdbscan as clustering_mod
    import lpmo_pipeline.analysis.prolif_ifp as prolif_mod
    import lpmo_pipeline.config as config_mod

    monkeypatch.setenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", str(runtime_config))

    try:
        config_mod.clear_runtime_paths_cache()
        config_mod = importlib.reload(config_mod)
        clustering_mod = importlib.reload(clustering_mod)
        prolif_mod = importlib.reload(prolif_mod)

        hdbscan_config = clustering_mod.load_hdbscan_config()
        prolif_config = prolif_mod.load_prolif_features_config()

        assert hdbscan_config.min_cluster_size == 7
        assert hdbscan_config.metric == "dice"
        assert hdbscan_config.cluster_selection_method == "leaf"
        assert hdbscan_config.locked is True
        assert prolif_config["active_interaction_types"] == [
          "ImplicitHBAcceptor",
          "ImplicitHBDonor",
          "VdWContact",
        ]
    finally:
        monkeypatch.delenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", raising=False)
        config_mod.clear_runtime_paths_cache()
        importlib.reload(config_mod)
        importlib.reload(clustering_mod)
        importlib.reload(prolif_mod)


def test_defaults_config_controls_chain_schema_and_target_mapping(monkeypatch, tmp_path: Path) -> None:
    runtime_config, project_root, _, _, _ = _write_runtime_paths_config(tmp_path)

    (project_root / "configs" / "defaults.yaml").write_text(
        """
chain_schema:
  protein: "P"
  glycans: ["L", "M"]
  metal: "Z"
target_metadata:
  substrate_by_prefix:
    XYZ: xylan
""".lstrip()
    )

    import lpmo_pipeline.analysis.analysis_orchestrator as orchestrator_mod
    import lpmo_pipeline.analysis.mdanalysis_metrics as geometry_mod
    import lpmo_pipeline.config as config_mod
    import lpmo_pipeline.qc.active_site_proximity as proximity_mod
    import lpmo_pipeline.qc.custom_geometry_checks as geometry_checks_mod
    import lpmo_pipeline.qc.hard_qc_orchestrator as hard_qc_mod
    import lpmo_pipeline.qc.posebusters_runner as posebusters_mod

    monkeypatch.setenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", str(runtime_config))

    try:
        config_mod.clear_runtime_paths_cache()
        config_mod = importlib.reload(config_mod)
        orchestrator_mod = importlib.reload(orchestrator_mod)
        geometry_mod = importlib.reload(geometry_mod)
        proximity_mod = importlib.reload(proximity_mod)
        geometry_checks_mod = importlib.reload(geometry_checks_mod)
        hard_qc_mod = importlib.reload(hard_qc_mod)
        posebusters_mod = importlib.reload(posebusters_mod)

        assert posebusters_mod._PROTEIN_CHAIN_ID == "P"
        assert hard_qc_mod.DEFAULT_CU_CHAIN == "Z"
        assert hard_qc_mod.DEFAULT_GLYCAN_CHAINS == ("L", "M")
        assert proximity_mod.DEFAULT_CU_CHAIN == "Z"
        assert proximity_mod.DEFAULT_GLYCAN_CHAINS == ("L", "M")
        assert geometry_checks_mod.DEFAULT_CU_CHAIN == "Z"
        assert geometry_checks_mod.DEFAULT_GLYCAN_CHAINS == ("L", "M")
        assert geometry_mod.DEFAULT_PROTEIN_CHAIN == "P"
        assert geometry_mod.DEFAULT_GLYCAN_CHAINS == ("L", "M")
        assert geometry_mod.DEFAULT_CU_CHAIN == "Z"
        assert orchestrator_mod._parse_target_metadata("XYZ7") == ("xylan", 7)
    finally:
        monkeypatch.delenv("LPMO_PIPELINE_RUNTIME_PATHS_CONFIG", raising=False)
        config_mod.clear_runtime_paths_cache()
        importlib.reload(config_mod)
        importlib.reload(orchestrator_mod)
        importlib.reload(geometry_mod)
        importlib.reload(proximity_mod)
        importlib.reload(geometry_checks_mod)
        importlib.reload(hard_qc_mod)
        importlib.reload(posebusters_mod)