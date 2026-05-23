"""Shared configuration loader for runtime paths and package assets.

This module is the documented mechanism for resolving changeable runtime
paths in the analyse pipeline. External tool paths and package asset
locations live in YAML config so code does not need to hardcode them.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml


_RUNTIME_PATHS_ENV_VAR = "LPMO_PIPELINE_RUNTIME_PATHS_CONFIG"
_PIPELINE_RUN_CONFIG_ENV_VAR = "LPMO_PIPELINE_RUN_CONFIG"
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_RUNTIME_PATHS_CONFIG = _PROJECT_ROOT / "configs" / "runtime_paths.yaml"


@dataclass(frozen=True)
class ExternalToolConfig:
    """Configured external tool executables and candidate container paths."""

    apptainer_executable: str
    privateer_sif_candidates: tuple[Path, ...]
    posebusters_sif_candidates: tuple[Path, ...]


@dataclass(frozen=True)
class RuntimeSettingsConfig:
    """Configured runtime defaults for external tool wrappers."""

    privateer_default_mode: str
    python_executable: str


@dataclass(frozen=True)
class PipelineAssetConfig:
    """Configured paths to versioned in-repo config and schema assets."""

    thresholds_config: Path
    prolif_features_config: Path
    defaults_config: Path
    geometry_rules_config: Path
    residue_rules_config: Path
    cv_hierarchy_config: Path
    qc_report_schema: Path


@dataclass(frozen=True)
class RuntimePathsConfig:
    """Resolved runtime path configuration for the analyse package."""

    config_path: Path
    project_root: Path
    external_tools: ExternalToolConfig
    runtime_settings: RuntimeSettingsConfig
    pipeline_assets: PipelineAssetConfig


def _expand_string(value: str) -> str:
    return os.path.expanduser(os.path.expandvars(value))


def _resolve_path(value: str | Path, *, base_dir: Path) -> Path:
    path = Path(_expand_string(str(value)))
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def _resolve_path_list(values: Any, *, base_dir: Path) -> tuple[Path, ...]:
    if values is None:
        return ()
    if not isinstance(values, list):
        raise TypeError(f"Expected a list of paths, got {type(values).__name__}")
    return tuple(_resolve_path(value, base_dir=base_dir) for value in values)


def _resolve_command(value: str, *, base_dir: Path) -> str:
    expanded = _expand_string(value)
    if "/" in expanded:
        return str(_resolve_path(expanded, base_dir=base_dir))
    return expanded


def _require_mapping(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key) or {}
    if not isinstance(value, dict):
        raise TypeError(f"Expected mapping for '{key}', got {type(value).__name__}")
    return value


def load_runtime_paths_config(config_path: str | Path | None = None) -> RuntimePathsConfig:
    """Load resolved runtime paths from YAML.

    Resolution order:
    1. explicit ``config_path``
    2. ``LPMO_PIPELINE_RUNTIME_PATHS_CONFIG`` environment variable
    3. bundled ``configs/runtime_paths.yaml``
    """

    cache_key = str(config_path) if config_path is not None else None
    return _load_runtime_paths_config_cached(cache_key)


def load_defaults_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load defaults.yaml through the shared runtime path mechanism."""

    cache_key = str(config_path) if config_path is not None else None
    return _load_defaults_config_cached(cache_key)


def load_pipeline_run_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load hybrid pipeline run config YAML.

    Resolution order:
    1. explicit ``config_path``
    2. ``LPMO_PIPELINE_RUN_CONFIG`` environment variable

    The config must define a top-level mapping. Optional sections are:
    - ``base``: shared defaults across command groups
    - ``commands``: per-command overrides, keyed by command name
      (for example ``run``, ``predictive``, ``cbm_paired``,
      ``family_enrichment``, ``tune``).
    """

    cache_key = str(config_path) if config_path is not None else None
    return _load_pipeline_run_config_cached(cache_key)


def resolve_pipeline_command_config(
    pipeline_run_config: dict[str, Any],
    command_name: str,
) -> dict[str, Any]:
    """Resolve merged command config from ``base`` + ``commands.<name>``.

    ``command_name`` accepts CLI names (for example ``cbm-paired``) and will
    be normalized to underscore keys (for example ``cbm_paired``).
    """

    if not isinstance(pipeline_run_config, dict):
        raise TypeError(
            "Pipeline run config must be a mapping, got "
            f"{type(pipeline_run_config).__name__}"
        )

    base = pipeline_run_config.get("base") or {}
    commands = pipeline_run_config.get("commands") or {}
    if not isinstance(base, dict):
        raise TypeError(
            "Pipeline run config field 'base' must be a mapping"
        )
    if not isinstance(commands, dict):
        raise TypeError(
            "Pipeline run config field 'commands' must be a mapping"
        )

    normalized_name = command_name.replace("-", "_").strip()
    command_section = commands.get(normalized_name) or {}
    if not isinstance(command_section, dict):
        raise TypeError(
            "Pipeline run config field "
            f"'commands.{normalized_name}' must be a mapping"
        )
    return _deep_merge_mappings(base, command_section)


@lru_cache(maxsize=4)
def _load_runtime_paths_config_cached(config_path: str | None) -> RuntimePathsConfig:
    raw_config_path = config_path or os.environ.get(_RUNTIME_PATHS_ENV_VAR) or str(
        _DEFAULT_RUNTIME_PATHS_CONFIG
    )
    resolved_config_path = _resolve_path(raw_config_path, base_dir=_PROJECT_ROOT)
    if not resolved_config_path.exists():
        raise FileNotFoundError(
            "Runtime paths config not found: "
            f"{resolved_config_path}. Set {_RUNTIME_PATHS_ENV_VAR} or create configs/runtime_paths.yaml"
        )

    payload = yaml.safe_load(resolved_config_path.read_text()) or {}
    runtime_paths = _require_mapping(payload, "runtime_paths")
    config_dir = resolved_config_path.parent
    project_root_value = runtime_paths.get("project_root") or ".."
    project_root = _resolve_path(project_root_value, base_dir=config_dir)

    external_tools = _require_mapping(runtime_paths, "external_tools")
    runtime_settings = _require_mapping(runtime_paths, "runtime_settings")
    pipeline_assets = _require_mapping(runtime_paths, "pipeline_assets")

    return RuntimePathsConfig(
        config_path=resolved_config_path,
        project_root=project_root,
        external_tools=ExternalToolConfig(
            apptainer_executable=_resolve_command(
                str(external_tools.get("apptainer_executable") or "apptainer"),
                base_dir=project_root,
            ),
            privateer_sif_candidates=_resolve_path_list(
                external_tools.get("privateer_sif_candidates"),
                base_dir=project_root,
            ),
            posebusters_sif_candidates=_resolve_path_list(
                external_tools.get("posebusters_sif_candidates"),
                base_dir=project_root,
            ),
        ),
        runtime_settings=RuntimeSettingsConfig(
            privateer_default_mode=str(
                runtime_settings.get("privateer_default_mode") or "ccp4i2"
            ),
            python_executable=_resolve_command(
                str(runtime_settings.get("python_executable") or "python"),
                base_dir=project_root,
            ),
        ),
        pipeline_assets=PipelineAssetConfig(
            thresholds_config=_resolve_path(
                pipeline_assets.get("thresholds_config") or "configs/thresholds.yaml",
                base_dir=project_root,
            ),
            prolif_features_config=_resolve_path(
                pipeline_assets.get("prolif_features_config") or "configs/prolif_features.yaml",
                base_dir=project_root,
            ),
            defaults_config=_resolve_path(
                pipeline_assets.get("defaults_config") or "configs/defaults.yaml",
                base_dir=project_root,
            ),
            geometry_rules_config=_resolve_path(
                pipeline_assets.get("geometry_rules_config") or "configs/geometry_rules.yaml",
                base_dir=project_root,
            ),
            residue_rules_config=_resolve_path(
                pipeline_assets.get("residue_rules_config") or "configs/residue_rules.yaml",
                base_dir=project_root,
            ),
            cv_hierarchy_config=_resolve_path(
                pipeline_assets.get("cv_hierarchy_config") or "configs/cv_hierarchy.yaml",
                base_dir=project_root,
            ),
            qc_report_schema=_resolve_path(
                pipeline_assets.get("qc_report_schema") or "schemas/qc_report_schema.json",
                base_dir=project_root,
            ),
        ),
    )


def clear_runtime_paths_cache() -> None:
    """Clear the cached runtime path config.

    Primarily useful for tests that swap the config file via environment
    variables and then reload modules.
    """

    _load_runtime_paths_config_cached.cache_clear()
    _load_defaults_config_cached.cache_clear()
    _load_pipeline_run_config_cached.cache_clear()


@lru_cache(maxsize=4)
def _load_defaults_config_cached(config_path: str | None) -> dict[str, Any]:
    path = (
        Path(config_path).resolve()
        if config_path is not None
        else load_runtime_paths_config().pipeline_assets.defaults_config
    )
    return yaml.safe_load(path.read_text()) or {}


@lru_cache(maxsize=4)
def _load_pipeline_run_config_cached(config_path: str | None) -> dict[str, Any]:
    raw_config_path = config_path or os.environ.get(_PIPELINE_RUN_CONFIG_ENV_VAR)
    if raw_config_path is None:
        raise FileNotFoundError(
            "Pipeline run config path is not set. "
            f"Pass config_path or set {_PIPELINE_RUN_CONFIG_ENV_VAR}."
        )

    resolved_config_path = _resolve_path(raw_config_path, base_dir=_PROJECT_ROOT)
    if not resolved_config_path.exists():
        raise FileNotFoundError(
            "Pipeline run config not found: "
            f"{resolved_config_path}. Pass config_path or set {_PIPELINE_RUN_CONFIG_ENV_VAR}."
        )

    payload = yaml.safe_load(resolved_config_path.read_text()) or {}
    if not isinstance(payload, dict):
        raise TypeError(
            "Pipeline run config root must be a mapping"
        )
    return payload


def _deep_merge_mappings(
    base: dict[str, Any],
    override: dict[str, Any],
) -> dict[str, Any]:
    merged: dict[str, Any] = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge_mappings(
                merged[key],
                value,
            )
        else:
            merged[key] = value
    return merged
