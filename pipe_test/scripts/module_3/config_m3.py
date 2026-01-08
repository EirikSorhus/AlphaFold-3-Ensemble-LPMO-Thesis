# scripts/module_3/config_m3.py
"""
Configuration and defaults for Module 3: Domain annotation and sequence variant generation.

The configuration is intentionally lightweight and only relies on the standard
library. If a YAML configuration file is present it will be used to override
the defaults, otherwise the built-in defaults are applied.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    yaml = None


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DATA_SEQUENCES_DIR = DATA_DIR / "sequences"
DATA_METADATA_DIR = DATA_DIR / "metadata"
DATA_DOMAINS_DIR = DATA_DIR / "domains"
RAW_DOMTBL_DIR = DATA_DOMAINS_DIR / "dbcan"  

# Module 2 outputs used as inputs here
# This config tries to find the most recent ones, or you can override via --fasta/--metadata
def _find_latest_m2_fasta() -> Path:
    """Find the most recent module 2 output FASTA file."""
    pattern = DATA_SEQUENCES_DIR / "lpmo_all_*_raw.fasta"
    files = list(DATA_SEQUENCES_DIR.glob("lpmo_all_*_raw.fasta"))
    if files:
        # Return the most recently modified
        return sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)[0]
    # Fallback to old filename
    return DATA_SEQUENCES_DIR / "lpmo_all_raw.fasta"

def _find_latest_m2_metadata() -> Path:
    """Find the most recent module 2 metadata CSV file."""
    pattern = DATA_METADATA_DIR / "m2_sequence_3d_metadata_*.csv"
    files = list(DATA_METADATA_DIR.glob("m2_sequence_3d_metadata_*.csv"))
    if files:
        # Return the most recently modified
        return sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)[0]
    # Fallback to old filename
    return DATA_METADATA_DIR / "m2_sequence_3d_metadata.csv"

M2_FASTA = _find_latest_m2_fasta()
M2_METADATA = _find_latest_m2_metadata()

# Module 3 outputs
DOMAINS_PARSED = DATA_DOMAINS_DIR / "m3_domains_parsed.csv"
DOMAIN_METADATA = DATA_METADATA_DIR / "m3_domain_metadata.csv"
RUN_METADATA = DATA_METADATA_DIR / "m3_run_metadata.json"
FULL_LENGTH_FASTA = DATA_SEQUENCES_DIR / "lpmo_full_length.fasta"
CATALYTIC_FASTA = DATA_SEQUENCES_DIR / "lpmo_catalytic_domain.fasta"

# Default location for a user-provided configuration file
DEFAULT_CONFIG_FILE = PROJECT_ROOT / "config_m3.yaml"


# -----------------------------------------------------------------------------
# Dataclass with defaults
# -----------------------------------------------------------------------------

@dataclass
class Module3Config:
    # HMM filtering
    evalue_cutoff: float = 1e-15
    coverage_cutoff: float = 0.35
    min_domain_length: int = 60
    max_domain_overlap_fraction: float = 0.25
    prefer_family_from_metadata: bool = True

    # HMM sources (new structure: data/domains/dbcan/)
    dbcan_hmm: Path = PROJECT_ROOT / "data" / "domains" / "dbcan" / "dbCAN-HMMdb-V14.hmm"
    dbcan_sub_hmm: Path = PROJECT_ROOT / "data" / "domains" / "dbcan" / "dbCAN_sub.hmm"
    cbm_hmm: Path = None  # Not used anymore - CBM domains are in dbCAN
    use_subfamily_hmms: bool = True
    hmmscan_binary: str = "hmmscan"
    hmmscan_cpu: int = 2
    reuse_existing_domtbl: bool = True

    # Policy for sequence variants
    generate_catalytic_for_cbm_only: bool = True
    generate_catalytic_if_no_cbm: bool = False
    fallback_to_full_length_if_no_domain: bool = True
    min_catalytic_domain_length: int = 150

    # Known LPMO families (used when prefer_family_from_metadata is False or metadata missing)
    lpmo_families: tuple = ("AA9", "AA10", "AA11", "AA13")

    # Optional notes recorded in run metadata
    notes: str = ""


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def ensure_output_directories() -> None:
    """Create output directories used by Module 3."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    DATA_SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)
    DATA_METADATA_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DOMAINS_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DOMTBL_DIR.mkdir(parents=True, exist_ok=True)


def _coerce_value(field_name: str, value: Any) -> Any:
    """Coerce loaded config values to correct type where possible."""
    path_fields = {"dbcan_hmm", "dbcan_sub_hmm", "cbm_hmm"}
    if field_name in path_fields and isinstance(value, str):
        return Path(value)
    return value


def _load_from_mapping(config: Module3Config, mapping: Dict[str, Any]) -> Module3Config:
    """Update a config dataclass from a mapping, ignoring unknown keys."""
    for f in fields(config):
        if f.name in mapping:
            setattr(config, f.name, _coerce_value(f.name, mapping[f.name]))
    return config


def load_config(config_file: Optional[Path] = None) -> Module3Config:
    """
    Load Module 3 configuration from YAML/JSON (if available) with sensible defaults.

    Args:
        config_file: Optional path to config file. If None, DEFAULT_CONFIG_FILE is used.
    """
    config = Module3Config()
    path = config_file or DEFAULT_CONFIG_FILE

    if not path.exists():
        return config

    try:
        if yaml:
            with open(path, "r", encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle)
        else:
            with open(path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[CONFIG] WARNING: Failed to read {path}: {exc}")
        return config

    if isinstance(loaded, dict):
        _load_from_mapping(config, loaded)
    else:
        print(f"[CONFIG] WARNING: Unsupported config format in {path}, using defaults")

    return config


def config_to_dict(config: Module3Config) -> Dict[str, Any]:
    """Convert the config dataclass to a dict for serialization."""
    serializable = asdict(config)
    for key in ("dbcan_hmm", "dbcan_sub_hmm", "cbm_hmm"):
        value = serializable.get(key)
        if isinstance(value, Path):
            serializable[key] = str(value)
    return serializable


__all__ = [
    "Module3Config",
    "PROJECT_ROOT",
    "DATA_DIR",
    "DATA_SEQUENCES_DIR",
    "DATA_METADATA_DIR",
    "DATA_DOMAINS_DIR",
    "RAW_DOMTBL_DIR",
    "M2_FASTA",
    "M2_METADATA",
    "DOMAINS_PARSED",
    "DOMAIN_METADATA",
    "RUN_METADATA",
    "FULL_LENGTH_FASTA",
    "CATALYTIC_FASTA",
    "DEFAULT_CONFIG_FILE",
    "load_config",
    "ensure_output_directories",
    "config_to_dict",
]
