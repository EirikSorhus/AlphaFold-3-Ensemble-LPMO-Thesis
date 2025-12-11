"""
Configuration for SignalPeptid module (SignalP6 + trimming).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict, fields
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    yaml = None

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DATA_SEQUENCES_DIR = DATA_DIR / "sequences"
DATA_METADATA_DIR = DATA_DIR / "metadata"
DATA_SIGNALP_DIR = DATA_DIR / "signalpeptide"
DEFAULT_CONFIG_FILE = PROJECT_ROOT / "config_signalp.yaml"


def _find_latest_m2_fasta() -> Path:
    """Find the latest Module 2 output FASTA (raw)."""
    files = sorted(DATA_SEQUENCES_DIR.glob("lpmo_all_*_raw.fasta"), key=lambda p: p.stat().st_mtime, reverse=True)
    if files:
        return files[0]
    return DATA_SEQUENCES_DIR / "lpmo_all_raw.fasta"


def _find_latest_m2_metadata() -> Path:
    """Find the latest Module 2 metadata CSV."""
    files = sorted(DATA_METADATA_DIR.glob("m2_sequence_3d_metadata_*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    if files:
        return files[0]
    return DATA_METADATA_DIR / "m2_sequence_3d_metadata.csv"


@dataclass
class SignalPConfig:
    signalp6_path: str = "/cluster/home/eisorhus/.local/bin/signalp6"
    signalp6_output_dir: Path = DATA_SIGNALP_DIR
    signalp6_mode: str = "fast"
    signalp6_format: str = "txt"
    signalp6_max_seqs_per_run: int = 5000

    min_signal_peptide_length: int = 10
    max_signal_peptide_length: int = 70
    treat_missing_as_no_sp: bool = True

    input_fasta: Path = _find_latest_m2_fasta()
    input_metadata: Path = _find_latest_m2_metadata()
    output_fasta: Path = DATA_SEQUENCES_DIR / "lpmo_mature.fasta"
    output_metadata: Path = DATA_METADATA_DIR / "m2_sequence_3d_metadata_with_sp.csv"
    parsed_signalp_tsv: Path = DATA_SIGNALP_DIR / "signalp6_parsed.tsv"
    run_metadata_json: Path = DATA_METADATA_DIR / "signalpeptide_run_metadata.json"


_PATH_FIELDS = {"signalp6_output_dir", "input_fasta", "input_metadata", "output_fasta", "output_metadata", "parsed_signalp_tsv", "run_metadata_json"}


def _coerce_value(field_name: str, value: Any) -> Any:
    if field_name in _PATH_FIELDS and isinstance(value, str):
        return Path(value)
    return value


def _load_from_mapping(config: SignalPConfig, mapping: Dict[str, Any]) -> SignalPConfig:
    for f in fields(config):
        if f.name in mapping:
            setattr(config, f.name, _coerce_value(f.name, mapping[f.name]))
    return config


def load_config(config_file: Optional[Path] = None) -> SignalPConfig:
    path = config_file or DEFAULT_CONFIG_FILE
    config = SignalPConfig()

    if not path.exists():
        return config

    try:
        if yaml:
            with open(path, "r", encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle)
        else:
            with open(path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
    except Exception:
        return config

    if isinstance(loaded, dict):
        _load_from_mapping(config, loaded)
    return config


def config_to_dict(config: SignalPConfig) -> Dict[str, Any]:
    serializable = asdict(config)
    for key in _PATH_FIELDS:
        value = serializable.get(key)
        if isinstance(value, Path):
            serializable[key] = str(value)
    return serializable


def ensure_output_dirs(config: SignalPConfig) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    DATA_SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)
    DATA_METADATA_DIR.mkdir(parents=True, exist_ok=True)
    DATA_SIGNALP_DIR.mkdir(parents=True, exist_ok=True)
    config.signalp6_output_dir.mkdir(parents=True, exist_ok=True)


__all__ = [
    "SignalPConfig",
    "load_config",
    "config_to_dict",
    "ensure_output_dirs",
    "PROJECT_ROOT",
    "DATA_SEQUENCES_DIR",
    "DATA_METADATA_DIR",
    "DATA_SIGNALP_DIR",
]
