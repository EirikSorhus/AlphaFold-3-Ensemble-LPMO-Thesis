# scripts/module_3/run_module3.py
"""
Module 3 entrypoint: Domain annotation and sequence variant generation.

Example:
    python -m scripts.module_3.run_module3 --config config_m3.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from scripts.module_3.config_m3 import (
    Module3Config,
    M2_FASTA,
    M2_METADATA,
    load_config,
)
from scripts.module_3.domain_annotation import run_domain_pipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Module 3: Domain annotation and sequence variant preparation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=Path, default=None, help="Path to Module 3 YAML/JSON config file")
    parser.add_argument("--fasta", type=Path, default=M2_FASTA, help="Input FASTA (from Module 2)")
    parser.add_argument("--metadata", type=Path, default=M2_METADATA, help="Module 2 metadata CSV")
    parser.add_argument("--hmmscan-binary", type=str, default=None, help="Override hmmscan binary path")
    parser.add_argument("--no-subfamily", action="store_true", help="Disable dbCAN subfamily HMMs")
    parser.add_argument("--force-rerun", action="store_true", help="Ignore existing domtblout files")
    parser.add_argument("--cpu", type=int, default=None, help="Number of CPUs for hmmscan")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    config: Module3Config = load_config(args.config)
    if args.hmmscan_binary:
        config.hmmscan_binary = args.hmmscan_binary
    if args.no_subfamily:
        config.use_subfamily_hmms = False
    if args.force_rerun:
        config.reuse_existing_domtbl = False
    if args.cpu:
        config.hmmscan_cpu = args.cpu

    try:
        run_domain_pipeline(config, fasta_path=args.fasta, metadata_path=args.metadata)
    except FileNotFoundError as exc:
        print(f"[MODULE3] ERROR: {exc}")
        return 1
    except Exception as exc:  # pragma: no cover - defensive
        print(f"[MODULE3] ERROR: Unexpected failure: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
