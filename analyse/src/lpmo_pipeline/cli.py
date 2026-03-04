"""
LPMO Pipeline: CLI Entry Point
Responsibility: Main command-line interface with two modes: tune, run (production)

Usage:
  python -m lpmo_pipeline tune --config configs/tuning_af3.yaml --output results/tuning_af3/
  python -m lpmo_pipeline run --mode production --config configs/production.yaml --output results/del_a/
"""

import argparse
import sys
from pathlib import Path
from typing import Optional
import json
from datetime import datetime

from lpmo_pipeline.tuning.tune_orchestrator import run_tuning
from lpmo_pipeline.utils.manifest import ManifestBuilder, ToolVersionFetcher


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="LPMO Structure Prediction → Analysis Pipeline (v2.1)"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Pipeline mode")
    
    # TUNE subcommand
    tune_parser = subparsers.add_parser("tune", help="Run parameter tuning")
    tune_parser.add_argument("--model", required=True, choices=["AF3", "RF3", "Boltz2"],
                             help="Prediction model to tune")
    tune_parser.add_argument("--config", type=Path, required=True,
                             help="Tuning config YAML (parameter grid + dataset)")
    tune_parser.add_argument("--output", type=Path, required=True,
                             help="Output directory for tuning results")
    tune_parser.add_argument("--n-jobs", type=int, default=1,
                             help="Number of parallel jobs")
    
    # RUN subcommand
    run_parser = subparsers.add_parser("run", help="Run production analysis")
    run_parser.add_argument("--mode", required=True, choices=["production"],
                            help="Execution mode (fixed to production)")
    run_parser.add_argument("--config", type=Path, required=True,
                            help="Main config YAML (locked from tuning)")
    run_parser.add_argument("--output", type=Path, required=True,
                            help="Output directory for analysis results")
    run_parser.add_argument("--del", required=True, choices=["del_a", "del_b"],
                            help="Analysis branch (all proteins or CBM subset)")
    run_parser.add_argument("--n-jobs", type=int, default=1,
                            help="Number of parallel jobs")
    
    args = parser.parse_args()
    
    if args.command == "tune":
        return cmd_tune(args)
    elif args.command == "run":
        return cmd_run(args)
    else:
        parser.print_help()
        return 1


def cmd_tune(args):
    """Execute tuning subcommand."""
    print(f"[TUNE] Starting {args.model} tuning...")
    print(f"  Config: {args.config}")
    print(f"  Output: {args.output}")
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Initialize manifest
    manifest_builder = ManifestBuilder(mode="tune", version="2.1")
    manifest_builder.set_git_info()
    manifest_builder.set_tool_versions(ToolVersionFetcher.get_versions())
    
    # Load config
    try:
        with open(args.config) as f:
            import yaml
            config = yaml.safe_load(f)
        manifest_builder.set_config_hash(config)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1
    
    # Run tuning (functional API)
    try:
        result = run_tuning(
            model=args.model,
            config=config,
            output_dir=args.output,
            n_jobs=args.n_jobs,
        )
        
        # Write manifest
        manifest_builder.record_gate("tuning_completed", True)
        manifest_path = args.output / "run_manifest.json"
        manifest_builder.write(manifest_path)
        
        print("[TUNE] Tuning completed successfully")
        print(f"  Best params: {args.output / 'best_params.yaml'}")
        print(f"  Summary: {args.output / 'tuning_summary.json'}")
        
        return 0
    
    except Exception as e:
        print(f"[ERROR] Tuning failed: {e}", file=sys.stderr)
        manifest_builder.record_gate("tuning_completed", False)
        manifest_path = args.output / "run_manifest.json"
        manifest_builder.write(manifest_path)
        return 1


def cmd_run(args):
    """Execute production mode subcommand."""
    print(f"[RUN] Starting {args.del} analysis...")
    print(f"  Config: {args.config}")
    print(f"  Output: {args.output}")
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Initialize manifest
    manifest_builder = ManifestBuilder(mode="production", version="2.1")
    manifest_builder.set_git_info()
    manifest_builder.set_tool_versions(ToolVersionFetcher.get_versions())
    
    # Load config
    try:
        with open(args.config) as f:
            import yaml
            config = yaml.safe_load(f)
        manifest_builder.set_config_hash(config)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1
    
    # Link to tuning results (if provided in config)
    tuning_ref = config.get("tuning_reference")
    if tuning_ref:
        manifest_builder.set_tuning_reference(tuning_ref)
    
    print("[RUN] Production mode (implementation placeholder)")
    print("  Steps: ingest → normalize → glycan_prep → placer → qc → analysis → reporting")
    
    # Write manifest
    manifest_builder.record_gate("production_pipeline_started", True)
    manifest_path = args.output / "run_manifest.json"
    manifest_builder.write(manifest_path)
    
    print(f"[RUN] Manifest written to {manifest_path}")
    
    # TODO: Implement full pipeline orchestration
    print("[RUN] (Full implementation pending)")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
