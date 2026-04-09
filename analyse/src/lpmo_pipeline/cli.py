"""
LPMO Pipeline: CLI Entry Point
Responsibility: Main command-line interface with production run mode and optional post-analysis tuning mode.

Usage:
  python -m lpmo_pipeline run --mode production --config configs/production.yaml --output results/del_a/
    python -m lpmo_pipeline tune --config configs/tuning_af3.yaml --output results/tuning_af3/
"""

import argparse
import json
import sys
from pathlib import Path

from lpmo_pipeline.io.discovery import discover_work_root
from lpmo_pipeline.tuning.tune_orchestrator import run_tuning
from lpmo_pipeline.utils.manifest import ManifestBuilder, ToolVersionFetcher


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="LPMO Structure Prediction → Analysis Pipeline (v2.1)"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Pipeline mode")
    
    # TUNE subcommand (optional post-analysis)
    tune_parser = subparsers.add_parser("tune", help="Run optional post-analysis parameter tuning")
    tune_parser.add_argument("--model", required=True, choices=["AF3", "RF3", "Boltz2"],
                             help="Prediction model to tune")
    tune_parser.add_argument("--config", type=Path, required=True,
                             help="Tuning config YAML (parameter grid + dataset)")
    tune_parser.add_argument("--output", type=Path, required=True,
                             help="Output directory for tuning results")
    tune_parser.add_argument("--n-jobs", type=int, default=1,
                             help="Number of parallel jobs")
    
    # DISCOVER subcommand
    disc_parser = subparsers.add_parser(
        "discover",
        help="Discover AF3/RF3 prediction artifacts in a work root",
    )
    disc_parser.add_argument(
        "--work-root", type=Path, required=True,
        help="Path to structure_pipeline work/ directory",
    )
    disc_parser.add_argument(
        "--output", type=Path, default=None,
        help="Write discovery manifest JSON to this file (default: stdout)",
    )
    disc_parser.add_argument(
        "--summary-only", action="store_true",
        help="Print only the summary, not the full manifest",
    )
    disc_parser.add_argument(
        "--af3-only", action="store_true",
        help="Only discover AF3 predictions (skip RF3)",
    )
    disc_parser.add_argument(
        "--latest-only", action="store_true",
        help="Only discover the run pointed to by the 'latest' symlink",
    )

    # RUN subcommand
    run_parser = subparsers.add_parser("run", help="Run production analysis")
    run_parser.add_argument("--mode", required=True, choices=["production"],
                            help="Execution mode (fixed to production)")
    run_parser.add_argument("--config", type=Path, required=True,
                            help="Main config YAML (locked from tuning)")
    run_parser.add_argument("--output", type=Path, required=True,
                            help="Output directory for analysis results")
    run_parser.add_argument("--del", dest="del_branch", required=True, choices=["del_a", "del_b"],
                            help="Analysis branch (all proteins or CBM subset)")
    run_parser.add_argument("--n-jobs", type=int, default=1,
                            help="Number of parallel jobs")
    
    args = parser.parse_args()
    
    if args.command == "tune":
        return cmd_tune(args)
    elif args.command == "run":
        return cmd_run(args)
    elif args.command == "discover":
        return cmd_discover(args)
    else:
        parser.print_help()
        return 1


def cmd_tune(args):
    """Execute tuning subcommand."""
    print(f"[TUNE] Starting optional post-analysis {args.model} tuning...")
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
    print(f"[RUN] Starting {args.del_branch} analysis...")
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
    print("  Steps: ingest -> normalize -> glycan_prep -> placer -> pre_qc_proximity -> qc -> analysis -> reporting")
    
    # Write manifest
    manifest_builder.record_gate("production_pipeline_started", True)
    manifest_path = args.output / "run_manifest.json"
    manifest_builder.write(manifest_path)
    
    print(f"[RUN] Manifest written to {manifest_path}")
    
    # TODO: Implement full pipeline orchestration
    print("[RUN] (Full implementation pending)")
    
    return 0


def cmd_discover(args):
    """Execute discovery subcommand: scan work root for AF3/RF3 artifacts."""
    work_root = args.work_root.resolve()
    if not work_root.is_dir():
        print(f"[ERROR] Not a directory: {work_root}", file=sys.stderr)
        return 1

    manifest = discover_work_root(
        work_root,
        af3_only=args.af3_only,
        latest_only=args.latest_only,
    )

    if args.summary_only:
        output_text = json.dumps(manifest.summary, indent=2)
    else:
        output_text = manifest.to_json()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_text + "\n")
        print(f"[DISCOVER] Manifest written to {args.output}")
        # Always print summary to stderr for visibility
        print(
            json.dumps(manifest.summary, indent=2),
            file=sys.stderr,
        )
    else:
        print(output_text)

    if manifest.errors:
        print(
            f"[DISCOVER] {len(manifest.errors)} warning(s) — see 'errors' in manifest",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
