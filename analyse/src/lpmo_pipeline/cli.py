"""
LPMO Pipeline: CLI Entry Point
Responsibility: Main command-line interface with production run mode and optional post-analysis tuning mode.

Usage:
  python -m lpmo_pipeline run --mode production --config configs/production.yaml --output results/del_a/
    python -m lpmo_pipeline tune --config configs/tuning_af3.yaml --output results/tuning_af3/
    python -m lpmo_pipeline predictive --condition-table results/condition_table.tsv --protein-metadata metadata/protein_metadata.tsv --output results/
"""

import argparse
import json
import sys
from pathlib import Path

from lpmo_pipeline.analysis.analysis_orchestrator import run_analysis_core
from lpmo_pipeline.analysis.family_enrichment_postprocess import run_family_enrichment_postprocess
from lpmo_pipeline.analysis.predictive_postprocess import run_predictive_postprocess
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

    # PREDICTIVE subcommand
    predictive_parser = subparsers.add_parser(
        "predictive",
        help="Run predictive postprocess over summary tables",
    )
    predictive_parser.add_argument(
        "--condition-table", type=Path, required=True,
        help="Path to condition_table.tsv",
    )
    predictive_parser.add_argument(
        "--protein-metadata", type=Path, required=True,
        help="Path to protein_metadata.tsv",
    )
    predictive_parser.add_argument(
        "--output", type=Path, required=True,
        help="Output directory for predictive analysis artifacts",
    )
    predictive_parser.add_argument(
        "--task", choices=["all", "c1_c4", "substrate"], default="all",
        help="Which predictive task set to run",
    )
    predictive_parser.add_argument(
        "--n-folds", type=int, default=5,
        help="Target number of grouped CV folds",
    )
    predictive_parser.add_argument(
        "--random-state", type=int, default=42,
        help="Random seed for grouped CV and logistic regression",
    )

    family_parser = subparsers.add_parser(
        "family-enrichment",
        help="Run optional AA9/AA10 family residue enrichment postprocess",
    )
    family_parser.add_argument(
        "--protein-condition-residue-scores", type=Path, required=True,
        help="Path to protein_condition_residue_scores.tsv",
    )
    family_parser.add_argument(
        "--protein-residue-regio-delta", type=Path, required=True,
        help="Path to protein_residue_regio_delta.tsv",
    )
    family_parser.add_argument(
        "--protein-metadata", type=Path, required=True,
        help="Path to protein metadata TSV",
    )
    family_parser.add_argument(
        "--core-fasta", type=Path, required=True,
        help="Path to the deduplicated catalytic-core FASTA",
    )
    family_parser.add_argument(
        "--output", type=Path, required=True,
        help="Output directory for family enrichment artifacts",
    )
    family_parser.add_argument(
        "--alignment-dir", type=Path, default=None,
        help="Optional directory with precomputed {AA9,AA10}.aligned.fasta files",
    )
    family_parser.add_argument(
        "--families", nargs="+", default=["AA9", "AA10"],
        help="Family labels to include (default: AA9 AA10)",
    )
    family_parser.add_argument(
        "--mafft-executable", default="mafft",
        help="MAFFT executable used when precomputed alignments are absent",
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
    elif args.command == "predictive":
        return cmd_predictive(args)
    elif args.command == "family-enrichment":
        return cmd_family_enrichment(args)
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

    manifest_builder.record_gate("production_pipeline_started", True)
    manifest_path = args.output / "run_manifest.json"

    try:
        result = run_analysis_core(
            config=config,
            output_dir=args.output,
            del_variant=args.del_branch,
            n_jobs=args.n_jobs,
        )
        manifest_builder.record_gate("analysis_core_completed", result.success)
        manifest_builder.record_gate("hard_qc_completed", result.qc_report_path is not None)
        manifest_builder.record_gate(
            "geometry_stage_completed",
            result.pose_geometry_tsv_path is not None,
        )
        manifest_builder.record_gate(
            "ifp_stage_completed",
            result.pose_ifp_table_tsv_path is not None,
        )
        manifest_builder.record_gate(
            "clustering_stage_completed",
            result.condition_cluster_summary_tsv_path is not None,
        )
        manifest_builder.record_gate(
            "cluster_annotation_stage_completed",
            result.cluster_annotation_stage_completed,
        )
        manifest_builder.record_gate(
            "crystal_anchoring_stage_completed",
            result.crystal_anchoring_stage_completed,
        )
        manifest_builder.write(manifest_path)

        if not result.success:
            print("[ERROR] Production analysis finished without any prepared poses", file=sys.stderr)
            return 1

        print("[RUN] Analysis core completed")
        if result.pose_manifest_tsv_path is not None:
            print(f"  Pose manifest: {result.pose_manifest_tsv_path}")
        if result.pose_confidence_tsv_path is not None:
            print(f"  Pose confidence: {result.pose_confidence_tsv_path}")
        if result.structure_index_tsv_path is not None:
            print(f"  Structure index: {result.structure_index_tsv_path}")
        if result.qc_attrition_tsv_path is not None:
            print(f"  QC attrition: {result.qc_attrition_tsv_path}")
        if result.qc_report_path is not None:
            print(f"  QC report: {result.qc_report_path}")
        if result.pose_geometry_tsv_path is not None:
            print(f"  Pose geometry: {result.pose_geometry_tsv_path}")
        if result.pose_ifp_table_tsv_path is not None:
            print(f"  Pose IFP table: {result.pose_ifp_table_tsv_path}")
        if result.cluster_assignments_tsv_path is not None:
            print(f"  Cluster assignments: {result.cluster_assignments_tsv_path}")
        if result.medoid_manifest_tsv_path is not None:
            print(f"  Medoid manifest: {result.medoid_manifest_tsv_path}")
        if result.condition_cluster_summary_tsv_path is not None:
            print(f"  Condition cluster summary: {result.condition_cluster_summary_tsv_path}")
        if result.cluster_table_tsv_path is not None:
            print(f"  Cluster table: {result.cluster_table_tsv_path}")
        if result.cluster_ifp_signature_tsv_path is not None:
            print(f"  Cluster IFP signature: {result.cluster_ifp_signature_tsv_path}")
        if result.cluster_residue_signature_tsv_path is not None:
            print(f"  Cluster residue signature: {result.cluster_residue_signature_tsv_path}")
        if result.cluster_signatures_json_path is not None:
            print(f"  Cluster signatures JSON: {result.cluster_signatures_json_path}")
        if result.condition_table_tsv_path is not None:
            print(f"  Condition table: {result.condition_table_tsv_path}")
        if result.protein_summary_table_tsv_path is not None:
            print(f"  Protein summary table: {result.protein_summary_table_tsv_path}")
        if result.crystal_anchor_tsv_path is not None:
            print(f"  Crystal anchor table: {result.crystal_anchor_tsv_path}")
        if result.metrics_csv_path is not None:
            print(f"  Metrics CSV: {result.metrics_csv_path}")
        if result.summary_json_path is not None:
            print(f"  Summary JSON: {result.summary_json_path}")
        if result.report_html_path is not None:
            print(f"  HTML report: {result.report_html_path}")
        print(f"  Run summary: {result.summary_path}")
        print(f"  Manifest: {manifest_path}")
        return 0

    except Exception as e:
        print(f"[ERROR] Production run failed: {e}", file=sys.stderr)
        manifest_builder.record_gate("analysis_core_completed", False)
        manifest_builder.record_gate("crystal_anchoring_stage_completed", False)
        manifest_builder.write(manifest_path)
        return 1


def cmd_predictive(args):
    """Execute predictive postprocess subcommand."""
    print("[PREDICTIVE] Starting predictive postprocess...")
    print(f"  Condition table: {args.condition_table}")
    print(f"  Protein metadata: {args.protein_metadata}")
    print(f"  Output: {args.output}")

    args.output.mkdir(parents=True, exist_ok=True)

    try:
        result = run_predictive_postprocess(
            condition_table_path=args.condition_table,
            protein_metadata_path=args.protein_metadata,
            output_dir=args.output,
            task=args.task,
            n_folds=args.n_folds,
            random_state=args.random_state,
        )
    except Exception as exc:
        print(f"[ERROR] Predictive postprocess failed: {exc}", file=sys.stderr)
        return 1

    print("[PREDICTIVE] Predictive postprocess completed")
    print(f"  Summary JSON: {result.summary_path}")
    for name, path in sorted(result.modeling_table_paths.items()):
        print(f"  Modeling table ({name}): {path}")
    for name, path in sorted(result.metrics_paths.items()):
        print(f"  Metrics ({name}): {path}")
    for name, path in sorted(result.predictions_paths.items()):
        print(f"  Predictions ({name}): {path}")
    return 0


def cmd_family_enrichment(args):
    """Execute optional family enrichment postprocess."""
    print("[FAMILY] Starting optional family residue enrichment postprocess...")
    print(f"  Protein-condition residue scores: {args.protein_condition_residue_scores}")
    print(f"  Protein residue regio delta: {args.protein_residue_regio_delta}")
    print(f"  Protein metadata: {args.protein_metadata}")
    print(f"  Core FASTA: {args.core_fasta}")
    print(f"  Output: {args.output}")

    args.output.mkdir(parents=True, exist_ok=True)

    try:
        result = run_family_enrichment_postprocess(
            protein_condition_residue_scores_path=args.protein_condition_residue_scores,
            protein_residue_regio_delta_path=args.protein_residue_regio_delta,
            protein_metadata_path=args.protein_metadata,
            core_fasta_path=args.core_fasta,
            output_dir=args.output,
            alignment_dir=args.alignment_dir,
            families=tuple(args.families),
            mafft_executable=args.mafft_executable,
        )
    except Exception as exc:
        print(f"[ERROR] Family enrichment postprocess failed: {exc}", file=sys.stderr)
        return 1

    print("[FAMILY] Family enrichment postprocess completed")
    print(f"  Family aligned residue table: {result.family_aligned_residue_table_path}")
    print(f"  Family residue enrichment: {result.family_residue_enrichment_path}")
    print(f"  Alignment manifest: {result.alignment_manifest_path}")
    print(f"  Summary: {result.summary_path}")
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
