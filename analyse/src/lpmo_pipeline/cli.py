"""
LPMO Pipeline: CLI Entry Point
Responsibility: Main command-line interface with production run mode and optional post-analysis tuning mode.

Usage:
  python -m lpmo_pipeline run --mode production --config configs/production.yaml --output results/del_a/
  python -m lpmo_pipeline tune --config configs/tuning_af3.yaml --output results/tuning_af3/
  python -m lpmo_pipeline predictive --condition-table results/condition_table.tsv --protein-metadata metadata/protein_metadata.tsv --output results/
  python -m lpmo_pipeline cbm-paired --condition-table results/condition_table.tsv --cluster-table results/cluster_table.tsv --output results/
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from lpmo_pipeline.analysis.analysis_orchestrator import run_analysis_core
from lpmo_pipeline.analysis.cbm_comparison import run_cbm_paired_analysis
from lpmo_pipeline.analysis.family_enrichment_postprocess import run_family_enrichment_postprocess
from lpmo_pipeline.analysis.predictive_postprocess import run_predictive_postprocess
from lpmo_pipeline.config import load_pipeline_run_config, resolve_pipeline_command_config
from lpmo_pipeline.io.discovery import discover_work_root
from lpmo_pipeline.tuning.tune_orchestrator import run_tuning
from lpmo_pipeline.utils.manifest import ManifestBuilder, ToolVersionFetcher


def _command_config_from_args(args: argparse.Namespace, command_name: str) -> dict[str, Any]:
    run_config_path = getattr(args, "run_config", None)
    if run_config_path is None:
        return {}
    run_config = load_pipeline_run_config(run_config_path)
    return resolve_pipeline_command_config(run_config, command_name)


def _resolve_required(
    args: argparse.Namespace,
    command_config: dict[str, Any],
    arg_key: str,
    config_key: str,
    *,
    command_name: str,
) -> Any:
    arg_value = getattr(args, arg_key, None)
    if arg_value is not None:
        return arg_value
    config_value = command_config.get(config_key)
    if config_value is not None:
        return config_value
    raise ValueError(
        f"Missing required '{config_key}' for '{command_name}'. "
        f"Pass --{arg_key.replace('_', '-')} or set commands.{command_name.replace('-', '_')}.{config_key} in run config"
    )


def _resolve_optional(
    args: argparse.Namespace,
    command_config: dict[str, Any],
    arg_key: str,
    config_key: str,
) -> Any:
    arg_value = getattr(args, arg_key, None)
    if arg_value is not None:
        return arg_value
    return command_config.get(config_key)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="LPMO Structure Prediction → Analysis Pipeline (v2.1)"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Pipeline mode")
    
    # TUNE subcommand (optional post-analysis)
    tune_parser = subparsers.add_parser("tune", help="Run optional post-analysis parameter tuning")
    tune_parser.add_argument("--run-config", type=Path, default=None,
                             help="Optional hybrid run config YAML for config-first execution")
    tune_parser.add_argument("--model", required=False, choices=["AF3", "RF3", "Boltz2"],
                             help="Prediction model to tune")
    tune_parser.add_argument("--config", type=Path, required=False,
                             help="Tuning config YAML (parameter grid + dataset)")
    tune_parser.add_argument("--output", type=Path, required=False,
                             help="Output directory for tuning results")
    tune_parser.add_argument("--n-jobs", type=int, default=None,
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
        "--run-config", type=Path, default=None,
        help="Optional hybrid run config YAML for config-first execution",
    )
    predictive_parser.add_argument(
        "--condition-table", type=Path, required=False,
        help="Path to condition_table.tsv",
    )
    predictive_parser.add_argument(
        "--protein-metadata", type=Path, required=False,
        help="Path to protein_metadata.tsv",
    )
    predictive_parser.add_argument(
        "--output", type=Path, required=False,
        help="Output directory for predictive analysis artifacts",
    )
    predictive_parser.add_argument(
        "--task", choices=["all", "c1_c4", "substrate"], default=None,
        help="Which predictive task set to run",
    )
    predictive_parser.add_argument(
        "--n-folds", type=int, default=None,
        help="Target number of grouped CV folds",
    )
    predictive_parser.add_argument(
        "--random-state", type=int, default=None,
        help="Random seed for grouped CV and logistic regression",
    )

    cbm_parser = subparsers.add_parser(
        "cbm-paired",
        help="Run condition-level CBM full-length vs domain-only paired analysis",
    )
    cbm_parser.add_argument(
        "--run-config", type=Path, default=None,
        help="Optional hybrid run config YAML for config-first execution",
    )
    cbm_parser.add_argument(
        "--condition-table", type=Path, required=False,
        help="Path to condition_table.tsv",
    )
    cbm_parser.add_argument(
        "--cluster-table", type=Path,
        help="Optional path to cluster_table.tsv for C1/C4 compatible fractions",
    )
    cbm_parser.add_argument(
        "--cluster-residue-signature-table", type=Path,
        help="Optional path to cluster_residue_signature.tsv for non-core/bridge fractions",
    )
    cbm_parser.add_argument(
        "--protein-metadata", type=Path,
        help="Optional path to protein metadata TSV",
    )
    cbm_parser.add_argument(
        "--output", type=Path, required=False,
        help="Output directory for CBM paired-analysis artifacts",
    )
    cbm_parser.add_argument(
        "--random-state", type=int, default=None,
        help="Random seed for bootstrap confidence intervals",
    )

    family_parser = subparsers.add_parser(
        "family-enrichment",
        help="Run optional AA9/AA10 family residue enrichment postprocess",
    )
    family_parser.add_argument(
        "--run-config", type=Path, default=None,
        help="Optional hybrid run config YAML for config-first execution",
    )
    family_parser.add_argument(
        "--protein-condition-residue-scores", type=Path, required=False,
        help="Path to protein_condition_residue_scores.tsv",
    )
    family_parser.add_argument(
        "--protein-residue-regio-delta", type=Path, required=False,
        help="Path to protein_residue_regio_delta.tsv",
    )
    family_parser.add_argument(
        "--protein-metadata", type=Path, required=False,
        help="Path to protein metadata TSV",
    )
    family_parser.add_argument(
        "--core-fasta", type=Path, required=False,
        help="Path to the deduplicated catalytic-core FASTA",
    )
    family_parser.add_argument(
        "--output", type=Path, required=False,
        help="Output directory for family enrichment artifacts",
    )
    family_parser.add_argument(
        "--alignment-dir", type=Path, default=None,
        help="Optional directory with precomputed {AA9,AA10}.aligned.fasta files",
    )
    family_parser.add_argument(
        "--families", nargs="+", default=None,
        help="Family labels to include (default: AA9 AA10)",
    )
    family_parser.add_argument(
        "--substrates", nargs="+", default=None,
        help="Predicted substrate classes to evaluate for family substrate enrichment (default: cellulose chitin)",
    )
    family_parser.add_argument(
        "--mafft-executable", default=None,
        help="MAFFT executable used when precomputed alignments are absent",
    )

    # RUN subcommand
    run_parser = subparsers.add_parser("run", help="Run production analysis")
    run_parser.add_argument("--run-config", type=Path, default=None,
                            help="Optional hybrid run config YAML for config-first execution")
    run_parser.add_argument("--mode", required=False, choices=["production"],
                            help="Execution mode (fixed to production)")
    run_parser.add_argument("--config", type=Path, required=False,
                            help="Main config YAML (locked from tuning)")
    run_parser.add_argument("--output", type=Path, required=False,
                            help="Output directory for analysis results")
    run_parser.add_argument("--del", dest="del_branch", required=False, choices=["del_a", "del_b"],
                            help="Analysis branch (all proteins or CBM subset)")
    run_parser.add_argument("--n-jobs", type=int, default=None,
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
    elif args.command == "cbm-paired":
        return cmd_cbm_paired(args)
    elif args.command == "family-enrichment":
        return cmd_family_enrichment(args)
    else:
        parser.print_help()
        return 1


def cmd_tune(args):
    """Execute tuning subcommand."""
    try:
        command_config = _command_config_from_args(args, "tune")
        model = str(_resolve_required(args, command_config, "model", "model", command_name="tune"))
        config_path = Path(_resolve_required(args, command_config, "config", "config", command_name="tune"))
        output_dir = Path(_resolve_required(args, command_config, "output", "output", command_name="tune"))
        n_jobs = int(_resolve_optional(args, command_config, "n_jobs", "n_jobs") or 1)
    except Exception as exc:
        print(f"[ERROR] Tune argument resolution failed: {exc}", file=sys.stderr)
        return 1

    print(f"[TUNE] Starting optional post-analysis {model} tuning...")
    print(f"  Config: {config_path}")
    print(f"  Output: {output_dir}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize manifest
    manifest_builder = ManifestBuilder(mode="tune", version="2.1")
    manifest_builder.set_git_info()
    manifest_builder.set_tool_versions(ToolVersionFetcher.get_versions())
    
    # Load config
    try:
        with open(config_path) as f:
            import yaml
            config = yaml.safe_load(f)
        manifest_builder.set_config_hash(config)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1
    
    # Run tuning (functional API)
    try:
        result = run_tuning(
            model=model,
            tuning_config_path=config_path,
            test_cases=[],
            output_dir=output_dir,
            pipeline_config=config,
            max_parallel=n_jobs,
        )
        
        # Write manifest
        manifest_builder.record_gate("tuning_completed", True)
        manifest_path = output_dir / "run_manifest.json"
        manifest_builder.write(manifest_path)
        
        print("[TUNE] Tuning completed successfully")
        print(f"  Best params: {output_dir / 'best_params.yaml'}")
        print(f"  Summary: {output_dir / 'tuning_summary.json'}")
        
        return 0
    
    except Exception as e:
        print(f"[ERROR] Tuning failed: {e}", file=sys.stderr)
        manifest_builder.record_gate("tuning_completed", False)
        manifest_path = output_dir / "run_manifest.json"
        manifest_builder.write(manifest_path)
        return 1


def cmd_run(args):
    """Execute production mode subcommand."""
    try:
        command_config = _command_config_from_args(args, "run")
        mode = str(_resolve_required(args, command_config, "mode", "mode", command_name="run"))
        if mode != "production":
            raise ValueError(f"Unsupported mode '{mode}'. Only 'production' is valid")
        config_path = Path(_resolve_required(args, command_config, "config", "config", command_name="run"))
        output_dir = Path(_resolve_required(args, command_config, "output", "output", command_name="run"))
        del_branch = str(_resolve_required(args, command_config, "del_branch", "del", command_name="run"))
        if del_branch not in {"del_a", "del_b"}:
            raise ValueError("run.del must be one of: del_a, del_b")
        n_jobs = int(_resolve_optional(args, command_config, "n_jobs", "n_jobs") or 1)
    except Exception as exc:
        print(f"[ERROR] Run argument resolution failed: {exc}", file=sys.stderr)
        return 1

    print(f"[RUN] Starting {del_branch} analysis...")
    print(f"  Config: {config_path}")
    print(f"  Output: {output_dir}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize manifest
    manifest_builder = ManifestBuilder(mode="production", version="2.1")
    manifest_builder.set_git_info()
    manifest_builder.set_tool_versions(ToolVersionFetcher.get_versions())
    
    # Load config
    try:
        with open(config_path) as f:
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
    manifest_path = output_dir / "run_manifest.json"

    try:
        result = run_analysis_core(
            config=config,
            output_dir=output_dir,
            del_variant=del_branch,
            n_jobs=n_jobs,
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
        if ((config.get("production") or {}).get("predictive") or {}).get("enabled", False):
            manifest_builder.record_gate(
                "predictive_stage_completed",
                result.predictive_summary_path is not None,
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
        if result.predictive_summary_path is not None:
            print(f"  Predictive summary: {result.predictive_summary_path}")
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
    try:
        command_config = _command_config_from_args(args, "predictive")
        condition_table = Path(
            _resolve_required(
                args,
                command_config,
                "condition_table",
                "condition_table",
                command_name="predictive",
            )
        )
        protein_metadata = Path(
            _resolve_required(
                args,
                command_config,
                "protein_metadata",
                "protein_metadata",
                command_name="predictive",
            )
        )
        output_dir = Path(_resolve_required(args, command_config, "output", "output", command_name="predictive"))
        task = str(_resolve_optional(args, command_config, "task", "task") or "all")
        n_folds = int(_resolve_optional(args, command_config, "n_folds", "n_folds") or 5)
        random_state = int(_resolve_optional(args, command_config, "random_state", "random_state") or 42)
    except Exception as exc:
        print(f"[ERROR] Predictive argument resolution failed: {exc}", file=sys.stderr)
        return 1

    print("[PREDICTIVE] Starting predictive postprocess...")
    print(f"  Condition table: {condition_table}")
    print(f"  Protein metadata: {protein_metadata}")
    print(f"  Output: {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_predictive_postprocess(
            condition_table_path=condition_table,
            protein_metadata_path=protein_metadata,
            output_dir=output_dir,
            task=task,
            n_folds=n_folds,
            random_state=random_state,
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


def cmd_cbm_paired(args):
    """Execute CBM paired postprocess subcommand."""
    try:
        command_config = _command_config_from_args(args, "cbm-paired")
        condition_table = Path(
            _resolve_required(
                args,
                command_config,
                "condition_table",
                "condition_table",
                command_name="cbm-paired",
            )
        )
        cluster_table = _resolve_optional(args, command_config, "cluster_table", "cluster_table")
        cluster_residue_signature_table = _resolve_optional(
            args,
            command_config,
            "cluster_residue_signature_table",
            "cluster_residue_signature_table",
        )
        protein_metadata = _resolve_optional(args, command_config, "protein_metadata", "protein_metadata")
        output_dir = Path(
            _resolve_required(
                args,
                command_config,
                "output",
                "output",
                command_name="cbm-paired",
            )
        )
        random_state = int(_resolve_optional(args, command_config, "random_state", "random_state") or 42)
    except Exception as exc:
        print(f"[ERROR] CBM argument resolution failed: {exc}", file=sys.stderr)
        return 1

    cluster_table_path = Path(cluster_table) if cluster_table is not None else None
    cluster_residue_signature_table_path = (
        Path(cluster_residue_signature_table)
        if cluster_residue_signature_table is not None
        else None
    )
    protein_metadata_path = Path(protein_metadata) if protein_metadata is not None else None

    print("[CBM] Starting CBM paired postprocess...")
    print(f"  Condition table: {condition_table}")
    print(f"  Cluster table: {cluster_table_path}")
    print(f"  Cluster residue signature table: {cluster_residue_signature_table_path}")
    print(f"  Protein metadata: {protein_metadata_path}")
    print(f"  Output: {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_cbm_paired_analysis(
            condition_table_path=condition_table,
            cluster_table_path=cluster_table_path,
            cluster_residue_signature_table_path=cluster_residue_signature_table_path,
            protein_metadata_path=protein_metadata_path,
            output_dir=output_dir,
            random_state=random_state,
        )
    except Exception as exc:
        print(f"[ERROR] CBM paired postprocess failed: {exc}", file=sys.stderr)
        return 1

    print("[CBM] CBM paired postprocess completed")
    print(f"  Summary JSON: {result.summary_path}")
    for name, path in sorted(result.table_paths.items()):
        print(f"  Table ({name}): {path}")
    return 0


def cmd_family_enrichment(args):
    """Execute optional family enrichment postprocess."""
    try:
        command_config = _command_config_from_args(args, "family-enrichment")
        protein_condition_residue_scores = Path(
            _resolve_required(
                args,
                command_config,
                "protein_condition_residue_scores",
                "protein_condition_residue_scores",
                command_name="family-enrichment",
            )
        )
        protein_residue_regio_delta = Path(
            _resolve_required(
                args,
                command_config,
                "protein_residue_regio_delta",
                "protein_residue_regio_delta",
                command_name="family-enrichment",
            )
        )
        protein_metadata = Path(
            _resolve_required(
                args,
                command_config,
                "protein_metadata",
                "protein_metadata",
                command_name="family-enrichment",
            )
        )
        core_fasta = Path(
            _resolve_required(
                args,
                command_config,
                "core_fasta",
                "core_fasta",
                command_name="family-enrichment",
            )
        )
        output_dir = Path(
            _resolve_required(
                args,
                command_config,
                "output",
                "output",
                command_name="family-enrichment",
            )
        )
        alignment_dir = _resolve_optional(args, command_config, "alignment_dir", "alignment_dir")
        families = _resolve_optional(args, command_config, "families", "families") or ["AA9", "AA10"]
        substrates = _resolve_optional(args, command_config, "substrates", "substrates") or [
            "cellulose",
            "chitin",
        ]
        mafft_executable = str(
            _resolve_optional(args, command_config, "mafft_executable", "mafft_executable") or "mafft"
        )
    except Exception as exc:
        print(f"[ERROR] Family enrichment argument resolution failed: {exc}", file=sys.stderr)
        return 1

    alignment_path = Path(alignment_dir) if alignment_dir is not None else None

    print("[FAMILY] Starting optional family residue enrichment postprocess...")
    print(f"  Protein-condition residue scores: {protein_condition_residue_scores}")
    print(f"  Protein residue regio delta: {protein_residue_regio_delta}")
    print(f"  Protein metadata: {protein_metadata}")
    print(f"  Core FASTA: {core_fasta}")
    print(f"  Output: {output_dir}")
    print(f"  Substrate enrichment targets: {', '.join(str(value) for value in substrates)}")

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_family_enrichment_postprocess(
            protein_condition_residue_scores_path=protein_condition_residue_scores,
            protein_residue_regio_delta_path=protein_residue_regio_delta,
            protein_metadata_path=protein_metadata,
            core_fasta_path=core_fasta,
            output_dir=output_dir,
            alignment_dir=alignment_path,
            families=tuple(str(value) for value in families),
            substrate_classes=tuple(str(value) for value in substrates),
            mafft_executable=mafft_executable,
        )
    except Exception as exc:
        print(f"[ERROR] Family enrichment postprocess failed: {exc}", file=sys.stderr)
        return 1

    print("[FAMILY] Family enrichment postprocess completed")
    print(f"  Family aligned residue table: {result.family_aligned_residue_table_path}")
    print(f"  Family residue enrichment: {result.family_residue_enrichment_path}")
    print(f"  Family substrate residue enrichment: {result.family_substrate_residue_enrichment_path}")
    print(f"  Family wrong-ligand residue enrichment: {result.family_wrong_ligand_residue_enrichment_path}")
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
