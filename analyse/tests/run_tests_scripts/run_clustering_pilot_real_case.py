#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

from lpmo_pipeline.analysis.clustering_pilot_real_case import (
    prepare_real_case_pilot_run,
    summarize_pilot_execution,
)
from lpmo_pipeline.cli import cmd_run


DEFAULT_SELECTION_MANIFEST_PATH = Path(__file__).with_name(
    "clustering_pilot_real_case_selection.yaml"
)
DEFAULT_ANALYSE_ENV_PYTHON = Path(
    "/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare or execute a resumable real-case clustering-pilot run from a "
            "protein-target manifest."
        ),
    )
    parser.add_argument("--run-dir", required=True, help="Output directory for checkpoints and artifacts.")
    parser.add_argument(
        "--selection-manifest",
        default=str(DEFAULT_SELECTION_MANIFEST_PATH),
        help="YAML manifest describing the real-case protein-target selections.",
    )
    parser.add_argument("--run-id", default="", help="Optional run identifier written into the config.")
    parser.add_argument(
        "--del-branch",
        default="del_a",
        choices=["del_a", "del_b"],
        help="Analysis branch to run through the production CLI.",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually run the production pipeline. By default only prepare checkpoints and config.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel workers for production execution.",
    )
    parser.add_argument(
        "--force-step",
        action="append",
        default=[],
        choices=["discovery", "staging", "config", "run", "summary"],
        help="Force one step to rebuild even if its checkpoint exists.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    run_dir = Path(args.run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.selection_manifest).resolve()
    run_id = args.run_id or run_dir.name
    summary_path = run_dir / "clustering_pilot_real_case_summary.json"

    summary: dict[str, object] = {
        "run_dir": str(run_dir),
        "selection_manifest_path": str(manifest_path),
        "run_id": run_id,
        "execute": bool(args.execute),
        "n_jobs": max(1, int(args.n_jobs)),
    }

    try:
        prepared = prepare_real_case_pilot_run(
            run_dir=run_dir,
            manifest_path=manifest_path,
            run_id=run_id,
            del_branch=args.del_branch,
            python_executable=(
                DEFAULT_ANALYSE_ENV_PYTHON
                if DEFAULT_ANALYSE_ENV_PYTHON.exists()
                else Path(sys.executable).resolve()
            ),
            script_path=Path(__file__).resolve(),
            force_steps=args.force_step,
            n_jobs=max(1, int(args.n_jobs)),
        )
        summary.update(prepared)

        if not prepared.get("ready_to_run", False):
            summary["run_step_status"] = "blocked_missing_selection_data"
            summary_path.write_text(json.dumps(summary, indent=2))
            print(json.dumps(summary, indent=2))
            return 1

        production_output = Path(str(prepared["production_output"]))
        analysis_summary_path = production_output / "analysis_core_summary.json"

        if args.execute:
            should_rerun = (
                "run" in set(args.force_step)
                or bool(prepared.get("manifest_changed"))
                or bool(prepared.get("config_changed"))
                or not analysis_summary_path.exists()
            )
            if should_rerun:
                exit_code = cmd_run(
                    argparse.Namespace(
                        mode="production",
                        config=Path(str(prepared["config_path"])),
                        output=production_output,
                        del_branch=args.del_branch,
                        n_jobs=max(1, int(args.n_jobs)),
                    )
                )
                summary["cli_exit_code"] = exit_code
                summary["run_step_status"] = "executed"
                if exit_code != 0:
                    summary_path.write_text(json.dumps(summary, indent=2))
                    print(json.dumps(summary, indent=2))
                    return exit_code
            else:
                summary["cli_exit_code"] = 0
                summary["run_step_status"] = "reused_existing_output"

            execution_summary, summary_reused = summarize_pilot_execution(
                run_dir=run_dir,
                production_output=production_output,
                force=("summary" in set(args.force_step) or summary["run_step_status"] == "executed"),
            )
            summary.update(execution_summary)
            summary["summary_step_status"] = "reused" if summary_reused else "written"
            summary_path.write_text(json.dumps(summary, indent=2))
            print(json.dumps(summary, indent=2))
            return 0 if summary.get("pilot_outputs_ready", False) else 1

        summary["run_step_status"] = "prepare_only"
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 0
    except Exception as exc:
        summary["error"] = f"{exc.__class__.__name__}: {exc}"
        summary["traceback"] = traceback.format_exc()
        summary_path.write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
