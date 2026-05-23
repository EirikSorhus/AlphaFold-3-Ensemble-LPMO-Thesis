# src/lpmo_pipeline/report/build_report_html.py
"""
Responsibility: Generate a human-readable HTML report from summary + metrics.
Input:  summary.json payload, metrics.csv, optional tuning summary.
Output: report.html with descriptive QC, clustering, geometry and crystal tables.
"""
from __future__ import annotations

import csv
from collections import Counter
from html import escape
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def _read_metrics_rows(metrics_csv_path: Path) -> list[dict[str, str]]:
    if not metrics_csv_path.exists():
        return []
    with open(metrics_csv_path, newline="") as handle:
        return list(csv.DictReader(handle))


def _non_empty(value: Any) -> bool:
    return value not in {None, "", "None", "nan", "NaN"}


def _to_float(value: Any) -> float | None:
    if not _non_empty(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fmt(value: Any, *, digits: int = 3) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    if isinstance(value, (list, tuple)):
        return " - ".join(_fmt(item, digits=digits) for item in value)
    return str(value)


def _pct(value: Any) -> str:
    numeric = _to_float(value)
    if numeric is None:
        return "N/A"
    return f"{numeric:.1%}"


def _numeric_range(rows: list[dict[str, str]], column: str) -> str:
    values = [_to_float(row.get(column)) for row in rows]
    values = [value for value in values if value is not None]
    if not values:
        return "N/A"
    return f"{min(values):.3f} - {max(values):.3f}"


def _counter_rows(counter: Counter[str], *, empty_label: str = "No rows") -> list[tuple[str, Any]]:
    if not counter:
        return [(empty_label, 0)]
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))


def _table(headers: list[str], rows: list[list[Any]], *, class_name: str = "") -> str:
    class_attr = f' class="{escape(class_name)}"' if class_name else ""
    head = "".join(f"<th>{escape(str(header))}</th>" for header in headers)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{escape(_fmt(cell))}</td>" for cell in row)
        body_rows.append(f"<tr>{cells}</tr>")
    return f"<table{class_attr}><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _key_value_table(items: list[tuple[str, Any]]) -> str:
    return _table(["Metric", "Value"], [[key, value] for key, value in items], class_name="kv")


def _render_tuning_section(tuning_data: dict[str, Any] | None) -> str:
    if not tuning_data:
        return ""
    rows = []
    for row in tuning_data.get("table_data", []):
        marker = " *" if row.get("is_best") else ""
        rows.append(
            [
                row.get("phase", ""),
                f"{row.get('param_hash', '')}{marker}",
                _fmt(row.get("qc_pass_rate"), digits=2),
                _fmt(row.get("outlier_rate"), digits=3),
                _fmt(row.get("mean_n_clusters"), digits=1),
                _fmt(row.get("crystal_sim"), digits=3),
            ]
        )
    return (
        "<section><h2>Tuning Results</h2>"
        f"<p>{escape(str(tuning_data.get('decision_rationale', '')))}</p>"
        + _table(["Phase", "Params", "QC pass", "Outlier", "Clusters", "Crystal"], rows)
        + "</section>"
    )


def build_report_html(
    summary: dict[str, Any],
    metrics_csv_path: Path,
    output_path: Path,
    tuning_data: dict[str, Any] | None = None,
) -> None:
    """Generate an HTML report from the current production summary surfaces."""
    metrics_rows = _read_metrics_rows(metrics_csv_path)
    dataset = summary.get("dataset_stats", {})
    qc = summary.get("qc_stats", {})
    clusters = summary.get("cluster_stats", {})
    geometry = summary.get("geometry_stats", {})
    crystal = summary.get("crystal_stats", {})

    qc_rows = [
        ("Total QC poses", qc.get("total_poses_qc", 0)),
        ("Passed", qc.get("passed", 0)),
        ("Flagged", qc.get("flagged", 0)),
        ("Dropped", qc.get("dropped", 0)),
        ("Pass rate", _pct(qc.get("pass_rate", 0.0))),
    ]
    dataset_rows = [
        ("Proteins", dataset.get("n_proteins", 0)),
        ("Ligands", dataset.get("n_ligands", 0)),
        ("Models", dataset.get("n_models", 0)),
        ("Total poses", dataset.get("n_total_poses", 0)),
        ("Rows in metrics.csv", len(metrics_rows)),
    ]
    cluster_rows = [
        ("Mean clusters per condition", _fmt(clusters.get("mean_n_clusters", 0.0), digits=2)),
        ("Mean noise/outlier rate", _pct(clusters.get("mean_outlier_rate", 0.0))),
        ("Condition count", clusters.get("n_protein_ligand_combinations", 0)),
    ]
    geometry_rows = [
        ("Cu-C1 range from summary", geometry.get("cu_c1_range") or "N/A"),
        ("Cu-C4 range from summary", geometry.get("cu_c4_range") or "N/A"),
        ("Cu-C1 range from metrics.csv", _numeric_range(metrics_rows, "min_cu_c1")),
        ("Cu-C4 range from metrics.csv", _numeric_range(metrics_rows, "min_cu_c4")),
        ("Geometry measurements in summary", geometry.get("n_measurements", 0)),
    ]
    crystal_rows = [
        ("Crystal comparison entries", crystal.get("n_comparisons", 0)),
        ("Mean best IFP Tanimoto", _fmt(crystal.get("mean_best_tanimoto", 0.0), digits=3)),
        ("Pocket RMSD range from metrics.csv", _numeric_range(metrics_rows, "pocket_rmsd_vs_crystal")),
        ("IFP similarity range from metrics.csv", _numeric_range(metrics_rows, "ifp_similarity_crystal")),
    ]

    qc_status_counts = Counter(row.get("qc_status", "unknown") or "unknown" for row in metrics_rows)
    cluster_counts = Counter(row.get("cluster_id", "") or "unassigned" for row in metrics_rows)
    contact_values = [_to_float(row.get("n_ifp_contacts")) for row in metrics_rows]
    contact_values = [value for value in contact_values if value is not None]
    contact_rows = [
        ("Mean IFP contacts", sum(contact_values) / len(contact_values) if contact_values else None),
        ("Max IFP contacts", max(contact_values) if contact_values else None),
        ("Rows with crystal similarity", sum(1 for row in metrics_rows if _non_empty(row.get("ifp_similarity_crystal")))),
    ]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>LPMO Pipeline Report - {escape(str(summary.get("run_id", "")))}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; line-height: 1.45; color: #1f2933; }}
    h1 {{ margin-bottom: 0.2rem; }}
    h2 {{ margin-top: 2rem; border-bottom: 1px solid #d5d9df; padding-bottom: 0.35rem; }}
    table {{ border-collapse: collapse; margin: 0.75rem 0 1.25rem; width: min(100%, 920px); }}
    th, td {{ border: 1px solid #d5d9df; padding: 0.45rem 0.65rem; text-align: left; }}
    th {{ background: #eef2f6; }}
    .meta {{ color: #52616f; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem 2rem; }}
    .note {{ background: #f7f9fb; border-left: 4px solid #527da3; padding: 0.75rem 1rem; max-width: 920px; }}
  </style>
</head>
<body>
  <h1>LPMO Structure Prediction Pipeline Report</h1>
  <p class="meta"><strong>Run:</strong> {escape(str(summary.get("run_id", "")))} |
  <strong>Mode:</strong> {escape(str(summary.get("mode", "")))} |
  <strong>DEL:</strong> {escape(str(summary.get("del_variant", "")))} |
  <strong>Generated:</strong> {escape(str(summary.get("timestamp", "")))}</p>

  <section class="grid">
    <div><h2>Dataset</h2>{_key_value_table(dataset_rows)}</div>
    <div><h2>QC Summary</h2>{_key_value_table(qc_rows)}</div>
  </section>

  <section class="grid">
    <div><h2>Cluster Landscape</h2>{_key_value_table(cluster_rows)}</div>
    <div><h2>IFP Contacts</h2>{_key_value_table(contact_rows)}</div>
  </section>

  <section>
    <h2>Geometry</h2>
    <p class="note">Downstream geometry is descriptive. Hard-QC gates are applied earlier; unavailable
    downstream geometry remains reportable as missing/not computable rather than removing retained poses.</p>
    {_key_value_table(geometry_rows)}
  </section>

  <section>
    <h2>Crystal Comparison</h2>
    {_key_value_table(crystal_rows)}
  </section>

  <section class="grid">
    <div><h2>QC Statuses In metrics.csv</h2>{_table(["QC status", "Rows"], _counter_rows(qc_status_counts))}</div>
    <div><h2>Cluster IDs In metrics.csv</h2>{_table(["Cluster ID", "Rows"], _counter_rows(cluster_counts))}</div>
  </section>

  {_render_tuning_section(tuning_data)}

  <hr>
  <p class="meta">Generated by LPMO Pipeline {escape(str(summary.get("pipeline_version", "")))}</p>
</body>
</html>
"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    logger.info("Wrote HTML report to %s", output_path)
