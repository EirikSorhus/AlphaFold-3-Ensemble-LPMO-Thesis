#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from lpmo_pipeline.analysis.clustering_pilot_real_case import (
    build_selection_manifest_from_overview,
    load_pilot_protein_overview_tsv,
    write_selection_manifest,
)


DEFAULT_OVERVIEW_TSV_PATH = (
    Path(__file__).resolve().parents[2] / "input_data" / "clustering_pilot_protein_overview.tsv"
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build domain-only and full-length pilot selection manifests from the pilot protein overview TSV.",
    )
    parser.add_argument(
        "--overview-tsv",
        default=str(DEFAULT_OVERVIEW_TSV_PATH),
        help="Curated TSV overview of pilot proteins.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where the generated manifests and summary should be written.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    overview_tsv = Path(args.overview_tsv).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    overview_rows = load_pilot_protein_overview_tsv(overview_tsv)

    shared_clustering_pilot = {
        "minimum_clusterable_n": 10,
        "agglomerative_linkage": "average",
        "agglomerative_distance_threshold": 0.5,
        "agglomerative_min_cluster_size": 10,
    }
    domain_manifest = build_selection_manifest_from_overview(
        overview_rows,
        construct_type="domain_only",
        clustering_pilot={
            **shared_clustering_pilot,
            "label": "real_case_selected_proteins_domain_only",
        },
    )
    full_manifest = build_selection_manifest_from_overview(
        overview_rows,
        construct_type="full_length",
        clustering_pilot={
            **shared_clustering_pilot,
            "label": "real_case_selected_proteins_full_length",
        },
    )

    domain_manifest_path = output_dir / "clustering_pilot_domain_only_selection.yaml"
    full_manifest_path = output_dir / "clustering_pilot_full_length_selection.yaml"
    write_selection_manifest(domain_manifest, domain_manifest_path)
    write_selection_manifest(full_manifest, full_manifest_path)

    summary = {
        "overview_tsv": str(overview_tsv),
        "n_overview_rows": len(overview_rows),
        "domain_only_manifest_path": str(domain_manifest_path),
        "full_length_manifest_path": str(full_manifest_path),
        "n_domain_only_selections": len(domain_manifest.selections),
        "n_full_length_selections": len(full_manifest.selections),
        "domain_only_proteins": [selection.protein_id for selection in domain_manifest.selections],
        "full_length_proteins": [selection.protein_id for selection in full_manifest.selections],
    }
    summary_path = output_dir / "clustering_pilot_manifest_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())