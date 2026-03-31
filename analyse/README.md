# LPMO Pipeline

Analyse-pipeline for LPMO enzyme–oligosaccharide structure predictions
from AlphaFold 3 only (RosettaFold 3 / RF3 og Boltz-2 er ekskludert fra analyse).

## Overview

This pipeline takes predicted structures (mmCIF), normalizes them,
runs quality control (PoseBusters, Privateer, custom Cu-geometry),
refines with PLACER, clusters interaction fingerprints (HDBSCAN),
and produces a comprehensive analysis report.

See [MASTERPLAN.md](MASTERPLAN.md) for full specification.
For pseudocode-only workflow, see [PSEUDOKODE_CLUSTER_FIRST.md](PSEUDOKODE_CLUSTER_FIRST.md).
Module-by-module pseudocode is in [PSEUDOKODE_MODUL_FOR_MODUL.md](PSEUDOKODE_MODUL_FOR_MODUL.md).
Statistics pseudocode is in [PSEUDOKODE_STATISTIKK.md](PSEUDOKODE_STATISTIKK.md).

## Execution Model

- Structure prediction (AF3 only) is run as a separate stage from analysis. RF3/RosettaFold 3 and Boltz-2 are excluded.
- Prediction outputs are generated once per chosen parameter setup and reused by analysis.
- Main analyses are run first; parameter tuning is optional and done after baseline analysis.
- Main analysis is cluster-primary (not enzyme-aggregated) for descriptive and predictive modeling.
- Operative sub-analyses are substrate x DP combinations for DP4, DP6, DP8.
- A pre-QC active-site proximity gate runs before PoseBusters and Privateer.
- Detailed geometric planarity requirements are currently missing and must be specified before final reporting.
- R is preferred for descriptive/predictive statistics where practical.

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Production mode (run full pipeline)
lpmo-pipeline run --config configs/defaults.yaml --systems-csv systems.csv

# Optional post-analysis tuning (comparative reruns)
lpmo-pipeline tune --config configs/tuning_af3.yaml

# EC metadata -> activity labels (R helper script)
Rscript scripts/ec_activity_mapping.R \
  --input metadata/enzyme_metadata.csv \
  --output results/activity_mapping.tsv
```

## Project Structure

```
configs/       – YAML configs (defaults, thresholds, tuning grids)
schemas/       – JSON schemas and contract docs for I/O validation
src/lpmo_pipeline/
  io/          – mmCIF ingest, normalization, protonation, CCD lookup
  mapping/     – Cross-model atom mapping and renaming
  qc/          – PoseBusters, Privateer, Cu-geometry, QC report
                and pre-QC active-site proximity gate
  placer/      – PLACER refinement wrapper and ensemble ranking
  analysis/    – MDAnalysis metrics, ProLIF IFP, HDBSCAN clustering,
                 cluster signatures, crystal anchoring, CBM analysis,
                 activity mapping, predictive models
  tuning/      – Parameter grid, sweep runner, orchestrator, summary
  report/      – Summary JSON, metrics CSV, HTML report builder
  utils/       – Data models, logging, manifest, hashing, paths, exceptions
tests/         – Unit and contract tests with fixtures
scripts/       – Helper scripts (including R-based EC activity mapping)
```

## Key Invariants

- **CCD monosaccharides only** for Privateer input
- **Active-site proximity pre-gate** before PoseBusters/Privateer
- **PLACER before scoring** — no scoring on raw predictions
- **atom_mapping_coverage = 100%** — all atoms must be mapped
- **Chain convention**: protein=A, glycans=B-D, metal=E
- **Cu-His distance**: 1.9–2.6 Å (hard gate)
- **Cluster rows are primary** for main descriptive and predictive analyses
- **HDBSCAN params fixed per analysis run** (no in-run p-hacking)

## Implementation Order

See [IMPLEMENTATION_PLAYBOOK.md](IMPLEMENTATION_PLAYBOOK.md) for the
20-step prioritized implementation plan with stop-points.

## Open Questions

See [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) for unresolved decisions.
