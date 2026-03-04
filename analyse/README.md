# LPMO Pipeline

Analyse-pipeline for LPMO enzyme–oligosaccharide structure predictions
from AlphaFold 3, RoseTTAFold 3, and Boltz-2.

## Overview

This pipeline takes predicted structures (mmCIF), normalizes them,
runs quality control (PoseBusters, Privateer, custom Cu-geometry),
refines with PLACER, clusters interaction fingerprints (HDBSCAN),
and produces a comprehensive analysis report.

See [MASTERPLAN.md](MASTERPLAN.md) for full specification.

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Tuning mode (find optimal parameters per predictor)
lpmo-pipeline tune --config configs/tuning_af3.yaml

# Production mode (run full pipeline)
lpmo-pipeline run --config configs/defaults.yaml --systems-csv systems.csv
```

## Project Structure

```
configs/       – YAML configs (defaults, thresholds, tuning grids)
schemas/       – JSON schemas and contract docs for I/O validation
src/lpmo_pipeline/
  io/          – mmCIF ingest, normalization, protonation, CCD lookup
  mapping/     – Cross-model atom mapping and renaming
  qc/          – PoseBusters, Privateer, Cu-geometry, QC report
  placer/      – PLACER refinement wrapper and ensemble ranking
  analysis/    – MDAnalysis metrics, ProLIF IFP, HDBSCAN clustering,
                 cluster signatures, crystal anchoring, CBM analysis,
                 activity mapping, predictive models
  tuning/      – Parameter grid, sweep runner, orchestrator, summary
  report/      – Summary JSON, metrics CSV, HTML report builder
  utils/       – Data models, logging, manifest, hashing, paths, exceptions
tests/         – Unit and contract tests with fixtures
```

## Key Invariants

- **CCD monosaccharides only** for Privateer input
- **PLACER before scoring** — no scoring on raw predictions
- **atom_mapping_coverage = 100%** — all atoms must be mapped
- **Chain convention**: protein=A, glycans=B-D, metal=E
- **Cu-His distance**: 1.9–2.6 Å (hard gate)
- **HDBSCAN params locked** after tuning — no re-tuning in production

## Implementation Order

See [IMPLEMENTATION_PLAYBOOK.md](IMPLEMENTATION_PLAYBOOK.md) for the
20-step prioritized implementation plan with stop-points.

## Open Questions

See [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) for unresolved decisions.
