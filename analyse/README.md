# LPMO Pipeline

Analyse-pipeline for LPMO enzyme–oligosaccharide structure predictions
from AlphaFold 3 only (RosettaFold 3 / RF3 og Boltz-2 er ekskludert fra analyse).

## Overview

This pipeline takes predicted structures (mmCIF), normalizes them,
runs quality control (PoseBusters, Privateer, custom Cu-geometry),
generates interaction fingerprints (ProLIF) and residue-level contact tables,
and clusters binding modes (HDBSCAN).
Performs residue importance analysis and produces summary reports.

Primary governing document: [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) (v1.0, primary source from 2026-04-21).

Current implementation status: steps 1-7b are verified. Stage 8 pre-QC
active-site proximity is now wired into hard QC, and targeted QC tests pass
in the existing `analyse_env` environment. Stage 10 Privateer wrapper parsing
and integration scaffolding are implemented and unit-tested. Stage 12 now has
schema-backed QC report tests passing on synthetic three-pose batches.
Real-data QC verification for stages 8-12 is still pending.

**⚠️ Known open conflict: see [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) item 11 (geometric planarity thresholds not yet operationalized).**

See [MASTERPLAN.md](MASTERPLAN.md) for full integrated specification.

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
- **atom_mapping_coverage = 100%** — all atoms must be mapped
- **Chain convention**: protein=A, glycans=B-D, metal=E
- **Cu-His distance**: 1.9–2.6 Å hard gate evaluated on selected coordinating nitrogens only (`His1:N`, `His1:ND1`, and one non-His1 histidine N; `His1:NE2` excluded)
- **Cluster rows are primary** for main descriptive and predictive analyses
- **HDBSCAN params fixed per analysis run** (no in-run p-hacking)

## Implementation Order

See [IMPLEMENTATION_PLAYBOOK.md](IMPLEMENTATION_PLAYBOOK.md) for the
prioritized implementation plan with stop-points.

## Open Questions

See [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) for unresolved decisions.

---

## Documentation Map

### Active documentation (governing)

| File | Role |
|------|------|
| [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) | **Primary governing plan** (v1.0) — pipeline stages, design decisions, output tables, geometry procedures |
| [MASTERPLAN.md](MASTERPLAN.md) | Integrated secondary plan — stage summaries, failure policy, RQ→output mapping, EC label mapping |
| [IMPLEMENTATION_PLAYBOOK.md](IMPLEMENTATION_PLAYBOOK.md) | Step-by-step implementation guide with status tracking and stop-points |
| [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) | Open decisions, resolved conflicts, and avklaringer |
| [copilot.md](copilot.md) | AI assistant (Copilot) guide — hard rules, architectural policy, AI decision boundaries |
| [ATTRIBUTION.md](ATTRIBUTION.md) | Third-party code attribution |
| [DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md](DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md) | Prioritized open problems and manual verification steps |

### Legacy / archived documentation

These files are preserved as historical context. They are **not governing**. Each file has an ARCHIVED header explaining its status and what has been superseded.

| File | Notes |
|------|-------|
| [PSEUDOKODE_CLUSTER_FIRST.md](PSEUDOKODE_CLUSTER_FIRST.md) | Legacy pseudocode — predates PLACER removal and AF3-only simplification |
| [PSEUDOKODE_MODUL_FOR_MODUL.md](PSEUDOKODE_MODUL_FOR_MODUL.md) | Legacy module pseudocode — contains obsolete PLACER package section |
| [PSEUDOKODE_STATISTIKK.md](PSEUDOKODE_STATISTIKK.md) | Legacy statistics pseudocode — references obsolete single `pose_table.tsv` and `cross_model_support` |
| [plan_analyse.txt](plan_analyse.txt) | Original analysis plan (not updated) — research questions RQ1-RQ4 still informative as background |
| [plan_implementation_spec.txt](plan_implementation_spec.txt) | Old implementation spec — PLACER steps 3-4 are obsolete |
