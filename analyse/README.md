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
schema-backed QC report tests passing on synthetic three-pose batches. Stages
8-12 hard QC are now also verified on real AF3 data via
`tests/run_tests_scripts/test_hard_qc_real_cifs.sh` (job 612234). Stage 14
downstream geometry is implemented, wired into the current production
analysis-core slice, and verified on real AF3 data via
`tests/run_tests_scripts/test_analysis_core_real_cifs.sh` (job 613251). Stage 13
ProLIF/IFP is now verified as a standalone real-data slice via
`tests/run_tests_scripts/test_prolif_real_cifs.sh` (job 617120); that slice
preserves ligand monosaccharides as separate residues in the ProLIF feature
space and writes the expanded nine-interaction count fields. The current
production `run` path still stops before ProLIF/IFP generation,
clustering, and crystal anchoring.

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

# Production mode (current analysis-core slice: QC + geometry + reports)
lpmo-pipeline run \
  --mode production \
  --config configs/production.analysis_core.example.yaml \
  --output results/del_a \
  --del del_a

# Optional post-analysis tuning (comparative reruns)
lpmo-pipeline tune \
  --model AF3 \
  --config configs/tuning_af3.yaml \
  --output results/tuning_af3

# Real-data smoke test for the current production entry path
sbatch tests/run_tests_scripts/test_analysis_core_real_cifs.sh

# Standalone real-data ProLIF validation
sbatch tests/run_tests_scripts/test_prolif_real_cifs.sh

# EC metadata -> activity labels (R helper script)
Rscript scripts/ec_activity_mapping.R \
  --input metadata/enzyme_metadata.csv \
  --output results/activity_mapping.tsv
```

## Production Config

The current production entry path is the analysis-core slice. It runs discovery,
normalization, hard QC, downstream geometry, and report generation. It stops
before ProLIF/IFP generation, clustering, and crystal anchoring.

Verified 2026-05-03 on 3 staged real AF3 CIFs (`analysis_core_real_cifs_613251`):
1 `pass`, 1 `soft_flag`, 1 `hard_fail`; 2 poses were analyzed downstream, and
no production `geometry_debug.pdb` files were written.

Use [configs/production.analysis_core.example.yaml](configs/production.analysis_core.example.yaml)
as the starting point.

Current keys:

- `production.work_root`: absolute path to the `structure_pipeline/work` directory
- `production.run_id`: optional identifier written into the production outputs
- `production.af3_only`: limit discovery to AF3 runs
- `production.latest_only`: prefer the `latest` symlink under each target/model
- `production.max_cases`: optional cap for smaller validation runs
- `production.include_targets`: optional target whitelist such as `NAG4`, `STA4`, `STA6`

Notes:

- The output directory is still controlled by the CLI `--output` argument, not the YAML file.
- The DEL branch is still controlled by the CLI `--del` argument.
- Production output must not contain `geometry_debug.pdb`; that file remains test-only.

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
