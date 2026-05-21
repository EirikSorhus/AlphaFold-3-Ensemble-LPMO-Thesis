# LPMO Pipeline

Analyse-pipeline for LPMO enzyme–oligosaccharide structure predictions
from AlphaFold 3 only (RosettaFold 3 / RF3 og Boltz-2 er ekskludert fra analyse).

## Overview

This pipeline takes predicted structures (mmCIF), normalizes them,
runs quality control (PoseBusters, Privateer, custom Cu-geometry),
generates interaction fingerprints (ProLIF) and residue-level contact tables,
and clusters binding modes through the shared production clustering interface
(primary pilot-selected method: agglomerative Jaccard with `distance_threshold=0.55`
and `min_cluster_size=3`; predefined agglomerative/HDBSCAN sensitivity paths retained).
Performs residue importance analysis and produces summary reports.

Primary governing document: [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) (v1.1, primary source; clustering method decision recorded in [DECISIONS.md](DECISIONS.md) and [clustering_pilot_plan.yaml](clustering_pilot_plan.yaml)).

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
space and writes the expanded nine-interaction count fields. Stage 13b
`pose_residue_contact_table.tsv` is now wired into the production `run` path
from the same ligand-resolved ProLIF features, with focused pytest coverage for
the extractor and orchestrator integration, and is now also verified on real
AF3 data via `tests/run_tests_scripts/test_residue_contact_real_cifs.sh`
(job 619027). Stage 14b convergence metrics are now implemented and wired into
the current production path with focused pytest coverage, and the standalone
convergence slice is now verified on real AF3 data via
`tests/run_tests_scripts/test_convergence_real_cifs.sh` (job 619047); integrated
real-data verification of that production path still remains.
Crystal anchoring is now also verified as a standalone real-data slice via
`tests/run_tests_scripts/test_crystal_anchoring_real_cifs.sh`; that slice
resolves crystal references from `input_data/pdb_structure_data.csv`, verifies
both `5ACI` and `7PXW` for `A0A0S2GKZ1`, prepares chemistry-preserving crystal
subsets, normalizes ligand-bound crystal references before protonation, and
writes compact `ifp_result.json` summaries plus detailed `pose_ifp_table.tsv`
and `ifp_matrix.csv` artifacts for the representative pose and each prepared
crystal reference. Crystal subsets use the same canonical chain scheme as AF3
poses (`A` protein, `B/C/D` glycans, `E` Cu). Cu stays in `complex_H.pdb` for
geometry/RMSD, but is excluded from `ligand_for_prolif.mol2` so ProLIF receives
the glycan ligand only. During ProLIF loading, ligand atom residue identities are
resolved against the sibling `ligand_only_for_prolif.pdb` when available, so
multiple glycan chains with repeated residue numbers remain distinct
(`BGC1.B`, `BGC1.C`, etc.) instead of collapsing to a single chain. In
production, crystal anchoring compares every retained
cluster medoid; conditions with no retained clusters may use the top-level AF3
model CIF only after hard QC passes, and that fallback is not counted in normal
pose/clustering denominators. Current pocket RMSD uses a local gemmi/numpy Kabsch
alignment on shared pocket C-alpha atoms, with the pocket defined as protein
residues within 5 A of ligand or Cu. Apo crystal references can reuse a
sequence-projected pocket from the representative/medoid pose, with
residue-name normalization for variants such as `HIC -> HIS`. Crystal-side IFP
Tanimoto is reported only when the prepared crystal IFP passes the same non-vdW
contact-eligibility rule as pose clustering; pocket RMSD and crystal geometry
remain reportable when the crystal IFP is VdW-only, zero-contact, or otherwise
low-specific-contact. After the crystal-prep hardening, those non-comparable IFP
classes should be interpreted as data/contact signal outcomes rather than known
Cu/ligand-export contamination. Ligand-bound crystal references also emit C1/C4
geometry for comparison against pose geometries. The same crystal-anchoring slice is now
also wired into `run_analysis_core`, writes an explicit
`crystal_anchoring_stage_completed` gate in `run_manifest.json`, and has passed
focused pytest plus production smoke validation. Remaining gap: a production
real-data run that reaches non-empty medoid-vs-crystal comparisons in the
integrated path still needs to be exercised explicitly.

Stage 16 cluster annotation now also writes `cluster_table.tsv` as a flat,
backward-compatible TSV export of the retained non-noise `cluster_summaries`
already emitted in `cluster_signatures.json`. The table keeps the original
identifier/type/geometry columns and appends condition metadata, cluster size,
C1/C4 computable/plausible/highly-plausible fractions, pose-confidence
aggregates, and convergence aggregates. This export is wired through
`run_analysis_core` and the CLI artifact listing, and focused pytest covering
`tests/test_cluster_signatures.py`, `tests/test_analysis_orchestrator.py`, and
`tests/test_cli_run.py` passes in `analyse_env`.

`analysis/activity_mapping.py` now has a separate postprocess builder for
`predictive_cluster_table.tsv` from enriched `cluster_table.tsv`, protein
metadata, and EC/activity labels. It is intentionally not part of the core
`lpmo-pipeline run` path.

`analysis/predictive_models.py` currently provides a tested provisional scaffold:
leakage-safe grouped folds on `protein_id` with a 5-fold default when enough
protein groups exist, plus two narrow baseline runners. These are not final
model choices. The predictive-analysis plans are explicitly marked for revision
before final model implementation. The current C1/C4 and substrate plans have
been revised into smaller, less-detailed activity-prediction plans that lock
compact predictor sets and keep the implementation scope narrow.

`analysis/cbm_comparison.py` now has a condition-level CBM paired-analysis
builder for `cbm_construct_condition_summary.tsv` and
`cbm_paired_comparison_table.tsv`, using `condition_table.tsv`,
`cluster_table.tsv`, and protein metadata. Pose-level CBM dual-IFP remains a
later side analysis.

Stage 16b residue importance is implemented and wired into the production
analysis-core path against the Stage 16 signature contract
(`cluster_residue_signature.tsv`, `cluster_ifp_signature.tsv`,
`cluster_signatures.json`). Conditions with observed contact residues but no
retained non-noise clusters emit explicit zero-valued residue rows instead of
header-only Stage 16b residue tables. Loop fraction is not part of the current
contract; patch summaries use residue class, hydrogen-bond contacts, and
explicit region flags where present.

**⚠️ Known open conflict: see [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) item 11 (geometric planarity thresholds not yet operationalized).**

See [MASTERPLAN.md](MASTERPLAN.md) for full integrated specification.

## Data Availability

Precomputed AF3 structures are fully available under the structure_pipeline at these locations:

- **Domain-only constructs:** `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core`
- **Full-length constructs:** `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length`

These directories contain the organized AF3 prediction artifacts (mmCIF files and confidence metrics) that feed directly into the analysis pipeline. See [Production Config](#production-config) for how to reference these in your configuration.

## Execution Model

- Structure prediction (AF3 only) is run as a separate stage from analysis. RF3/RosettaFold 3 and Boltz-2 are excluded.
- Prediction outputs are generated once per chosen parameter setup and reused by analysis.
- Main analyses are run first; parameter tuning is optional and done after baseline analysis.
- Main analysis is cluster-primary (not enzyme-aggregated) for descriptive and predictive modeling.
- Operative sub-analyses are substrate x DP combinations for DP4, DP6, DP8.
- Primary IFP clustering method is locked to agglomerative Jaccard with `linkage=average`, `distance_threshold=0.55`, and `min_cluster_size=3`, based on the completed clustering pilot.
- A pre-QC active-site proximity gate runs before PoseBusters and Privateer.
- PoseBusters is run in built-in `dock` mode on auto-split ligand/protein inputs. The pipeline does not override PoseBusters' intermolecular-distance defaults, so `protein-ligand_maximum_distance` uses `max_distance=5.0 Å`, while `minimum_distance_to_protein` is the separate no-clashes check.
- Production `--n-jobs` is now used for independent per-pose prepare work, hard-QC/Privateer dispatch, and ProLIF batch work. The full clustering-pilot wrapper defaults it from `SLURM_CPUS_PER_TASK` so requested cores are actually used.
- Routine normalization keeps `normalize_report.json` as the normal audit surface. `atom_map.tsv` and `rename_log.json` are now debug/failure artifacts only, because writing them for every successful pose produced heavy I/O without adding useful signal.
- Detailed geometric planarity requirements are currently missing and must be specified before final reporting.
- R is preferred for descriptive/predictive statistics where practical.

## Method And Downstream Plans

One YAML plan records the completed clustering-method pilot decision, and three downstream
YAML plans define the intended statistical analyses after the core AF3
structure-analysis outputs have been generated:

- [clustering_pilot_plan.yaml](clustering_pilot_plan.yaml): completed pilot plan and decision record for the selected global IFP clustering method.
- [c1_c4_predictive_analysis_plan_simplified.yaml](c1_c4_predictive_analysis_plan_simplified.yaml): revised compact C1/C4 activity-prediction plan with small-effective-n constraints and a locked limited predictor set.
- [substrate_activity_prediction_plan.yaml](substrate_activity_prediction_plan.yaml): revised compact substrate-activity prediction plan with small-effective-n constraints and a locked limited predictor set.
- [cbm_full_length_vs_domain_only_analysis_plan.yaml](cbm_full_length_vs_domain_only_analysis_plan.yaml): paired full-length vs domain-only CBM analysis for proteins with both construct types.

These files are implementation plans, not runtime configuration for
`lpmo-pipeline run`. The two activity-prediction plans were simplified on
2026-05-21 and no longer try to be exhaustive implementation blueprints. They
refine the Stage 14 predictive-analysis and Stage 15 CBM paired-analysis
specifications in [MASTERPLAN.md](MASTERPLAN.md).

## Configuration Policy

- All changeable pipeline values must live in config files under `configs/`, not as duplicated literals in scripts.
- This includes thresholds, limits, numeric defaults, flags, model/runtime choices, selections, filenames, folder names, absolute paths, external tool paths, schema/config asset locations, and similar operational settings.
- Shared resolution of runtime paths and bundled config/schema assets happens through `lpmo_pipeline.config`.
- External tool/container paths and bundled asset paths are centralized in [configs/runtime_paths.yaml](configs/runtime_paths.yaml).
- Scripts may load a config value once locally inside a module or function for readability/performance, but the source of truth must still be the config files.
- Stable algorithmic details may remain inline only when they are not realistic user/runtime configuration and moving them would reduce clarity.

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Production mode (current analysis-core slice: QC + geometry + ProLIF + clustering + cluster annotation + residue importance + crystal anchoring + reports)
lpmo-pipeline run \
  --mode production \
  --config configs/production.analysis_core.example.yaml \
  --output results/del_a \
  --del del_a \
  --n-jobs 4

# Optional post-analysis tuning (comparative reruns)
lpmo-pipeline tune \
  --model AF3 \
  --config configs/tuning_af3.yaml \
  --output results/tuning_af3

# Optional AA9/AA10 family residue enrichment postprocess
lpmo-pipeline family-enrichment \
  --protein-condition-residue-scores results/del_a/protein_condition_residue_scores.tsv \
  --protein-residue-regio-delta results/del_a/protein_residue_regio_delta.tsv \
  --protein-metadata input_data/metadata_final_ec_fixed.tsv \
  --core-fasta input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta \
  --output results/del_a

# Sbatch-backed validation for the optional family enrichment layer
sbatch tests/run_tests_scripts/test_family_enrichment_validation.sh

# Real-data smoke test for the current production entry path
sbatch tests/run_tests_scripts/test_analysis_core_real_cifs.sh

# Standalone real-data ProLIF validation
sbatch tests/run_tests_scripts/test_prolif_real_cifs.sh

# Standalone real-data crystal anchoring validation
sbatch tests/run_tests_scripts/test_crystal_anchoring_real_cifs.sh

# EC metadata -> activity labels (R helper script)
Rscript scripts/ec_activity_mapping.R \
  --input metadata/enzyme_metadata.csv \
  --output results/activity_mapping.tsv
```

## Production Config

The current production entry path is the analysis-core slice. It runs discovery,
normalization, hard QC, downstream geometry, ProLIF/IFP, residue-contact
extraction, convergence metrics, condition-wise clustering, Stage 7 cluster
annotation, residue importance, crystal anchoring, and report generation. It
now writes the implemented pose/QC TSV surfaces (`pose_manifest.tsv`,
`pose_confidence.tsv`, `structure_index.tsv`, `qc_attrition_table.tsv`) plus
raw clustering/convergence/IFP tables, `cluster_table.tsv`,
`cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`,
`cluster_signatures.json`, `condition_table.tsv`, `protein_summary_table.tsv`,
`crystal_anchor_table.tsv`,
`crystal_geometry_table.tsv`, and `crystal_ifp_diagnostic_summary.tsv`. The crystal-anchoring stage gate in this
production path is smoke-validated, but a real-data production run that reaches
actual medoid-vs-crystal comparisons still remains.

The optional family enrichment layer is not part of this production contract.
It is a separate postprocess over Stage 16b outputs and currently targets only
AA9/AA10 domain-only rows through `lpmo-pipeline family-enrichment`.

Verified 2026-05-03 on 3 staged real AF3 CIFs (`analysis_core_real_cifs_613251`):
1 `pass`, 1 `soft_flag`, 1 `hard_fail`; 2 poses were analyzed downstream, and
no production `geometry_debug.pdb` files were written.

Use [configs/production.analysis_core.example.yaml](configs/production.analysis_core.example.yaml)
as the starting point.

**Data input configuration:**

The pipeline discovers AF3-predicted structures from the structure_pipeline work directories.
Configure which data to use via the `work_roots` keys:

- `work_roots.domain_only`: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core` — domain-only LPMO constructs
- `work_roots.full_length`: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length` — full-length LPMO constructs

Select which to run via `construct_type: "domain_only"` or `construct_type: "full_length"` in your config.

**Other configuration keys:**

- `production.run_id`: optional identifier written into the production outputs
- `production.af3_only`: limit discovery to AF3 runs
- `production.latest_only`: prefer the `latest` symlink under each target/model
- `production.max_cases`: optional cap for smaller validation runs
- `production.include_targets`: optional target whitelist such as `NAG4`, `STA4`, `STA6`
- `production.n_jobs`: optional YAML fallback for parallel workers; the CLI `--n-jobs` value takes precedence.

**Notes:**

- The output directory is still controlled by the CLI `--output` argument, not the YAML file.
- The DEL branch is still controlled by the CLI `--del` argument.
- `--n-jobs` controls parallel workers for independent pose preparation, hard QC/Privateer dispatch, and ProLIF batch computation. Scaling is useful but not expected to be perfectly linear because process startup, container startup, filesystem I/O, and unequal per-pose runtime still contribute fixed overhead.
- Production output must not contain `geometry_debug.pdb`; that file remains test-only.
- Standalone crystal-anchoring harnesses may auto-select a best-ranked AF3 pose for test coverage. The production path compares all retained cluster medoids; if a condition has no retained clusters, it may compare the top-level AF3 model CIF only after that fallback passes hard QC.
- The optional `family-enrichment` command is intentionally decoupled from `run`. It consumes existing Stage 16b TSVs plus metadata and the deduplicated catalytic-core FASTA, filters to AA9/AA10 domain-only rows, and writes outputs under `08_family_residue_enrichment/`. If no precomputed family alignment is provided and `mafft` is unavailable, the command records a clean skip reason instead of failing the production outputs.
- Shared external tool paths and default asset references are controlled separately through [configs/runtime_paths.yaml](configs/runtime_paths.yaml).

## Config Files

- [configs/runtime_paths.yaml](configs/runtime_paths.yaml): shared runtime paths for external tools plus central references to bundled config/schema assets.
- [configs/production.analysis_core.example.yaml](configs/production.analysis_core.example.yaml): current production run entry config for discovery scope and run selection.
- [configs/defaults.yaml](configs/defaults.yaml): global defaults such as chain schema, confidence policy, atom-mapping defaults, and target-prefix substrate mapping.
- [configs/thresholds.yaml](configs/thresholds.yaml): hard/soft QC thresholds, downstream geometry thresholds, convergence settings, and clustering inclusion settings.
- [clustering_pilot_plan.yaml](clustering_pilot_plan.yaml): completed method-selection plan and decision record for contact eligibility, IFP feature audit, agglomerative/HDBSCAN comparison, bootstrap stability, and the final clustering choice.
- [configs/geometry_rules.yaml](configs/geometry_rules.yaml): Cu/his-brace identification, virtual oxyl/H placement rules, and geometry output contracts.
- [configs/prolif_features.yaml](configs/prolif_features.yaml): ProLIF interaction set, cutoffs, selections, feature naming, and residue-contact settings.
- [configs/residue_rules.yaml](configs/residue_rules.yaml): residue/atom normalization rules, CCD mappings, and region tagging rules.
- [configs/cv_hierarchy.yaml](configs/cv_hierarchy.yaml): grouped cross-validation hierarchy and stratification settings for later predictive analyses.
- [configs/tuning_af3.yaml](configs/tuning_af3.yaml), [configs/tuning_boltz2.yaml](configs/tuning_boltz2.yaml), [configs/tuning_rf3.yaml](configs/tuning_rf3.yaml): model-specific tuning grids.
- [human_readability_code_walkthrough/FAMILY_ENRICHMENT_POSTPROCESS.md](human_readability_code_walkthrough/FAMILY_ENRICHMENT_POSTPROCESS.md): detailed walkthrough of the optional AA9/AA10 family enrichment postprocess.

## Project Structure

```
configs/       – YAML configs (defaults, thresholds, tuning grids)
schemas/       – JSON schemas and contract docs for I/O validation
src/lpmo_pipeline/
  io/          – mmCIF ingest, normalization, protonation, CCD lookup
  mapping/     – Cross-model atom mapping and renaming
  qc/          – PoseBusters, Privateer, Cu-geometry, QC report
                and pre-QC active-site proximity gate
  analysis/    – MDAnalysis metrics, ProLIF IFP, clustering pilot/HDBSCAN/agglomerative paths,
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
- **Clustering primary method selected from pilot**: agglomerative Jaccard on contact-eligible IFP rows, `linkage=average`, `distance_threshold=0.55`, `min_cluster_size=3`. The pilot selection evidence used `main_contact_eligible_ifp_matrix.csv`. Sensitivity settings are agglomerative `distance_threshold=0.45`, `min_cluster_size=3`, and HDBSCAN Jaccard `min_cluster_size=3`, `min_samples=null`.

## Implementation Order

See [IMPLEMENTATION_PLAYBOOK.md](IMPLEMENTATION_PLAYBOOK.md) for the
prioritized implementation plan with stop-points.

## Implemented Decision Logic

See [DECISIONS.md](DECISIONS.md) for a code-derived description of the
runtime decision rules that are currently active in the implementation,
including contact eligibility, clustering feature selection, cluster typing,
the flat Stage 7 `cluster_table.tsv` export, residue signature construction,
and Stage 16b residue aggregation.

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
| [DECISIONS.md](DECISIONS.md) | Implementation-grounded decision log — documents the rules that are actually active in code |
| [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md) | Open decisions, resolved conflicts, and avklaringer |
| [c1_c4_predictive_analysis_plan_simplified.yaml](c1_c4_predictive_analysis_plan_simplified.yaml) | Revised compact plan for exploratory C1/C4 activity prediction under small effective n |
| [substrate_activity_prediction_plan.yaml](substrate_activity_prediction_plan.yaml) | Revised compact plan for exploratory substrate activity prediction under small effective n |
| [cbm_full_length_vs_domain_only_analysis_plan.yaml](cbm_full_length_vs_domain_only_analysis_plan.yaml) | Detailed implementation plan for paired full-length vs domain-only CBM analysis |
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
