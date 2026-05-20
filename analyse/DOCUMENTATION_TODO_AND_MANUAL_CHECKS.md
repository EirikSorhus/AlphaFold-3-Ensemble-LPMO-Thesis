# Documentation TODO and Manual Checks

Version: 1.0  
Created: 2026-04-21  
Purpose: Prioritized open problems with proposed solutions, and a manual-check section with step-by-step instructions.

---

## Data Availability Status

**RESOLVED 2026-05-07**: All precomputed AF3 structures are fully available and ready for analysis.

**Data locations:**
- **Domain-only constructs:** `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core`
- **Full-length constructs:** `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length`

These directories contain the complete set of AF3 prediction artifacts (mmCIF files and confidence metrics) organized by target, substrate, and construct type. Configure which dataset to analyze via the `work_roots.domain_only` and `work_roots.full_length` keys in your production config.

---

## Detailed Downstream Analysis Plans

The following YAML files have been added as detailed implementation plans for
downstream analyses that run after the core AF3 structure-analysis pipeline
produces condition-, cluster-, residue-, and geometry-level outputs:

- `c1_c4_predictive_analysis_plan_simplified.yaml` — detailed exploratory
  predictive-analysis plan for protein-level C1 and C4 regioactivity. It
  defines the modeling row unit, required input tables, target labels,
  repeated ligand-context handling, feature engineering, nested protein-grouped
  cross-validation, and output interpretation limits.
- `substrate_activity_prediction_plan.yaml` — detailed exploratory
  predictive-analysis plan for substrate activity across chitin, cellulose, and
  starch. It defines protein-substrate activity labels, condition-level
  modeling rows, cluster aggregation, feature blocks, grouped CV, models, and
  reporting outputs.
- `cbm_full_length_vs_domain_only_analysis_plan.yaml` — detailed paired-analysis
  plan for comparing full-length and domain-only constructs in CBM-containing
  proteins. It defines matched construct-pair inclusion rules, CBM/linker
  region requirements, condition summaries, primary paired endpoints,
  statistics by sample size, stratified reporting, and expected CBM outputs.

These files should be treated as active downstream implementation plans that
refine `MASTERPLAN.md` Stage 14 and Stage 15. They are not runtime
configuration files for the current `lpmo-pipeline run` command.

---

## Part 1 — Prioritized Open Problems

Problems are ordered from highest to lowest priority.

---

### P0 — Missing-glycan normalization false-success path (RESOLVED 2026-05-08)

**Priority:** RESOLVED

**Description:**
Real-case debugging of `contact_eligibility_real_cifs_646218` showed that the
problematic STA8 inputs from AF3 run `336981` already arrive without any parsed
glycan chain/residues in the source CIF (`gemmi.read_structure(...)` sees only
chains `A` and `B`, not the expected branched chain `C`). By contrast, newer
STA8 outputs from run `408007` contain the expected branched glycan chain and
GLC residues.

The bug in normalization was that
`"No glycan residues found in normalized chains B..D"` was treated as a soft
flag and normalization returned success. This allowed malformed source CIFs to
slip into later QC and fail there as downstream `"No ligand atoms found"`
cases instead of stopping at normalization.

**Implemented fix:**
1. `normalize_mmcif.py` now hard-fails when no glycan residues remain after chain remap.
2. `normalize_report.json` now records explicit `failure_reason = "missing_glycan_chain"` for this path.
3. Regression coverage was added for source CIFs that contain protein + Cu but no glycan chain, while preserving the separate `no _struct_conn` soft-flag test case.
4. `io/discovery.py` now honors the intended newest AF3 run when `latest_only=True`, even when the `latest` symlink points across `work` / `work_core`; if that symlink is absent or unusable, discovery falls back to the highest numeric run ID instead of scanning all historical runs.

**Verification:**
- `pytest tests/test_normalize.py -q` → 15 passed (2026-05-08)
- `pytest tests/test_discovery.py -q` → 80 passed (2026-05-08)

**Audit update 2026-05-08:**
- A work-root audit across all 9 AF3 targets under `structure_pipeline/work_core` showed the same historical pattern, not just STA8.
- The earliest `336973`–`336981` runs are systematically ligandless at source level for their respective targets (`CEL4`, `CEL6`, `CEL8`, `NAG4`, `NAG6`, `NAG8`, `STA4`, `STA6`, `STA8`): sampled CIFs contain only chains `A` and `B`, entity types `polymer` + `non-polymer`, and no branched glycan entity/residues.
- The current latest runs (`407999`–`408007`, plus `406216` for the one-off CEL4 rerun) contain the expected branched glycan metadata and glycan residues.
- Audit artifact: `tests/tests_results/latest_run_audit_20260508/work_core_af3_run_audit.json`.

---

### P0a — Clustering-pilot runtime and summary/artifact volume (RESOLVED 2026-05-19)

**Priority:** RESOLVED for the current production path and clustering-pilot Slurm orchestration.

**Description:**
The full clustering pilot timed out before reaching the method-comparison stage because per-pose preparation and QC/IFP work were mostly serialized, despite Slurm cores being requested. The wrapper also wrote several large, near-duplicate JSON summaries/checkpoints, and normalization wrote atom-mapping debug files for every successful pose.

**Implemented fix:**
1. `lpmo-pipeline run`, the real-case pilot runner, and `run_clustering_pilot_full.sh` now pass `n_jobs` into the production pipeline. The wrapper defaults to `SLURM_CPUS_PER_TASK`.
2. Independent per-pose preparation now uses process workers; hard QC/Privateer dispatch and ProLIF batch work use bounded worker pools.
3. Routine successful normalization writes only `normalize_report.json`; `atom_map.tsv` and `rename_log.json` are emitted only for incomplete atom mapping or CCD-validation failure.
4. Pilot prepare summaries are compacted so the main summary stores counts/sample records instead of duplicating full discovery/staging payloads. Detailed checkpoint files remain the place for deep debugging.
5. `gemmi_compat.py` now groups CIF remap work by loop, so each affected mmCIF loop is rewritten once even when several tags in that loop need chain remapping.
6. `submit_clustering_pilot_staged.sh` is now the primary full-pilot launcher. It builds selection manifests locally, splits them into protein-level shard manifests, submits domain-only and full-length Slurm arrays, and collects one compact aggregate shard summary after the arrays finish.

**Why:**
This reduces wall time by using the cores already requested, lowers metadata/file-system overhead, and makes pilot output easier to inspect after long Slurm runs.

**Verification:**
- Focused pytest slice covering normalization, orchestration, CLI, hard QC, Privateer, ProLIF, and pilot runner: 57 passed (2026-05-19).
- Small real-pose timing probe on 4 poses, not the full pilot: `n_jobs=1` took 56.23 s; `n_jobs=2` took 28.94 s. This is about 1.9x for this tiny prepare-heavy probe, but full-run scaling is not expected to be perfectly linear.
- After the CIF remap optimization, a 2-pose Slurm check on the first two domain-only poses reduced normalization from 37.48 s total to 15.30 s total. This was a targeted sanity check, not a new full timing campaign.
- Static wrapper checks passed for the Slurm-array launcher (`bash -n` and `git diff --check`). No full pilot was submitted during implementation.

**Remaining action:**
Use the staged Slurm-array launcher for the next full pilot. Separate resource profiles for later non-pilot production stages can still be added if the full analysis outgrows the current shard-level execution model.

---

### P1 — Geometry plausibility thresholds — RESOLVED 2026-04-21

**Priority:** RESOLVED

**Description:**  
All geometry plausibility threshold keys in `configs/thresholds.yaml` section `geometry_plausibility` are now set with literature-backed values and `locked: true` has been set. MASTERPLAN.md and the plan document have been updated to reflect the locked state.

Locked values:
- `oxyl_h_min_a: 1.5`, `oxyl_h_max_a: 4.0` (plausibility window)
- `oxyl_h_optimum_a: 2.1` (DFT-derived HAT optimum)
- `oxyl_h_score_normalization: 1.9`
- `oxyl_h_highly_plausible_min_a: 1.8`, `oxyl_h_highly_plausible_max_a: 2.5`
- `locked: true`

**Remaining action:** Record the lock state and config hash in `run_manifest.json` at the start of the full analysis run.

---

### P2 — atom mapping modules (steps 5–6) exist but are not verified (PIPELINE BLOCKER) — RESOLVED 2026-04-22

**Priority:** RESOLVED

**Description:**  
Verified 2026-04-22 via `test_mapping_contracts.sh` on real AF3 data (A0A0A1ED04_NAG8). Both `mapping/cross_model_atom_mapping.py` and `mapping/rename_atoms.py` have been validated:
- Atom mapping coverage = 100% achieved.
- No unmapped atoms in structure.
- Round-trip rename validation (forward→reverse) passed successfully.
- `atom_map.tsv` and `rename_log.json` generated with correct schemas.

**Update 2026-05-19:** The mapping modules can still generate these files, but routine successful normalization no longer writes `atom_map.tsv` and `rename_log.json` for every pose. They are kept as debug/failure artifacts only, because the normal audit signal is already captured in `normalize_report.json` and per-pose mapping files created excessive output.

---

### P3 — Step 7 protonation/export and step 7b CIF→PDB (RESOLVED 2026-04-23)

**Priority:** RESOLVED

**Description:**
Step 7 and 7b were rewritten and verified 2026-04-23 with PDBFixer/OpenMM as primary backend:

- `io/cif_to_pdb.py`: `_convert_auto()` tries PDBFixer first; falls back to gemmi.
  New report field `backend_fallback_reason` records why fallback was used (empty = PDBFixer succeeded).
- `io/protonate_export.py`: protonation priority is PDBFixer/OpenMM → reduce → obabel → BLOCKER.
  No silent copy-fallback. New report field `complex_h_backend` replaces old `used_reduce`/`used_obabel_for_complex_h`.
  Ligand extraction validates that extracted chains are not protein-only.
  Final connectivity solution:
  - name-based glykan-linking før protonering: C1(i)→O4(i+1) for NAG/BGC/GLC
  - eksplisitt CONECT rewrite i `complex_H.pdb` med komplett glykan-konnektivitet
    (intra-residue + inter-residue link)
- Test script `test_protonation_contracts.sh` now verifies:
  - `complex_h_backend` ≠ "none"
  - `complex_H.pdb` atom count > `for_posebusters.pdb` atom count (hydrogens added)
  - MOL2 has `@<TRIPOS>ATOM` and `@<TRIPOS>BOND` sections
  - MOL2 is free of protein residue labels

Verification runs:
- sbatch jobs: 572928, 572982, 572996, 573050
- `complex_h_backend=pdbfixer`, hydrogen atom delta > 0
- `ligand_for_prolif.mol2` passes TRIPOS + ligand-only checks
- manual PyMOL spot-check confirms expected glykan-bonding appearance

---

### P4 — CCD lookup verification on real AF3 data not yet done (RESOLVED)

**Priority:** HIGH — step 4 is in progress; remaining work is verification. This is the current implementation step.

**Description:**  
IMPLEMENTATION_PLAYBOOK.md step 4: CCD lookup is wired into `normalize_mmcif.py` but not yet verified on real AF3 data. Until verified, normalization is unconfirmed.

**Proposed solution:**  
1. Run `pytest tests/test_io_contracts.py::TestCCDLookup -v` on the target system with a real AF3 output file.
2. Inspect `normalize_report.json` for the test run; verify that all glycan CCD codes are recognized (recognition_rate = 1.0 for valid input).
3. Run a negative test: feed a structure with an invalid glycan CCD code and confirm hard fail is reported correctly.
4. Once verified, mark step 4 as ✅ in IMPLEMENTATION_PLAYBOOK.md and proceed to step 5.

---

### P6 — Crystal reference set needs curation/coverage review

**Priority:** MEDIUM — current reference set exists, but family coverage/curation is still a secondary data task.

**Description:**  
OPEN_QUESTIONS.md item 5 is no longer about missing operational inputs: the current pipeline resolves crystal references from `input_data/pdb_structure_data.csv` plus files under `crystal_structures/`. The remaining work is to review whether that set should be expanded or curated further per family and per holo/apo coverage.

**Proposed solution:**  
1. Review the existing `input_data/pdb_structure_data.csv` against the proteins/families currently in scope.
2. Identify gaps: missing holo coverage, missing apo coverage, or families with only a single unresolved reference type.
3. If the set is expanded, document the chosen PDB codes in `metadata/crystal_reference_list.tsv` with columns: `protein_id`, `pdb_code`, `chain`, `has_ligand`, `family`, `regio_label`, `notes`.
4. Update OPEN_QUESTIONS.md item 5 when the coverage review is finished and mark P6 as RESOLVED.

---

### P7 — Family-specific substrate-recognition/pocket residues not defined per family

**Priority:** MEDIUM — the current heuristic works operationally, but literature-based residue definitions may still improve interpretation.

**Description:**  
OPEN_QUESTIONS.md item 13: the current operational crystal-anchoring pocket uses protein residues within 5 A of ligand or Cu. Apo crystal references can inherit a sequence-projected pocket from the representative/medoid pose, with residue-name normalization for variants such as `HIC -> HIS`. Literature-based residue definitions per family are still preferred if they should replace or refine this heuristic.

**Proposed solution:**  
1. For each LPMO family in the dataset, perform a literature search to identify known substrate-binding residues.
2. Document in `metadata/crystal_reference_list.tsv` or a separate `metadata/alignment_residue_definitions.tsv` with columns: `family`, `residue_number`, `residue_name`, `source`, `notes`.
3. If literature data is unavailable for a family, keep the current operational fallback: residues within 5 A of ligand or Cu in holo references, and sequence-projected pocket residues for apo references. Flag this as `source = proximity_fallback` in the metadata.
4. Document which method was used per system in `crystal_anchor_table.tsv` (the `alignment_method` column).

---

### P8 — Privateer SIF execution path and output contract (RESOLVED 2026-04-30)

**Priority:** RESOLVED

**Description:**  
Privateer is now integrated via SIF instead of host PATH installation. The resolved execution path is `apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif ...` with explicit bind mounts. The resolved output contract is `-mode ccp4i2` plus parsing of `validation_data-privateer`, not JSON stdout.

Implemented and verified points:
- `qc/privateer_runner.py` builds bind-aware SIF invocations and parses `validation_data-privateer`.
- `get_privateer_version()` queries the SIF with `-list`, and `utils/manifest.py` now records that version instead of calling `privateer -V` on the host.
- hard QC uses the runner through `hard_qc_orchestrator.py`, including batched Privateer dispatch for eligible poses; this dispatch now respects the production worker count.
- successful runs now keep `validation_data-privateer` by default; raw stdout/stderr are retained only on failure or when `debug_output=True`.
- targeted pytest coverage passes for the runner, dry-run path, QC gates, QC report, manifest path, and hard-QC orchestration.

**Remaining action:**  
Operational verification now belongs to the real-pose hard-QC check in the implementation playbook, not to CLI/version/output-format discovery.

---

### P8a — Native Gemmi remapping for normalized/privateer inputs (RESOLVED 2026-05-01)

**Priority:** RESOLVED

**Description:**
Native Gemmi importability alone was not sufficient: earlier remap logic only mutated the fallback compatibility objects correctly. On native Gemmi, loop-backed tags were not being rewritten the same way, which could leave regenerated `normalized.cif` and `privateer_input.cif` with stale chain assignments even though parsing itself succeeded.

Implemented and verified points:
- native-Gemmi-compatible mutation paths are now used for remapped scalar and loop-backed mmCIF values.
- the target analysis environment now reports functional native Gemmi (`gemmi_is_functional() = True`, version `0.7.5`).
- a fresh real-case parse of `privateer_input.cif` shows the expected chain layout `['A', 'E', 'B']` after regeneration.
- downstream real-case Privateer runs again recognize the expected glycan residues from regenerated input.

**Remaining action:**
No separate blocker remains; final confirmation now belongs to the planned 3-pose hard-QC manual run.

---

### P8b — PoseBusters metal-containing input PDB (RESOLVED 2026-05-01)

**Priority:** RESOLVED

**Description:**
PoseBusters could return no usable results on the real exported `for_posebusters.pdb` when a Cu atom was present. The resolved contract is now:
- `for_posebusters.pdb` is a PoseBusters-specific export with incompatible metal atom records removed.
- metal-only `TER` / `LINK` / `CONECT` references are removed at the same time.
- non-metal `CONECT` records are retained.
- `complex_H.pdb` is still generated from a full-complex temporary export, so Cu remains present in downstream analysis artifacts.

Implemented and verified points:
- `io/cif_to_pdb.py` now strips incompatible metal atoms from the public PoseBusters PDB and records the stripped element/count in `cif_to_pdb_report.json`.
- `io/protonate_export.py` now protonates a full-complex temporary export so `complex_H.pdb` retains Cu.
- focused regression coverage now checks both Cu removal and retention of non-metal connectivity.
- on the real one-case validation rerun, the regenerated `for_posebusters.pdb` yields PoseBusters CSV output, while the raw full-complex export and `complex_H.pdb` still retain Cu.

**Remaining action:**
No separate blocker remains; the same export contract was reconfirmed in
the integrated 3-pose hard-QC rerun `hard_qc_real_cifs_612234` (job 612234).

---

### P8c — PoseBusters combined-input contract and real hard-QC rerun (RESOLVED 2026-05-03)

**Priority:** RESOLVED

**Description:**
After the Cu-stripping fix, PoseBusters still produced false real-case hard
fails because the pipeline treated a combined AF3 protein+glycan
`for_posebusters.pdb` as standalone `mol_pred` in `mol` mode. The resolved
contract is now:
- combined AF3 PoseBusters PDBs are auto-split into ligand-only `mol_pred`
  and protein `mol_cond`
- PoseBusters is invoked in documented `dock` mode for that case
- the hard-QC sbatch harness now includes `tests/test_posebusters_runner.py`
  in its focused regression phase

Implemented and verified points:
- `qc/posebusters_runner.py` now auto-prepares ligand/protein inputs from
  combined AF3 exports before invoking PoseBusters.
- focused regression `tests/test_posebusters_runner.py` passes and is now run
  from `tests/run_tests_scripts/test_hard_qc_real_cifs.sh`.
- integrated 3-pose hard-QC rerun `hard_qc_real_cifs_612234` (job 612234)
  completed with schema-valid `qc_report.json` and no longer shows the prior
  false PoseBusters hard fails.
- rerun outcome on real data: 1 `pass`, 1 `soft_flag`, 1 `hard_fail`.
- the only remaining PoseBusters hard-fail code on real data is
  `minimum_distance_to_protein`; the same STA4 case also retains a real
  Privateer anomer failure.

**Remaining action:**
Severity policy for PoseBusters fail codes is still open; see
`OPEN_QUESTIONS.md` item 18.

---

### P8d — Downstream dropped-vs-flagged filtering audit not completed

**Priority:** MEDIUM — hard/soft QC semantics are implemented in QC itself, but downstream consumers still need confirmation.

**Description:**
Hard QC now drops poses from downstream analysis while preserving all computed metrics and QC artifacts. Soft failures remain analyzable and should only be flagged. This contract is implemented in the QC orchestration/reporting path, but downstream consumers outside the QC module have not yet been re-audited against the updated behavior.

Update 2026-05-03:
The current analysis-core production slice has now been checked at the actual QC → downstream-analysis boundary via `tests/run_tests_scripts/test_analysis_core_real_cifs.sh` (job 613251, run `analysis_core_real_cifs_613251`). In that run:
- the dropped STA4 pose was excluded from downstream geometry (`analysis_status = skipped_dropped`)
- the soft-flagged STA6 pose remained analyzable and was kept in downstream outputs
- the passed NAG4 pose remained analyzable as expected
- `metrics.csv`, `summary.json`, and `report.html` reflected the same pass/flag/drop split

Update 2026-05-08:
- `hard_qc_orchestrator.py` now short-circuits on hard distance-gate failures: active-site proximity runs first, Cu-His/Cu-substrate geometry runs next, and poses that fail there do not continue to PoseBusters or Privateer.
- The focused real-CIF contact-eligibility harness now supports disabling PoseBusters and Privateer entirely before the ProLIF/clustering slice when requested.

Update 2026-05-16:
The production analysis-core path now writes the implemented pose/QC TSV
surfaces that can be generated before downstream cluster annotation:
- `pose_manifest.tsv`
- `pose_confidence.tsv`
- `structure_index.tsv`
- `qc_attrition_table.tsv`

Dropped-pose numeric metrics still persist primarily in `qc_report.json` and
`analysis_core_summary.json`; the new TSVs provide the tabular pose/QC index
and attrition surface, not a full replacement for the schema-backed QC report.

**Proposed solution:**
1. Decide whether `qc_report.json` + `analysis_core_summary.json` remain the intended persistent surfaces for dropped-pose numeric metrics.
2. If a fully tabular dropped-pose metric surface is needed, extend `pose_manifest.tsv` or add a dedicated QC metrics TSV without re-admitting dropped poses to downstream analysis.
3. Recheck the same contract again after real-data clustering output inspection.

---

### P8e — Downstream geometry analysis branch implemented and real-data checked (RESOLVED 2026-05-03)

**Priority:** RESOLVED

**Description:**
Step 14 is no longer just a stub. `analysis/mdanalysis_metrics.py` now produces the
planned `pose_geometry.tsv` row contract and has been validated in focused tests plus
a dedicated real-CIF harness.

Implemented and verified points:
- `PoseGeometryMetrics.to_row()` and `write_pose_geometry_tsv()` now emit the planned row contract.
- The module reuses `qc/custom_geometry_checks.check_geometry()` only for Cu / His-brace identification and `his_brace_angle_deg` context.
- Downstream `Cu_C1_distance` and `Cu_C4_distance` are recalculated against the repositioned Cu used by the oxyl model; they are not copied from hard-QC `min_cu_c1` / `min_cu_c4`.
- `tests/run_tests_scripts/test_geometry_real_cifs.sh` and `tests/run_tests_scripts/run_geometry_real_cifs.py` now validate the branch on real normalized AF3 CIFs and write `pose_geometry.tsv`.
- `analysis/analysis_orchestrator.py` and `lpmo-pipeline run` now write `pose_geometry.tsv` from the current production analysis-core slice for `passed` and `flagged` poses.
- `geometry_debug.pdb` and `geometry_metrics.json` are explicitly test-only debug artifacts for manual inspection and are not intended as production pipeline outputs.
- The HC1 misplacement seen in the STA4 debug export was fixed by including cross-residue glycosidic neighbors from `structure.connections` / `_struct_conn` during virtual-H placement.

**Remaining action:**
Extend the current production analysis-core slice forward into ProLIF/IFP, clustering, and RMSD-backed fields while keeping debug PDB export behind test/debug-only entry points.

---

### P8f — Analysis-core production orchestration wired and real-data checked (RESOLVED 2026-05-03)

**Priority:** RESOLVED

**Description:**
The production `run` entry path is no longer just a placeholder. The current production slice now covers discovery, normalization, hard QC, downstream geometry, ProLIF/IFP, residue contacts, convergence, condition-wise clustering, cluster signatures, residue importance, crystal anchoring, and report generation.

Implemented and verified points:
- `analysis/analysis_orchestrator.py` now provides a real production control path through the implemented analysis-core stages.
- `cli.py` now routes `lpmo-pipeline run` into that analysis-core orchestrator instead of the old placeholder path.
- `configs/production.analysis_core.example.yaml` now documents the current production config surface.
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py` and `tests/run_tests_scripts/test_analysis_core_real_cifs.sh` now exercise the real production CLI path on staged real AF3 CIFs.
- real-data smoke run `analysis_core_real_cifs_613251` (job 613251) staged 3 CIFs and produced `qc_report.json`, `pose_geometry.tsv`, `metrics.csv`, `summary.json`, `report.html`, and `run_manifest.json`.
- the smoke-run outcome was 1 `pass`, 1 `soft_flag`, 1 `hard_fail`; 2 poses were analyzed downstream.
- no production `geometry_debug.pdb` files were written in that run.

**Remaining action:**
Keep extending the same production path into the higher-level descriptive, predictive, CBM, and final report-table layers instead of creating separate side paths.

---

### P9 — Clustering primary method locked from pilot

**Status:** RESOLVED 2026-05-20 for method selection and production Stage 6 wiring.

**Description:**  
The clustering pilot and parameter-sensitivity run selected agglomerative
Jaccard as the global primary method for full analysis:
`linkage=average`, `distance_threshold=0.55`, `min_cluster_size=3`.
Sensitivity settings are agglomerative Jaccard `distance_threshold=0.45`,
`min_cluster_size=3`, and HDBSCAN Jaccard `min_cluster_size=3`,
`min_samples=null`.

Primary evidence: in the 94 formally clusterable pilot conditions,
agglomerative `0.55/min3` recovered clusters in 84 conditions, had median noise
fraction 0.52, and produced non-noise clusters with mean size 4.77 poses
(SD 2.80; median 4; range 3-21).

**Implemented solution:**  
1. The production full-analysis Stage 6 path now uses the selected
   agglomerative primary method rather than the older HDBSCAN default.
2. The `0.45/min3` agglomerative and `HDBSCAN min3` outputs remain predefined
   sensitivity analyses.
3. Parameter sweeps should not be re-run during the main analysis.

---

### P10 — R model choice for cluster-level regioselectivity not finalized

**Priority:** MEDIUM — needed for predictive modeling stage (Stage 14). Non-blocking for descriptive stages.

**Description:**  
OPEN_QUESTIONS.md item 9: primary model is `glmnet` (penalized logistic) with `glmer` (mixed effects) as sensitivity check. This is still a default, not a confirmed decision.

Update 2026-05-16:
Two detailed predictive implementation plans now exist:
- `c1_c4_predictive_analysis_plan_simplified.yaml` for C1/C4 regioactivity prediction.
- `substrate_activity_prediction_plan.yaml` for substrate activity prediction.

These plans specify elastic-net logistic regression, grouped CV by
`protein_id`, compact condition-level feature sets, and strict exploratory
interpretation. The remaining decision is whether to keep the originally
preferred R/glmnet implementation, use the sklearn implementation described in
the YAML plans, or keep both with one marked as the primary implementation.

**Proposed solution:**  
1. Confirm the primary implementation backend for the YAML predictive plans: R/glmnet, sklearn elastic-net logistic regression, or both with one primary.
2. Confirm whether `glmer` remains a sensitivity check, or whether the YAML-defined nested grouped CV models replace it for the first implementation.
3. Update OPEN_QUESTIONS.md item 9 to mark the model/backend choice as resolved.
4. Implement `scripts/ec_activity_mapping.R` first (step 19), then build the predictive modeling tables and models according to `c1_c4_predictive_analysis_plan_simplified.yaml` and `substrate_activity_prediction_plan.yaml`.

---

### P10b — CBM paired-analysis detailed plan added but not implemented

**Priority:** MEDIUM — needed for Stage 15 CBM side analysis after the core condition/cluster tables exist.

**Description:**
`cbm_full_length_vs_domain_only_analysis_plan.yaml` has been added as the
detailed implementation plan for paired comparison of full-length and
domain-only constructs. It specifies the paired row unit
`protein_id x substrate_class x dp`, required domain/CBM/linker region labels,
condition-level summaries, primary endpoints such as `bridge_fraction`,
`catalytic_domain_ifp_jaccard_distance`, `delta_C4_minus_C1_geometry_bias`,
`delta_qc_pass_fraction`, and `delta_cluster_entropy`, plus small-n reporting
rules and expected output tables.

**Proposed solution:**
1. Ensure upstream outputs needed by the plan exist: `condition_table.tsv`, `cluster_table.tsv`, `cluster_residue_signature.tsv`, `protein_condition_residue_scores.tsv`, and the construct-specific condition summaries.
2. Define or import residue region annotations for catalytic domain, CBM, linker, and other regions before computing CBM metrics.
3. Implement `analysis/cbm_comparison.py` and/or related scripts against `cbm_full_length_vs_domain_only_analysis_plan.yaml`.
4. Validate the paired table on a small set of proteins with both construct types before running the full CBM side analysis.

---

### P10c — Planned downstream TSVs still depend on unimplemented annotation layers

**Priority:** MEDIUM — these are required before descriptive, predictive, residue-importance, and CBM analyses can be run from the production outputs.

**Description:**
The current analysis-core production path now writes all TSV surfaces that are
directly supported by implemented stages: pose manifest/confidence/index,
QC attrition, pose geometry, pose IFP, pose residue contacts, convergence,
raw clustering, medoids, condition cluster summary, cluster IFP/residue
signatures, residue-importance tables, crystal anchor table, crystal geometry
table, crystal IFP diagnostic summary, and metrics.
For short smoke runs with observed contacts but zero retained non-noise
clusters, the Stage 16b residue tables now emit explicit zero-valued residue
rows instead of remaining header-only. Cluster-dependent tests can use
`tests/fixtures/clustering_stage_outputs/`, which contains the selected pilot
agglomerative Jaccard outputs (`distance_threshold=0.55`,
`min_cluster_size=3`), instead of depending on the short one-pose real-data
smoke to produce meaningful clusters. The following planned TSVs are still not
implemented because their upstream analysis layers are not implemented yet:

- `cluster_table.tsv`
- `family_aligned_residue_table.tsv`
- `family_residue_enrichment.tsv`
- `condition_table.tsv`
- `protein_summary_table.tsv`
- predictive modeling tables under `10_predictive/modeling_tables/`
- CBM paired-analysis outputs such as `cbm_construct_condition_summary.tsv` and `cbm_paired_comparison_table.tsv`

**Proposed solution:**
1. Promote pilot clustering outputs for at least one multi-pose real-data condition with retained clusters into a reusable fixture or smoke input, then validate non-zero cluster signatures and non-zero Stage 16b residue weights against it.
2. Implement `cluster_table.tsv` plus the later `condition_table.tsv` / `protein_summary_table.tsv` summary layers on top of the now-stable cluster-signature and residue-importance outputs.
3. Implement family-aligned residue enrichment only after the within-protein residue outputs are validated on clustered real data.
4. Implement predictive modeling and CBM paired analysis after `condition_table.tsv` exists.

---

### P11 — `enzyme_id` → `protein_id` refactor not completed in existing code

**Priority:** LOW — not a blocker for implementation, but creates inconsistency between old code and new table schemas.

**Description:**  
MASTERPLAN.md section 5 and OPEN_QUESTIONS.md item 16 (resolved 2026-04-21): canonical term is `protein_id`. Existing code still uses `enzyme_id` in some places.

**Proposed solution:**  
At next code refactor (not as a blocker), use project-wide symbol rename to replace all occurrences of `enzyme_id` with `protein_id` in Python source files under `src/`. Verify tests pass after rename.

---

### P12 — EC 1.14.99.- non-AA17 substrate label is provisional

**Priority:** LOW — a provisional label exists and does not block analysis.

**Description:**  
OPEN_QUESTIONS.md item 10: proteins with EC `1.14.99.-` that are not AA17 are assigned `substrate_class = xylan_or_other`, `regio_class = unknown`, `activity = xylan_like_oxidative`. This is a provisional mapping.

**Proposed solution:**  
Review the actual proteins in the dataset that fall into this category and check if more specific mappings are available from the CAZy database or primary literature. Update `scripts/ec_activity_mapping.R` if more specific labels are warranted. Document any change in the mapping script with a comment.

---

### P13 — `_chem_comp_bond` requirement appears legacy for AF3

**Priority:** LOW — consistency/documentation issue; currently non-blocking in code.

**Description:**  
AF3 mmCIF inputs generally do not include `_chem_comp_bond`. The normalization code now treats missing `_chem_comp_bond` as informational (not warning/failure), but multiple docs/config fields still describe it as required.

Potential risk:
- If a downstream AF3 step later depends on explicit bond-table data and no alternative source is used, that downstream step should fail with a clear error. This should not be silently assumed by normalization.

**Proposed solution:**
1. Confirm whether any downstream AF3 path truly requires `_chem_comp_bond`.
2. If not required, remove legacy “required” wording from config/spec/docs and eventually remove the presence check.
3. If required, implement and verify an explicit fallback source for bond information, then keep/upgrade the check where appropriate.

---

## Part 2 — Manual Checks

These checks must be performed by a human. They cannot be automated. Each check includes exact step-by-step instructions.

### Pending checks only

Completed checks have been removed from this section. Keep only checks that are still required.


### MC14 — Define and verify crystal anchoring reference metadata

**Why:** P6 and P7 remain open; crystal checks are weak until reference structures and residue definitions are explicit.

**Steps:**

1. Create/update `metadata/crystal_reference_list.tsv` with at least:
   - `protein_id`, `pdb_code`, `chain`, `has_ligand`, `family`, `regio_label`, `notes`

2. Create/update residue-definition metadata:
   - Either in `metadata/crystal_reference_list.tsv` notes
   - Or in `metadata/alignment_residue_definitions.tsv`

3. For each family, label residue source as `literature` or `proximity_fallback`.

4. Confirm crystal anchoring output includes the chosen alignment method per system.

---

*End of document.*
