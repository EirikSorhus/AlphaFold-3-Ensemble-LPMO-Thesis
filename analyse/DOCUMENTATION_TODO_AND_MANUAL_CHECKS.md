# Documentation TODO and Manual Checks

Version: 1.0  
Created: 2026-04-21  
Purpose: Prioritized open problems with proposed solutions, and a manual-check section with step-by-step instructions.

---

## Part 1 — Prioritized Open Problems

Problems are ordered from highest to lowest priority.

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

### P6 — Crystal reference structures not specified

**Priority:** MEDIUM — blocks crystal sanity-check stage (Stage 11). Non-blocking for main pipeline steps 1–10.

**Description:**  
OPEN_QUESTIONS.md item 5: which PDB codes to use as crystal references is not decided. The plan requires ligand-bound crystals for AA9 at minimum (4EIS, 5ACF mentioned as defaults), but other families and the full reference list are not specified.

**Proposed solution:**  
1. For each LPMO family in the dataset, identify available PDB structures with bound oligosaccharide ligands (search PDB with family annotation and filter for saccharide ligands).
2. Minimum: at least one C1-crystallized and one C4-crystallized structure per family if available.
3. Document the chosen PDB codes in `metadata/crystal_reference_list.tsv` with columns: `protein_id`, `pdb_code`, `chain`, `has_ligand`, `family`, `regio_label`, `notes`.
4. Update OPEN_QUESTIONS.md item 5 when the list is finalized and mark P6 in DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md as RESOLVED.

---

### P7 — Substrate-recognition residues for PyMOL pair_fit not defined per family

**Priority:** MEDIUM — required for optimal crystal anchoring alignment. The fallback (proximity-based) is available but literature-based is preferred.

**Description:**  
OPEN_QUESTIONS.md item 13: the substrate-recognition surface residues used in PyMOL `pair_fit` alignment must be identified per LPMO family. Literature-based selection is preferred; proximity-based is the fallback.

**Proposed solution:**  
1. For each LPMO family in the dataset, perform a literature search to identify known substrate-binding residues.
2. Document in `metadata/crystal_reference_list.tsv` or a separate `metadata/alignment_residue_definitions.tsv` with columns: `family`, `residue_number`, `residue_name`, `source`, `notes`.
3. If literature data is unavailable for a family, use proximity cutoff (5 Å from ligand in the predicted structure), and flag this as `source = proximity_fallback` in the metadata.
4. Document which method was used per system in `crystal_anchor_table.tsv` (the `alignment_method` column).

---

### P8 — Privateer SIF execution path and output contract (RESOLVED 2026-04-30)

**Priority:** RESOLVED

**Description:**  
Privateer is now integrated via SIF instead of host PATH installation. The resolved execution path is `apptainer run --cleanenv /cluster/projects/nn1003k/prog/privateer/privateer.sif ...` with explicit bind mounts. The resolved output contract is `-mode ccp4i2` plus parsing of `validation_data-privateer`, not JSON stdout.

Implemented and verified points:
- `qc/privateer_runner.py` builds bind-aware SIF invocations and parses `validation_data-privateer`.
- `get_privateer_version()` queries the SIF with `-list`, and `utils/manifest.py` now records that version instead of calling `privateer -V` on the host.
- hard QC uses the runner through `hard_qc_orchestrator.py`, including batched Privateer dispatch for eligible poses.
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

**Proposed solution:**
1. Trace each consumer of `qc_report.json`, pose summary tables, and later-stage analysis inputs.
2. Verify dropped poses are excluded only at the downstream-analysis boundary, not earlier when artifacts are generated.
3. Verify soft-flagged poses remain included and visibly marked.
4. Add one integration-level check at the final filtering boundary once the consuming surface is identified.

---

### P9 — HDBSCAN parameters not locked in thresholds.yaml

**Priority:** MEDIUM — parameters must be locked before the full analysis run to prevent in-run parameter searching (anti p-hacking rule from MASTERPLAN.md stage 6).

**Description:**  
MASTERPLAN.md Stage 6 requires HDBSCAN `min_cluster_size` and `min_samples` to be stored in `configs/thresholds.yaml` and locked after tuning. The values are currently unspecified (tuning not yet run). The minimum cluster occupancy (≥ 0.05) must also be confirmed.

**Proposed solution:**  
1. Run HDBSCAN on a representative subset of real AF3 data (e.g., 2–3 protein–ligand conditions) with a small parameter sweep to identify stable parameter choices.
2. Commit the chosen parameters to `configs/thresholds.yaml` under `clustering:`.
3. Document the selection rationale in a comment in `thresholds.yaml`.
4. Do not re-run the sweep during the main analysis.

---

### P10 — R model choice for cluster-level regioselectivity not finalized

**Priority:** MEDIUM — needed for predictive modeling stage (Stage 14). Non-blocking for descriptive stages.

**Description:**  
OPEN_QUESTIONS.md item 9: primary model is `glmnet` (penalized logistic) with `glmer` (mixed effects) as sensitivity check. This is still a default, not a confirmed decision.

**Proposed solution:**  
1. Confirm the primary model as penalized logistic regression (`glmnet`) with grouped CV at protein level, as specified in MASTERPLAN.md stage 14.
2. Confirm `glmer` as the sensitivity check.
3. Update OPEN_QUESTIONS.md item 9 to mark as resolved.
4. Implement `scripts/ec_activity_mapping.R` first (step 19), then `analysis/predictive_models.py` + R scripts (step 20).

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
