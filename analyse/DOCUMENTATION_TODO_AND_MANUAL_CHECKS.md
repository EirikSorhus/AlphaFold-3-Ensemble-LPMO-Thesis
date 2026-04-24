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

### P8 — Privateer CLI version not autodetected (BLOCKED — Privateer not installed)

**Priority:** MEDIUM — Privateer v1 and v2 produce different JSON output formats. Incorrect parsing will produce silent errors.

**Description:**  
OPEN_QUESTIONS.md item 2: Privateer v1 vs v2 has different JSON output. **BLOCKER: Privateer is not installed in `analyse_env`.** Job 574004 (sbatch run on 2026-04-23 17:39:46) confirmed: `[ERROR] privateer is not available in analyse_env PATH`. 

Privateer must be installed before real-CIF verification can proceed.

**Proposed solution:**  
1. **FIRST:** Install Privateer in `analyse_env`. Check with cluster admins or project manager for installation path/method.
   - Likely source: Singularity container or conda package (if available in project channels).
   - Confirm installation by running `privateer --version` in analyse_env.
2. Once installed, in `qc/privateer_runner.py`, implement version autodetection: run `privateer --version` and parse the version string at startup.
3. Implement two JSON-parsing branches: one for v1 format, one for v2 format.
4. Log the detected version to `run_manifest.json`.
5. Add a contract test that runs Privateer on a known structure and checks that the output is parsed correctly for the detected version.
6. Add and run a manual real-CIF probe script under `tests/run_tests_scripts/` that stores `privateer --version`, help text, attempted JSON command templates, and raw stdout/stderr for several AF3 CIFs. Use this output to confirm what the JSON actually looks like before locking parser branches.

Manual verification target files for the first probe run (three distinct CIFs are currently referenced in context; pass any additional CIF paths as script arguments):
- `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG4/af3/latest/Q7SCE9_NAG4/seed-4_sample-1/Q7SCE9_NAG4_seed-4_sample-1_model.cif`
- `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA6/af3/latest/Q59930_STA6/seed-2_sample-0/Q59930_STA6_seed-2_sample-0_model.cif`
- `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/STA4/af3/latest/A0A0S2GKZ1_STA4/seed-2_sample-2/A0A0S2GKZ1_STA4_seed-2_sample-2_model.cif`

Outcome to record after the first run:
- which command variant produced JSON, if any
- whether stdout parsed as JSON
- top-level JSON keys or top-level type
- whether failures were tool-syntax failures, parse failures, or structure/content failures

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


### MC13 — Confirm Privateer version handling is implemented and tested

**Why:** P8 remains open; output parsing can silently fail if version format differs.

**Steps:**

1. Run `privateer --version` in the analysis environment and capture the exact version string.

2. Confirm `qc/privateer_runner.py` detects version and selects correct parser branch.

3. Run one contract test on known structure and verify parsed fields match expected schema.

4. Ensure detected version is written to `run_manifest.json`.

---

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
