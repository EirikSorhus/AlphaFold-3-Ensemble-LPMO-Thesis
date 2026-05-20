# LPMO Structure Prediction -> Analysis Pipeline - Integrated Master Plan (v2.3)

## 0. Priority Order

This file follows the active priority stack for implementation decisions:
1. `AF3_LPMO_pipeline_detailed_plan.md` — **primary and governing source** (updated 2026-04-21)
2. Latest user comments in the working thread
3. `MASTERPLAN.md` (this file) and `IMPLEMENTATION_PLAYBOOK.md`
4. Detailed downstream implementation plans:
  - `clustering_pilot_plan.yaml` — completed pilot plan and decision record for the selected IFP clustering method
   - `c1_c4_predictive_analysis_plan_simplified.yaml` — exploratory C1/C4 regioactivity predictive modeling
   - `substrate_activity_prediction_plan.yaml` — exploratory substrate activity prediction from AF3 ligand-condition summaries
   - `cbm_full_length_vs_domain_only_analysis_plan.yaml` — paired full-length vs domain-only CBM analysis
5. `plan_implementation_spec.txt` — **archived legacy spec** (PLACER steps obsolete; use for historical context only)
6. `plan_analyse.txt` — **archived legacy reference only** (research questions RQ1-RQ4 and AF3 tuning rationale as background)

**⚠️ CRITICAL: Large changes to order of analyses, which analyses are run, or core pipeline structure must NOT be made by AI without explicit user approval first. Any proposed changes in these areas require user review and confirmation.**

**⚠️ KNOWN RESOLVED DECISIONS (2026-04-21):**
- **A (RESOLVED)**: AF3 runs completed with `num_diffusion_samples=5` → 75 poses per system. All documents updated (IMPLEMENTATION_PLAYBOOK.md and OPEN_QUESTIONS.md item 14 updated 2026-04-21).
- **B (RESOLVED)**: PLACER is removed from the pipeline entirely.
- **C (RESOLVED)**: Use `protein_id` / "protein" as canonical term.
- **D (RESOLVED)**: Five separate pose-level tables (no merged `pose_table.tsv`).

## 1. Scope And Invariants

- Models: AF3 (**RF3/RosettaFold 3 og Boltz-2 er ekskludert fra analyse**)
- **Data availability**: Precomputed AF3 structures fully available at:
  - Domain-only: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core`
  - Full-length: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length`
- AF3 fixed parameters (main analysis): `num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5` (75 poses per protein–ligand condition; runs completed)
- Operational analyses: 9 independent sub-analyses
  - chitin_DP4, chitin_DP6, chitin_DP8
  - cellulose_DP4, cellulose_DP6, cellulose_DP8
  - starch_DP4, starch_DP6, starch_DP8
- Master format: mmCIF (ModelCIF-compatible)
- Prediction stage and analysis stage are separate; analysis consumes precomputed prediction artifacts.
- Primary analysis unit after pose generation is cluster.
- Main analyses must not aggregate cluster rows to one enzyme row.
- Enzyme-level summaries are secondary sensitivity analyses only.
- IFP clustering uses IFP features only; geometry is linked after clustering.
- Primary IFP clustering method is locked from the completed pilot: agglomerative Jaccard on contact-eligible IFP rows with `linkage=average`, `distance_threshold=0.55`, and `min_cluster_size=3`. HDBSCAN remains a predefined sensitivity path only.
- Atom names are not assumed consistent across models; mapping key is (element, CCD, local bond graph, 3D proximity).
- Chain schema: protein=A, glycans=B..D, metal=E.
- Normalization must hard-fail if no glycan residues remain after chain remap; protein+metal-only source CIFs are invalid analysis inputs and must not continue downstream as soft warnings.
- Routine successful normalization should keep one compact audit report (`normalize_report.json`). Per-pose atom-map/rename debug files are reserved for incomplete mapping or validation failure to avoid unnecessary output volume.
- `_chem_comp_bond` must be complete for all `comp_id`, and `_struct_conn` must include glycosidic + Cu coordination links.
- Preserve all computed numeric metrics; do not drop distance/angle/support fields from output tables.
- Independent per-pose stages may run with bounded workers from `n_jobs`. The clustering pilot now has a Slurm-array wrapper that splits domain-only and full-length selections into independent protein-level shards, so the pilot can use multiple jobs/nodes instead of one long 24 h allocation. Future production-scale runs can extend this pattern to stage-specific resource profiles when needed.
- Geometry plausibility thresholds are operationalized and locked in `configs/thresholds.yaml` (`geometry_plausibility.locked = true`, 2026-04-21).
- Statistical/descriptive analysis should prefer R where possible (Python wrappers can orchestrate).

## 2. Required Output Tables

Tables are grouped by analysis level. The new detailed plan (v1.0) introduces a split structure for
pose-level outputs. **See OPEN_QUESTIONS.md item 17 for the decision on merged vs split pose table.**

### 2.1 Pose-level tables

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| pose_manifest.tsv | TSV | pose | Pose metadata, paths, AF3 confidence fields, parsing status |
| pose_confidence.tsv | TSV | pose | AF3 confidence metrics only (ipTM, pTM, pLDDT, etc.) |
| pose_ifp_table.tsv | TSV | pose | IFP vector, feature names, per-type contact counts |
| pose_residue_contact_table.tsv | TSV | pose × residue × interaction_type | Residue-level contact table (new in v1.0) |
| pose_geometry.tsv | TSV | pose | Geometry metrics and status labels |

### 2.2 Cluster-level tables

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| cluster_table.tsv | TSV | cluster | Main descriptive unit; includes geometry annotation |
| cluster_assignments.tsv | TSV | pose | Cluster membership, noise flag |
| medoid_manifest.tsv | TSV | cluster | Medoid pose ID and structure path |
| cluster_ifp_signature.tsv | TSV | cluster × feature | IFP feature frequency per cluster |
| cluster_residue_signature.tsv | TSV | cluster × residue | Residue contact frequency per cluster (new in v1.0) |
| predictive_cluster_table.tsv | TSV | cluster | Main predictive modeling table |

### 2.3 Condition / protein-level tables

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| condition_table.tsv | TSV | protein–ligand condition | Central repeated-measures table (new in v1.0) |
| condition_cluster_summary.tsv | TSV | condition | Cluster count, entropy, noise fraction |
| protein_summary_table.tsv | TSV | protein | Secondary sensitivity table only (was enzyme_summary_table.tsv) |
| qc_attrition_table.tsv | TSV | condition | Attrition by reason code |

### 2.4 Interpretation / side-analysis tables

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| protein_condition_residue_scores.tsv | TSV | protein × condition × residue | Occupancy-weighted contact scores (new in v1.0) |
| protein_residue_regio_delta.tsv | TSV | protein × residue | C1 vs C4 delta contact scores (new in v1.0) |
| family_aligned_residue_table.tsv | TSV | family × alignment position | Within-family alignment (optional, new in v1.0) |
| family_residue_enrichment.tsv | TSV | family × residue | Family-level residue enrichment (optional, new in v1.0) |
| condition_patch_summary.tsv | TSV | condition | Residue-property patch summaries (new in v1.0) |
| protein_patch_summary.tsv | TSV | protein | Protein-level patch summaries (new in v1.0) |
| crystal_anchor_table.tsv | TSV | medoid/fallback × crystal reference | Crystal sanity check; includes pocket RMSD, IFP comparability, and crystal geometry summary fields |
| crystal_geometry_table.tsv | TSV | prepared crystal reference | Crystal-side C1/C4 geometry for ligand-bound references |
| crystal_ifp_diagnostic_summary.tsv | TSV | diagnostic scope | Crystal IFP eligibility/exclusion counts and percentages |
| cbm_comparison_table.tsv | TSV | protein–subanalysis | Domain-only vs full-length paired analyses |

### 2.5 Run artifacts

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| qc_report.json | JSON | pose | Pre-QC proximity + PB + Privateer + geometry gates |
| clusters.json | JSON | cluster | Occupancy, support, medoid, outlier stats |
| cluster_signatures.json | JSON | cluster | Median/IQR geometry + IFP enrichments |
| metrics.csv | CSV | mixed | Flat export for plotting and QA |
| summary.json | JSON | run | Dataset and gate summaries |
| run_manifest.json | JSON | run | Hashes, tool versions, seeds, gate outcomes |

## 3. Core Pipeline (Per Sub-Analysis)

The pipeline follows the stage numbering in `AF3_LPMO_pipeline_detailed_plan.md` (v1.0).
Implementation order follows `IMPLEMENTATION_PLAYBOOK.md`.

### Stage 1 - AF3 Parsing and Pose Manifest
- Input: precomputed predictions from AF3 only (RF3 and Boltz-2 excluded).
- Parse AF3 outputs into `pose_manifest.tsv` and `pose_confidence.tsv`.
- Store ipTM, pTM, mean_pLDDT, ligand_interface_confidence, pae_summary per pose.
- ipTM must be retained as descriptive metric but not used as primary validity criterion.
- Compute and store per-condition ipTM mean/SD/median/IQR (separately for all poses and geometry-computable poses only).
- AF3 parameters (RESOLVED 2026-04-21): `num_seeds=15`, `num_diffusion_samples=5` → 75 poses per protein–ligand condition.

### Stage 2 - Pre-QC / Hard QC (was Step 2–3)
- Pre-QC active-site proximity gate: `min_cu_ligand_distance <= 10.0 Å` (hard gate before PoseBusters/Privateer; see `thresholds.yaml: hard_gates.active_site_proximity_max_a`).
- PoseBusters currently runs in built-in `dock` mode on auto-split ligand/protein inputs. The pipeline does not override PoseBusters' `intermolecular_distance` defaults, so the far-away check uses `max_distance=5.0 Å` and `minimum_distance_to_protein` is the renamed clash/no-clashes check rather than that distance threshold.
- Hard fail criteria: severe PoseBusters fail, severe Privateer fail, Cu missing, ligand missing/broken/unparsable, structure corrupt.
- Soft flags: retained as metadata; do NOT exclude before IFP.
- Required gate: atom mapping coverage = 100%.
- Keep all computed distances in output tables even for dropped poses.

### Stage 3 - IFP Generation (was Step 4 / Branch A)
- Use AF3 structures as-is. Do NOT reposition Cu. Do NOT add virtual oxyl/H.
- Generate binary ProLIF vectors with fixed interaction types per sub-analysis.
- Current validated standalone audit set uses all 9 interactions in `configs/prolif_features.yaml`: `HBDonor`, `HBAcceptor`, `Hydrophobic`, `PiStacking`, `Anionic`, `Cationic`, `CationPi`, `PiCation`, `VdWContact`.
- Main clustering features are selected after the pilot feature audit. Default main-clustering includes are Hbond donor, Hbond acceptor, and aromatic/stacking; hydrophobic and cation-pi are conditional; van der Waals/close-contact features are descriptive only by default.
- Ligand handling must keep each monosaccharide as a separate ligand residue in the feature space. Current flattened feature naming contract: `ligand_residue|protein_residue|interaction`.
- `pose_ifp_table.tsv` is currently validated on real data as a standalone slice; `pose_residue_contact_table.tsv` remains planned and should reuse the same ligand-resolved residue labels.

### Stage 4 - Geometry Branch (was Step 5 / Branch B)
- Identify histidine-brace atoms, build brace plane.
- Reposition Cu for geometry branch only (not fed back to IFP or main analysis).
- Place virtual oxyl. Place virtual H for C1 and C4.
- Compute: `oxyl_H_C1_distance`, `oxyl_H_C4_distance`, `Cu_C1_distance`, `Cu_C4_distance`, `attack_angle_C1`, `attack_angle_C4`, `sugar_face_orientation`, `ring_normal_vs_brace_normal`.
- Geometry status: `geometry_not_computable` / `geometry_computable_implausible` / `geometry_plausible` / `geometry_highly_plausible`.
- Plausibility thresholds (operational, specified and locked in `thresholds.yaml: geometry_plausibility`): oxyl-H window 1.5–4.0 Å, reference optimum ~2.1 Å, tighter window 1.8–2.5 Å.
- Output: `pose_geometry.tsv`.

### Stage 5 - Convergence Metrics (new in v1.0)
- Per protein–ligand condition: compute ligand RMSD to a fixed reference pose.
- Current production implementation uses the top-ranked QC-passing pose within the condition as the fixed convergence reference. Cluster medoids are used later for crystal anchoring and structural figures, not as an implicit filter on convergence.
- Per pose: `ligand_rmsd_to_reference`, `convergent_flag` (RMSD < `convergent_rmsd_max_a`; see `thresholds.yaml: convergence.convergent_rmsd_max_a`).
- Per condition: `convergence_fraction`, `median_ligand_rmsd`, `iqr_ligand_rmsd`.

### Stage 6 - Clustering (was Step 6)
- Input: QC-passing poses with successful IFP generation and enough non-vdW contact signal.
- One clustering per protein–ligand condition (AF3-only simplification).
- Method locked from the 2026-05-20 clustering pilot decision:
  agglomerative Jaccard on contact-eligible IFP rows,
  `linkage=average`, `distance_threshold=0.55`, `min_cluster_size=3`.
  The pilot selection evidence used `main_contact_eligible_ifp_matrix.csv`.
- Contact eligibility is applied after IFP; main pilot rule is
  `min_non_vdw_interactions >= 2` and `min_non_vdw_contact_residues >= 1`,
  with lenient/strict sensitivity rules.
- Conditions with `n_contact_eligible < 10` are reported as
  `insufficient_clusterable_signal` rather than formally clustered.
- Use one global primary method and fixed parameters; do not choose method
  separately per condition.
- Sensitivity settings: agglomerative Jaccard with `distance_threshold=0.45`,
  `min_cluster_size=3`, and HDBSCAN Jaccard with `min_cluster_size=3`,
  `min_samples=null`.
- Noise/low-support groups are reported separately and retained in raw output.
- Minimum cluster occupancy for main summaries: `>= 0.05` (confirmed; see `thresholds.yaml: cluster_inclusion.min_occupancy`).
- Output: `cluster_assignments.tsv`, `medoid_manifest.tsv`, `condition_cluster_summary.tsv`.

### Stage 7 - Cluster Annotation (was Step 7–8)
- Attach geometry/QC/support features to each cluster (median/IQR summaries).
- Assign cluster type: `C1_compatible`, `C4_compatible`, `mixed_compatible`, `non_plausible`, `uncertain`.
- Thresholds for cluster type stored in `configs/thresholds.yaml`.
- Build `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, and
  `cluster_signatures.json` from existing Stage 6 cluster membership and
  medoids; do not re-cluster or re-select medoids during annotation.

### Stage 8 - Residue Importance Analysis (new in v1.0)
- Consumes Stage 7 cluster annotation outputs (`cluster_residue_signature.tsv`,
  `cluster_ifp_signature.tsv`, `cluster_signatures.json`) rather than
  recomputing residue/cluster membership from raw Stage 6 tables.
- **Within-protein**: occupancy-weighted residue contact scores per protein–condition.
  - `residue_contact_score = sum(cluster_occupancy × residue_contact_frequency_in_cluster)`
  - Output: `protein_condition_residue_scores.tsv`
- **C1 vs C4 delta**: compare residues enriched in C1-compatible vs C4-compatible clusters within each protein.
  - Output: `protein_residue_regio_delta.tsv`
- **Within-family** (optional): align proteins within family, map residues to alignment columns.
  - Output: `family_aligned_residue_table.tsv`, `family_residue_enrichment.tsv`
- **Cross-dataset patch-level**: aromatic/polar/charged/hydrogen-bond contact density and region flags where explicitly annotated.
  - Output: `condition_patch_summary.tsv`, `protein_patch_summary.tsv`
- **Important**: raw residue numbers are NOT directly comparable across unrelated proteins. Split into within-protein, within-family, and property-level analyses.

### Stage 9 - Protein–Ligand-Condition Summaries (new in v1.0)
- One row per protein–ligand condition in `condition_table.tsv`.
- This is the central table for repeated-measures comparisons across ligand conditions.
- Repeated measures: ligand type and DP as repeated conditions within protein.

### Stage 10 - Protein-Level Summaries (secondary)
- Build `protein_summary_table.tsv` (was `enzyme_summary_table.tsv`).
- Canonical term: **protein** / `protein_id` (RESOLVED 2026-04-21, OPEN_QUESTIONS.md item 16). Existing code using `enzyme_id` is updated at next refactor; not a blocker.
- Secondary sensitivity table only; do NOT use as main analysis table.

### Stage 11 - Crystal Sanity-Check (was Step 12)
- Input: every retained cluster medoid, not all poses.
- If a condition has no retained clusters, use the top-level AF3 model CIF as a crystal-anchoring fallback only when it passes hard QC; this fallback is not included in normal pose, clustering, or medoid denominators.
- Alignment: current operational implementation uses local Kabsch alignment on shared pocket C-alpha atoms; pocket residues come from ligand/Cu proximity in holo references or sequence-projected representative pockets for apo references.
- Crystal IFP similarity is comparable only when the prepared crystal IFP passes the same non-vdW contact-eligibility rule used for pose clustering. VdW-only, zero-contact, and low-specific-contact crystal IFPs keep RMSD/geometry outputs but leave IFP Tanimoto non-comparable.
- Ligand-bound crystal references get C1/C4 geometry computed with the same downstream geometry fields used for poses.
- Output: `crystal_anchor_table.tsv`, `crystal_geometry_table.tsv`, `crystal_ifp_diagnostic_summary.tsv`.
- Crystal mismatch does NOT invalidate a cluster; context and plausibility only.

### Stage 12 - Removed

~~Optional PLACER Sensitivity Branch~~ — **PLACER removed entirely (decision 2026-04-21).**

### Stage 13 - Descriptive Analysis (was Step 9)
- Primary descriptive unit: cluster.
- Required outputs: QC attrition plot, cluster count distribution, top-cluster occupancy distribution, cluster entropy, geometry plausibility fractions, ipTM vs QC/geometry/occupancy, residue contact heatmaps, medoid structural figures, and clustering rate versus crystal-structure coverage.
- Cross-condition comparisons: chitin vs cellulose vs starch; DP4 vs DP6 vs DP8; domain-only vs full-length.
- Crystal-coverage descriptive figure: plot clustering rate against the number of available crystal structures per protein as the primary view, with a binary `has_crystal_reference` yes/no summary as a secondary panel or grouped overlay.

### Stage 14 - Exploratory Predictive Analysis (was Step 10)
- Keep simple: at most 1–2 tasks (e.g., C1 vs C4, chitin vs cellulose preference).
- Primary rows: cluster.
- Grouped CV at protein level (all clusters from same protein in one fold).
- Report: balanced accuracy, macro F1, AUROC where applicable.
- Results are exploratory; do NOT overinterpret as causal biology.
- Detailed implementation plans:
  - `c1_c4_predictive_analysis_plan_simplified.yaml` specifies the C1/C4 regioactivity modeling table, labels, compact feature set, nested grouped CV, and binary C1/C4 model outputs.
  - `substrate_activity_prediction_plan.yaml` specifies the substrate activity modeling table, protein-substrate activity labels, grouped CV, feature aggregation, and reporting for chitin/cellulose/starch prediction.

### Stage 15 - CBM Paired Analysis (was Step 13)
- Only proteins with both domain-only and full-length constructs.
- Paired comparisons within each protein–ligand condition.
- Output: `cbm_comparison_table.tsv`.
- Statistics: paired Wilcoxon signed-rank; report effect sizes and direction.
- Detailed implementation plan: `cbm_full_length_vs_domain_only_analysis_plan.yaml` specifies the matched construct-pair design, CBM/linker region definitions, condition summaries, primary endpoints, paired statistics, and required CBM analysis outputs.

## 4. Activity Label Mapping From EC

Maintain explicit EC-to-activity mapping from metadata (implemented in `scripts/ec_activity_mapping.R`):

- `1.14.99.54` -> cellulose, C1-hydroxylating
- `1.14.99.56` -> cellulose, C4-dehydrogenating
- `1.14.99.53` -> chitin, C1/C4 mixed
- `1.14.99.55` -> starch, C1-hydroxylating
- `1.14.99.-` -> special check:
  - AA17 family: homogalacturonan, C4 oxidation
  - otherwise: xylan/other oxidative activity (provisional label)

## 5. Entity Identifiers

Canonical term: **protein** / `protein_id` (resolved 2026-04-21; see OPEN_QUESTIONS.md item 16). Use `protein_id` in all new code and table schemas. Existing code using `enzyme_id` will be updated at next refactor.

Stable ID formats (from `AF3_LPMO_pipeline_detailed_plan.md`):

- **Pose ID**: `{protein_id}__{construct_type}__{substrate_class}_DP{dp}__seed{seed}__sample{sample}`
- **Condition ID**: `{protein_id}__{construct_type}__{substrate_class}_DP{dp}`
- **Cluster ID**: `{condition_id}__cluster{n}`

## 6. Predictive Cluster Table Requirements

Minimum required columns:

- analysis_id
- protein_id (⚠️ was enzyme_id — see OPEN_QUESTIONS.md item 16)
- family
- cbm_status
- substrate_class
- DP
- cluster_id
- occupancy
- qc_factor
- support_factor
- cluster_weight
- convergence_support (note: cross_model_support is not applicable for AF3-only pipeline)
- invalid_geometry_fraction
- plausible_geometry_fraction
- scorable_geometry_fraction
- median_oxyl_H_C1
- median_oxyl_H_C4
- median_Cu_C1
- median_Cu_C4
- median_attack_angle_C1
- median_attack_angle_C4
- face_orientation_summary
- ifp_features_selected
- experimental_regio_label
- experimental_substrate_label (was experimental_ligand_specificity_label)

Suggested weight definition:

- cluster_weight = occupancy * qc_factor * support_factor

## 7. Optional Post-Analysis Tuning

Tuning is optional and runs after baseline analysis as a comparative/sensitivity workflow.

- It must not block the first main analysis delivery.
- Only AF3 is used; no cross-model comparison in tuning.
- If tuning is executed later, tuned settings are applied in follow-up reruns and documented in the manifest.

## 8. Failure Policy

| Event | Hard? | Action | Keep numeric metrics? |
|---|---|---|---|
| atom_mapping < 100% | YES | skip | yes, store coverage + reason; write atom-map/rename debug artifacts |
| pre-QC ligand too far from active site | YES | drop before PB/Privateer | yes, store distances |
| glycan not CCD-valid for Privateer prep | YES | skip/drop | yes |
| critical PoseBusters error | YES | drop pose | yes |
| severe Privateer fail | YES | drop pose | yes |
| Cu-His outside 1.9-2.6 A (selected brace N only: His1:N, His1:ND1, third non-His1 histidine N) | YES | drop pose | yes |
| soft PB warning | NO | keep + flag | yes |
| clustering noise / low-support group | NO | keep with method-specific label | yes |
| low crystal similarity | NO | keep + flag | yes |

## 9. Reproducibility

Every run writes `run_manifest.json` including:

- pipeline_version, timestamp, git_commit, config_hash
- tool_versions (Gemmi, PoseBusters, Privateer, MDAnalysis, ProLIF, selected clustering backend, R)
- worker/resource settings such as `n_jobs` where they affect runtime behavior
- seeds and input checksums (SHA-256)
- gate outcomes including pre-QC active-site proximity
- references to prediction artifacts reused by analysis

## 10. RQ To Output Mapping

| RQ | Primary analysis unit | Key metrics | Artifact |
|---|---|---|---|
| RQ1 C1/C4 regioselectivity | cluster | Cu-C1/C4, oxyl-H, attack angles, occupancy | cluster_table.tsv, predictive_cluster_table.tsv |
| RQ2 substrate specificity | cluster | substrate/DP-specific occupancy + IFP signatures | cluster_table.tsv, predictive_cluster_table.tsv |
| RQ3 CBM effect | paired protein–subanalysis | occupancy shifts, support, CBM proximity | cbm_comparison_table.tsv |
| RQ4 crystal anchoring | cluster/protein–ligand condition | pocket RMSD (current operational implementation: optimized local C-alpha alignment on ligand/Cu-proximal pocket via gemmi/numpy Kabsch), ligand/proximal RMSD, contact-eligible IFP similarity, crystal-side C1/C4 geometry comparison | crystal_anchor_table.tsv, crystal_geometry_table.tsv, crystal_ifp_diagnostic_summary.tsv |
| RQ4b crystal coverage vs clustering | protein / protein–ligand condition | clustering rate by number of available crystal structures, plus binary has/no-has crystal-reference comparison | condition_table.tsv, crystal_anchor_table.tsv, crystal_ifp_diagnostic_summary.tsv, figure_manifest.tsv |

## 11. Implementation Status

Pipeline status has moved from pseudocode-only to step-by-step implementation.
Execution follows the IMPLEMENTATION_PLAYBOOK.md priority order.
Each step is implemented, tested, and verified before the next begins.

## 12. Third-Party Code

Adapted patterns from the following MIT-licensed projects (see ATTRIBUTION.md):
- **PoseBench** (BioinfoMachineLearning/PoseBench): CIF→PDB conversion via PDBFixer, PoseBusters Python API usage, mol_table pattern.
- **benchmarking-af3** (lyulab/benchmarking-af3): Pocket residue identification via proximity cutoff, residue mapping between structures, metrics aggregation patterns.

## 13. Alignment Strategy (Crystal Anchoring)

Alignment of predicted structures against crystal references currently uses a
local pocket C-alpha superposition implemented with gemmi/numpy Kabsch:

1. **Alignment atoms**: Shared protein C-alpha atoms from the selected pocket.
2. **Pocket residues**: Current operational heuristic is protein residues within 5 A of ligand or Cu in holo references. Apo references can reuse a sequence-projected pocket from the medoid/representative pose, with residue-name normalization for variants such as `HIC -> HIS`.
3. **Literature refinement**: Family-specific substrate-recognition residues from literature are still preferred if they should replace or refine the current heuristic. Document which method was used per system.
4. **Measurement**: After optimal local superposition, report pocket RMSD, per-residue deviations, and ligand/proximal RMSD where available. Always label results as optimized local alignment to distinguish from global alignment.
5. **Note**: Only AF3 is used; cross-model differences are not applicable. PyMOL `pair_fit` parity/hardening remains possible future work, but is not the current operational backend.
