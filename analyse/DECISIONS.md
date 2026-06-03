# Implemented Decision Logic

This file documents decision logic that is active in the current codebase.
It is intentionally implementation-grounded:

- Source of truth for this file is the runtime code under `src/` plus the config files that the code actually loads.
- Tests are used as confirmation of intended behavior.
- If a plan document and the current implementation differ, this file follows the implementation.

## ProLIF Implicit-H Migration (2026-05-22)

ProLIF production code now uses
`/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env` and the
implicit hydrogen-bond interactions from ProLIF. The raw IFP surface is limited
to `ImplicitHBAcceptor`, `ImplicitHBDonor`, and `VdWContact`.

ProLIF inputs are non-protonated PDB/PDB artifacts written under
`analysis_export/`: `complex_for_prolif.pdb` and
`ligand_only_for_prolif.pdb`. `complex_H.pdb` and `ligand_for_prolif.mol2`
remain legacy protonation artifacts and are not ProLIF inputs.

Main production clustering excludes `VdWContact` and uses only
`ImplicitHBAcceptor`/`ImplicitHBDonor`; VdW remains available for raw IFP,
residue-contact diagnostics, and contact-eligibility accounting. Hydrophobic
interaction predictors/fractions were removed, and the H-bond fraction now
counts only the two implicit H-bond types.

## Scope And Traceability

The main decision surfaces currently live in these code paths:

- `src/lpmo_pipeline/analysis/prolif_ifp.py` for raw IFP feature construction and contact eligibility.
- `src/lpmo_pipeline/analysis/clustering_pilot.py` for pilot-only feature prevalence audits and main-matrix selection.
- `src/lpmo_pipeline/analysis/clustering_agglomerative.py` for the current production clustering path.
- `src/lpmo_pipeline/analysis/clustering_hdbscan.py` for the retained HDBSCAN sensitivity path and shared clustering helpers.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py` for which poses and matrices actually move downstream.
- `src/lpmo_pipeline/analysis/cluster_signatures.py` for cluster type assignment, the flat `cluster_table.tsv` export, and Stage 16 signature tables.
- `src/lpmo_pipeline/analysis/residue_contact_extraction.py` for the residue-contact surface consumed by Stage 16.
- `src/lpmo_pipeline/analysis/residue_importance.py` for Stage 16b residue and patch aggregation.
- `src/lpmo_pipeline/qc/posebusters_runner.py` for the active PoseBusters input contract and failure classification.

## 1. Raw IFP Feature Space

### Interaction types used to build raw IFPs

The raw ProLIF feature surface is controlled by `load_prolif_features_config()` in `src/lpmo_pipeline/analysis/prolif_ifp.py`.

- The code first reads `configs/prolif_features.yaml`.
- If `active_interaction_types` is present there, that list becomes the active raw IFP interaction set.
- If the config is missing or empty, the fallback default is the three-type set hard-coded in `prolif_ifp.py`:
  - `ImplicitHBAcceptor`
  - `ImplicitHBDonor`
  - `VdWContact`

Current `configs/prolif_features.yaml` explicitly activates those three interaction types. ProLIF uses non-protonated PDB/PDB inputs from `analysis_export/`.

### Feature naming contract

The flattened feature name is parsed as:

`ligand_residue|protein_residue|interaction`

Examples of current downstream assumptions:

- `parse_ifp_feature_name()` in `prolif_ifp.py` requires exactly three `|`-separated components.
- `build_pose_residue_contact_rows()` in `residue_contact_extraction.py` uses the same three-part format directly.
- Protein residue labels must match the current `RESNAME<number>.<chain>` shape, for example `ASN10.A`.

This means Stage 13, Stage 13b, Stage 16, and Stage 16b all currently depend on the same flattened feature-name contract.

## 2. Contact Eligibility Before Clustering

Formal clustering is not gated on total contact count alone. The active rule is `evaluate_contact_eligibility()` in `src/lpmo_pipeline/analysis/prolif_ifp.py`.

### Active thresholds

`load_contact_eligibility_rule()` reads `contact_eligibility.main_rule` from `configs/prolif_features.yaml`.

Current active thresholds are:

- `min_non_vdw_interactions = 2`
- `min_non_vdw_contact_residues = 1`

The config also contains sensitivity variants, but the active production rule remains the 2/1 gate above. A single non-VdW interaction against one residue is not currently treated as sufficient specific-contact signal for clustering or pose-vs-crystal IFP comparison.

### Active exclusion classes

The code assigns one of three exclusion classes before clustering:

- `null_ifp`
  - Assigned when `n_total_contacts == 0`.
- `vdw_only`
  - Assigned when the pose has contacts, but every active contact is `VdWContact`.
- `low_specific_contact`
  - Assigned when the pose has some non-vdW signal, but either:
    - `n_non_vdw_interactions < min_non_vdw_interactions`, or
    - `n_non_vdw_contact_residues < min_non_vdw_contact_residues`.

A pose is contact-eligible only if none of those exclusion classes applies.

### How the orchestrator applies the gate

`run_analysis_core()` in `src/lpmo_pipeline/analysis/analysis_orchestrator.py` applies the rule as follows:

- Contact eligibility is evaluated for poses with `IFPResult.status` in `{ "ok", "zero_contacts" }`.
- Only poses with `status == "ok"` and `eligible == True` are added to `clusterable_indices` and sent into Stage 6 clustering.
- Poses with `zero_contacts`, `vdw_only`, or `low_specific_contact` are still recorded in per-pose and per-condition summaries, but they do not appear in the actual clustering input matrix.

As a result, the condition-level clustering summary fractions are currently calculated over the IFP-success surface that received an eligibility classification, not over all discovered poses and not over raw QC-pass counts alone.

## 3. Which IFP Features Are Actually Used For Clustering

This is currently split between the production Stage 6 path and the pilot-only comparison path.

### Current production Stage 6 behavior

The current production clustering path is `AgglomerativeJaccardClusterer` in
`src/lpmo_pipeline/analysis/clustering_agglomerative.py`, wired from
`run_analysis_core()` in
`src/lpmo_pipeline/analysis/analysis_orchestrator.py`.

The important implementation detail is:

- Production Stage 6 first evaluates contact eligibility on the full raw
  ProLIF result.
- For the agglomerative clustering input, it then builds a condition-local
  `main_clustering_ifp_matrix.csv` from the contact-eligible poses.
- That matrix uses the same main-feature policy as the pilot path:
  `ImplicitHBAcceptor` and `ImplicitHBDonor` by default, with explicit
  `extra_include_interaction_types` support and `VdWContact` excluded.
- Formal production clustering is allowed only when
  `n_contact_eligible >= minimum_clusterable_n` and the condition-local main
  matrix has at least one feature. Otherwise the condition is reported with
  `insufficient_clusterable_signal`, `empty_main_matrix`, or
  `no_clusterable_poses` in `condition_cluster_summary.tsv` and
  `condition_table.tsv`.
- The raw all-interaction IFP matrix is still written and retained for audit,
  contact summaries, and downstream descriptive signatures.

In practice, residue-level feature columns can differ by protein-condition
because the observed contacting residues differ, but the production interaction
type policy is fixed unless config explicitly opts in to extra interaction
types.

2026-05-21 real-data smoke validation:

- `analysis_core_real_cifs_1152258` exercised a no-contact-eligible condition
  and recorded `clustering_status = no_clusterable_poses` with no main matrix,
  as expected.
- `analysis_core_real_cifs_1152273` exercised a contact-eligible but
  insufficient-support condition. It wrote
  `ifp_matrices/A0A0S2GKZ1__domain_only__chitin_DP8/main_clustering_ifp_matrix.csv`
  with one row and two main features, propagated
  `clustering_status = insufficient_clusterable_signal` to
  `condition_cluster_summary.tsv` and `condition_table.tsv`, and did not emit
  cluster assignments or medoids.

### Pilot-only main feature selection

The pilot path in `src/lpmo_pipeline/analysis/clustering_pilot.py` builds two matrices for the same selected pose subset:

- `raw_contact_eligible_ifp_matrix.csv`
  - all raw features for contact-eligible poses.
- `main_contact_eligible_ifp_matrix.csv`
  - pilot-filtered features for the same pose subset.

The selection logic in `select_main_clustering_features()` is currently:

- Default included interaction types:
  - `ImplicitHBAcceptor`
  - `ImplicitHBDonor`
- Default excluded interaction types:
  - `VdWContact`
- Any interaction type not in the include set is dropped even if it is prevalent.
- Extra types can be admitted through `extra_include_interaction_types` in the production config.

### Rarity filter rule in the pilot

The pilot rarity filter is more specific than the plan text suggests. A feature is dropped only when both conditions are true:

- `pose_prevalence < rare_feature_pose_prevalence_lt`
- `n_conditions_with_feature < rare_feature_condition_prevalence_lt_n_conditions`

Because the code uses `and`, a feature survives if it is rare by pose prevalence but appears across enough conditions, or if it is concentrated in one condition but not below the condition-count cutoff.

### When formal pilot clustering is allowed

`summarize_condition_matrices()` currently marks formal clustering as allowed only when both conditions are met:

- `n_selected_poses >= minimum_clusterable_n`
- `main_feature_count > 0`

Otherwise the pilot records one of these statuses instead of running the method comparison on that condition:

- `insufficient_clusterable_signal`
- `empty_main_matrix`

## 4. Current Clustering Method Behavior

### Pilot-selected primary method for full analysis

Decision date: 2026-05-23.

The refreshed staged rerun and follow-up parameter-sensitivity review selected
one global primary clustering method for the full analysis:

- Primary method: `hdbscan_jaccard`
- Primary distance metric: binary Jaccard on contact-eligible IFP rows
- Selection evidence matrix: refreshed pilot `main_contact_eligible_ifp_matrix.csv`
- Primary `min_cluster_size`: `5`
- Primary `min_samples`: `null`
- Primary `cluster_selection_method`: `eom`
- Lenient HDBSCAN sensitivity setting:
  - `min_cluster_size = 3`
  - `min_samples = null`
  - `cluster_selection_method = eom`
- Orthogonal agglomerative sensitivity setting:
  - `distance_threshold = 0.55`
  - `min_cluster_size = 5`
  - `linkage = average`
- Conservative negative-control setting:
  - `min_cluster_size = 10`
  - `min_samples = null`
  - `cluster_selection_method = eom`

This decision is based on the refreshed staged pilot results under
`tests/tests_results/clustering_pilot_full_20260522_224240/` and the associated
parameter-sensitivity outputs under its
`clustering_parameter_sensitivity/` subtree.

The updated decision used the refreshed summary set and selected HDBSCAN
`min_cluster_size=5` because it reduced excessive subdivision relative to
HDBSCAN `min_cluster_size=3` while maintaining a lower noise burden than the
orthogonal agglomerative `0.55/min5` comparison. The decision criterion is now
an explicit balance between cluster recovery, cluster granularity, and noise,
not simply maximizing recovered cluster count.

Implementation note: some production code paths may still be wired to the older
agglomerative primary method. That is implementation backlog, not the current
governing Stage 6 decision.

### Production method

The governing production target is now `HdbscanJaccardClusterer` with:

- `min_cluster_size = 5`
- `min_samples = null`
- `metric = jaccard`
- `cluster_selection_method = eom`

If the standard production path is still instantiating
`AgglomerativeJaccardClusterer` in
`src/lpmo_pipeline/analysis/analysis_orchestrator.py`, that wiring is stale
relative to the governing documentation and should be treated as backlog.

The effective intended primary Stage 6 behavior is therefore:

- Distance metric is binary Jaccard.
- HDBSCAN clustering is fit on contact-eligible IFP rows with the selected
  `min_cluster_size=5` policy.
- Noise labels remain part of the reporting contract rather than being forced
  into retained clusters.
- Exact medoids are still selected by minimum summed within-cluster Jaccard
  distance for retained clusters.

Sensitivity comparison methods remain available as documented:

- HDBSCAN `min_cluster_size = 3`, `min_samples = null`
- Agglomerative Jaccard `distance_threshold = 0.55`, `min_cluster_size = 5`, `linkage = average`
- HDBSCAN `min_cluster_size = 10`, `min_samples = null`

### Current PoseBusters hard-QC behavior

The active PoseBusters runtime path is `run_posebusters_single()` in
`src/lpmo_pipeline/qc/posebusters_runner.py`.

- Combined AF3 exports are auto-split into ligand-only `mol_pred` and protein
  `mol_cond` inputs before PoseBusters is called.
- The runner uses built-in PoseBusters `dock` mode and does not inject a custom
  PoseBusters config.
- In the installed PoseBusters build, the active
  `posebusters.modules.intermolecular_distance.check_intermolecular_distance()`
  defaults are:
  - `radius_type = vdw`
  - `radius_scale = 1.0`
  - `clash_cutoff = 0.75`
  - `ignore_types = {hydrogens, organic_cofactors, inorganic_cofactors, waters}`
  - `max_distance = 5.0`
  - `search_distance = 6.0`
- `protein-ligand_maximum_distance` is the far-away check (`not_too_far_away`).
- `minimum_distance_to_protein` is the renamed `no_clashes` result, not the
  distance-cutoff field.
- If an enabled PoseBusters invocation raises an unexpected exception in hard
  QC, the orchestrator fails the pose closed with a
  `posebusters_runner_error` critical error instead of treating the missing
  result as pass.
- Current explicit hard-fail PoseBusters codes are:
  - `sanitization`
  - `all_atoms_connected`
  - `no_radicals`
  - `internal_steric_clash`
  - `bond_lengths`
  - `bond_angles`
  - `tetrahedral_chirality`
  - `protein-ligand_maximum_distance`
  - `minimum_distance_to_protein`
  - `volume_overlap_with_protein`
- Current explicit soft-warning PoseBusters codes are:
  - `aromatic_ring_flatness`
  - `double_bond_flatness`
  - `double_bond_stereochemistry`
  - `non-aromatic_ring_non-flatness`
  - `internal_energy`
  - `inchi_convertible`
- Any failing PoseBusters code outside those soft-warning names is still treated
  conservatively as critical by the runner.

### Current Privateer hard-QC runner-error behavior

Privateer is dispatched as a batch from `run_hard_qc()` after active-site
proximity and Cu-His geometry pass. If the batch backend raises unexpectedly,
or if it returns no result for an eligible pose, the orchestrator now creates an
explicit failing `PrivateerResult` (`privateer_batch_runner_error` or
`privateer_missing_result`). These are hard-fail verdicts in `qc_report.json`;
disabled Privateer (`run_privateer=False`) remains an intentional skip, not a
failure.

### Pilot agglomerative method

The comparison path in `src/lpmo_pipeline/analysis/clustering_agglomerative.py` behaves differently in a few important places:

- It uses precomputed Jaccard distance with linkage and distance threshold from pilot options.
- Clusters smaller than `min_cluster_size` are relabeled to noise after fitting.
- Non-noise labels are then renumbered sequentially.
- If the condition has a single pose or an all-zero Jaccard matrix, it is returned as one cluster.
- Medoids are still chosen by minimum summed within-cluster Jaccard distance, matching the HDBSCAN path.

### Condition summary metrics

`build_condition_cluster_summary()` in `src/lpmo_pipeline/analysis/clustering_hdbscan.py` currently uses these conventions:

- `n_ifp_clustered` means the number of rows that actually entered clustering, not all QC-pass poses.
- `n_noise` and `noise_fraction` are reported separately.
- `top_cluster_occupancy`, `cluster_entropy`, and `occupancy_gini` are computed over non-noise cluster sizes only.

## 5. Stage 16 Cluster Signature Logic

The current production Stage 16 path uses `build_cluster_signature_tables()` in `src/lpmo_pipeline/analysis/cluster_signatures.py`.

This is the active path wired by `src/lpmo_pipeline/analysis/analysis_orchestrator.py`.
The older `compute_cluster_signatures()` helper in the same file is not the function used to build the current Stage 16 outputs.

### Which clusters are included

Stage 16 membership is reconstructed from `cluster_assignments.tsv`-style rows.

- Rows with `cluster_id == -1` are excluded.
- Rows with `cluster_member_flag == False` are excluded.
- Every remaining non-noise cluster is annotated.

Current implementation note:

- There is no active `cluster_inclusion.min_occupancy` contract in
  `configs/thresholds.yaml`. The old plan-level key was removed from the active
  config because Stage 6 and Stage 16 do not read it.
- In practice, all non-noise clusters in `cluster_assignments.tsv` flow into
  Stage 16 regardless of occupancy.

### How occupancy is currently defined

The current denominator for cluster occupancy is not QC-pass poses and not all IFP-success poses.

- `total_clustered_by_condition` is counted from all rows present in `cluster_assignment_rows` for the condition.
- That count includes non-noise rows and noise rows that entered clustering.
- Therefore:

`occupancy = n_poses_in_cluster / n_rows_that_reached_clustering_for_that_condition`

This means Stage 16 occupancy is relative to the contact-eligible clustering input, including noise in the denominator.

### Cluster type classification

`_classify_cluster_type()` assigns `cluster_type` from occupancy plus geometry-plausibility fractions.

Only the geometry statuses below count as plausible in `_plausible_fraction()`:

- `geometry_plausible`
- `geometry_highly_plausible`

Current thresholds are loaded from `cluster_type` in `configs/thresholds.yaml`:

- `c1_compatible_min_plausible_fraction = 0.50`
- `c4_compatible_min_plausible_fraction = 0.50`
- `mixed_compatible_min_plausible_fraction = 0.40`
- `non_plausible_max_plausible_fraction = 0.10`
- `uncertain_max_occupancy_for_classification = 0.05`

The active decision order is:

1. `uncertain` if `occupancy <= uncertain_max_occupancy`
2. `non_plausible` if both C1 and C4 plausible fractions are `<= non_plausible_max`
3. `mixed_compatible` if:
   - both fractions are above their respective C1/C4 minima, or
   - the smaller of the two fractions is at least `mixed_min`
4. `C1_compatible` if C1 meets its minimum and is greater than C4
5. `C4_compatible` if C4 meets its minimum and is greater than C1
6. `uncertain` otherwise

### Current flat cluster-table export behavior

`write_cluster_table_tsv()` in `src/lpmo_pipeline/analysis/cluster_signatures.py`
currently writes `cluster_table.tsv` as a flat TSV export of the same
`cluster_summaries` rows that are also serialized into `cluster_signatures.json`.

- It is not a second aggregation pass; it reuses the already built Stage 16 summary rows.
- It includes only retained non-noise clusters because `cluster_summaries` are derived from non-noise membership.
- Conditions with no retained clusters therefore produce a header-only `cluster_table.tsv` rather than synthetic null rows.
- Current exported columns are:
  - identifiers: `condition_id`, `cluster_id`, `protein_id`, `ligand_id`, `medoid_pose_id`
  - summary labels: `cluster_type`, `n_poses`, `occupancy`, `c1_plausible_fraction`, `c4_plausible_fraction`
  - geometry summaries: median/IQR pairs for every field in `_GEOMETRY_STAT_FIELDS`
  - appended downstream fields: `construct_type`, `substrate_class`, `dp`, `cluster_size`, C1/C4 computable/plausible/highly-plausible fractions, pose-confidence mean/median/medoid fields, and convergence mean/median/medoid fields
- Stage 16b residue-importance code does not consume `cluster_table.tsv`; it still works from the richer Stage 16 signature/summary surfaces already passed in memory.

### Current Stage 16 IFP signature behavior

`_ifp_frequency_rows()` builds `cluster_ifp_signature.tsv` from `IFPResult` rows belonging to cluster members.

- Only member poses with `status == "ok"` and aligned `feature_names`/`flat_bitvector` are counted.
- `contact_frequency` is normalized by `n_ifp_poses`, which is the number of cluster members with usable IFP rows.
- The builder initializes counts for every feature present in the member pose layouts, not only for active bits.

As a result, the current implementation can emit explicit zero-frequency Stage 16 IFP rows for features that exist in the aligned cluster feature layout but are inactive across the cluster.

### Current Stage 16 residue signature behavior

`_residue_frequency_rows()` builds `cluster_residue_signature.tsv` from the Stage 13b residue-contact table.

Current grouping and aggregation rules are:

- Grouping key is `(residue_chain, residue_number, residue_name, condition_id)`.
- `contact_frequency` is `n_member_poses_with_contact / n_cluster_members`.
- `interaction_type` is the sorted comma-joined union of active interaction types observed for that residue inside the cluster.
- `ligand_residue_label` is the sorted comma-joined union of active ligand residue labels observed for that residue inside the cluster.
- Region flags are OR-reduced across member rows.

Because Stage 13b writes one row per feature bit, including zero bits, the current Stage 16 residue signature builder can also emit explicit zero-frequency rows. In those rows, `interaction_type` and `ligand_residue_label` fall back to `none` when the residue never has an active contact inside that cluster.

## 6. Stage 13b Residue-Contact Extraction Rules

`build_pose_residue_contact_rows()` in `src/lpmo_pipeline/analysis/residue_contact_extraction.py` currently follows a very direct rule:

- One emitted row per existing ProLIF feature bit.
- No additional residue aggregation happens at Stage 13b.
- `contact_present` is just the raw bit value for that exact feature.

This design is important because Stage 16 residue signatures are derived from this surface. The current pipeline is deliberately preserving the IFP feature layout at residue-contact level instead of rebuilding a separate residue abstraction first.

## 7. Stage 16b Residue Importance Rules

The current Stage 16b implementation is `compute_residue_importance_outputs()` in `src/lpmo_pipeline/analysis/residue_importance.py`.

### Protein x condition x residue score

The current aggregation is:

`residue_contact_score += cluster_occupancy * cluster_residue_contact_frequency`

Two additional weighted variants are then accumulated from the same residue row:

- `c1_weighted_residue_score += occupancy * residue_frequency * c1_plausible_fraction`
- `c4_weighted_residue_score += occupancy * residue_frequency * c4_plausible_fraction`

Current implementation details:

- `cluster_support_count` counts contributing cluster rows, not poses.
- `total_cluster_occupancy_with_contact` sums the occupancies of clusters contributing that residue row.
- `max_cluster_residue_frequency` stores the maximum within-cluster residue frequency seen for that residue.

### Zero-row backfill when no valid clusters exist

The current code explicitly preserves observed contact residues even when a condition has no valid Stage 16 clusters.

- `observed_residue_contact_rows` is used as a fallback source of residue identities.
- If a condition has observed contact residues but no retained cluster summaries, Stage 16b emits explicit zero-valued residue rows for those residues.

This is the current implementation reason that no-valid-cluster conditions stay visible downstream instead of collapsing to header-only files.

### Protein-level regio delta

`protein_residue_regio_delta.tsv` is currently built with these rules:

- `n_conditions_with_contact` counts how many condition-level residue rows exist for that protein-residue.
- `n_conditions_with_valid_clusters` is counted per protein from `condition_patch_summary.any_valid_cluster`.
- Mean C1 and C4 weighted scores are divided by `n_conditions_with_valid_clusters`, not by `n_conditions_with_contact`.

So the current protein-level delta averages over valid-cluster conditions only, while still reporting how often the residue was observed at condition level.

### Patch summaries

The current patch logic is split across two different weighted surfaces.

For `condition_patch_summary.tsv`:

- Feature-class fractions are computed from `cluster_ifp_signature_rows`.
- Each IFP feature contributes `occupancy * contact_frequency` mass.
- Only positive weighted mass contributes.
- Residue classes are currently defined in code as:
  - aromatic: `PHE`, `TRP`, `TYR`, `HIS`
  - polar: `SER`, `THR`, `ASN`, `GLN`, `CYS`
  - charged: `ASP`, `GLU`, `LYS`, `ARG`
- Hydrophobic interaction contribution is removed from the active patch summaries.
- H-bond contribution is interaction-driven and currently counts `ImplicitHBAcceptor` and `ImplicitHBDonor`.

For the region fractions in the same condition-level table:

- The code uses weighted residue contact mass from `cluster_residue_signature_rows`, not IFP feature mass.
- `catalytic_surface_contact_fraction` and `non_core_contact_fraction` are the active interpretation fields.
- `non_core_contact_fraction` collapses all available non-core signal into one bucket; downstream CBM analysis no longer distinguishes CBM versus linker contact attribution.
- Legacy `cbm_contact_fraction` and `linker_contact_fraction` may still be carried in intermediate tables for compatibility, but they are not the active biological interpretation contract.

### Protein patch summary

`protein_patch_summary.tsv` is currently a simple mean over all condition patch rows for a protein.

- Conditions without valid clusters remain in the average with zero-valued fractions.
- `n_conditions_with_valid_clusters` is tracked separately rather than being used as the averaging denominator.

## 8. Geometry Branch Decision Rules

The active downstream geometry row is built by `compute_pose_metrics_from_structure()` in `src/lpmo_pipeline/analysis/mdanalysis_metrics.py`.

### How the geometry row is currently built

- The geometry branch reuses `check_geometry()` from hard QC to identify Cu and the histidine brace.
- The downstream geometry branch is descriptive/flagging, not a hard QC gate. If `compute_pose_metrics_from_structure()` raises for a pose that already passed or was flagged by hard QC, the orchestrator keeps the pose in analysis export/IFP, writes a `geometry_not_computable` row, and records `geometry_metrics_error` in the pose manifest/case metadata.
- If Cu is not found, the function returns early and the pose keeps the default geometry state:
  - most geometry distances stay `None`
  - both `geometry_status_C1` and `geometry_status_C4` stay `geometry_not_computable`
- The proximal sugar is selected across the configured glycan chains by minimum heavy-atom distance to the original Cu position.
- `brace_integrity_flag` is `True` only when:
  - exactly three brace positions were selected, and
  - every Cu-His measurement is within `his_brace_max_search_a`
- The repositioned Cu used by the downstream geometry branch is the centroid of the selected brace positions.

### Virtual atom placement rules

- Virtual oxyl placement uses the vector from the N-terminal brace nitrogen to the repositioned Cu.
- That vector is projected into the brace plane and normalized.
- The virtual oxyl is then placed at `cu_oxyl_bond_length_a` from the repositioned Cu.
- Virtual H placement for C1 or C4 uses inferred bonded heavy neighbors.
- The neighbor search is not same-residue only:
  - same-residue heavy atoms within the local bond-length heuristic are included
  - heavy neighbors connected through `structure.connections` / `_struct_conn` are also included

That cross-residue connection handling is what currently lets glycosidic neighbors influence virtual H placement.

### Current geometry score and status logic

Thresholds are loaded from `geometry_plausibility` and `geometry_rules` in `configs/thresholds.yaml`.

For each target carbon separately, the active rules are:

- If the target atom is missing, the repositioned Cu is missing, the virtual oxyl is missing, or virtual H placement fails:
  - status is `geometry_not_computable`
- Otherwise:
  - `Cu_C1_distance` or `Cu_C4_distance` is measured from the repositioned Cu
  - `oxyl_H_*_distance` is measured from the virtual oxyl to the virtual H
  - `Cu_oxyl_H_*_angle` is measured as `repositioned_Cu -> virtual_oxyl -> virtual_H`, with the oxyl as the vertex

The older `attack_angle_*` / `oxyl-C-H` definition is not part of the current
output contract.

The oxyl-H score is currently:

`max(1.0 - abs(distance - oxyl_h_optimum_a) / oxyl_h_score_normalization, 0.0)`

but only when the distance is within `[oxyl_h_min_a, oxyl_h_max_a]`.
Outside that window the score is `NaN`.

The active status assignment order is:

1. `geometry_not_computable` if the distance is missing
2. `geometry_computable_implausible` if the score is `NaN`
3. `geometry_highly_plausible` if the distance is within `[oxyl_h_highly_plausible_min_a, oxyl_h_highly_plausible_max_a]`
4. `geometry_plausible` otherwise

### Sugar-face orientation

The current `sugar_face_orientation` is based on the sign of the dot product between:

- the sugar ring normal, and
- the vector from ring centroid to repositioned Cu

Current labels are:

- `ambiguous` when the value is near zero
- `ring_normal_toward_cu` when the dot product is positive
- `ring_normal_away_from_cu` when the dot product is negative

### Current geometry RMSD ownership

`core_rmsd_vs_reference` and `pocket_rmsd_vs_crystal` are present in the row contract, but this geometry module does not currently compute them.

## 9. Convergence Reference And Flag Rules

The current convergence path is implemented in `src/lpmo_pipeline/analysis/convergence_metrics.py`.

### Reference pose selection

`select_reference_pose()` currently applies this decision order:

- Prefer poses with `seed == 1`
- If no seed-1 pose exists, use the full pose set
- Within the remaining candidates, sort by:
  - highest `ranking_score`
  - then highest `mean_plddt`
  - then lowest `seed`
  - then lowest `sample`
  - then `pose_id` for deterministic tie-breaking

This is implemented by minimizing `_reference_sort_key()`, where ranking score and pLDDT are negated so that larger values win.

### Per-pose RMSD calculation

`compute_condition_convergence()` currently:

- requires exactly one `condition_id` per call
- aligns each pose to the reference with the configured protein alignment selection
- then computes ligand RMSD on the configured ligand selection after alignment

Default selections loaded from `configs/thresholds.yaml` are currently:

- alignment: `protein and name CA`
- ligand: `chainID B C D and not name H*`

If the protein alignment selection yields zero atoms, or the mobile and reference protein selections have different atom counts, convergence calculation errors out rather than silently degrading.

### Convergent flag and low-convergence flag

The current threshold rules are strict less-than comparisons:

- `convergent_flag = ligand_rmsd_to_reference < convergent_rmsd_max_a`
- `low_convergence_flag = convergence_fraction < low_convergence_flag_threshold`

With the current defaults from `configs/thresholds.yaml`:

- `convergent_rmsd_max_a = 2.0`
- `low_convergence_flag_threshold = 0.3`

### Condition-level summary

`build_condition_convergence_summary()` currently uses:

- `convergence_fraction` as the mean of boolean `convergent_flag`
- `median_ligand_rmsd` as the median of per-pose ligand RMSD values
- `iqr_ligand_rmsd` as `p75 - p25`

The summary also requires that every row agree on the same `condition_id` and the same `reference_pose_id`.

## 10. Crystal Anchoring Decision Rules

The active production crystal-anchoring path is the representative screen built from:

- `_choose_crystal_representatives()` in `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, and
- `run_crystal_reference_screen()` in `src/lpmo_pipeline/analysis/crystal_anchoring.py`

This is the current production path. The older `run_crystal_anchoring()` medoid-vs-crystal API still exists in the codebase, but it is not what `run_analysis_core()` calls.

### Representative pose selection in production

`_choose_crystal_representatives()` currently uses this decision order:

- If clustering produced retained medoids:
  - compare every retained cluster medoid
  - write `representative_role = cluster_medoid`
  - write `medoid_pose_id` only for true medoid rows
- If no retained medoids exist:
  - discover the top-level AF3 model CIF outside the `seed-*_sample-*` folders
  - prepare and hard-QC that model through the same normalization/PoseBusters/Privateer/geometry hard-QC path
  - run crystal anchoring only if the fallback verdict is not `dropped`
  - write `representative_role = af3_top_model_fallback` and leave `medoid_pose_id` blank

The top-level AF3 fallback is not included in generated-pose, QC, IFP, clustering, or medoid denominators.

### What gets written back to the production summaries

When crystal anchoring succeeds for a representative pose, the orchestrator currently reduces the comparison set as follows:

- `best_tanimoto` is the maximum non-null `ifp_tanimoto` across comparisons
- `best_pocket_rmsd` is the minimum non-null `pocket_rmsd` across comparisons

For medoid representatives, the medoid pose gets those values written back into the production metrics/case records. Fallback rows are reported only in crystal outputs and summary metadata because the top-level AF3 model is not part of the normal pose table.

In `crystal_anchor_table.tsv`, the current row mapping is:

- `representative_pose_id` = compared medoid or fallback pose
- `representative_role` = `cluster_medoid` or `af3_top_model_fallback`
- `medoid_pose_id` = true medoid pose ID, blank for fallback rows
- `contact_overlap_score` = `ifp_tanimoto`
- `same_binding_region_flag` = `pocket_rmsd_below_threshold`

So the table no longer overloads `medoid_pose_id` for non-medoid fallback representatives.

### Crystal reference loading and selection

`load_crystal_reference_records()` currently treats `input_data/pdb_structure_data.csv` as authoritative for protein-to-PDB membership.

- The columns `Oligo_Activity` and `Comment` are intentionally ignored.
- For each referenced PDB entry, the loader prefers:
  - `<PDB>_<protein_id>_ligand.cif`
  - and only falls back to `<PDB>_<protein_id>_no_ligand.cif` when the ligand file is missing
- If neither file exists on disk, the record is skipped.

### Selecting the protein copy and ligand site

`select_crystal_reference_site()` currently uses these rules:

- Preferred protein chain comes from the PDB field hint if present, otherwise chain `A`
- Ligand chains are resolved first from `_pdbx_branch_scheme` and mapped onto actual model chain IDs
- If branch-scheme mapping is unavailable, the code falls back to carbohydrate-like `_pdbx_nonpoly_scheme` chains
- Copper chains are resolved per selected protein site:
  - first from `_struct_conn` metal coordination rows involving the selected protein copy
  - then from `_pdbx_nonpoly_scheme` rows whose deposited strand/auth chain matches the selected protein copy
  - finally from nearest Cu distance as a fallback

Ligand-bound chain choice is then decided by heavy-atom distance:

- `preferred_chain_has_ligand` is `True` only if at least one ligand chain lies within `ligand_distance_cutoff_a` of the preferred protein chain
- If a ligand is expected and the preferred chain is not ligand-bound under that cutoff:
  - fallback candidates are all protein chains with at least one ligand chain within cutoff
  - the selected chain is the fallback candidate with the smallest ligand distance

The current default site-selection cutoff is `6.0 A`.

### Subset writing and chain remapping

`_write_selected_reference_cif()` does not keep the deposited assembly unchanged.
It writes a filtered subset mmCIF containing only:

- the selected protein chain
- the selected ligand chains
- the selected copper chains

The subset writer then remaps chain IDs onto the pipeline's default chain scheme:

- selected protein chain -> `A`
- selected ligand chains -> `B`, `C`, `D`, then later available chain IDs if needed
- selected site Cu -> `E`

The filtered file rewrites the relevant `_atom_site`, `_entity`, `_struct_asym`,
`_chem_comp`, `_chem_comp_bond`, `_pdbx_branch_scheme`,
`_pdbx_nonpoly_scheme`, `_struct_conn_type`, and `_struct_conn` loops to match
the remapped subset. Atom rows are filtered by role, not chain alone: protein
residues are kept only for the selected protein chain, carbohydrate residues
only for selected ligand chains, and Cu only for the selected site Cu. Water,
chloride, common buffer/solvent residues, and off-site external ligands are
excluded from the selected crystal subset.

### Prepared crystal states

`prepare_crystal_reference()` currently has an explicit no-ligand branch:

- If the selected crystal site has no ligand chains after selection:
  - return `status = "prepared_no_ligand"`
  - keep the subset/normalized reference
  - skip analysis export and crystal IFP generation
- If ligand chains are present:
  - normalize the selected subset with `NormalizeMMCIFRunner`
  - export non-protonated `analysis_export/complex_for_prolif.pdb`
  - export non-protonated `analysis_export/ligand_only_for_prolif.pdb`
  - compute a crystal IFP from those artifacts

The representative pose follows a similar normalize -> analysis_export -> IFP path before screening.

Crystal analysis export keeps Cu in `complex_for_prolif.pdb` for geometry and
pocket RMSD, but the ProLIF ligand export is glycan-only. If Cu or another
non-ligand residue is present in a copied glycan chain, `io/analysis_export.py`
skips it from `ligand_only_for_prolif.pdb` and records a warning and
skipped-residue count in `analysis_export_report.json`.

Current ProLIF loading reads the ligand directly from `ligand_only_for_prolif.pdb`,
so PDB residue names, numbers, and chain IDs define the ligand residue identity.
This keeps contacts from `BGC1.B` and `BGC1.C` distinct in the IFP feature space.

### IFP comparison rule in the current screen path

`run_crystal_reference_screen()` currently compares pose-vs-crystal IFPs with `compute_feature_aligned_tanimoto()`.

That means the active comparison rule is:

- align both bitvectors onto the union of feature names
- fill missing features with zero
- compute Tanimoto after union-feature alignment

So the current production crystal screen does not assume that the pose and crystal bitvectors already share the same raw column layout.

The prepared crystal IFP must also pass the same contact-eligibility rule used for clustering:

- at least 2 non-VdW interactions
- at least 1 non-VdW contact residue

There is no crystal-specific lenient override in the active screen path. A 1/1 crystal IFP remains non-comparable because one non-VdW contact alone is not treated as enough specific-contact evidence for IFP similarity scoring.

If a crystal IFP is `vdw_only`, `null_ifp`, or `low_specific_contact`, pocket RMSD and crystal geometry remain reportable, but `ifp_tanimoto` is left non-comparable and `crystal_ifp_exclusion_class` records the reason.

`crystal_ifp_diagnostic_summary.tsv` reports these eligibility/exclusion counts both by unique crystal reference and by medoid/fallback comparison row.

After the crystal-reference prep hardening, those exclusion classes should no
longer be attributed to known Cu contamination or missing crystal-side
normalization without further evidence. They may still occur because the
prepared deposited ligand genuinely has too little non-vdW ProLIF signal under
the shared contact-eligibility rule. The current working interpretation is that
crystal IFPs that remain empty or non-comparable are more likely to reflect
real biological-structural contact sparsity or non-specificity in the deposited
reference than a known bug in the crystal preparation or comparison path.

### Crystal-side C1/C4 geometry

Ligand-bound crystal references now run the downstream C1/C4 geometry branch after crystal subsetting and non-protonated analysis export. The geometry code uses the remapped ligand/Cu chains from the selected crystal subset and writes:

- full rows in `crystal_geometry_table.tsv`
- summary fields in `crystal_anchor_table.tsv`: `crystal_Cu_C1_distance`, `crystal_Cu_C4_distance`, `crystal_geometry_status_C1`, and `crystal_geometry_status_C4`

Apo references keep `prepared_no_ligand` status and do not produce C1/C4 crystal geometry.

### Pocket definition and pocket RMSD

`identify_pocket_residues_by_proximity()` currently defines the pocket as protein residues whose heavy atoms lie within the cutoff of:

- ligand heavy atoms, or
- Cu atoms

The current default cutoff in that function is `5.0 A`.

Pocket RMSD is then handled differently for ligand-bound and apo references.

For ligand-bound crystal references:

- identify crystal pocket residues from the prepared crystal complex
- map those crystal residue numbers onto the representative-pose sequence with `_map_residue_number_pairs_by_sequence()`
- invert the mapped pairs so RMSD is computed as representative residue number vs crystal residue number

For apo references:

- identify representative-pose pocket residues from the representative pose
- project those residue numbers onto the apo reference sequence with `_map_residue_number_pairs_by_sequence()`

### Sequence mapping and residue-name normalization

`_map_residue_number_pairs_by_sequence()` currently:

- collects protein residue identities from both structures
- normalizes modified residue names through the alias map before matching
  - for example `HIC -> HIS` and `MSE -> MET`
- if the normalized residue-name series are identical, maps indices one-to-one
- otherwise uses `SequenceMatcher` over the normalized residue-name sequences
- only retains mapped pairs whose normalized residue names still match exactly

### Pocket RMSD computation and threshold flag

`_compute_pocket_rmsd_from_residue_pairs()` currently:

- keeps only residue pairs where both sides have a C-alpha position
- returns `None` unless at least two shared residue pairs remain
- performs a Kabsch superposition on the paired C-alpha coordinates
- reports the resulting RMSD after removing rigid-body transforms

The current threshold flag is:

- `pocket_rmsd_below_threshold = pocket_rmsd < 2.5`

### Finalized interpretation policy (2026-05-23)

The current release freezes the crystal-anchoring interpretation contract as:

- IFP similarity is interpreted only when `ifp_comparison_eligible=True`.
- `ifp_tanimoto >= 0.50` is treated as moderate crystal-overlap support.
- `0.30 <= ifp_tanimoto < 0.50` is treated as weak/noisy support.
- `ifp_tanimoto < 0.30` is treated as low overlap.
- `local_pocket_rmsd < 2.5 A` remains a binding-region support flag and is not used as a standalone biological proof.

VdW-dominance / low-specific-contact interpretation is frozen to the
diagnostic summary in `crystal_ifp_diagnostic_summary.tsv`:

- if `contact_eligible_fraction < 0.50`, crystal-IFP interpretation is downgraded to geometry-first/IFP-limited,
- if `n_vdw_only > n_contact_eligible`, report the condition as VdW-dominated and avoid strong substrate-specific IFP claims.

### Finalized PoseBusters severity policy (2026-05-23)

The active hard-QC policy keeps PoseBusters dock defaults unchanged and freezes:

- `minimum_distance_to_protein` as hard-fail,
- explicit soft warning classes as soft-flag only,
- unknown PoseBusters failures as conservative hard-fail.

## 11. Current Code/Plan Clarifications And Remaining Gaps

These are the main places where the current implementation needs explicit clarification or still differs from a naive reading of the plan documents:

- Production Stage 6 now clusters on a condition-local main IFP matrix rather
  than the full raw contact-eligible matrix; broad real-data inspection of this
  production output is still pending.
- The active Stage 6/16 contract is that `cluster_inclusion.min_occupancy` is
  not used. All retained non-noise clusters are included in summaries and
  signatures; there is no pending occupancy gate to implement unless the code
  contract is changed explicitly later.
- Stage 16 documentation should follow `build_cluster_signature_tables()`, not the older `compute_cluster_signatures()` helper.
- Stage 13b and Stage 16 deliberately preserve explicit zero rows in some situations instead of reporting only positive contacts.
- `core_rmsd_vs_reference` and `pocket_rmsd_vs_crystal` are part of the geometry row contract, but the geometry module itself does not currently compute them.
- The current production crystal-anchoring path is all-medoid/fallback representative screening via `run_crystal_reference_screen()`.
  - Legacy constants and behavior in `run_crystal_anchoring()` are not the same thing as the active production path.

When updating plans or README text, these implementation-grounded rules should be treated as the current behavior until the code changes.

## 12. Optional Family Residue Enrichment Postprocess

The current within-family residue comparison is implemented as a standalone
optional postprocess, not as part of the required production `run` path.

### Execution boundary

- The production `run_analysis_core()` path still ends at the existing Stage 16b
  outputs plus summary/reporting layers.
- The optional family layer is entered only through the dedicated CLI command
  `lpmo-pipeline family-enrichment`, which calls
  `run_family_enrichment_postprocess()` in
  `src/lpmo_pipeline/analysis/family_enrichment_postprocess.py`.
- No analysis-core stage, summary table, or report depends on the family outputs.

### Input contract

The postprocess consumes these inputs directly:

- `protein_condition_residue_scores.tsv`
- `protein_residue_regio_delta.tsv`
- `input_data/metadata_final_ec_fixed.tsv`
- `input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta`

The current implementation does not read sequences from the production metadata
TSVs themselves. The catalytic-core FASTA is the canonical sequence source.

### Family and construct scope

- The active family scope is limited to `AA9` and `AA10`.
- The active construct scope is limited to `domain_only` rows, determined from
  the second `condition_id` token.
- Full-length rows are intentionally ignored in this first version because the
  family comparison is based on catalytic-core sequence alignment and assumes
  domain-local residue numbering.

### Catalytic-sequence mapping and deduplication

- FASTA headers are parsed from the `UniProtIDs|...|...` accession field in the
  deduplicated catalytic-core FASTA.
- Metadata rows are resolved against that FASTA using the primary UniProt ID,
  `CoAccessions`, and `Cat_CoAccessions`.
- The deduplication and aggregation unit is `Cat_Seq_Group` when available,
  with fallback to the protein identifier only if no sequence-group field is
  present.
- Family-level enrichment statistics are computed over unique catalytic
  sequence groups, while `family_aligned_residue_table.tsv` still keeps the
  per-protein provenance rows.

### Alignment source rules

- If `--alignment-dir` contains `<family>.aligned.fasta`, that precomputed
  alignment is used as the source of truth.
- Otherwise the postprocess attempts a MAFFT L-INS-i style alignment using
  `mafft --localpair --maxiterate 1000` on one representative sequence per
  catalytic sequence group.
- The alignment IDs are the catalytic sequence-group identifiers, not raw
  protein IDs.

### Skip and failure semantics

This feature is explicitly soft-failing. A family is skipped instead of failing
the overall postprocess when any of these conditions holds:

- no resolved AA9/AA10 family members after metadata and construct filtering
- fewer than 2 catalytic sequence groups remain in the family
- no precomputed alignment is available and `mafft` is missing
- MAFFT returns a non-zero exit code
- the aligned FASTA omits one of the expected sequence-group IDs
- the ungapped aligned sequence does not match the canonical catalytic-core
  sequence for that group
- no Stage 16b residues can be mapped to alignment columns

The postprocess summary JSON records processed versus skipped families and the
reason for each skipped family.

### Output contract

The postprocess writes under `08_family_residue_enrichment/`:

- `family_aligned_residue_table.tsv`
- `family_residue_enrichment.tsv`
- `family_alignment_manifest.tsv`
- `family_enrichment_summary.json`

`family_aligned_residue_table.tsv` now includes `family_aggregation_id`,
`construct_type`, `catalytic_core_residue_index`, and `alignment_residue` as
explicit provenance fields. `family_residue_enrichment.tsv` keeps raw protein
counts and also adds catalytic sequence-group counts so duplicate catalytic
cores are not overcounted in the summary means.
