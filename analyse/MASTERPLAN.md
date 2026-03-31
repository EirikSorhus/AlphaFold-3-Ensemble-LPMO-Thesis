# LPMO Structure Prediction -> Analysis Pipeline - Integrated Master Plan (v2.2)

## 0. Priority Order

This file follows the active priority stack for implementation decisions:
1. Latest user comments in the working thread
2. `plan_implementation_spec.txt`
3. `plan_analyse.txt` (legacy reference only)

**⚠️ CRITICAL: Large changes to order of analyses, which analyses are run, or core pipeline structure must NOT be made by AI without explicit user approval first. Any proposed changes in these areas require user review and confirmation.**

## 1. Scope And Invariants

- Models: AF3 (**RF3/RosettaFold 3 og Boltz-2 er ekskludert fra analyse**)
- AF3 fixed parameters (main analysis): `num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5`
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
- PLACER is mandatory before main QC and analysis.
- Atom names are not assumed consistent across models; mapping key is (element, CCD, local bond graph, 3D proximity).
- Chain schema: protein=A, glycans=B..D, metal=E.
- `_chem_comp_bond` must be complete for all `comp_id`, and `_struct_conn` must include glycosidic + Cu coordination links.
- Preserve all computed numeric metrics; do not drop distance/angle/support fields from output tables.
- Geometric planarity criteria are still not operationalized; final thresholds remain an open decision.
- Statistical/descriptive analysis should prefer R where possible (Python wrappers can orchestrate).

## 2. Required Output Tables

| Artifact | Format | Primary unit | Notes |
|---|---|---|---|
| pose_table.tsv | TSV | pose | Includes all QC, IFP, and geometry calculations |
| cluster_table.tsv | TSV | cluster | Main descriptive unit |
| predictive_cluster_table.tsv | TSV | cluster | Main predictive modeling table |
| enzyme_summary_table.tsv | TSV | enzyme | Secondary sensitivity table only |
| crystal_anchor_table.tsv | TSV | cluster/enzyme-ligand | Optional sanity check |
| cbm_comparison_table.tsv | TSV | enzyme-subanalysis | Domain-only vs full-length paired analyses |
| qc_report.json | JSON | pose | Includes pre-QC proximity + PB + Privateer + geometry gates |
| clusters.json | JSON | cluster | Occupancy, support, medoid, outlier stats |
| cluster_signatures.json | JSON | cluster | Median/IQR geometry + IFP enrichments |
| metrics.csv | CSV | mixed | Flat export for plotting and QA |
| summary.json | JSON | run | Dataset and gate summaries |
| run_manifest.json | JSON | run | Hashes, tool versions, seeds, gate outcomes |

## 3. Core Pipeline (Per Sub-Analysis)

### Step 1 - Ingest And Normalize
- Input: precomputed predictions from AF3 only (RF3 and Boltz-2 excluded).
- Validate mmCIF parse, categories, atom count.
- Normalize chains and connectivity.
- Required gate: atom mapping coverage = 100%.

### Step 2 - PLACER Refinement (Mandatory)
- Input: normalized structures.
- Output: PLACER ensemble + ranking metadata.
- Required gate: at least one refined pose.

### Step 3 - Pre-QC Active-Site Proximity (Before PoseBusters/Privateer)
- Goal: ensure ligand is close enough to active site before expensive chemistry checks.
- Compute and store at pose level:
  - nearest ligand atom to Cu
  - minimum Cu-ligand distance
  - minimum Cu-C1 and Cu-C4 where identifiable
- Hard pre-gate: `min_cu_ligand_distance <= 10.0 A`.
- Keep all computed distances in `pose_table.tsv` even for dropped poses (with status metadata).

### Step 4 - Hard QC
- Tools: PoseBusters, Privateer, Cu-His geometry gate.
- Hard fail examples:
  - critical PoseBusters error
  - severe Privateer fail (recognition/anomer)
  - Cu-His outside 1.9-2.6 A
- Soft flags are retained as metadata.

### Step 5 - Branch A: IFP
- Use post-PLACER structures.
- Do not reposition Cu and do not add virtual oxyl/H in IFP branch.
- Generate binary ProLIF vectors with consistent feature set.

### Step 6 - Branch B: Geometry
- Identify histidine-brace atoms.
- Build brace plane.
- Reposition Cu (geometry branch only).
- Place virtual oxyl.
- Place virtual H for C1 and C4.
- Compute pose-level geometry metrics:
  - oxyl_H_C1_distance
  - oxyl_H_C4_distance
  - Cu_C1_distance
  - Cu_C4_distance
  - attack_angle_C1
  - attack_angle_C4
  - sugar_face_orientation
  - optional ring_normal_vs_brace_normal
- Suggested plausibility window for oxyl-H: 1.5-4.0 A, reference optimum around 2.1 A.

### Step 7 - Clustering (IFP Only)
- 7A within-model clustering by (enzyme, substrate_class, DP, model_platform).
- 7B cross-model clustering by (enzyme, substrate_class, DP).
- Input for clustering: IFP only.
- Output: occupancy, outlier_rate, cluster_size, medoid pose, support across model platforms.

### Step 8 - Cluster Annotation
- Attach geometry/QC/support features to each final cluster.
- Use median/IQR summaries for geometry.
- Do not rely only on medoid geometry for predictive features.

### Step 9 - Descriptive Analysis (Cluster-Level)
- Primary descriptive unit: cluster.
- Recommended outputs: occupancy profiles, cluster counts, geometry distributions, IFP enrichments, family and CBM stratification.
- Recommended methods: proportions, medians/IQR, odds ratios, Fisher/chi-square, FDR correction.

### Step 10 - Predictive Analysis (Cluster Rows)
- Primary table: `predictive_cluster_table.tsv`.
- Each row is one cluster with occupancy/support/QC/geometry/IFP features.
- Main target: regioselectivity (C1 vs C4, optional mixed-class extension).
- Grouped CV at enzyme level; all clusters from same enzyme stay in one fold.
- Do not collapse to one mean geometry per enzyme as the primary model input.

### Step 11 - Optional Enzyme-Level Summaries (Secondary)
- Build occupancy-weighted enzyme summaries only for sensitivity analyses.
- Keep this explicitly separate from primary predictive modeling.

### Step 12 - Optional Crystal Anchoring
- Use only as sanity check (not primary validation criterion).
- **Alignment method**: PyMOL `pair_fit` for optimal local superposition.
  - Align on: (1) histidine-brace residues near Cu, (2) residues near Cu-site, (3) surface residues involved in substrate recognition.
  - Finding substrate-recognition residues: prefer literature-based selection (best); fallback is proximity-based (all residues within a cutoff of ligand in predicted structure). Note: proximity-based selection can give different residues per prediction depending on ligand placement.
  - Be explicit that RMSD is measured after optimized local alignment (not global).
- Metrics: local pocket RMSD (post pair_fit), ligand RMSD or proximal sugar RMSD, IFP similarity where comparable.

### Step 13 - CBM Paired Analysis
- Compare catalytic-domain-only vs full-length for CBM enzymes.
- Keep clustering domain-focused; treat CBM contacts as annotation/secondary features.

### Step 14 - Reporting
- Produce `summary.json`, `metrics.csv`, `report.html`.

## 4. Activity Label Mapping From EC

Maintain explicit EC-to-activity mapping from metadata (implemented in `scripts/ec_activity_mapping.R`):

- `1.14.99.54` -> cellulose, C1-hydroxylating
- `1.14.99.56` -> cellulose, C4-dehydrogenating
- `1.14.99.53` -> chitin, C1/C4 mixed
- `1.14.99.55` -> starch, C1-hydroxylating
- `1.14.99.-` -> special check:
  - AA17 family: homogalacturonan, C4 oxidation
  - otherwise: xylan/other oxidative activity (provisional label)

## 5. Predictive Cluster Table Requirements

Minimum required columns:

- analysis_id
- enzyme_id
- family
- cbm_status
- substrate_class
- DP
- cluster_id
- occupancy
- qc_factor
- support_factor
- cluster_weight
- cross_model_support
- convergence_support
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
- experimental_ligand_specificity_label

Suggested weight definition:

- cluster_weight = occupancy * qc_factor * support_factor

## 6. Optional Post-Analysis Tuning

Tuning is optional and runs after baseline analysis as a comparative/sensitivity workflow.

- It must not block the first main analysis delivery.
- Only AF3 is used; no cross-model comparison in tuning.
- If tuning is executed later, tuned settings are applied in follow-up reruns and documented in the manifest.

## 7. Failure Policy

| Event | Hard? | Action | Keep numeric metrics? |
|---|---|---|---|
| atom_mapping < 100% | YES | skip | yes, store coverage + reason |
| pre-QC ligand too far from active site | YES | drop before PB/Privateer | yes, store distances |
| glycan not CCD-valid for Privateer prep | YES | skip/drop | yes |
| PLACER returns 0 poses | YES | skip | yes, store run metadata |
| critical PoseBusters error | YES | drop pose | yes |
| severe Privateer fail | YES | drop pose | yes |
| Cu-His outside 1.9-2.6 A | YES | drop pose | yes |
| soft PB warning | NO | keep + flag | yes |
| HDBSCAN outlier | NO | keep with outlier label | yes |
| low crystal similarity | NO | keep + flag | yes |

## 8. Reproducibility

Every run writes `run_manifest.json` including:

- pipeline_version, timestamp, git_commit, config_hash
- tool_versions (Gemmi, PLACER, PoseBusters, Privateer, MDAnalysis, ProLIF, HDBSCAN, R)
- seeds and input checksums (SHA-256)
- gate outcomes including pre-QC active-site proximity
- references to prediction artifacts reused by analysis

## 9. RQ To Output Mapping

| RQ | Primary analysis unit | Key metrics | Artifact |
|---|---|---|---|
| RQ1 C1/C4 regioselectivity | cluster | Cu-C1/C4, oxyl-H, attack angles, occupancy | cluster_table.tsv, predictive_cluster_table.tsv |
| RQ2 substrate specificity | cluster | substrate/DP-specific occupancy + IFP signatures | cluster_table.tsv, predictive_cluster_table.tsv |
| RQ3 CBM effect | paired enzyme-subanalysis | occupancy shifts, support, CBM proximity | cbm_comparison_table.tsv |
| RQ4 crystal anchoring | cluster/enzyme-ligand | pocket RMSD (optimized local alignment via PyMOL pair_fit), ligand/proximal RMSD, IFP similarity | crystal_anchor_table.tsv |

## 10. Implementation Status

Pipeline status has moved from pseudocode-only to step-by-step implementation.
Execution follows the IMPLEMENTATION_PLAYBOOK.md priority order.
Each step is implemented, tested, and verified before the next begins.

## 11. Third-Party Code

Adapted patterns from the following MIT-licensed projects (see ATTRIBUTION.md):
- **PoseBench** (BioinfoMachineLearning/PoseBench): CIF→PDB conversion via PDBFixer, PoseBusters Python API usage, mol_table pattern.
- **benchmarking-af3** (lyulab/benchmarking-af3): Pocket residue identification via proximity cutoff, residue mapping between structures, metrics aggregation patterns.

## 12. Alignment Strategy (Crystal Anchoring)

Alignment of predicted structures against crystal references uses PyMOL `pair_fit`:

1. **Alignment atoms**: Histidine-brace Cα/Nε2 + Cu-coordinating residues + substrate-recognition surface residues.
2. **Substrate-recognition residues**: Identified via literature (preferred) or proximity to ligand in predicted structure (fallback). Document which method was used per system.
3. **Measurement**: After optimal local superposition, report pocket RMSD, per-residue deviations, and ligand RMSD. Always label results as "optimized local alignment" to distinguish from global alignment.
4. **Note**: Only AF3 is used; cross-model differences are not applicable.
