# AF3-LPMO Analysis Pipeline: Detailed Practical Plan

Version: 1.1
Status: Recommended main-analysis design; clustering method selected from pilot
Scope: AF3-only pipeline for LPMO–oligosaccharide complexes with pose-level QC, IFP-based clustering, geometry annotation, residue-level interpretation, and limited exploratory prediction.

Clustering update: the 2026-05-20 pilot decision selected agglomerative Jaccard
as primary method (`distance_threshold=0.55`, `min_cluster_size=3`) with
agglomerative `0.45/min3` and HDBSCAN `min3` retained as sensitivity settings.

Runtime update 2026-05-19: production execution now passes `n_jobs` into independent per-pose preparation, hard-QC/Privateer dispatch, and ProLIF batch work. Routine successful normalization keeps only `normalize_report.json`; atom-map/rename debug files are written only for unexpected mapping or validation failures. CIF loop remapping in normalization now rewrites each affected mmCIF loop once instead of once per tag, which removes the previous dominant `gemmi_compat.py` bottleneck. The full clustering pilot should be launched through the staged Slurm-array wrapper, which splits domain-only and full-length selections into independent protein-level shards and runs them across multiple jobs/nodes. The legacy single-job wrapper is retained only as a stable fallback.

\---

## 1\. Purpose

This document specifies a practical, implementable, AI-friendly analysis pipeline for AF3-predicted LPMO–oligosaccharide complexes.

The pipeline is designed to answer five linked goals:

1. Describe which binding modes AF3 produces for each protein–ligand condition.
2. Quantify how stable or multimodal those binding modes are across seeds and samples.
3. Identify contact residues and contact patterns associated with ligand class, ligand length, and putative C1/C4 compatibility.
4. Summarize how geometry-based plausibility behaves across poses and clusters without using that geometry to pre-filter away most of the model output.
5. Evaluate what this kind of AF3 analysis can realistically be used for, including where ipTM and similar confidence metrics appear to have limited practical value.

This is not a strong-validation pipeline. It is a structured computational interpretation pipeline.

\---

## 2\. High-level design decisions

### 2.1 Locked decisions

* Main structure source: **AF3 only**.
* **Data availability**: Precomputed AF3 structures fully available at:
  - Domain-only: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core`
  - Full-length: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length`
* Main pose ensemble per protein–ligand condition: **15 seeds x 5 samples = 75 poses** (AF3 `num_diffusion_samples=5`, runs completed).
* Main clustering input: **ProLIF binary IFP only**.
* Primary clustering method: **locked from pilot** to agglomerative Jaccard clustering on contact-eligible IFP rows with `linkage=average`, `distance_threshold=0.55`, and `min_cluster_size=3`. HDBSCAN `min_cluster_size=3`, `min_samples=null`, `cluster_selection_method=eom` is retained as a sensitivity path.
* Main descriptive unit: **cluster**.
* Main biological repeated-measures unit: **protein**.
* Main inferential caution: **poses and seeds are not independent biological replicates**.
* Main QC philosophy: **remove clearly invalid structures early; do not hard-filter primarily on downstream mechanistic geometry before IFP**.
* Main geometry philosophy: **compute after IFP; use as annotation, plausibility, and interpretation**.
* Main medoid philosophy: **use for visualization, figure choice, manual inspection, and crystal sanity-check; do not use as the only analytical representation**.

### 2.2 Optional / weak-priority components

* Crystal comparison is a sanity-check, not a primary validation criterion.
* Predictive modeling is exploratory and should remain small and simple.

### 2.3 Current clustering decision state

The pipeline previously treated HDBSCAN on binary IFPs as the recommended
default. That was superseded by the clustering pilot and parameter-sensitivity
run completed on 2026-05-20.

Current decision:

* use one global primary method for all protein-ligand conditions, not a separate method per condition
* primary method: agglomerative Jaccard clustering on contact-eligible IFP rows
* primary parameters: `linkage=average`, `distance_threshold=0.55`, `min_cluster_size=3`
* pilot selection evidence used `main_contact_eligible_ifp_matrix.csv`
* conservative agglomerative sensitivity: `distance_threshold=0.45`, `min_cluster_size=3`
* HDBSCAN sensitivity: `min_cluster_size=3`, `min_samples=null`, `cluster_selection_method=eom`
* conditions below the contact-eligible/formal clustering threshold remain reported as insufficient signal rather than interpreted as binding-mode clusters

Pilot rationale:

* In the 94 formally clusterable pilot conditions, agglomerative `0.55/min3`
  recovered clusters in 84 conditions, versus 79 for `0.45/min3` and 53 for
  HDBSCAN `min3`.
* The selected primary setting reduced median noise fraction to 0.52, compared
  with 0.68 for the conservative `0.45/min3` setting.
* Mean non-noise cluster size was 4.77 poses (SD 2.80; median 4; range 3-21),
  which improved coverage without collapsing the pilot into a few large
  clusters.

\---

## 3\. Dependency structure and statistical hierarchy

### 3.1 Nested dependence

Dependence exists at several levels:

* sample within seed
* seed within protein–ligand condition
* ligand condition within protein
* constructs within protein (domain-only vs full-length for CBM cases)
* proteins within family / homologous groups
* identical proteins sharing the same MSA or near-identical sequence background

### 3.2 Consequences

* **Pose** is a technical sampling unit, not an independent biological replicate.
* **Cluster** is a compressed description of repeated pose output from one protein–ligand condition.
* **Protein–ligand condition** is the smallest unit that can reasonably summarize one biological condition.
* **Protein** is the natural grouping variable for grouped CV and repeated-measures comparisons across ligand conditions.
* Ligand type and DP should be treated as repeated conditions within protein, not as globally independent groups.

### 3.3 Practical rule

Never report pose count as effective biological n.

\---

## 4\. Main questions the pipeline must support

### 4.1 Binding mode questions

* Does one protein–ligand condition yield one dominant IFP-defined mode or several competing modes?
* Are modes stable across seeds and samples?
* Are some ligand classes or DP values associated with more multimodal or unstable behavior?

### 4.2 Regioselectivity questions

* Do QC-passing clusters show geometry more compatible with C1 or C4 attack?
* Does the same protein shift between C1-like and C4-like geometry across ligand classes or DP values?
* Is there at least one robust plausible C1-like or C4-like cluster per protein–ligand condition?

### 4.3 Residue interpretation questions

* Which residues recurrently contact ligand in stable clusters?
* Which residues are enriched in conditions associated with C1-like vs C4-like geometry?
* Which residue classes or patch types are enriched for chitin vs cellulose vs starch conditions?
* For CBM-containing proteins, which CBM or linker residues gain contact occupancy in full-length predictions?

### 4.4 Confidence / model-utility questions

* How informative is AF3-reported ipTM relative to actual pose plausibility, clustering stability, or cluster geometry?
* Does high ipTM identify biologically useful modes, or does it mainly vary weakly / unhelpfully in this setting?

\---

## 5\. Scope of sub-analyses

Run each ligand condition as an independent sub-analysis in IFP space.

Recommended primary set:

* chitin\_DP4
* chitin\_DP6
* chitin\_DP8
* cellulose\_DP4
* cellulose\_DP6
* cellulose\_DP8
* starch\_DP4
* starch\_DP6
* starch\_DP8

Important rule:

* Do not cluster different substrate classes together.
* Do not cluster different DP values together.
* Cross-condition comparisons happen **after** each condition has been summarized.

\---

## 6\. Required inputs

### 6.1 Protein-level metadata

Required columns in `protein\_metadata.tsv`:

* `protein\_id`
* `source\_db`
* `family\_label`
* `organism`
* `sequence\_full\_length`
* `sequence\_catalytic\_domain`
* `has\_cbm`
* `cbm\_type`
* `construct\_type` (`domain\_only`, `full\_length`)
* `experimental\_regio\_label` (`C1`, `C4`, `mixed`, `unknown`)
* `experimental\_substrate\_label` (`chitin`, `cellulose`, `starch`, `mixed`, `unknown`)
* `experimental\_oligo\_activity\_label`
* `has\_crystal\_reference`
* `crystal\_reference\_path`
* `notes`

### 6.2 Ligand-level metadata

Required columns in `ligand\_metadata.tsv`:

* `ligand\_id`
* `substrate\_class`
* `dp`
* `glycan\_name`
* `input\_cif\_path`
* `proximal\_residue\_definition`
* `c1\_atom\_name`
* `c4\_atom\_name`
* `ring\_atom\_names`
* `reducing\_end\_definition`
* `notes`

### 6.3 Run manifest

Required columns in `run\_manifest.tsv`:

* `analysis\_id`
* `protein\_id`
* `construct\_type`
* `ligand\_id`
* `substrate\_class`
* `dp`
* `seed`
* `sample\_index`
* `expected\_output\_dir`
* `run\_status`

### 6.4 Static reference definitions

Store in versioned config or YAML:

* Cu atom naming conventions
* his-brace identification rules
* ligand proximal sugar identification rules
* allowed residue renaming / atom normalization rules
* QC thresholds
* clustering parameters
* convergence thresholds
* geometry plausibility thresholds
* residue importance thresholds

\---

## 7\. Directory structure

Recommended project root, but is to be adapted to current structure:

```text
project\_root/
  config/
    config.yaml
    thresholds.yaml
    prolif\_features.yaml
    geometry\_rules.yaml
    residue\_rules.yaml
  metadata/
    protein\_metadata.tsv
    ligand\_metadata.tsv
    run\_manifest.tsv
  01\_inputs/
    proteins/
    ligands/
    af3\_json/
  02\_af3\_raw/
    {analysis\_id}/
      seed\_{seed}/
        sample\_{sample}/
  03\_parsed/
    pose\_manifest.tsv
    pose\_confidence.tsv
    structure\_index.tsv
  04\_qc/
    posebusters/
    privateer/
    qc\_tables/
  05\_ifp/
    fingerprints/
    residue\_contacts/
  06\_geometry/
    pose\_geometry.tsv
    geometry\_debug/
  07\_clustering/
    distance\_matrices/
    cluster\_assignments.tsv
    medoid\_manifest.tsv
  08\_cluster\_annotation/
    cluster\_table.tsv
    residue\_importance\_table.tsv
  09\_descriptive/
    plots/
    summary\_tables/
  10\_predictive/
    modeling\_tables/
    cv\_results/
  11\_crystal/
    crystal\_anchor\_table.tsv
    crystal\_geometry\_table.tsv
    crystal\_ifp\_diagnostic\_summary.tsv
  12\_reports/
    report\_ready\_tables/
    figure\_manifest.tsv
```

\---

## 8\. Entity identifiers

Every object must have a stable ID.

### 8.1 Pose ID

```text
{protein\_id}\_\_{construct\_type}\_\_{substrate\_class}\_DP{dp}\_\_seed{seed}\_\_sample{sample}
```

### 8.2 Condition ID

```text
{protein\_id}\_\_{construct\_type}\_\_{substrate\_class}\_DP{dp}
```

### 8.3 Cluster ID

```text
{condition\_id}\_\_cluster{n}
```

### 8.4 Medoid ID

Same as pose ID of the selected medoid.

\---

## 9\. Stage 1: AF3 generation and parsing

### 9.1 Inputs

For each run:

* protein construct sequence
* ligand CIF/mmCIF
* Cu included explicitly if present in the modeling setup
* AF3 JSON input
* seed definition
* sample control / parsing metadata

### 9.2 Outputs to parse from AF3

Per pose, store at least:

* `pose\_id`
* `protein\_id`
* `construct\_type`
* `substrate\_class`
* `dp`
* `seed`
* `sample\_index`
* `raw\_structure\_path`
* `confidence\_json\_path`
* `ranking\_score` if present
* `iptm` if AF3 reports it for that run
* `ptm` if present
* `mean\_plddt`
* `ligand\_interface\_confidence` if present
* `pae\_summary` if present
* `run\_status`

### 9.3 ipTM handling

ipTM must be retained as a descriptive metric, but not treated as a main validity criterion.

Rules:

* If AF3 reports `ipTM`, store it directly.
* If `ipTM` is absent for a given pose or configuration, store `NA`; do not derive an unofficial substitute.
* Compute and report:

  * cluster-level `iptm\_mean`, `iptm\_sd`, `iptm\_median`, `iptm\_iqr`
  * protein–ligand-condition `iptm\_mean`, `iptm\_sd`, `iptm\_median`, `iptm\_iqr`
  * protein-level `iptm\_mean`, `iptm\_sd`, `iptm\_median`, `iptm\_iqr`
* Also report the same metrics separately for:

  * all generated poses
  * QC-passing poses only
  * geometry-computable poses only

Purpose:

* Describe whether ipTM tracks anything useful.
* Show explicitly that ipTM may have limited practical interpretability in this setting.

### 9.4 Output tables after Stage 1

#### `pose\_manifest.tsv`

One row per pose.

Core columns:

* pose metadata
* structure paths
* AF3 confidence fields
* parsing success/fail fields

#### `pose\_confidence.tsv`

One row per pose with only AF3 confidence metrics.

Purpose:

* keeps downstream plotting simple
* lets ipTM be analyzed independently of structural parsing

\---

## 10\. Stage 2: Hard QC before IFP

### 10.1 Goal

Remove only poses that are clearly unusable for any meaningful structural analysis.

### 10.2 Hard-fail criteria

A pose is excluded before IFP if any of the following applies:

* severe PoseBusters fail
* severe Privateer fail
* ligand missing / broken / unparsable
* Cu missing
* protein–ligand structure cannot be parsed consistently
* atom naming prevents residue or ligand mapping
* ligand is far from the relevant binding region
* structure file corrupt or incomplete

Contact-poor poses are not automatically hard QC failures. Null IFP, van der Waals-only contact, and low-specific-contact cases are classified after IFP/contact eligibility so that condition-level summaries can distinguish structural failure from lack of binding-mode signal.

### 10.3 Soft QC flags

These do not exclude a pose before IFP, but must be recorded:

* mild stereochemical issue
* non-fatal clash issue
* uncertain sugar atom mapping
* unusual but parsable Cu environment
* very low local confidence near binding surface
* uncertain histidine brace identification candidate set

### 10.4 Required QC output fields

Add to `pose\_table.tsv`:

* `posebusters\_status`
* `posebusters\_fail\_codes`
* `privateer\_status`
* `privateer\_fail\_codes`
* `qc\_hard\_fail`
* `qc\_soft\_flag\_count`
* `qc\_soft\_flag\_codes`
* `qc\_stage1\_pass`

### 10.5 QC attrition table

Build `qc\_attrition\_table.tsv` with one row per protein–ligand condition:

* `n\_generated`
* `n\_stage1\_hard\_fail`
* `n\_stage1\_pass`
* `hard\_fail\_rate`
* fail counts by reason code

Purpose:

* quantify uneven attrition
* show where selection bias enters
* support later interpretation of missingness patterns

\---

## 11\. Stage 3: IFP generation

Documentation is available at: https://prolif.readthedocs.io/en/stable/source/modules/interaction-fingerprint.html

### 11.1 Input set

Use only poses with `qc\_stage1\_pass = True`.

### 11.2 Important rule

Use the AF3 structure as-is for IFP generation.

Do not:

* reposition Cu
* add virtual oxyl
* add virtual H on C1/C4

### 11.3 IFP settings

For each pose:

* make a copy of pose
* add hydrogens to copy if needed by the IFP engine
* use the same residue and ligand selection masks for every pose in the same sub-analysis
* generate binary IFP over a fixed set of interaction types
* preserve each monosaccharide as a separate ligand residue in the feature space

The broad audit interaction set, matching the current `configs/prolif_features.yaml` validation surface, includes:

* Hbond donor
* Hbond acceptor
* hydrophobic
* aromatic / stacking
* anionic
* cationic
* cation-pi
* pi-cation
* van der Waals contact

This broad set is retained for audit and descriptive tables. The pilot-defined
main clustering feature set is narrower unless the feature audit justifies
additional interaction types:

* default include for main clustering: Hbond donor, Hbond acceptor, aromatic / stacking
* conditional include: hydrophobic if it is residue- and ligand-unit-specific; cation-pi if it occurs in enough conditions to be informative
* default exclude for main clustering: van der Waals / close-contact features

Feature identity for clustering is:

```text
protein_residue_id + interaction_type + ligand_unit_id
```

The current flattened implementation represents the same concept as
`ligand_residue|protein_residue|interaction`. The important invariant is that
monosaccharide/ligand-unit identity is preserved. Clustering must not be based
only on total contact counts.

The pilot feature audit must compute prevalence by interaction type and by
individual feature before model comparison. Extremely rare features are removed
from the main clustering matrix only if both conditions hold:

* pose prevalence < 0.01
* feature is present in fewer than 2 pilot conditions

Rare features remain in raw contact/IFP tables. Common polar or aromatic
features are not removed just because they are common; they are removed only if
they are uninformative for binding-mode separation.

### 11.4 Per-pose outputs

`pose\_ifp\_table.tsv`:

* `pose\_id`
* `condition\_id`
* `ifp\_generation\_status`
* `ifp\_vector`
* `ifp\_feature\_names`
* feature names use `ligand_residue|protein_residue|interaction`
* `n\_total\_contacts`
* `n\_hbond\_donor`
* `n\_hbond\_acceptor`
* `n\_hydrophobic`
* `n\_aromatic`

Current standalone implementation note (validated 2026-05-04):

* also writes `n\_anionic`, `n\_cationic`, `n\_cation\_pi`, `n\_pi\_cation`, `n\_vdw\_contact`
* writes `ifp\_interaction\_counts` as a JSON map of all active interaction types
* treats collapsed `UNL`-style ligand labels as invalid for the monosaccharide-resolved branch

### 11.5 Residue-contact extraction

In parallel with IFP, build a residue-level contact table.

`pose\_residue\_contact\_table.tsv`:

One row per `pose\_id x residue\_id x interaction\_type` with:

* `pose\_id`
* `protein\_id`
* `condition\_id`
* `residue\_chain`
* `residue\_number`
* `residue\_name`
* `interaction\_type`
* `contact\_present`
* `ligand\_residue\_label`
* `distance\_if\_available`
* `is\_catalytic\_surface\_region`
* `is\_cbm\_region`
* `is\_linker\_region`

This table is essential for later residue interpretation.

Current status 2026-05-04: `pose\_residue\_contact\_table.tsv` remains planned;
the standalone ProLIF implementation has not yet emitted this table.

\---

## 12\. Stage 4: Geometry branch

### 12.1 Input set

Use the same `qc\_stage1\_pass` poses that went into IFP.

### 12.2 Philosophy

Geometry is computed after IFP and cluster assignment will later be linked to geometry. Geometry should not drive the main pre-IFP exclusion unless the pose is fundamentally uncomputable.

### 12.3 Steps

#### Step 4A. Histidine brace identification

Rule: identify the 3 copper-coordinating nitrogen atoms by proximity (atom-first, not residue-first).

Locked atom-selection constraints for LPMO His-brace:

* The selected 3 coordinating nitrogens must include `His1:N` and `His1:ND1`.
* `His1:NE2` is never allowed in the selected 3 coordinating nitrogens.
* The third coordinating nitrogen must come from a histidine with residue number other than 1.
* Cu-His hard-gate evaluation uses only these selected 3 nitrogens.

Procedure:

1. Require `His1:N` and `His1:ND1` from chain A residue 1.
2. Enumerate candidate histidine nitrogens in chain A excluding `His1:NE2`.
3. Compute Euclidean distance from each candidate N atom to the Cu position (chain E).
4. Select the closest valid third nitrogen from histidine residue != 1 within `geometry_rules.yaml: his_brace.max_search_dist_a`.
5. Build the coordinating set as exactly three atoms: `His1:N`, `His1:ND1`, and the selected third nitrogen.
6. If the third nitrogen is missing within cutoff, record `brace_identification_status = partial`; if mandatory His1 atoms are missing, `failed`.
7. Use `His1:N` as the N-terminal nitrogen for oxyl direction.

Brace integrity is **not a hard QC gate** in Stage 4. Poses where the brace is partially distorted continue to geometry computation with reduced `brace_confidence`. A pose is blocked from geometry computation only if brace identification fully fails.

Store:

* `brace_identification_status` (`ok` / `partial` / `failed`)
* `brace_atom_ids` — the 3 identified N atom positions
* `brace_n_found` — count of N atoms found within cutoff
* `brace_confidence` (`high` / `medium` / `low`)
* `brace_integrity_flag` — soft flag: True if all 3 N atoms found within cutoff

#### Step 4B. Brace plane construction

The brace plane is defined by the 3 N atoms identified in Step 4A.

Procedure: compute the unit normal vector to the plane spanned by the 3 N atom positions.

Store:

* `brace_plane_status`
* `brace_normal_vector`

#### Step 4C. Cu standardization for geometry only

Store:

* `cu\_reposition\_status`
* `repositioned\_cu\_coordinates`

#### Step 4D. Virtual oxyl placement

The oxyl intermediate occupies the equatorial position trans to the N-terminal amine in the square-planar Cu coordination environment.

Procedure:

1. Identify the N-terminal N atom among the 3 brace N atoms (Step 4A): the one belonging to residue 1, chain A.
2. Compute the unit vector from N-terminal N toward the repositioned Cu.
3. Project this vector onto the brace plane (Step 4B).
4. Place the virtual oxyl O at distance `cu_oxyl_bond_length_a` from repositioned Cu, along this projected vector in the direction **away** from the N-terminal N (i.e., trans across Cu).

All distances are stored in `geometry_rules.yaml: virtual_oxyl`.

Store:

* `oxyl_placement_status`
* `virtual_oxyl_coordinates`

#### Step 4E. Virtual H placement on C1 and C4

The abstractable hydrogen is placed opposite the heavy-atom neighbors of C1 (or C4), representing the H that would be abstracted in the LPMO reaction.

Procedure (same for C1 and C4):

1. For the target carbon (C1 or C4 of the proximal sugar), retrieve bonded heavy atoms from the CCD topology.
2. Compute the centroid of all bonded heavy atoms.
3. Place H along the vector centroid → C (pointing away from neighbors), at distance `c_h_bond_length_a` from C (see `geometry_rules.yaml: virtual_h.c_h_bond_length_a`).

Store:

* `c1_identification_status`
* `c4_identification_status`
* `virtual_H_C1_coordinates`
* `virtual_H_C4_coordinates`

#### Step 4F. Geometry metrics

Required fields:

* `oxyl\_H\_C1\_distance`
* `oxyl\_H\_C4\_distance`
* `Cu\_C1\_distance`
* `Cu\_C4\_distance`
* `attack\_angle\_C1`
* `attack\_angle\_C4`
* `sugar\_face\_orientation`
* `ring\_normal\_vs\_brace\_normal`
* `proximal\_sugar\_id`

### 12.4 Geometry scoring and status

#### Oxyl-H continuous score

For each target (C1 and C4), compute a continuous plausibility score from the oxyl–H distance. All threshold values live in `thresholds.yaml: geometry_plausibility`.

```
d = euclidean_distance(virtual_oxyl, virtual_H_Cx)

if d < oxyl_h_min_a or d > oxyl_h_max_a:
    score = NaN   # outside physically meaningful range
else:
    score = max(1.0 - abs(d - oxyl_h_optimum_a) / oxyl_h_score_normalization, 0.0)
```

The gate at `oxyl_h_min_a`–`oxyl_h_max_a` rejects poses where the substrate is either clashing with the active site or too far away for HAT. The optimum at `oxyl_h_optimum_a` (~2.1 Å) is derived from DFT studies of Cu-oxyl HAT barriers. Score approaches 1.0 near the optimum and 0.0 at the window boundaries.

Store per pose:

* `oxyl_H_score_C1` — continuous score [0, 1] or NaN
* `oxyl_H_score_C4` — continuous score [0, 1] or NaN

#### Status label assignment

Assign one status per target (C1 and C4), stored as `geometry_status_C1` and `geometry_status_C4`:

* `geometry_not_computable` — brace, Cu, C1, or C4 identification failed; no score computed
* `geometry_computable_implausible` — computed but oxyl–H distance outside plausibility window; score = NaN
* `geometry_plausible` — within the plausibility window; score > 0
* `geometry_highly_plausible` — within the tighter window around the optimum (see `thresholds.yaml: geometry_plausibility.oxyl_h_highly_plausible_*`)

Thresholds are operational, not mechanistic truth. They must be fixed in `thresholds.yaml` before the full analysis run.

#### Brace integrity (soft, descriptive only)

`brace_integrity_flag` is stored as a soft QC column in `pose_geometry.tsv`. It does not affect geometry status or cause pose exclusion. It is used in:

* descriptive summaries of how often AF3 produces intact brace geometry
* as an optional feature for exploratory predictive analysis

### 12.5 Geometry output table

`pose\_geometry.tsv`

One row per pose with all geometry fields and status labels.

### 12.6 How geometry feeds the next stage

Geometry is not used to define clusters.

Geometry is used later to:

* annotate clusters
* define cluster-level plausibility fractions
* identify C1-like vs C4-like compatible clusters
* generate descriptive plots
* create exploratory predictive features

\---

## 13\. Stage 5: Convergence / internal reproducibility metrics

### 13.1 Goal

Quantify whether AF3 repeatedly generates similar ligand placements for the same protein–ligand condition.

### 13.2 Method

Per protein–ligand condition:

1. **Superimpose** all QC-passing poses by protein Cα atoms onto a fixed reference pose (see alignment reference rule below).
2. **Compute ligand RMSD** of each pose versus the superimposed reference: Euclidean RMSD over all ligand heavy atoms.
3. **Flag each pose**: `convergent_flag = True` if ligand RMSD < `convergent_rmsd_max_a` (see `thresholds.yaml: convergence.convergent_rmsd_max_a`).
4. **Aggregate per condition**: `convergence_fraction = n_convergent_poses / n_qc_pass_poses`.

Alignment reference rule:

* **Per-pose computation**: use the top-ranked QC-passing pose within the condition (by `ranking_score`, or `mean_plddt` if ranking score is absent).
* **Reporting context**: cluster medoids are reported separately for clustering, structural figures, and crystal anchoring; convergence remains a descriptive per-condition metric and does not re-filter poses.

**Convergence metrics are descriptive only.** They are not used as a filter to drop poses. All QC-passing poses continue to IFP and geometry regardless of their convergence flag.

### 13.3 Output fields

Per pose:

* `ligand\_rmsd\_to\_reference`
* `convergent\_flag` — True if RMSD < `convergent_rmsd_max_a`

Per condition:

* `convergence\_fraction` — fraction of poses flagged as convergent
* `median\_ligand\_rmsd`
* `iqr\_ligand\_rmsd`
* `low\_convergence\_flag` — True if `convergence_fraction` < `low_convergence_flag_threshold` (see `thresholds.yaml`)

Purpose:

* Robustness meta-metric: a `convergence_fraction` near 1.0 means AF3 consistently places the ligand in the same orientation; geometric scores from the top-ranked pose are then more representative.
* Values below the flag threshold suggest the binding mode is ambiguous and geometric scores from any single pose should be interpreted cautiously.
* Not a direct measure of mechanistic correctness.

\---

## 14\. Stage 6: Clustering

### 14.1 Input set

Per condition:

* all hard-QC-passing poses with successful IFP generation
* contact-eligible poses only for formal binding-mode clustering

### 14.2 AF3-only simplification

Because the main pipeline is AF3-only:

* remove within-model vs cross-model clustering separation
* perform **one clustering per protein–ligand condition**

### 14.3 Contact eligibility before formal clustering

Contact eligibility is applied after IFP generation and decides which poses enter binding-mode clustering. It is not a replacement for hard QC.

Main rule for the pilot:

* `min_non_vdw_interactions >= 2`
* `min_non_vdw_contact_residues >= 1`

Sensitivity rules:

* lenient: `min_non_vdw_interactions >= 1` and `min_non_vdw_contact_residues >= 1`
* strict: `min_non_vdw_interactions >= 2` and `min_non_vdw_contact_residues >= 2`

Excluded poses must still be reported in condition summaries with explicit labels:

* `null_ifp`: no interactions
* `vdw_only`: van der Waals / close contacts only
* `low_specific_contact`: non-vdW signal below the contact-eligibility threshold

Formal clustering is skipped when `n_contact_eligible < 10`; report the condition as `insufficient_clusterable_signal`. This protects against interpreting cluster counts and noise estimates as meaningful when the input set is too small.

### 14.4 Clustering method pilot

The primary clustering method is not selected yet. The pilot compares:

**Agglomerative Jaccard**

* metric: Jaccard distance on binary IFP vectors
* main linkage: average
* sensitivity linkage: complete
* distance cutoffs: 0.35, 0.45, 0.55
* minimum main cluster size: 3
* smaller groups are labelled `low_support_cluster` for reporting rather than silently dropped

**HDBSCAN Jaccard**

* metric: precomputed Jaccard distance
* `min_cluster_size`: 3 and 4
* `min_samples`: 1 and 2
* cluster selection method: `eom`
* run only on contact-eligible poses

Optional secondary check:

* PAM / k-medoids may be used only as a medoid-stability check, not as the primary model-selection method.

### 14.5 Pilot evaluation metrics

Per condition, report at least:

* `n_qc_pass`
* `n_ifp_success`
* `n_contact_eligible`
* `contact_eligible_fraction`
* `null_ifp_fraction`
* `vdw_only_fraction`
* `low_specific_contact_fraction`
* `median_n_non_vdw_interactions`
* `median_n_non_vdw_contact_residues`

Per method and parameter set, report at least:

* `n_clusters`
* `top_cluster_occupancy`
* `cluster_entropy`
* `occupancy_gini`
* `noise_fraction` for HDBSCAN
* `low_support_pose_fraction` for agglomerative clustering
* median within-cluster Jaccard distance
* median between-cluster Jaccard distance
* seed-mixing score
* bootstrap stability metrics where available
* post-hoc geometry annotation stability

Seed-mixing flags clusters where one seed contributes > 70% of members. Geometry is evaluated only after clustering; if C1/C4/non-plausible annotations change strongly between clustering methods, the biological interpretation is method-sensitive.

Bootstrap/subsampling stability is part of the pilot when `n_contact_eligible >= 10`:

* 100 repeats
* 80% subsampling
* stable if bootstrap ARI >= 0.60, top-cluster occupancy SD <= 0.15, and cluster-count SD <= 1.0

### 14.6 Method selection rule

Select one global primary method for the full analysis. Do not choose agglomerative for some conditions and HDBSCAN for others.

Prefer agglomerative Jaccard if:

* >= 70% of pilot conditions have `n_contact_eligible >= 10`
* top-cluster occupancy and cluster count are stable across cutoffs
* bootstrap stability is acceptable in most pilot conditions
* HDBSCAN gives similar dominant modes or no clear improvement

Prefer HDBSCAN if:

* substantial residual noise remains after contact eligibility filtering
* dense clusters are stable across `min_cluster_size` / `min_samples`
* bootstrap stability is better than agglomerative clustering
* noise fraction is not extremely parameter-sensitive
* agglomerative clustering produces many low-support or seed-specific clusters

Prefer coarse contact-profile reporting if:

* both formal methods are unstable
* many pilot conditions have `n_contact_eligible < 10`
* cluster assignments are strongly parameter-sensitive

If agglomerative and HDBSCAN give the same qualitative picture, use agglomerative Jaccard as the default primary method and keep HDBSCAN as a sensitivity analysis because agglomerative clustering is simpler and more transparent for small per-condition pose sets.

The selected method, parameters, retained/excluded interaction types, contact-eligibility rule, non-interpretable conditions, and required sensitivity analyses must be recorded in `clustering_method_decision.md` before the full analysis.

### 14.7 Output fields

`cluster\_assignments.tsv`:

* `pose\_id`
* `condition\_id`
* `cluster\_id`
* `cluster\_member\_flag`
* `noise\_flag` for HDBSCAN outputs
* `low_support_cluster_flag` for agglomerative outputs
* `contact_eligibility_status`
* `distance\_to\_cluster\_representative` if available

### 14.8 Medoid selection

For each retained cluster, choose one medoid pose in IFP space.

Medoid uses:

* visualization
* figure selection
* manual inspection
* representative structure for crystal sanity-check

Medoid does **not** define cluster statistics by itself.

### 14.9 Cluster inclusion rules

Recommended main analysis cluster set after method selection:

* exclude HDBSCAN noise from cluster-based inference, but report noise fraction
* exclude agglomerative low-support clusters from main cluster summaries, but report low-support fraction
* retain clusters with at least a small minimum support, for example:

  * `cluster\_size >= 3`, or
  * `occupancy >= 0.05`

Important:

* this threshold must be fixed in advance
* small clusters should still be retained in raw outputs even if excluded from main summaries

### 14.10 Outputs after clustering

#### `medoid\_manifest.tsv`

* `cluster\_id`
* `medoid\_pose\_id`
* `medoid\_structure\_path`
* `medoid\_ifp\_distance\_sum`

#### `condition\_cluster\_summary.tsv`

Per condition:

* `n\_qc\_pass\_poses`
* `n\_ifp\_success`
* `n\_contact\_eligible`
* `contact_eligible_fraction`
* `n\_ifp\_clustered`
* `n\_noise` where applicable
* `noise\_fraction` where applicable
* `n_low_support` where applicable
* `low_support_fraction` where applicable
* `n\_clusters`
* `top\_cluster\_occupancy`
* `cluster\_entropy`
* `occupancy\_gini`

\---

## 15\. Stage 7: Cluster annotation

### 15.1 Cluster-level aggregation

Build one row per cluster in `cluster\_table.tsv`.

Current implementation status 2026-05-20:

* the production analysis-core path now writes a first `cluster\_table.tsv`
  as a flat export of retained non-noise Stage 7 `cluster_summaries`
* this current implementation reuses the same Stage 7 aggregation already
  written to `cluster\_signatures.json`
* the richer field set below remains the recommended target contract for later
  enrichment where those fields are not yet materialized in code

Required fields:

* identifiers:

  * `cluster\_id`
  * `condition\_id`
  * `protein\_id`
  * `construct\_type`
  * `substrate\_class`
  * `dp`
* size / support:

  * `cluster\_size`
  * `occupancy`
  * `noise\_excluded\_denominator`
* QC:

  * `qc\_soft\_flag\_fraction`
  * `n\_hard\_qc\_fail\_pre\_ifp` at condition level
* geometry:

  * `geometry\_computable\_fraction`
  * `geometry\_plausible\_fraction`
  * `geometry\_highly\_plausible\_fraction`
  * median/IQR of all geometry metrics
* AF3 confidence:

  * `iptm\_mean`
  * `iptm\_sd`
  * `iptm\_median`
  * `iptm\_iqr`
  * `mean\_plddt\_mean`
  * `mean\_plddt\_sd`
* stability:

  * `median\_ligand\_rmsd\_to\_reference`
  * `convergence\_fraction`
* representation:

  * `medoid\_pose\_id`

### 15.2 Cluster type annotation

Assign operational labels:

* `C1\_compatible`
* `C4\_compatible`
* `mixed\_compatible`
* `non\_plausible`
* `uncertain`

Operational logic should be defined in YAML using geometry thresholds.

### 15.3 Cluster-level IFP signature

For each cluster, compute:

* frequency of every IFP feature across member poses
* frequency of each contacting residue across member poses
* dominant interaction type per residue

Store:

* `cluster\_ifp\_signature.tsv`
* `cluster\_residue\_signature.tsv`

`cluster\_residue\_signature.tsv` columns:

* `cluster\_id`
* `residue\_chain`
* `residue\_number`
* `residue\_name`
* `contact\_frequency`
* `dominant\_interaction\_type`
* `contact\_frequency\_hbond`
* `contact\_frequency\_hydrophobic`
* `contact\_frequency\_aromatic`
* `is\_cbm\_region`
* `is\_linker\_region`
* `is\_catalytic\_surface\_region`

This table is the core bridge from clustering to residue interpretation.

\---

## 16\. Stage 8: Residue importance analysis

This stage addresses the requirement to say something about residues relevant to C1/C4 activity and ligand specificity.

### 16.1 Important caution

Raw residue numbers are not directly comparable across unrelated proteins.

Therefore residue interpretation must be split into:

1. **within-protein residue importance**
2. **within-family or structurally aligned residue comparison**
3. **cross-dataset residue-property / patch-level interpretation**, not naive pooled residue numbering

### 16.2 Within-protein residue importance

For each protein and condition, compute occupancy-weighted residue contact score:

```text
residue\_contact\_score = sum(cluster\_occupancy x residue\_contact\_frequency\_in\_cluster)
```

Generate:

* `protein\_condition\_residue\_scores.tsv`

Columns:

* `protein\_id`
* `condition\_id`
* `residue\_chain`
* `residue\_number`
* `residue\_name`
* `occupancy\_weighted\_contact\_score`
* `occupancy\_weighted\_hbond\_score`
* `occupancy\_weighted\_hydrophobic\_score`
* `occupancy\_weighted\_aromatic\_score`
* `top\_cluster\_contact\_flag`
* `geometry\_weighted\_C1\_contact\_score`
* `geometry\_weighted\_C4\_contact\_score`

Interpretation:

* identifies which residues dominate ligand engagement for that protein under that ligand condition

### 16.3 Residues linked to C1 vs C4 compatibility within protein

For each protein, compare residues enriched in:

* C1-compatible clusters vs non-C1-compatible clusters
* C4-compatible clusters vs non-C4-compatible clusters

Use cluster-level presence/absence weighted by occupancy.

Outputs:

* `protein\_residue\_regio\_delta.tsv`

Columns:

* `protein\_id`
* `residue\_id`
* `delta\_contact\_C1\_minus\_C4`
* `delta\_hbond\_C1\_minus\_C4`
* `delta\_aromatic\_C1\_minus\_C4`
* `supporting\_cluster\_count`

Interpretation:

* residues with strong positive delta are candidate C1-associated contact residues for that protein
* residues with strong negative delta are candidate C4-associated contact residues for that protein

### 16.4 Within-family residue comparison

If sequence alignment or structural alignment within family is feasible:

* align proteins within family
* map residues to alignment columns or local structural pocket positions
* test whether aligned positions or local pocket bins show different contact occupancy in proteins with experimental C1 vs C4 labels or different ligand specificity labels

Outputs:

* `family\_aligned\_residue\_table.tsv`
* `family\_residue\_enrichment.tsv`

This is optional but strongly recommended if family sample size supports it.

### 16.5 Cross-dataset residue-property interpretation

Because exact residue identities are not portable across all proteins, also compute patch/property-level summaries.

Examples:

* aromatic contact density near catalytic surface
* polar contact density
* charged residue contact density
* fraction of contacts contributed by loop regions
* fraction of contacts contributed by residues within distance shells from Cu

Outputs:

* `condition\_patch\_summary.tsv`
* `protein\_patch\_summary.tsv`

Interpretation:

* better for general statements about ligand specificity across proteins

### 16.6 Recommended final residue outputs for the thesis

At minimum include:

1. protein-specific heatmaps of top contacting residues for representative proteins
2. occupancy-weighted residue scores for each ligand class
3. delta-contact maps for C1-like vs C4-like clusters within selected proteins
4. family-level aligned residue summaries only where alignment is trustworthy
5. patch/property summaries across the broader dataset

\---

## 17\. Stage 9: Protein–ligand-condition summaries

Build `condition\_table.tsv` with one row per protein–ligand condition.

Required columns:

* identifiers:

  * `condition\_id`
  * `protein\_id`
  * `construct\_type`
  * `substrate\_class`
  * `dp`
* pose counts:

  * `n\_generated`
  * `n\_qc\_pass`
  * `n\_ifp\_success`
  * `n\_geometry\_computable`
  * `n\_geometry\_plausible`
* cluster metrics:

  * `n\_clusters`
  * `noise\_fraction` where applicable
  * `low_support_fraction` where applicable
  * `insufficient_clusterable_signal_flag`
  * `top\_cluster\_occupancy`
  * `cluster\_entropy`
  * `occupancy\_gini`
* geometry summaries:

  * occupancy-weighted medians for main geometry metrics
  * fraction occupancy in C1-compatible clusters
  * fraction occupancy in C4-compatible clusters
  * fraction occupancy in non-plausible clusters
* residue summaries:

  * top 5 residues by occupancy-weighted contact score
  * aromatic\_contact\_fraction
  * polar\_contact\_fraction
  * loop\_contact\_fraction
* AF3 confidence summaries:

  * `iptm\_mean\_all`
  * `iptm\_sd\_all`
  * `iptm\_mean\_qc\_pass`
  * `iptm\_sd\_qc\_pass`
  * `iptm\_mean\_clustered`
  * `iptm\_sd\_clustered`

Purpose:

* this is the central table for repeated-measures comparisons across ligand conditions

\---

## 18\. Stage 10: Protein-level summaries

Build `protein\_summary\_table.tsv`.

One row per protein.

Purpose:

* compact overview
* secondary summaries
* exploratory predictive inputs if needed

Recommended fields:

* `protein\_id`
* `family\_label`
* `has\_cbm`
* `n\_conditions\_attempted`
* `mean\_qc\_pass\_rate`
* `mean\_n\_clusters`
* `mean\_cluster\_entropy`
* `mean\_top\_cluster\_occupancy`
* `mean\_geometry\_plausible\_fraction`
* `mean\_iptm`
* `sd\_iptm`
* occupancy-weighted fraction of C1-compatible clusters across all conditions
* occupancy-weighted fraction of C4-compatible clusters across all conditions
* ligand-class-specific cluster entropy and contact-density summaries

Do not use this as the only main analysis table unless forced by time.

\---

## 19\. Stage 11: Crystal sanity-check

### 19.1 Input

Use every retained cluster medoid for crystal comparison, not all poses. If a
condition has no retained clusters, the top-level AF3 model CIF may be used as a
fallback only after it passes hard QC. This fallback is not included in normal
pose, clustering, or medoid denominators.

### 19.2 Why medoid is appropriate here

Crystal comparison is a representative sanity-check, so retained medoids are the
right compression target for clustered conditions. The hard-QC-passing top-level
AF3 fallback exists only to keep no-cluster conditions interpretable against
crystal references; it does not redefine clustering results.

### 19.3 Required outputs

`crystal\_anchor\_table.tsv`:

* `cluster\_id`
* `representative\_pose\_id`
* `representative\_role`
* `medoid\_pose\_id`
* `crystal\_reference\_id`
* `local\_pocket\_rmsd`
* `ligand\_rmsd\_if\_comparable`
* `proximal\_sugar\_rmsd`
* `contact\_overlap\_score` / `ifp\_tanimoto` only when the prepared crystal IFP passes the same non-vdW contact-eligibility rule used for pose clustering
* `crystal\_ifp\_contact\_eligible`
* `crystal\_ifp\_exclusion\_class`
* `ifp\_comparison\_eligible`
* `crystal\_Cu\_C1\_distance`
* `crystal\_Cu\_C4\_distance`
* `crystal\_geometry\_status\_C1`
* `crystal\_geometry\_status\_C4`
* `same\_binding\_region\_flag`
* `same\_general\_orientation\_flag`
* `notes`

`crystal\_geometry\_table.tsv` stores full C1/C4 geometry metrics for prepared
ligand-bound crystal references. `crystal\_ifp\_diagnostic\_summary.tsv` stores
eligibility/exclusion counts and percentages by unique crystal reference and by
medoid/fallback comparison row.

Current implementation detail (2026-05-21): ligand-bound crystal references are
first written as chemistry-preserving selected mmCIFs, remapped to the canonical
`A` protein / `B-C-D` glycan / `E` Cu scheme, then normalized before
protonation. The full crystal complex keeps Cu for geometry and pocket RMSD, but
the ProLIF ligand MOL2 is glycan-only and excludes Cu, water, ions, and common
buffer/solvent residues. ProLIF resolves ligand residue IDs from the companion
`ligand_only_for_prolif.pdb` when available, preserving distinct glycan chains
even when Obabel MOL2 substructure labels repeat as `BGC1`, `BGC2`, and so on.
A crystal IFP that remains `vdw_only`,
`null_ifp`, or `low_specific_contact` after this prep is treated as a
non-comparable contact-signal outcome, not as known prep contamination.

### 19.4 Interpretation rule

Crystal mismatch does not invalidate a cluster. It is used only for context and plausibility discussion.

\---

## 20\. Stage 12: Removed

~~Optional PLACER sensitivity branch~~ — **PLACER is removed from the pipeline entirely (decision 2026-04-21).**

\---

## 21\. Stage 13: Descriptive analysis

This is the main analysis layer.

### 21.1 Required descriptive outputs per sub-analysis

1. QC attrition plot
2. distribution of cluster counts per protein
3. top-cluster occupancy distribution
4. cluster entropy distribution
5. fraction of conditions with at least one plausible C1-compatible cluster
6. fraction of conditions with at least one plausible C4-compatible cluster
7. ipTM distribution vs QC pass / geometry plausibility / cluster occupancy
8. cluster-level residue contact heatmaps
9. selected medoid structure figures
10. clustering rate versus crystal-structure coverage per protein

### 21.2 Required descriptive outputs across conditions

For each protein:

* compare chitin vs cellulose vs starch
* compare DP4 vs DP6 vs DP8
* compare domain-only vs full-length when CBM exists
* compare clustering rate against crystal-reference coverage:
  * primary view: number of available crystal structures for the protein
  * secondary/fallback view: binary `has_crystal_reference` yes/no

Recommended tables:

* `descriptive\_condition\_comparison.tsv`
* `descriptive\_cluster\_distribution.tsv`

### 21.3 Residue-focused outputs

Must include:

* top residues by occupancy-weighted contact score for each representative protein-condition
* residues recurrently found in high-occupancy clusters
* residues enriched in C1-compatible vs C4-compatible clusters within selected proteins
* residue-class summaries for ligand specificity

\---

## 22\. Stage 14: Exploratory predictive analysis

### 22.1 Scope

Keep simple. Do not build a large model zoo.

Update 2026-05-21: the downstream activity-prediction plans were revised into
smaller, less-detailed planning documents. They should now be treated as
compact exploratory plans with small-effective-n constraints rather than full
implementation blueprints.

### 22.2 Acceptable prediction targets

Choose at most one or two main tasks, for example:

1. `experimental\_regio\_label`: C1 vs C4
2. `experimental\_substrate\_label`: e.g. chitin-preferring vs cellulose-preferring where labels are credible

### 22.3 Allowed rows

Recommended primary row type:

* cluster

Alternative if cluster count becomes unmanageable:

* protein–ligand condition

### 22.4 Allowed features

Use a limited set only:

* occupancy
* geometry plausibility fractions
* top-cluster occupancy
* cluster entropy
* selected contact-feature frequencies
* selected residue-property summaries

Do not use huge sparse raw IFP vectors unless feature selection is tightly controlled.

### 22.5 Validation

Must use grouped CV at protein level.

Rule:

* all conditions and clusters from the same protein stay in the same fold

### 22.6 Metrics

* balanced accuracy
* macro F1
* AUROC only if class setting makes sense

### 22.7 Interpretation rule

Results are exploratory. They test whether the summarized AF3 output contains recoverable signal, not whether the model has discovered causal biology.

\---

## 23\. Stage 15: CBM paired analysis

### 23.1 Input

Only proteins with both constructs available.

### 23.2 Paired comparisons

Within each protein and ligand condition, compare:

* domain-only vs full-length
* difference in top-cluster occupancy
* difference in cluster entropy
* difference in geometry plausibility fraction
* difference in ligand proximity to catalytic surface
* gain/loss of CBM-region contacts

### 23.3 Output table

`cbm\_comparison\_table.tsv`:

* `protein\_id`
* `substrate\_class`
* `dp`
* `domain\_only\_condition\_id`
* `full\_length\_condition\_id`
* `delta\_top\_cluster\_occupancy`
* `delta\_cluster\_entropy`
* `delta\_geometry\_plausible\_fraction`
* `delta\_C1\_compatible\_fraction`
* `delta\_C4\_compatible\_fraction`
* `delta\_cbm\_contact\_fraction`
* `cbm\_recruitment\_flag`

### 23.4 Statistics

If used:

* paired Wilcoxon signed-rank
* report effect sizes and direction, not only p-values

\---

## 24\. How outputs feed the next stage

This section is critical for implementation.

### 24.1 Flow of data

#### Stage 1 -> Stage 2

AF3 raw outputs become parsed pose-level metadata and confidence tables.

#### Stage 2 -> Stage 3 and Stage 4

Only QC-passing poses continue to IFP and geometry.

#### Stage 3 -> Stage 6

Pose-level binary IFP vectors are the direct input to clustering.

#### Stage 4 -> Stage 7

Pose-level geometry is joined back onto cluster assignments after clustering.

#### Stage 6 -> Stage 7

Cluster assignments and medoid selection define the cluster table and medoid manifest.

#### Stage 7 -> Stage 8

Cluster residue signatures and cluster types feed residue interpretation.

#### Stage 7 -> Stage 9

Cluster summaries are aggregated to protein–ligand condition summaries.

#### Stage 9 -> Stage 10 and Stage 14

Condition summaries feed protein summaries and exploratory prediction.

#### Stage 7 -> Stage 11

Retained cluster medoids feed crystal sanity-check. No-cluster conditions may
use a hard-QC-passing top-level AF3 model CIF fallback for crystal anchoring
only; that fallback does not enter clustering or generated-pose denominators.

### 24.2 Hard rule

No later stage should overwrite the main AF3 pose set. Later stages only annotate or compare.

\---

## 25\. Required main tables

### 25.1 Pose-level tables

* `pose_manifest.tsv`
* `pose_confidence.tsv`
* `pose_ifp_table.tsv`
* `pose_residue_contact_table.tsv`
* `pose_geometry.tsv`

### 25.2 Cluster-level tables

* `cluster\_assignments.tsv`
* `medoid\_manifest.tsv`
* `cluster\_table.tsv`
* `cluster\_ifp\_signature.tsv`
* `cluster\_residue\_signature.tsv`

### 25.3 Condition / protein-level tables

* `condition\_table.tsv`
* `condition\_cluster\_summary.tsv`
* `protein\_summary\_table.tsv`

### 25.4 Interpretation / side-analysis tables

* `protein\_condition\_residue\_scores.tsv`
* `protein\_residue\_regio\_delta.tsv`
* `qc\_attrition\_table.tsv`
* `crystal\_anchor\_table.tsv`
* `crystal\_geometry\_table.tsv`
* `crystal\_ifp\_diagnostic\_summary.tsv`
* `cbm\_comparison\_table.tsv`

\---

## 26\. Required figure classes

At minimum:

1. QC attrition waterfall per ligand condition
2. histogram / violin of number of clusters per condition
3. distribution of top-cluster occupancy
4. cluster entropy by substrate class and DP
5. geometry plausibility fractions by substrate class and DP
6. ipTM vs QC pass / geometry plausibility / cluster occupancy
7. residue contact heatmaps for representative proteins
8. medoid structural figures for selected clusters
9. crystal sanity overlays for representative medoids, with crystal-side C1/C4 geometry comparison tables/figures where ligand-bound crystal references exist
10. clustering rate versus crystal-structure coverage:
    primary x-axis is number of available crystal structures per protein;
    include a binary has/no-has crystal-reference summary as a secondary panel
    or grouped overlay

\---

## 27\. Sensitivity analyses

Recommended minimum set:

1. top-cluster-only summaries
2. only clusters above predefined occupancy threshold
3. only geometry-plausible clusters
4. medoid-only descriptive comparison

Interpretation rule:

If the same qualitative pattern appears across the main analysis and at least one or two sensitivity analyses, confidence increases.

\---

## 28\. Decision rules that must be predefined

Store these in `thresholds.yaml` before full analysis:

* severe PoseBusters fail definitions
* severe Privateer fail definitions
* cluster minimum support rule
* geometry plausibility thresholds
* convergence cutoff
* residue importance display threshold
* top-residue reporting count
* ipTM missingness handling rule

\---

## 29\. Recommended implementation order

1. Input manifest generation
2. AF3 output parser
3. QC wrappers and QC tables
4. ProLIF generation
5. Pose-level residue contact extraction
6. Geometry engine
7. Condition-wise clustering
8. Cluster annotation
9. Condition and protein summary builders
10. Residue-importance analysis
11. Crystal sanity-check
12. Descriptive plotting
13. Small exploratory predictive module

\---

## 30\. Non-negotiable rules

* AF3 is the only main prediction source.
* Pose-level observations are not treated as independent biological replicates.
* Main clustering uses IFP only.
* Geometry is linked after clustering.
* Medoid is for representation, not as the only analytical feature source.
* QC removes clearly invalid structures, but reactive geometry is not the main pre-IFP exclusion basis.
* ipTM is descriptive and comparative, not a primary validity criterion.
* Condition comparisons across ligand type and DP are treated as repeated measures within protein.
* Grouped CV, if used, is done at protein level.
* Crystal checks do not overwrite main AF3 results.

\---

## 31\. Minimal robust version if time is limited

If the pipeline must be reduced aggressively, keep:

1. AF3 parsing
2. hard QC
3. ProLIF IFP
4. condition-wise clustering
5. cluster annotation with medoid selection
6. condition summary table
7. residue contact scoring
8. geometry annotation
9. ipTM descriptive summaries
10. crystal sanity-check on retained medoids, plus hard-QC-passing top-level AF3 fallback for no-cluster conditions

Cut first:

* family-aligned residue analysis
* predictive modeling beyond one very small model
* complex CBM statistics

\---

## 32\. Final practical interpretation framework

This pipeline can support claims of the following kind:

* AF3 often produces one / several stable binding modes for a given protein–ligand condition.
* Certain ligand classes or DP values are associated with more stable or more multimodal output.
* Certain residues or residue patches recurrently contribute to stable contact patterns.
* Some clusters are more C1-like or C4-like under predefined operational geometry rules.
* Full-length models may alter contact distribution relative to domain-only models.
* ipTM may or may not track useful pose-level or cluster-level plausibility; this can be shown directly.

This pipeline should not support strong claims that:

* one cluster is the true biological state
* high ipTM proves correct ligand pose
* virtual oxyl geometry proves mechanism
* pose count equals independent sample size
