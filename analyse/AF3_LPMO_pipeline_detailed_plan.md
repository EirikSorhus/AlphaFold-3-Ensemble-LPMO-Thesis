# AF3-LPMO Analysis Pipeline: Detailed Practical Plan

Version: 1.0  
Status: Recommended main-analysis design  
Scope: AF3-only pipeline for LPMO–oligosaccharide complexes with pose-level QC, IFP-based clustering, geometry annotation, residue-level interpretation, and limited exploratory prediction.

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
* Main pose ensemble per protein–ligand condition: **15 seeds x 5 samples = 75 poses** (AF3 `num_diffusion_samples=5`, runs completed).
* Main clustering input: **ProLIF binary IFP only**.
* Main descriptive unit: **cluster**.
* Main biological repeated-measures unit: **protein**.
* Main inferential caution: **poses and seeds are not independent biological replicates**.
* Main QC philosophy: **remove clearly invalid structures early; do not hard-filter primarily on downstream mechanistic geometry before IFP**.
* Main geometry philosophy: **compute after IFP; use as annotation, plausibility, and interpretation**.
* Main medoid philosophy: **use for visualization, figure choice, manual inspection, and crystal sanity-check; do not use as the only analytical representation**.

### 2.2 Optional / weak-priority components

* Crystal comparison is a sanity-check, not a primary validation criterion.
* Predictive modeling is exploratory and should remain small and simple.

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
* ligand missing / broken / unparsable / long distance from Cu
* Cu missing
* protein–ligand structure cannot be parsed consistently
* atom naming prevents residue or ligand mapping
* structure file corrupt or incomplete

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

Current validated default interaction types (matching `configs/prolif_features.yaml`):

* Hbond donor
* Hbond acceptor
* hydrophobic
* aromatic / stacking
* anionic
* cationic
* cation-pi
* pi-cation
* van der Waals contact

This broad default is useful for early validation and inspection, but avoid
creating many fragile interaction classes that explode sparsity. The final set
of interaction types to keep can be decided later after reviewing the real-data
behavior. Keep the interaction set locked within one sub-analysis. If later
sparsity review motivates pruning, do it as a config change before the full run
and rerun the entire sub-analysis; do not mix feature spaces within the same
batch.

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

* **Pre-clustering / per-pose computation**: use the top-ranked QC-passing pose within the condition (by `ranking_score`, or `mean_plddt` if ranking score is absent).
* **Final reporting**: replace reference with the top-occupancy cluster medoid once clustering is complete.

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

* all QC-passing poses with successful IFP generation

### 14.2 AF3-only simplification

Because the main pipeline is AF3-only:

* remove within-model vs cross-model clustering separation
* perform **one clustering per protein–ligand condition**

### 14.3 Clustering method

Recommended default:

* HDBSCAN on binary IFP vectors
* Jaccard distance

Store parameters in config:

* `min\_cluster\_size`
* `min\_samples`
* `distance\_metric`

### 14.4 Output fields

`cluster\_assignments.tsv`:

* `pose\_id`
* `condition\_id`
* `cluster\_id`
* `cluster\_member\_flag`
* `noise\_flag`
* `distance\_to\_cluster\_representative` if available

### 14.5 Medoid selection

For each non-noise cluster, choose one medoid pose in IFP space.

Medoid uses:

* visualization
* figure selection
* manual inspection
* representative structure for crystal sanity-check

Medoid does **not** define cluster statistics by itself.

### 14.6 Cluster inclusion rules

Recommended main analysis cluster set:

* exclude HDBSCAN noise as its own category from cluster-based inference, but still report noise fraction
* retain clusters with at least a small minimum support, for example:

  * `cluster\_size >= 3`, or
  * `occupancy >= 0.05`

Important:

* this threshold must be fixed in advance
* small clusters should still be retained in raw outputs even if excluded from main summaries

### 14.7 Outputs after clustering

#### `medoid\_manifest.tsv`

* `cluster\_id`
* `medoid\_pose\_id`
* `medoid\_structure\_path`
* `medoid\_ifp\_distance\_sum`

#### `condition\_cluster\_summary.tsv`

Per condition:

* `n\_qc\_pass\_poses`
* `n\_ifp\_clustered`
* `n\_noise`
* `noise\_fraction`
* `n\_clusters`
* `top\_cluster\_occupancy`
* `cluster\_entropy`
* `occupancy\_gini`

\---

## 15\. Stage 7: Cluster annotation

### 15.1 Cluster-level aggregation

Build one row per cluster in `cluster\_table.tsv`.

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
  * `noise\_fraction`
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

Use cluster medoids for crystal comparison, not all poses.

### 19.2 Why medoid is appropriate here

Crystal comparison is a representative sanity-check, so medoid is the right compression target.

### 19.3 Required outputs

`crystal\_anchor\_table.tsv`:

* `cluster\_id`
* `medoid\_pose\_id`
* `crystal\_reference\_id`
* `local\_pocket\_rmsd`
* `ligand\_rmsd\_if\_comparable`
* `proximal\_sugar\_rmsd`
* `contact\_overlap\_score`
* `same\_binding\_region\_flag`
* `same\_general\_orientation\_flag`
* `notes`

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

### 21.2 Required descriptive outputs across conditions

For each protein:

* compare chitin vs cellulose vs starch
* compare DP4 vs DP6 vs DP8
* compare domain-only vs full-length when CBM exists

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

Cluster medoids feed crystal sanity-check.

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
9. crystal sanity overlays for representative medoids

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
10. crystal sanity-check on a limited medoid subset

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

