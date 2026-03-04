# LPMO Structure Prediction → Analysis Pipeline — Integrated Master Plan (v2.1)

## 1. SCOPE & INVARIANTS

- **Models:** AF3, RF3, Boltz-2
- **Targets:** LPMO protein + glykan ligands (DP2–DP8) + Cu
- **Master format:** mmCIF (ModelCIF-compatible)
- **Prediction code:** lives in `structure_pipeline/`; `analysis/` orchestrates only
- **Hard rules (non-negotiable):**
  1. Privateer input: ONLY valid CCD monosaccharides (NAG, BGC, MAN …); NO custom oligomer IDs
  2. Glykan DP = sequence of monosakkarid residues + linkers (explicit `_struct_conn`)
  3. PLACER mandatory before any final scoring/reporting
  4. Atom names NOT assumed same across models; use cross-model mapping `(element, CCD, bond_graph, 3D)`
  5. Chain schema: protein=A, glycans=B..D, metal=E
  6. `_chem_comp_bond` complete for ALL `comp_id`; `_struct_conn` for glycosidic + Cu-coord
  7. Confidence: atom-level `_atom_site.B_iso_or_equiv` (all); residue-level `_ma_qa_metric_local` (AF3 only); preserve raw; derive median from `B_iso` only if field missing
  8. HDBSCAN hyperparams LOCKED post-tuning (anti p-hack)
  9. Grouped CV (hierarchy: structure→protein→ligand→model→seed)
  10. Stratified CV (maintain class balance: C1/C4/mixed, polymer type)

---

## 2. I/O CONTRACTS

| Artifact | Format | Key columns / fields |
|---|---|---|
| `normalized.cif` | mmCIF | canonical chains, `_chem_comp_bond`, `_struct_conn`, confidence |
| `atom_map.tsv` | TSV | old_atom, new_atom, element, ccd, confidence, reason |
| `rename_log.json` | JSON | total_atoms, mapped, coverage, unmapped list |
| `qc_report.json` | JSON | per-pose: PB pass/fail, Privateer pass/fail, Cu-His dist, flags |
| `placer_scores.json` | JSON | per-pose: placer_score, clash_count, rank |
| `ifp_matrix.csv` | CSV | rows=poses, cols=residue×interaction_type, binary 0/1 |
| `clusters.json` | JSON | cluster_id → {poses, occupancy, medoid, geometry_stats} |
| `cluster_signatures.json` | JSON | per cluster: dominant interactions, geometry aggregate |
| `crystal_metrics.json` | JSON | Tanimoto, pocket_RMSD per protein×ligand |
| `activity_features.json` | JSON | per protein: C1/C4 occupancy, interaction freq, CBM metrics |
| `metrics.csv` | CSV flat | run_id, protein_id, ligand_id, model, pose_id, cluster_id, Cu_C1, Cu_C4, angle, ifp_sim, … |
| `summary.json` | JSON | dataset stats, QC stats, cluster stats, geometry stats |
| `report.html` | HTML | human-readable with embedded plots |
| `run_manifest.json` | JSON | git_commit, config_hash, tool_versions, seeds, checksums, gates |

---

## 3. PIPELINE STEPS

### STEP 1 — INGEST & QC
- **In:** raw prediction (PDB/mmCIF)
- **Out:** ingest_qc.log, pass/fail
- **Tools:** Gemmi
- **Checks:** parse OK, `_atom_site`/`_chem_comp`/`_entity` present, atom_count > 0
- **Gate:** Fail → log + skip

### STEP 2 — NORMALIZE MMCIF
- **In:** raw mmCIF
- **Out:** `normalized.cif`, `atom_map.tsv`, `rename_log.json`
- **Tools:** Gemmi, cross-model atom mapper
- **Procedure:** assign chains (A/B..D/E), map atoms (element+CCD+bond_graph+3D), ensure `_chem_comp_bond` + `_struct_conn`, derive confidence if missing
- **Gate:** `atom_mapping_coverage = 100%`; fail → skip + log

### STEP 3 — GLYCAN EXPANSION (Privateer prep)
- **In:** `normalized.cif`
- **Out:** `privateer_input.cif`
- **Checks:** every glykan `comp_id` ∈ CCD monosaccharides whitelist; glycosidic linkage in `_struct_conn`; Cu-His present
- **Gate:** 100% recognized monosacchs; fail → skip + log

### STEP 4 — PROTONATION & EXPORT
- **In:** `normalized.cif`
- **Out:** `for_posebusters.pdb`, `ligand_for_prolif.mol2`, `complex_H.pdb`
- **Tools:** Reduce (H network, verify Cu-His), OpenBabel (bond orders, Gasteiger charges)
- **Checks:** CONECT present, bond_orders present, partial_charges present

### STEP 5 — PLACER REFINEMENT (MANDATORY ← hard rule)
- **In:** `normalized.cif`
- **Out:** `placer_ensemble/` (50–200 poses), `placer_scores.json`
- **Tools:** PLACER (multi-mode GPU)
- **Rank by:** PLACER_score, clash_count, Cu_geom_sanity
- **Gate:** ≥1 pose; fail → skip + log

### STEP 6 — VALIDATION (QC)
- **In:** `privateer_input.cif`, `best_placer.pdb`
- **Out:** `qc_report.json`
- **Tools:** PoseBusters (complex mode), Privateer (sugar_id, anomer, ring_pucker, linkage), Cu-geometry checker
- **Hard-fail:** critical PB error → drop; Privateer anomer fail → drop; Cu–His > 2.6 Å → drop
- **Soft-flag:** minor PB warning → keep + mark
- **Gates:** `no_critical_posebusters_errors = true`; `privateer_recognized_sugars = 100%`

### STEP 7 — ANALYSIS (per pose)
- **In:** `best_placer.pdb`, `complex_H.pdb`, `ligand.mol2`
- **Out:** `ifp_matrix.csv`, `geometry_metrics.json`
- **Tools:** ProLIF (IFP), MDAnalysis (Cu–C1, Cu–C4, His-brace angle, planarity, SASA, H-bonds)

### STEP 8 — CLUSTERING
- **In:** `ifp_matrix.csv` per protein×ligand×model
- **Out:** `clusters.json`
- **Tools:** HDBSCAN (Jaccard on binary IFP)
- **Phase 1:** within-run (compress seed noise)
- **Phase 2:** cross-run per protein×ligand (medoid IFPs)
- **Params:** LOCKED from tuning (`min_cluster_size`, `metric='jaccard'`, `epsilon=0.0`, EOM)

### STEP 9 — 3D GEOMETRY PER CLUSTER
- **In:** cluster assignments + `geometry_metrics.json` + `ifp_matrix.csv`
- **Out:** `cluster_signatures.json`
- **Per cluster:** select medoid; aggregate Cu–C1/C4 (mean/std), His-brace angle, planarity; dominant interactions (top 3–5)

### STEP 10 — CRYSTAL ANCHORING (optional)
- **In:** predicted clusters, crystal.pdb (if exists)
- **Out:** `crystal_metrics.json`
- **If ligand-bound crystal:** Tanimoto(pred_IFP, crystal_IFP); threshold ~0.5 empirical
- **If apo crystal:** pocket RMSD (local); flag if < 2.5 Å

### STEP 11 — CLUSTER→ACTIVITY MAPPING
- **In:** cluster signatures (all ligands, models)
- **Out:** `activity_features.json`
- **Per protein:** C1-like occupancy (Cu–C1 < Cu–C4 + 0.5 Å), C4-like, interaction freq by substrate

### STEP 12 — DEL B VARIANT (CBM, parallel)
- Same as steps 1–11 but with full-length (LPMO+CBM) input
- **Extra:** dual IFP (IFP_LPMO, IFP_CBM), proximity metrics (CBM–ligand, CBM–Cu), paired stats vs DEL A

### STEP 13 — REPORTING
- **In:** all step outputs
- **Out:** `summary.json`, `metrics.csv`, `report.html`

---

## 4. TUNING MODE

**Dataset:** AA9 (4 prot × 2 lig), AA10–17 subsample; C1/C4/both; ±crystal

| Model | Phase | Vary | Hold | Metric |
|-------|-------|------|------|--------|
| AF3 | Refinement | `num_recycles ∈ {10,15,20}` | diffusion=5, seeds=10 | QC pass, outlier%, clusters |
| AF3 | Diversity | `seeds ∈ {20,50,100}` | recycles=chosen | cluster stability |
| RF3 | Refinement | `n_recycles ∈ {10,20,30}` | seed=42, diff=5, steps=200 | QC pass, outlier% |
| RF3 | Diversity | `seed ∈ {42,80,120}` | n_recycles=chosen | same |
| Boltz2 | Refine-1 | `recycling ∈ {3,6,10}` | diff=5, sampling=200 | QC, outlier% |
| Boltz2 | Refine-2 | `sampling ∈ {200,400,600}` | recycling=chosen | same |
| Boltz2 | Diversity | `diffusion ∈ {3,10,15}` | recycling+sampling | same |

**Always for Boltz2:** `use_potentials=False`, `step_scale=1.638`

**Decision rule:** best QC pass → lowest outlier_rate → stable clusters → crystal sanity ≥ baseline

**Output:** `best_{model}.yaml` (locked), `tuning_summary.json`, `tuning_report.html`

---

## 5. PRODUCTION MODE

- **Input:** locked params from tuning
- **DEL A:** all proteins × 21 ligands × 3 models
- **DEL B:** CBM proteins × 21 ligands × 3 models (parallel pipeline)
- **Anti p-hack:** HDBSCAN params FIXED; no parameter search; grouped CV by structure; stratified CV for class balance

---

## 6. DATAMODEL

- **Run ID:** `{model}_{protein_id}_{ligand_id}_{seed}_{timestamp}`
- **Pose ID:** `{run_id}_pose_{N}`
- **Cluster ID:** `{protein_id}_{ligand_id}_cluster_{N}`
- **Hierarchy:** Structure → Protein family → Protein → Ligand (type+DP) → Model → Seed/pose

---

## 7. FAILURE POLICY

| Event | Hard? | Action | Log file |
|---|---|---|---|
| atom_mapping < 100% | YES | skip | mapping_failures.json |
| glykan not CCD | YES | skip | privateer_prep_failures.json |
| PLACER 0 poses | YES | skip | placer_failures.json |
| critical PB error | YES | drop pose | qc_report.json |
| soft PB warning | NO | keep + mark | qc_report.json (warnings) |
| Privateer anomer fail | YES | drop | qc_report.json |
| Cu–His > 2.6 Å | YES | drop | geometry_failures.json |
| HDBSCAN outlier | NO | keep, label=-1 | cluster_report.json |
| crystal_sim < 0.3 | NO | flag | crystal_report.json |

**Dropped poses:** logged metadata only (not in final `metrics.csv`); in failures log.

---

## 8. REPRODUSERBARHET

Every run writes `run_manifest.json`:
- pipeline_version, timestamp, git_commit, config_hash
- tool_versions (Gemmi, PLACER, PB, Privateer, MDA, ProLIF, HDBSCAN)
- all seeds, input checksums (SHA-256)
- gates_passed: atom_mapping_100, privateer_100, no_critical_pb, cu_his_ok
- tuning_reference (if production)

All random seeds logged. CCD results cached. Atom mapping deterministic.

---

## RQ → Outputs

| RQ | Analysis | Key metric | Artifact |
|----|----------|-----------|----------|
| RQ1 (C1/C4 geometry) | cluster geometry signatures | Cu–C1/C4 dist + His angle | cluster_signatures.json |
| RQ2 (substrate spec) | occupancy per polymer type | per-substrate occupancy | activity_features.json |
| RQ3 (CBM) | DEL A vs DEL B | Jaccard(IFP), paired occupancy | cbm_comparison |
| RQ4 (crystal) | crystal anchoring | Tanimoto, pocket RMSD | crystal_metrics.json |
