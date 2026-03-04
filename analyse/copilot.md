# Copilot folder guide — LPMO_structure_pipeline_v2_1

**Purpose**
This folder contains a structure-prediction → normalization → refinement → validation → interaction-fingerprinting → clustering → geometry + reporting pipeline for LPMO enzyme–oligosaccharide complexes.

---

## Repository layout + code ownership (project map)

**Root layout (important paths)**

- `eirik/Masteroppgave/`
  - `ligands/`
    - Ligandfiler + noe filkonvertering. **Ikke del av pipeline-arbeidet her.**
  - `structure_pipeline/`
    - **Source-of-truth for å sende inn/kjøre prediksjonsmodellene**:
      - AF3
      - RF3
      - Boltz-2
    - Denne mappen inneholder modell-“runners”/adapters og skal være stedet som faktisk kaller eksterne verktøy.
  - `analysis/`
    - Analyse-/orkestreringskode og dokumentasjon. **Denne `copilot.md` ligger her.**
    - Koden her skal *ikke* re-implementere modellkall; den skal orkestrere pipeline og bygge jobber/batcher.
  - `pipe_test/`
    - Uryddig/eksperimentell. **Skal ikke brukes** (ikke importer, ikke les configs herfra).
- `eirik/conda/`
  - Conda environments (kjøremiljø), men ikke pipeline-logikk.

### Hard rule: where prediction code lives
All code that runs AF3 / RF3 / Boltz-2 must be invoked via:
- `eirik/Masteroppgave/structure_pipeline/`

**Policy**
- `analysis/` may generate inputs/configs and schedule jobs, but must call the prediction runners/adapters from `structure_pipeline/`.
- Do not duplicate AF3/RF3/Boltz-2 invocation logic inside `analysis/`.
- If `analysis/` needs a new capability (new flags, new output parsing), implement it in `structure_pipeline/` adapters and call it from `analysis/`.

---

**Model preference (VSCode Copilot)**
When generating code, prefer **Claude Opus 4.6** if available in the Copilot model selector. Otherwise use the best available reasoning-capable model.

---

## Source-of-truth priority (conflict policy)

1) **Analysis plan (attached in prompt / doc)** is the **highest priority** and must be followed if anything conflicts. :contentReference[oaicite:1]{index=1}  
2) The **program + IO contract** below (the spec `LPMO_structure_pipeline_v2_1`).  
3) Reasonable engineering defaults.

If you (Copilot) see ambiguous/contradicting instructions, **follow (1)** and document the decision in code comments or `docs/decisions.md`.

---

## High-level goals (from analysis plan)

We need a reproducible pipeline that:
- Generates complex structures with **AF3 / RF3 / Boltz2**
- Performs **hard QC** with **PoseBusters + Privateer**
- Runs **PLACER refinement** (mandatory, GPU multi-mode)
- Creates **ProLIF interaction fingerprints (IFP)**
- Clusters binding modes with **HDBSCAN** (Jaccard/Tanimoto on binary IFP)
- Measures mechanistic geometry with **MDAnalysis**
- Produces combined reporting artifacts:
  - `summary.json`, `metrics.csv`, `report.html`

The analysis plan includes tuning workflows (parameter sweeps) and a main analysis (Del A: all proteins; Del B: full-length + CBM subset). :contentReference[oaicite:2]{index=2}

---

## HARD RULES (must not be violated)

1) **Privateer input requires monosaccharides as valid CCD `comp_id`** (e.g. `NAG`, `BGC`, `MAN`) — not only custom oligomer codes.
2) **PLACER is mandatory** before final analysis. Run multi-mode ensemble on GPU; runtime is not considered limiting.
3) **Atom names differ across models**; mapping must be **topology + geometry based**, not name-based.
4) **mmCIF is master format**. PDB/MOL2/SDF are tool-specific derivatives.

---

## Canonical formats + confidence policy

- Canonical internal format: **mmCIF**
- Confidence fields:
  - per-atom (all models): `_atom_site.B_iso_or_equiv`
  - per-residue (AF3): `_ma_qa_metric_local`
- Policy: preserve raw confidence values; only derive missing metadata.

---

## Normalization requirements

**Chain schema**
- Protein: chain `A`
- Glycans: chains `B..D`
- Metal: chain `E`

**Required annotations in normalized mmCIF**
- complete `_chem_comp` / `_chem_comp_bond`
- `_struct_conn` for:
  - glycosidic linkages
  - metal coordination

**Cross-model mapping key**
- element
- residue CCD
- local bond graph
- 3D proximity

---

## Pipeline steps (must reflect analysis plan + IO contract)

### step_1_prediction
Input: target definitions  
Programs:
- AlphaFold3.predict
- RoseTTAFold3.predict
- Boltz2.predict  
Output:
- `raw_models/*.cif`
- `raw_confidence/*`

### step_2_normalize_mmcif
Input: `raw_models/*.cif`, CCD dictionary  
Programs:
- Gemmi.parse_validate
- Gemmi.normalize_edit
- custom atom-mapper (topology+geometry)  
Output:
- `normalized/*.cif`
- `mapping_logs/*`

### step_3_privateer_view_prep
Input: `normalized/*.cif`  
Programs:
- custom glycan-expander/renamer → ensure **CCD monosaccharides** for each residue  
Output:
- `privateer_input/*.cif`

### step_4_refinement (MANDATORY)
Input: `normalized/*.cif` (or protonated complex if required)  
Programs:
- PLACER.multi_mode_ensemble  
Output:
- `placer_models/*.pdb|*.cif`
- `placer_scores/*.csv|*.json`

### step_5_validation
Input:
- `privateer_input/*.cif`
- `placer_models/best_model.*`  
Programs:
- Privateer (structure-only or map-assisted if map exists)
- PoseBusters (complex pose checks)  
Output:
- `reports/privateer/*`
- `reports/posebusters/*`

### step_6_protonation_and_conversion
Input: selected refined model  
Programs:
- Reduce.add_optimize_hydrogens (protein-focused H network; verify Cu-coordinating His)
- OpenBabel convert/protonate/charges (for ProLIF: bond orders + partial charges)  
Output:
- `analysis_input/complex_H.pdb`
- `analysis_input/ligand.mol2`

### step_7_analysis
Input: `analysis_input/complex_H.pdb`, optional trajectory  
Programs:
- MDAnalysis (single structure or trajectory)
- ProLIF (single pose or trajectory)  
Output:
- `reports/mdanalysis/*.csv|*.json`
- `reports/prolif/*.csv|*.json|*.html`

### step_8_reporting
Input: all partial reports  
Programs:
- custom report assembler  
Output:
- `summary.json`
- `metrics.csv`
- `report.html`

---

## Analysis-plan specifics Copilot must implement (key points)

### Separation of concerns: tuning vs main analysis (mandatory design rule)

Implement **two distinct orchestration paths**:

1) **Tuning phase (parameter sweeps)**
- `analysis/` must generate *batches* for parameter sweeps (grid/random as specified).
- A tuning run produces:
  - a machine-readable job list (e.g. JSON/YAML) describing each sweep point
  - per-job outputs in a dedicated tuning output tree
  - a summarizer that aggregates QC + analysis metrics and selects best parameters
- The sweep runner must call **prediction runners** from `structure_pipeline/` (AF3/RF3/Boltz-2),
  then run normalization → PLACER → QC → light analysis as required by tuning.

2) **Main analysis pipeline**
- Uses the chosen tuned parameters (locked) and runs DEL A / DEL B.
- Must not search/retune parameters during main analysis (“anti p-hacking” rule already below).
- Organize outputs separately from tuning runs.

**Do not mix tuning logic and main-analysis logic in the same script without a clear subcommand boundary.**
Preferred approach: separate CLI subcommands and separate config namespaces.

### Tuning phase (parameter sweeps)
Implement scripts to run tuning on a defined tuning subset and compare settings using:
- PoseBusters: pass-rate + error types
- Privateer: ring/anomer/stereo validity
- MDAnalysis “light”: align on core, measure Cu–C1 and Cu–C4 (flag if within 7 Å), pocket RMSD vs crystal when available
- ProLIF IFP + HDBSCAN: outlier-rate, #clusters, occupancy stability
- Optional: compare IFP similarity vs crystal (moderate threshold to be decided empirically)

AF3 tuning:
- Refinement: `num_recycles ∈ {10,15,20}`, keep `num_diffusion_samples=5`, `num_seeds=10`
- Diversity: `num_seeds ∈ {20,50,100}` with chosen `num_recycles`

RF3 tuning:
- Refinement: `n_recycles ∈ {10,20,30}`, no early stop; keep defaults for diffusion batch size & steps; `one_model_per_file=TRUE`
- Diversity: `seed ∈ {42,80,120}`

Boltz2 tuning:
- Always: `use_potentials=False`, `step_scale=1.638`, do not tune affinity
- Refinement 1: `recycling_steps ∈ {3,6,10}` (hold `diffusion_samples=5`, `sampling_steps=200`)
- Refinement 2: `sampling_steps ∈ {200,400,600}` (hold chosen recycling, diffusion=5)
- Diversity: `diffusion_samples ∈ {3,10,15}` (hold chosen recycling + sampling)

Decision rule for picking best parameters:
- prioritize: higher QC pass (PoseBusters/Privateer) + lower outlier-rate + stable clusters
- crystal similarity is a sanity check: avoid systematic degradation

(See analysis plan for full details.) :contentReference[oaicite:3]{index=3}

### Main analysis: DEL A (all proteins)
All proteins × 21 ligands × 3 models, using tuned parameters.
A2: hard QC; drop hard-fails, keep soft flags as metadata (do not mix as “pass”).
A3: run PLACER and re-run PoseBusters post-PLACER.
A4: ProLIF IFP per pose (residue-level + interaction type).
A5: Cluster within run to compress sampling noise; then cross-run clustering per (protein, ligand) using cluster representatives (medoid) or aggregate IFP.
A6: For each cluster: pick medoid; compute geometry metrics with MDAnalysis:
- Cu–C1, Cu–C4
- angle/tilt vs His-brace (define robust rule)
- optional planar stacking metrics
A7: Crystal anchoring where available:
- if ligand-bound crystal: IFP similarity
- if apo: pocket RMSD + key residues near Cu
A8: Connect to activity (C1/C4, substrate specificity, family) with aggregated cluster occupancies and signatures.

### DEL B (full-length proteins with CBM)
Run parallel to DEL A but with CBM-specific analysis:
- MDAnalysis: distances CBM–ligand and CBM–Cu (and dual proximity criterion)
- ProLIF: tag interactions from LPMO vs CBM residues; produce IFP subsets
- Compare domain-only vs full-length with paired statistics / mixed models (avoid treating as independent)

### Predictive modeling (later stage)
Two modeling targets:
1) enzyme → activity class (C1/C4/mixed) using:
   - occupancy of geometric C1-like/C4-like clusters
   - signature interaction frequencies
   - CBM proximity metrics for CBM enzymes
   Use grouped + stratified CV, class weights, balanced accuracy.
2) (enzyme, ligand) → binding mode/orientation using DP, polymer type, IFP+geometry.

---

## Implementation expectations

### Language & structure
Prefer **Python** for pipeline orchestration and analysis, plus small **bash** entrypoints.
Organize into:
- `src/` Python package
- `scripts/` CLI entrypoints
- `configs/` YAML/JSON for parameter settings and datasets
- `data/` (optional) small test fixtures only (no large binaries)
- `reports/` outputs
- `docs/` pipeline + decisions

### CLI
Provide a top-level CLI (e.g. `python -m lpmo_pipeline ...`) with subcommands:
- `predict`
- `normalize`
- `privateer-prep`
- `placer`
- `validate`
- `prep-analysis`
- `analyze`
- `report`
- `tune` (parameter sweep runner + summarizer)

### Reproducibility & logging
- Every run writes:
  - `run_manifest.json` with parameters, tool versions, timestamps, input hashes
  - structured logs
- Use deterministic seeds where applicable.

### Tool adapters
Abstract external tools behind adapters:
- `tools/af3.py`, `tools/rf3.py`, `tools/boltz2.py`
- `tools/gemmi.py`, `tools/privateer.py`, `tools/placer.py`, `tools/posebusters.py`
- `tools/reduce.py`, `tools/openbabel.py`
Adapters should:
- validate inputs/outputs
- capture stdout/stderr
- write machine-readable metadata (JSON)

### Data integrity checks
- Always keep original raw outputs immutable.
- Normalization should be idempotent.
- mmCIF master: derivations must cite source file + transformation in metadata.

### Cluster policy (anti p-hacking)
- Lock HDBSCAN hyperparameters from tuning before main analysis.
- Do not “search” cluster params during main runs.

---

## Output contract (must be stable)

- `summary.json`: top-level summary of dataset, settings, QC stats, cluster stats, geometry stats
- `metrics.csv`: flat table for downstream plotting/statistics (one row per pose/cluster/condition as relevant)
- `report.html`: human-readable report linking to key plots/tables

---

## Notes on common pitfalls

- PoseBusters complex mode often needs ligand `HETATM` + `CONECT` for robust checks.
- ProLIF requires correct bond orders/atom types/charges → prefer MOL2/SDF via OpenBabel before fingerprints.
- Reduce may alter histidine flips; verify Cu-coordinating His after protonation.
- mmCIF bond/conn categories must be populated for downstream correctness; use Gemmi editing to ensure `_struct_conn` and CCD completeness.

---

## When unsure

Default behavior:
1) Follow analysis plan. :contentReference[oaicite:4]{index=4}  
2) Keep mmCIF canonical and preserve raw confidence values.
3) Make decisions explicit in code comments and in `docs/decisions.md`.