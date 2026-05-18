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
      - ~~RF3~~ (ekskludert fra analyse)
      - ~~Boltz-2~~ (ekskludert fra analyse)    - **Data output** (fully available):
      - `structure_pipeline/work_core` — domain-only construct predictions
      - `structure_pipeline/work_full_length` — full-length construct predictions    - Denne mappen inneholder modell-“runners”/adapters og skal være stedet som faktisk kaller eksterne verktøy.
  - `analysis/`
    - Analyse-/orkestreringskode og dokumentasjon. **Denne `copilot.md` ligger her.**
    - Koden her skal *ikke* re-implementere modellkall; den skal orkestrere pipeline og bygge jobber/batcher.
  - `pipe_test/`
    - Uryddig/eksperimentell. **Skal ikke brukes** (ikke importer, ikke les configs herfra).
- `eirik/conda/`
  - Conda environments (kjøremiljø), men ikke pipeline-logikk.

### Hard rule: where prediction code lives
All code that runs AF3 must be invoked via:
- `eirik/Masteroppgave/structure_pipeline/`

**Policy**
- `analysis/` may generate inputs/configs and schedule jobs, but must call the prediction runners/adapters from `structure_pipeline/`.
- Do not duplicate AF3 invocation logic inside `analysis/`.
- If `analysis/` needs a new capability (new flags, new output parsing), implement it in `structure_pipeline/` adapters and call it from `analysis/`.

---

**Model preference (VSCode Copilot)**
When generating code, prefer **Claude Opus 4.6** if available in the Copilot model selector. Otherwise use the best available reasoning-capable model.

---

## Source-of-truth priority (conflict policy)

1) **`AF3_LPMO_pipeline_detailed_plan.md`** is the **primary and governing source** (v1.0, updated 2026-04-21). Always defer to this document for pipeline design, stage order, output tables, and analysis decisions.
2) **Latest user comments in active thread** — override all documents when present.
3) `MASTERPLAN.md` and `IMPLEMENTATION_PLAYBOOK.md` — integrated secondary plans.
4) `plan_implementation_spec.txt` (legacy, archived) — historical context only.
5) `plan_analyse.txt` (legacy, archived) — historical context only.
6) Reasonable engineering defaults.

If you (Copilot) see ambiguous/contradicting instructions, **follow (2)** for user-stated decisions and **(1)** for pipeline design, and document the decision in code comments or `docs/decisions.md`.

---

## High-level goals (from analysis plan)

We need a reproducible pipeline that:
- Generates complex structures with **AF3 only** (RF3 og Boltz-2 er ekskludert)
- Runs **active-site proximity pre-QC** before chemistry validators
- Performs **hard QC** with **PoseBusters + Privateer**
- Creates **ProLIF interaction fingerprints (IFP)**
- Clusters binding modes with **HDBSCAN** (Jaccard/Tanimoto on binary IFP)
- Measures mechanistic geometry with **MDAnalysis**
- Produces combined reporting artifacts:
  - `summary.json`, `metrics.csv`, `report.html`

Main analysis is cluster-primary (Del A: all proteins; Del B: full-length + CBM subset).
Optional tuning workflows are run as post-analysis sensitivity/comparison.

---

## HARD RULES (must not be violated)

1) **Privateer input requires monosaccharides as valid CCD `comp_id`** (e.g. `NAG`, `BGC`, `MAN`) — not only custom oligomer codes.
2) **Atom names differ across models**; mapping must be **topology + geometry based**, not name-based.
3) **mmCIF is master format**. PDB/MOL2/SDF are tool-specific derivatives.

---

## Configuration Policy (mandatory)

- All values that can realistically change between runs, systems, environments, or analyses must live in config files under `configs/`.
- This includes thresholds, limits, numeric defaults, flags, model/runtime choices, filenames, folder names, absolute paths, external tool commands/containers, schema/config asset paths, selection strings, and other operational settings.
- Analysis code must resolve such values through the shared loader in `lpmo_pipeline.config`; do not duplicate path/threshold literals across scripts.
- If a value is used repeatedly inside one script or function, load it once locally from config and reuse the local variable.
- New scripts must follow the same rule. Backward-compatible fallbacks are acceptable only when they are clearly documented and do not replace config as the source of truth.

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
- AlphaFold3.predict  (`num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5`)
- ~~Boltz2~~ (ekskludert fra analyse)
- ~~RoseTTAFold3~~ (ekskludert fra analyse)  
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

### step_4_ifp_and_geometry
Input: `normalized/*.cif` (QC-passed poses, AF3 directly)
Programs:
- ProLIF (IFP — no Cu repositioning or virtual atoms)
- MDAnalysis (geometry metrics)
Output:
- `analysis/ifp/*.csv`
- `analysis/geometry/*.csv`

### step_5_clustering
Input: IFP vectors from step_4
Programs:
- HDBSCAN (Jaccard metric on binary IFP only)
Output:
- `analysis/clusters/*.json`
- Medoid pose references

### step_6_reporting
Input: all partial reports  
Programs:
- custom report assembler  
Output:
- `summary.json`
- `metrics.csv`
- `report.html`

---

## Analysis-plan specifics Copilot must implement (key points)

### Separation of concerns: main analysis vs optional tuning (mandatory design rule)

Implement **two distinct orchestration paths**:

1) **Main analysis pipeline (first priority)**
- Consumes precomputed prediction artifacts from `structure_pipeline/`.
- Runs DEL A / DEL B with cluster as primary analysis unit.
- Must not collapse cluster rows to enzyme rows for primary descriptive/predictive analyses.
- Must not search/retune parameters during the same analysis run.

2) **Optional post-analysis tuning (if time/resources)**
- Generates parameter sweep batches and summaries in a dedicated output tree.
- Used for comparison/sensitivity and later reruns, not as a blocker for first analysis delivery.

**Do not mix tuning logic and main-analysis logic in the same script without a clear subcommand boundary.**
Preferred approach: separate CLI subcommands and separate config namespaces.

### Optional tuning phase (parameter sweeps)
Implement scripts to run tuning on a defined tuning subset and compare settings using:
- PoseBusters: pass-rate + error types
- Privateer: ring/anomer/stereo validity
- MDAnalysis “light”: align on core, measure Cu–C1 and Cu–C4 (flag if within 7 Å), pocket RMSD vs crystal when available
- ProLIF IFP + HDBSCAN: outlier-rate, #clusters, occupancy stability
- Optional: compare IFP similarity vs crystal (moderate threshold to be decided empirically)

AF3 tuning:
- **Valgte (faste) parametre for hovedanalyse: `num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5`**
- Tuning grid (utforsket, ikke aktiv): `num_recycles ∈ {10,15,20}`, keep `num_diffusion_samples=5`, `num_seeds=10`
- Diversity sweep (utforsket, ikke aktiv): `num_seeds ∈ {20,50,100}` with chosen `num_recycles`

~~RF3 tuning~~ (RF3 er ekskludert — ikke relevant):
- ~~Refinement: `n_recycles ∈ {10,20,30}`, no early stop; keep defaults for diffusion batch size & steps; `one_model_per_file=TRUE`~~
- ~~Diversity: `seed ∈ {42,80,120}`~~

~~Boltz2 tuning~~ (Boltz-2 er ekskludert — ikke relevant):
- ~~Always: `use_potentials=False`, `step_scale=1.638`, do not tune affinity~~
- ~~Refinement 1: `recycling_steps ∈ {3,6,10}` (hold `diffusion_samples=5`, `sampling_steps=200`)~~
- ~~Refinement 2: `sampling_steps ∈ {200,400,600}` (hold chosen recycling, diffusion=5)~~
- ~~Diversity: `diffusion_samples ∈ {3,10,15}` (hold chosen recycling + sampling)~~

Decision rule for picking best parameters:
- prioritize: higher QC pass (PoseBusters/Privateer) + lower outlier-rate + stable clusters
- crystal similarity is a sanity check: avoid systematic degradation

(See analysis plan for full details.) :contentReference[oaicite:3]{index=3}

### Main analysis: DEL A (all proteins)
All proteins × 21 ligands × AF3 only, using fixed parameters (`num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5`).
A2: active-site proximity pre-QC followed by hard QC; drop hard-fails, keep soft flags as metadata.
A3: ProLIF IFP per pose (residue-level + interaction type).
A4: Cluster within run to compress sampling noise; then cross-run clustering per (protein, ligand) using cluster representatives (medoid) or aggregate IFP.
A5: For each cluster: pick medoid; compute geometry metrics with MDAnalysis:
- Cu–C1, Cu–C4
- angle/tilt vs His-brace (define robust rule)
- optional planar stacking metrics
A6: Crystal anchoring where available:
- if ligand-bound crystal: IFP similarity
- if apo: pocket RMSD + key residues near Cu
A8: Connect to activity (C1/C4, substrate specificity, family) with aggregated cluster occupancies and signatures.

Primary modeling/data table remains cluster-level; enzyme-level aggregation is secondary only.

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
Prefer **Python** for pipeline orchestration and tool adapters, plus small **bash** entrypoints.
Prefer **R** for descriptive/predictive statistical analysis where practical.
Organize into:
- `src/` Python package
- `scripts/` CLI entrypoints
- `configs/` YAML/JSON for parameter settings and datasets
- `lpmo_pipeline.config` as the common path/config resolution mechanism
- `data/` (optional) small test fixtures only (no large binaries)
- `reports/` outputs
- `docs/` pipeline + decisions

### CLI
Provide a top-level CLI (e.g. `python -m lpmo_pipeline ...`) with subcommands:
- `predict`
- `normalize`
- `privateer-prep`
- `pre-qc-proximity`
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
- `tools/af3.py` (ikke `tools/rf3.py` eller `tools/boltz2.py` — RF3 og Boltz-2 er ekskludert)
- `tools/gemmi.py`, `tools/privateer.py`, `tools/posebusters.py`
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
- Lock HDBSCAN hyperparameters for each analysis run (and document origin if imported from optional tuning).
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

## AI Tool Use Preferences

**Preferanse: bruk innebygde lese-/søkefunksjoner fremfor shell-kommandoer for enkle operasjoner.**

Når AI trenger å sjekke enkle ting, foretrekk AIs innebygde verktøy (f.eks. `read_file`, `file_search`, `list_dir`, `grep_search`) fremfor å kjøre `ls`, `cat`, `head`, `grep` o.l. i terminalen.

Eksempler der innebygde verktøy er å foretrekke:
- Sjekke om en fil eksisterer → bruk `file_search` / `list_dir`
- Lese innholdet i en fil → bruk `read_file`
- Søke etter tekst i filer → bruk `grep_search`

Dette er en preferanse, ikke et absolutt krav. Bruk terminal når transformasjon, kjøring eller operasjoner utenfor workspace kreves.

---

## AI Decision Boundaries — Decisions Requiring User Approval

**The following changes must NOT be made by AI without explicit user approval:**

1. **Changes to pipeline order** — sequence of analysis steps (step_1 → step_8)
2. **Changes to which analyses are run** — enabling/disabling entire sub-analyses (DEL A, DEL B, tuning phases)
3. **Changes to primary analysis unit** — current primary unit is **cluster** (not enzyme)
4. **Changes to output artifacts** — schema or structure of required output tables
5. **Changes to core validation gates** — PoseBusters, Privateer mandatory status
6. **Changes to HARD RULES** — including chain schema, atom mapping strategy, confidence field handling

**Before proposing any of the above**, AI must:
- Ask for explicit confirmation via code comment with reasoning
- Wait for user response in the active thread
- Document the decision in `docs/decisions.md` with user approval reference

**ai may autonomously make:**
- Bug fixes to existing code
- Implementation of approved logic in new modules
- Parameter tuning within pre-approved ranges
- Performance optimizations that preserve correctness
- Test additions and documentation improvements

---

## When unsure

Default behavior:
1) Follow analysis plan. :contentReference[oaicite:4]{index=4}  
2) Keep mmCIF canonical and preserve raw confidence values.
3) Make decisions explicit in code comments and in `docs/decisions.md`.
4) If a decision falls on the boundary between "autonomous" and "requires approval", ask the user.