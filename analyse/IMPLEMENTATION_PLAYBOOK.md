# IMPLEMENTATION PLAYBOOK

Prioritert rekkefølge for implementasjon. Hvert steg kan testes isolert
før neste startes. **Stopp-punkt** = manuell verifikasjon før du går videre.

**⚠️ VIKTIG: Rekkefølgen på stegene og hvilke analyser som skal kjøres skal ikke endres av AI uten godkjenning fra deg først.**

**Status 2026-04-21: Steg 1–3 ferdig og verifisert på ekte AF3-data. PLACER er fjernet fra analysen helt (beslutning 2026-04-21). Ny primær kilde: `AF3_LPMO_pipeline_detailed_plan.md` (v1.0).**

Hvert steg implementeres, testes og verifiseres steg for steg.

Forutsetninger (oppdatert 21.04.2026):
- Strukturprediksjon (kun AF3) kjores separat fra analysekoden og produseres kun en gang per valgt parameteroppsett. (**RF3/RosettaFold 3 og Boltz-2 er ekskludert fra analyse.**)
- AF3 faste parametre (hovedanalyse): `num_recycles=10`, `num_seeds=15`, **`num_diffusion_samples=5`** → **75 poser per protein–ligand betingelse. AVKLART 2026-04-21: AF3-kjøringene er fullført med 5 diffusion samples. Se OPEN_QUESTIONS.md item 14.**
- Denne playbooken dekker analysepipen som konsumerer ferdige prediksjonsartifakter.
- Hovedanalyse kjores forst. Tuning er valgfri etteranalyse/sensitivitet.
- Primarenhet etter posegenerering er cluster (ikke tidlig aggregering til protein/enzym).
- Beregnede numeriske metrikker beholdes i tabeller/logg og skal ikke slettes ved gating.
- Bruk R der det er praktisk for beskrivende/prediktiv statistikk.
- PoseBusters kjøres via SIF-container: `/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif` (primær), `/cluster/projects/nn1003k/prog/posebusters/posebusters.sif` (fallback).
- Conda env: `/cluster/work/projects/nn1003k/eirik/conda/analyse_env/`.
- CIF→PDB konvertering bruker PDBFixer (OpenMM), adaptert fra PoseBench (MIT).
- Crystal anchoring alignment bruker PyMOL `pair_fit` på histidine-brace, Cu-nære residuer og substrat-recognition surface residues.

---

## Statusoversikt

| Steg | Modul | Status |
|------|-------|--------|
| 1 | `io/discovery.py` | ✅ Ferdig. Verifisert med 9 targets, 156K samples, af3_only + latest_only filter |
| 2 | `io/mmcif_ingest.py` | ✅ Ferdig. Verifisert på ekte AF3-CIF (CEL6, NAG6, STA8) |
| 3 | `io/normalize_mmcif.py` | ✅ Ferdig. Chain mapping A→A, C→B, B→E. Verifisert på alle 3 substrattyper |
| 4 | `io/ccd_lookup.py` | 🔄 Pågår. Wired inn i normalize med CCD-gate + rapportering; gjenstår verifikasjon på ekte AF3-data |
| 5–6 | `mapping/` | ⚠️ Kode finnes (`cross_model_atom_mapping.py`, `rename_atoms.py`), men mangler verifikasjon/kontraktstester |
| 7 | `io/protonate_export.py` | ❌ Broken (feil gemmi API), må omskrives |
| 7b | `io/cif_to_pdb.py` | ❌ Ikke startet |
| 8 | `qc/active_site_proximity.py` | ⚠️ Kode finnes, ikke testet på ekte data |
| 9–12 | `qc/` (PoseBusters, Privateer, Cu-His, QC report) | ⚠️ Kode finnes, ikke testet |
| 13 | `analysis/prolif_ifp.py` + `pose_ifp_table.tsv` | ⚠️ Kode finnes, ikke verifisert på ekte data |
| 13b | `analysis/residue_contact_extraction.py` + `pose_residue_contact_table.tsv` | ❌ Ikke startet (ny i v1.0) |
| 14 | `analysis/mdanalysis_metrics.py` | ⚠️ Kode finnes, ikke verifisert på ekte data |
| 14b | `analysis/convergence_metrics.py` | ❌ Ikke startet (ny i v1.0) |
| 15 | `analysis/clustering_hdbscan.py` | ⚠️ Kode finnes, ikke verifisert på ekte data |
| 16 | `analysis/cluster_signatures.py` | ⚠️ Kode finnes, ikke verifisert på ekte data |
| 16b | `analysis/residue_importance.py` | ❌ Ikke startet (ny i v1.0) |
| 17 | `placer/run_placer.py` | ❌ FJERNET — PLACER er fjernet fra analysen (beslutning 2026-04-21) |
| 18–25 | Analyse, aktivitetsmapping, modeller, rapportering | ⚠️ Delvis påbegynt (flere moduler finnes), ikke ende-til-ende verifisert |

---

## Implementasjonssteg

1. ✅ **Installer avhengigheter** — `pip install -e ".[dev]"` fra pyproject.toml.
   Verifiser at gemmi, posebusters, hdbscan, prolif, mdanalysis importeres.

2. ✅ **`io/mmcif_ingest.py`** — Implementer `ingest()` med ekte gemmi-parsing.
   ⛔ STOPP: Kjør på én AF3-CIF og inspiser resultatet manuelt.

3. ✅ **`io/normalize_mmcif.py`** — Chain-rename (A/B-D/E), confidence-ekstraksjon.
   ⛔ STOPP: Verifiser `normalized.cif` med `gemmi validate`.

4. **`io/ccd_lookup.py`** — CCD-cache + monosakkarid-validering.
    Status 2026-04-09: wired inn i `normalize_mmcif.py` med hard fail ved ugyldig glykan-CCD og status i `normalize_report.json`.
   Test: `pytest tests/test_io_contracts.py::TestCCDLookup -v`.

5. **`mapping/cross_model_atom_mapping.py`** — 3-tier atom-matching.
   ⛔ STOPP: Inspiser `atom_map.tsv` for én kjent NAG-residue. Coverage=100%?

6. **`mapping/rename_atoms.py`** — Last og anvend atom_map.
   Test: round-trip test med atom_map.tsv.

7. **`io/protonate_export.py`** — Reduce + OpenBabel for protonering.
   Forutsetning: `io/cif_to_pdb.py` omgjører mmCIF→PDB via PDBFixer først.
   ⛔ STOPP: Inspiser for_posebusters.pdb (CONECT) og .mol2 (bond orders).

7b. **`io/cif_to_pdb.py`** — CIF→PDB konvertering via PDBFixer (adaptert fra PoseBench MIT).
    Bruk PDBFixer for å håndtere AF3 mmCIF-format korrekt.
    ⛔ STOPP: Verifiser at output PDB har CONECT records og korrekt bondgraf.

8. **`qc/active_site_proximity.py`** — pre-QC gate (ligand nær aktivt sete).
   Krav: beregn og logg `min_cu_ligand_distance`, `min_cu_c1`, `min_cu_c4` for alle poser.
    ⛔ STOPP: Verifiser at poser med stor avstand droppes for PB/Privateer, men at metrikker beholdes i pose-nivå tabellene (`pose_manifest.tsv` + `pose_geometry.tsv`).

9. **`qc/posebusters_runner.py`** — PoseBusters-wrapper.
   Kjøres på AF3-poser direkte.
   SIF: `/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif`.
   Test: `pytest tests/test_qc_gates.py::TestPoseBustersGate`.

10. **`qc/privateer_runner.py`** — Privateer-wrapper + 100%-recog gate.
    Test: `pytest tests/test_qc_gates.py::TestPrivateerGate`.

11. **`qc/custom_geometry_checks.py`** — Cu-His 1.9–2.6 Å gate.
    Test: `pytest tests/test_qc_gates.py::TestCuHisGate`.

12. **`qc/qc_report.py`** — Samle pre-QC + PB + Privateer + geometri → verdict.
    ⛔ STOPP: Kjør full QC på 3 poser, verifiser JSON-rapport mot skjema.

13. **`analysis/prolif_ifp.py`** — IFP-beregning med ProLIF.
    Bruker AF3-poser direkte.
    Output: `pose_ifp_table.tsv` med IFP-vektor, feature-navn og per-type kontakttelling.
    ⛔ STOPP: Visualiser IFP-matrise for 5 poser.

13b. **`analysis/residue_contact_extraction.py`** — Residue-nivå kontakttabell (ny i v1.0).
    Kjøres parallelt med IFP-generering.
    Output: `pose_residue_contact_table.tsv` (én rad per pose × residue × interaction_type).
    Kolonner: pose_id, protein_id, condition_id, residue_chain, residue_number, residue_name,
              interaction_type, contact_present, ligand_residue_label, distance_if_available,
              is_catalytic_surface_region, is_cbm_region, is_linker_region.
    ⛔ STOPP: Sjekk at residue-mapping er konsistent med IFP-features.

14. **`analysis/mdanalysis_metrics.py`** — 3D-metrikker per pose → `pose_geometry.tsv`.
    ⛔ STOPP: Sjekk Cu-C1, Cu-C4, RMSD mot forventa størrelsesorden.

14b. **`analysis/convergence_metrics.py`** — Konvergensmetrikker per betingelse (ny i v1.0).
    Per pose: `ligand_rmsd_to_reference`, `convergent_flag` (RMSD < 2.0 Å anbefalt).
    Per betingelse: `convergence_fraction`, `median_ligand_rmsd`, `iqr_ligand_rmsd`.
    Referansepose: bruk QC-passing seed-1 pose under pre-clustering; oppdater til
    top-occupancy cluster medoid etter clustering for endelig rapportering.

15. **`analysis/clustering_hdbscan.py`** — HDBSCAN med Jaccard-metrikk på IFP-only.
    Kjør innen-modell clustering per (enzym, substrat, DP).
    Test: `pytest tests/test_clustering.py`.
    ⛔ STOPP: Inspiser klynge-output for ≥2 systemer.

16. **`analysis/cluster_signatures.py`** — clusterannotering med median/IQR.
    Krav: geometri brukes som annotering etter clustering, ikke som clusterinput.
    Output: `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`.
    Cluster type-etiketter: C1_compatible, C4_compatible, mixed_compatible, non_plausible, uncertain.
    Terskler defineres i `configs/thresholds.yaml` og låses før full analyse.

16b. **`analysis/residue_importance.py`** — Residueimportansanalyse (ny i v1.0).
    Steg 1 (within-protein): occupancy-weighted residue contact score per protein × betingelse.
      formula: residue_contact_score = sum(cluster_occupancy × residue_contact_frequency_in_cluster)
      Output: `protein_condition_residue_scores.tsv`
    Steg 2 (C1 vs C4 delta): delta-kontaktscore for C1-compatible vs C4-compatible clusters.
      Output: `protein_residue_regio_delta.tsv`
    Steg 3 (valgfritt, within-family): kartlegg residuer til alignmentkolonner per familie.
      Output: `family_aligned_residue_table.tsv`, `family_residue_enrichment.tsv`
    Steg 4 (cross-dataset, patch): aromatisk/polær/ladet kontakttetthet, loop-fraksjon, Cu-distanseskall.
      Output: `condition_patch_summary.tsv`, `protein_patch_summary.tsv`
    ⛔ STOPP: Verifiser at residue-score er konsistent med cluster_residue_signature.tsv.

17. ~~**`placer/run_placer.py`**~~ — **FJERNET. PLACER er fjernet fra analysen helt (beslutning 2026-04-21).** Steg 17 er avviklet. Se statusoversikt.

18. **`analysis/activity_mapping.py`** — bygg `predictive_cluster_table` fra cluster-rader.
    Krav: ingen primær aggregering cluster→protein. Bruk `protein_id` som kolonnenavn.

19. **`scripts/ec_activity_mapping.R`** — metadata EC# → aktivitet/substrat/regio.
    Krav: implementer eksplisitte regler for 1.14.99.53/54/55/56 og 1.14.99.- med AA17-spesialtilfelle.

20. **`analysis/predictive_models.py` + R scripts** — cluster-baserte modeller.
    Krav: grouped CV på enzymnivå, cluster-rader beholdes, avhengighet modelleres.

21. **`analysis/crystal_anchoring.py`, `cbm_variant.py`, `cbm_comparison.py`**
    som sekundæranalyser.
    Crystal anchoring bruker PyMOL `pair_fit` for optimal lokal superposisjon.
    ⛔ STOPP: Verifiser cluster_signatures/predictive_cluster_table mot skjema.

22. **`report/`** — build_summary_json, build_metrics_csv, build_report_html.
    Verifiser summary.json mot schema, metrics.csv kolonner, HTML visuelt.

23. **`cli.py`** — Integrer alle steg. Kjør `lpmo-pipeline run --config ...`.
    ⛔ STOPP: End-to-end test på 1 system (AF3). Sammenlign med forventet output.

24. **`tuning/` (valgfritt etteranalyse)** — implementer sweep_runner + tune_orchestrator.
    Kjør sammenlignende tuning etter baseline leveranse, dokumenter eventuell parameteroppdatering.
