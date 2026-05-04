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
| 4 | `io/ccd_lookup.py` | ✅ Verifisert via test_io_contracts.sh. CCD-gate + rapportering fungerer på ekte AF3-data |
| 5–6 | `mapping/` | ✅ Verifisert via test_mapping_contracts.sh. Coverage=100%, round-trip rename validert |
| 7 | `io/protonate_export.py` | ✅ Verifisert 2026-04-23. PDBFixer/OpenMM primær-backend, ingen stille kopi-fallback. Name-based glycan-linking (C1→O4 for NAG/BGC/GLC) før protonering, og eksplisitt CONECT rewrite for komplett ligand-konnektivitet i `complex_H.pdb`. Verifisert via sbatch (572928, 572982, 572996, 573050). |
| 7b | `io/cif_to_pdb.py` | ✅ Verifisert 2026-04-23. PDBFixer-backend er primær (AF3-riktig), gemmi fallback med `backend_fallback_reason`. Verifisert via test_protonation_contracts.sh (backend=pdbfixer). |
| 8 | `qc/active_site_proximity.py` | ⚠️ Integrert og real-data-verifisert i hard-QC/analysis-core (2026-05-03), men endelig pose-manifestflate for metrikker fra droppede poser er fortsatt uavklart |
| 9–12 | `qc/` (PoseBusters, Privateer, Cu-His, QC report) | ✅ Reell 3-pose hard-QC-kjøring er verifisert 2026-05-03 via `tests/run_tests_scripts/test_hard_qc_real_cifs.sh` (jobb 612234, run `hard_qc_real_cifs_612234`). PoseBusters auto-splitter kombinerte AF3 protein+glykan-input til ligand+protein og kjører i `dock`-modus. Resultat: 1 pass, 1 soft_flag, 1 hard_fail; eneste gjenstående PB-hardfail er `minimum_distance_to_protein`, og samme STA4-case har fortsatt en reell Privateer anomer-feil. |
| 13 | `analysis/prolif_ifp.py` + `pose_ifp_table.tsv` | ✅ Standalone slice verifisert på ekte data 2026-05-04 via `test_prolif_real_cifs.sh` (jobb 617120); produksjonsintegrasjon gjenstår |
| 13b | `analysis/residue_contact_extraction.py` + `pose_residue_contact_table.tsv` | ❌ Ikke startet (ny i v1.0) |
| 14 | `analysis/mdanalysis_metrics.py` | ⚠️ Implementert, koblet inn i analysis-core produksjonssti, og verifisert på ekte data; RMSD-felter og videre geometry-hardening gjenstår |
| 14b | `analysis/convergence_metrics.py` | ❌ Ikke startet (ny i v1.0) |
| 15 | `analysis/clustering_hdbscan.py` | ⚠️ Refaktorert 2026-05-04 til AF3-only condition-wise HDBSCAN med eksakte Jaccard-medoider og focused pytest; real-data verifikasjon og produksjonsintegrasjon gjenstår |
| 16 | `analysis/cluster_signatures.py` | ⚠️ Kode finnes, ikke verifisert på ekte data |
| 16b | `analysis/residue_importance.py` | ❌ Ikke startet (ny i v1.0) |
| 17 | `placer/run_placer.py` | ❌ FJERNET — PLACER er fjernet fra analysen (beslutning 2026-04-21) |
| 18–25 | Analyse, aktivitetsmapping, modeller, rapportering | ⚠️ Delvis påbegynt. `report/` + `cli.py` analysis-core produksjonssti er verifisert på ekte data (jobb 613251), men ProLIF/clustering/crystal anchoring og senere analyser er ikke ende-til-ende integrert |

---

## Implementasjonssteg

1. ✅ **Installer avhengigheter** — `pip install -e ".[dev]"` fra pyproject.toml.
   Verifiser at gemmi, posebusters, hdbscan, prolif, mdanalysis importeres.
    Status 2026-05-04: tidlig funksjonell smoke-test er lagt til via
    `tests/run_tests_scripts/test_analysis_tools_smoke.sh` +
    `tests/run_tests_scripts/run_analysis_tools_smoke.py`.
    Denne sjekker i `analyse_env` at:
    - gemmi, MDAnalysis, ProLIF og HDBSCAN ikke bare importeres, men bestar en liten funksjonell probe
    - `apptainer` finnes på PATH
    - PoseBusters SIF finnes og er kjoreklar
    - Privateer-versjon kan hentes via eksisterende SIF-helper

2. ✅ **`io/mmcif_ingest.py`** — Implementer `ingest()` med ekte gemmi-parsing.
   ⛔ STOPP: Kjør på én AF3-CIF og inspiser resultatet manuelt.

3. ✅ **`io/normalize_mmcif.py`** — Chain-rename (A/B-D/E), confidence-ekstraksjon.
   ⛔ STOPP: Verifiser `normalized.cif` med `gemmi validate`.

4. ✅ **`io/ccd_lookup.py`** — CCD-cache + monosakkarid-validering.
    Status 2026-04-22: Verifisert via test_io_contracts.sh på ekte AF3-data. Hard fail ved ugyldig glykan-CCD fungerer korrekt.

5. ✅ **`mapping/cross_model_atom_mapping.py`** — 3-tier atom-matching.
   Status 2026-04-22: Verifisert via test_mapping_contracts.sh. Coverage=100%, atom_map.tsv generert korrekt.

6. ✅ **`mapping/rename_atoms.py`** — Last og anvend atom_map.
   Status 2026-04-22: Verifisert via test_mapping_contracts.sh. Round-trip rename (forward/reverse) validert.

7. ✅ **`io/protonate_export.py`** — Protonering og eksport til QC/analyse-artefakter.
    Status 2026-04-23: verifisert via sbatch-kontrakter (572928, 572982, 572996, 573050).
    Protonerings-backend-prioritet:
      1. PDBFixer/OpenMM (primær) — `addMissingHydrogens(pH=7.0)`, preserverer chain IDs
      2. reduce (sekundær, AmberTools)
      3. obabel (tertiær)
      4. Blocker i rapport — INGEN stille kopi-fallback lenger
    Ligandekstraksjon validerer at glycan-kjeder ikke er protein-residuenames.
    Rapport-felt: `complex_h_backend` erstatter `used_reduce`/`used_obabel_for_complex_h`.
    For glykaner (NAG/BGC/GLC) legges C1(i)→O4(i+1)-koblinger til før protonering,
    og `complex_H.pdb` får eksplisitte CONECT-linjer skrevet fra komplett bond-graf
    (inkludert name-based intra-residue glycan-bonds + inter-residue C1→O4).
    Verifisert resultat:
    - `complex_h_backend = pdbfixer`
    - `complex_H.pdb` atom-count > `for_posebusters.pdb`
    - `ligand_for_prolif.mol2` har `@<TRIPOS>ATOM` og `@<TRIPOS>BOND`
    - `ligand_for_prolif.mol2` inneholder ikke protein-residuenavn

7b. ✅ **`io/cif_to_pdb.py`** — CIF→PDB konvertering.
    Status 2026-04-23: verifisert med PDBFixer som primær-backend.
     Backend-valg: `_convert_auto()` → prøver PDBFixer, faller tilbake til gemmi med grunn.
     `CIFToPDBReport` har nytt felt `backend_fallback_reason` som alltid settes.
     Manuelt verifisert i PyMOL (2026-04-22): rå AF3 vs normalized vs PDB er identiske.
    Sbatch-verifisert: `cif_to_pdb_report.json` rapporterer `backend=pdbfixer`.

8. **`qc/active_site_proximity.py`** — pre-QC gate (ligand nær aktivt sete).
   Krav: beregn og logg `min_cu_ligand_distance`, `min_cu_c1`, `min_cu_c4` for alle poser.
    Status 2026-04-23: integrert i `hard_qc_orchestrator.py`. Targeted tester bekrefter at poser med stor avstand droppes for PB/Privateer, og at pre-QC-metrikker beholdes i verdict-data.
    Status 2026-05-03: verifisert på ekte data i den nye analysis-core produksjonsstien via `tests/run_tests_scripts/test_analysis_core_real_cifs.sh` (jobb 613251, run `analysis_core_real_cifs_613251`). I denne kjøringen ble STA4-caset droppet før downstream geometri (`analysis_status = skipped_dropped`), STA6-caset beholdt som `flagged` og analysert videre, og pre-QC-metrikker ble bevart i `qc_report.json` og `analysis_core_summary.json`.
    ⛔ STOPP: Avklar om `qc_report.json` + `analysis_core_summary.json` er den endelige flaten for metrikker fra droppede poser, eller implementer en eksplisitt `pose_manifest.tsv` som beholder disse uten å slippe droppede poser inn i downstream geometri/IFP.

9. **`qc/posebusters_runner.py`** — PoseBusters-wrapper.
   Kjøres på AF3-poser direkte.
    Merk: "fallback-kandidater" betyr kun valg mellom to identiske SIF-paths
    (ikke at PoseBusters-gaten hoppes over).
     SIF (fallback-kandidater):
    `/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif`,
    `/cluster/projects/nn1003k/prog/posebusters/posebusters.sif`.
    Status 2026-05-03: wrapperen er verifisert både med targeted pytest og i
    full hard-QC real-run. Kombinerte AF3 protein+glykan `for_posebusters.pdb`
    auto-splittes til ligand + protein før PoseBusters, og real-run
    `hard_qc_real_cifs_612234` kjørte `dock`-modus for alle 3 caser. De tidligere
    falske real-case hard-failene `all_atoms_connected` og
    `internal_steric_clash` er borte; eneste gjenstående PB-hardfail på ekte data
    er `minimum_distance_to_protein` i STA4-caset.
   Test: `pytest tests/test_qc_gates.py::TestPoseBustersGate`.

10. **`qc/privateer_runner.py`** — Privateer-wrapper + 100%-recog gate.
    Status 2026-04-30: wrapperen kjøres via SIF (`/cluster/projects/nn1003k/prog/privateer/privateer.sif`) med `apptainer run --cleanenv` og eksplisitte bind mounts. Avklart parse-kilde er `validation_data-privateer` fra `-mode ccp4i2`, ikke JSON stdout. Dry-run, batch-kjøring, SIF-basert versjonsdeteksjon, og filtrert artefakt-retensjon er implementert. Targeted pytest passer i `analyse_env`. Verifikasjon i full QC på ekte poser gjenstår.
    Test: `pytest tests/test_qc_gates.py::TestPrivateerGate`.

11. **`qc/custom_geometry_checks.py`** — Cu-His 1.9–2.6 Å gate.
    Status 2026-04-23: gate-logikk og verdict-integrasjon dekkes av targeted pytest og passer i `analyse_env`.
    Test: `pytest tests/test_qc_gates.py::TestCuHisGate`.

12. **`qc/qc_report.py`** — Samle pre-QC + PB + Privateer + geometri → verdict.
    Status 2026-05-03: verdict-aggregasjon er verifisert i full 3-pose real-run
    via `tests/run_tests_scripts/test_hard_qc_real_cifs.sh`. Run
    `hard_qc_real_cifs_612234` skrev schema-valid `qc_report.json` med
    sluttstatus 1 `pass`, 1 `soft_flag`, 1 `hard_fail`.
    ✅ STOPP-PUNKT verifisert 2026-05-03: full QC på 3 poser kjørt, og JSON-rapport
    validert mot skjema.

13. **`analysis/prolif_ifp.py`** — IFP-beregning med ProLIF.
    Bruker AF3-poser direkte.
    Output: `pose_ifp_table.tsv` med IFP-vektor, feature-navn og per-type kontakttelling.
        Status 2026-05-04: standalone ProLIF-slice er verifisert på ekte data via
        `tests/run_tests_scripts/test_prolif_real_cifs.sh` (jobb 617120, run `prolif_real_cifs_617120`).
        Implementert nå:
        - `ligand_for_prolif.mol2` lastes residue-aware slik at monosakkarider beholdes
            som egne ligandresiduer (`NAG1.B`, `GLC2.B`, ...) i stedet for å kollapse til `UNL1`
        - alle 9 interaksjonstyper eksponert av installert ProLIF-build brukes som
            låst standardsett (`HBDonor`, `HBAcceptor`, `Hydrophobic`, `PiStacking`,
            `Anionic`, `Cationic`, `CationPi`, `PiCation`, `VdWContact`)
        - `ifp_feature_names` skrives som `ligand_residue|protein_residue|interaction`
        - `pose_ifp_table.tsv` skriver utvidede per-type tellekolonner + `ifp_interaction_counts`
        - real-run `prolif_real_cifs_617120` ga 3/3 `normalize_ok`, 3/3 `protonate_ok`,
            3/3 ProLIF `status = ok`, og ingen kollapsede `UNL`-ligandetiketter
        Gjenstår:
        - koble samme slice inn i analysis-core / `lpmo-pipeline run`
        - implementere `pose_residue_contact_table.tsv`
        - oppdatere downstream signatur/parsing til det nye feature-navneskjemaet
        ✅ STOPP-PUNKT verifisert 2026-05-04: `pose_ifp_table.tsv` og `ifp_matrix.csv`
        skrevet med ligand-residue-oppløste features på ekte AF3-data (jobb 617120).

13b. **`analysis/residue_contact_extraction.py`** — Residue-nivå kontakttabell (ny i v1.0).
    Kjøres parallelt med IFP-generering.
    Output: `pose_residue_contact_table.tsv` (én rad per pose × residue × interaction_type).
    Kolonner: pose_id, protein_id, condition_id, residue_chain, residue_number, residue_name,
              interaction_type, contact_present, ligand_residue_label, distance_if_available,
              is_catalytic_surface_region, is_cbm_region, is_linker_region.
    ⛔ STOPP: Sjekk at residue-mapping er konsistent med IFP-features.

14. **`analysis/mdanalysis_metrics.py`** — 3D-metrikker per pose → `pose_geometry.tsv`.
    Status 2026-05-03: downstream geometri-grenen er nå implementert, koblet inn
    i analysis-core produksjonsstien og verifisert med focused pytest
    (`tests/test_mdanalysis_metrics.py`), dedikert real-data-harness
    (`tests/run_tests_scripts/test_geometry_real_cifs.sh` +
    `tests/run_tests_scripts/run_geometry_real_cifs.py`), og produksjonsnær
    real-data smoke run (`tests/run_tests_scripts/test_analysis_core_real_cifs.sh`,
    jobb 613251, run `analysis_core_real_cifs_613251`).
    Implementert nå:
    - `PoseGeometryMetrics.to_row()` + `write_pose_geometry_tsv()` skriver den planlagte `pose_geometry.tsv`-kontrakten.
    - `qc/custom_geometry_checks.check_geometry()` gjenbrukes for Cu-/His-brace-identifikasjon og `his_brace_angle_deg`.
    - Downstream `Cu_C1_distance` og `Cu_C4_distance` beregnes deretter pa nytt mot reposisjonert Cu i oxyl-modellen; de kopieres ikke fra hard-QC-feltene `min_cu_c1` / `min_cu_c4`.
    - `geometry_metrics.json` og `geometry_debug.pdb` finnes bare som test-/debugartefakter for manuell inspeksjon og skal ikke vaere standard output i produksjonspipen.
    - Virtuell H-plassering inkluderer naa kryss-residue glycosidiske naboer via `structure.connections` / `_struct_conn`, som fikset HC1-feilplasseringen i STA4-debugcaset.
    - `analysis/analysis_orchestrator.py` + `cli.py` skriver nå `pose_geometry.tsv` i produksjonsoutput for `passed` og `flagged` poser, mens droppede poser stoppes ved analysegrensen.
    - Den integrerte produksjonsstien skriver ikke `geometry_debug.pdb`.
    Gjenstar for full ferdigstillelse:
    - erstatte dagens `_struct_conn`-baserte naboheuristikk med full CCD-topologi for virtuell H-plassering
    - implementere RMSD-feltene som fortsatt er placeholders
    - utvide den samme produksjonsstien videre inn i ProLIF/IFP, clustering og crystal anchoring
    ✅ STOPP-PUNKT verifisert 2026-05-03: integrert analysis-core-kjøring skrev `pose_geometry.tsv` fra analysepipen uten å eksportere `geometry_debug.pdb` i produksjon (jobb 613251).

14b. **`analysis/convergence_metrics.py`** — Konvergensmetrikker per betingelse (ny i v1.0).
    Per pose: `ligand_rmsd_to_reference`, `convergent_flag` (RMSD < 2.0 Å anbefalt).
    Per betingelse: `convergence_fraction`, `median_ligand_rmsd`, `iqr_ligand_rmsd`.
    Referansepose: bruk QC-passing seed-1 pose under pre-clustering; oppdater til
    top-occupancy cluster medoid etter clustering for endelig rapportering.

15. **`analysis/clustering_hdbscan.py`** — HDBSCAN med Jaccard-metrikk på IFP-only.
    Kjør innen-modell clustering per (enzym, substrat, DP).
    Status 2026-05-04: modul refaktorert mot gjeldende AF3-only plan slik at den
    klynger én protein-ligand-betingelse om gangen, laster HDBSCAN-parametre fra
    `configs/thresholds.yaml`, handterer degenerert identisk IFP-matrise som én
    cluster i stedet for falsk all-noise, returnerer all-noise når antall poser er
    lavere enn `min_cluster_size`, og velger medoid med minimum summert Jaccard-avstand
    innen cluster. Nye hjelpefunksjoner bygger rå-rader for `cluster_assignments.tsv`,
    `medoid_manifest.tsv` og `condition_cluster_summary.tsv`. Focused pytest passer,
    men real-data harness og produksjonskobling gjenstår.
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
    Status 2026-05-03: `build_summary_json`, `build_metrics_csv` og `build_report_html` brukes nå fra analysis-core produksjonsstien. Real-data smoke run `analysis_core_real_cifs_613251` skrev `summary.json`, `metrics.csv` og `report.html`; det droppede STA4-caset ble ekskludert fra `metrics.csv`, mens `passed`/`flagged` poser ble beholdt.
    ⛔ STOPP: Når ProLIF/clustering er koblet inn, verifiser schema og HTML visuelt på nytt med de rikere feltene.

23. **`cli.py`** — Integrer alle steg. Kjør `lpmo-pipeline run --config ...`.
    Status 2026-05-03: `lpmo-pipeline run` er nå koblet til analysis-core produksjonsstien (`analysis/analysis_orchestrator.py`) for discovery → normalize → hard QC → downstream geometri → rapportering. Real-data smoke run via `tests/run_tests_scripts/test_analysis_core_real_cifs.sh` (jobb 613251) staged 3 CIFs gjennom discovery og produserte 1 `pass`, 1 `soft_flag`, 1 `hard_fail`; 2 poser ble analysert videre.
    ⛔ STOPP: Full end-to-end hovedanalyse gjenstår etter at ProLIF/IFP, clustering og crystal anchoring er integrert.

24. **`tuning/` (valgfritt etteranalyse)** — implementer sweep_runner + tune_orchestrator.
    Kjør sammenlignende tuning etter baseline leveranse, dokumenter eventuell parameteroppdatering.
