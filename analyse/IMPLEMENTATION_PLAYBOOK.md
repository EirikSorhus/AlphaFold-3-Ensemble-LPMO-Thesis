# IMPLEMENTATION PLAYBOOK

Prioritert rekkefølge for implementasjon. Hvert steg kan testes isolert
før neste startes. **Stopp-punkt** = manuell verifikasjon før du går videre.

**⚠️ VIKTIG: Rekkefølgen på stegene og hvilke analyser som skal kjøres skal ikke endres av AI uten godkjenning fra deg først.**

**Status 2026-04-21: Steg 1–3 ferdig og verifisert på ekte AF3-data. PLACER er fjernet fra analysen helt (beslutning 2026-04-21). Ny primær kilde: `AF3_LPMO_pipeline_detailed_plan.md` (v1.0).**

Hvert steg implementeres, testes og verifiseres steg for steg.

Forutsetninger (oppdatert 21.04.2026):
- Strukturprediksjon (kun AF3) kjores separat fra analysekoden og produseres kun en gang per valgt parameteroppsett. (**RF3/RosettaFold 3 og Boltz-2 er ekskludert fra analyse.**)
- AF3 faste parametre (hovedanalyse): `num_recycles=10`, `num_seeds=15`, **`num_diffusion_samples=5`** → **75 poser per protein–ligand betingelse. AVKLART 2026-04-21: AF3-kjøringene er fullført med 5 diffusion samples. Se OPEN_QUESTIONS.md item 14.**
- **Data tilgjengelighet**: Precomputede AF3-strukturer er fullt tilgjengelige under:
  - Domain-only: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_core`
  - Full-length: `/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work_full_length`
- Denne playbooken dekker analysepipen som konsumerer ferdige prediksjonsartifakter.
- Hovedanalyse kjores forst. Tuning er valgfri etteranalyse/sensitivitet.
- Primarenhet etter posegenerering er cluster (ikke tidlig aggregering til protein/enzym).
- Beregnede numeriske metrikker beholdes i tabeller/logg og skal ikke slettes ved gating.
- Bruk R der det er praktisk for beskrivende/prediktiv statistikk.
- Endringsbare paths, verktøyvalg, terskler og andre runtime-innstillinger skal ligge i config-filer under `configs/` og lastes via `lpmo_pipeline.config`.
- Eksterne verktøy-paths og sentrale config/schema-referanser styres fra `configs/runtime_paths.yaml`.
- Conda env: `/cluster/work/projects/nn1003k/eirik/conda/analyse_env/`.
- CIF→PDB konvertering bruker PDBFixer (OpenMM), adaptert fra PoseBench (MIT).
- Crystal anchoring pocket-RMSD bruker per i dag lokal gemmi/numpy Kabsch-superposisjon på delte pocket C-alpha-atomer. Operativ pocket-definisjon er proteinrester innen 5 A fra ligand eller Cu; apo-krystaller kan få pocket via sekvensprojeksjon fra representant/medoid med normalisering av residunavn som `HIC -> HIS`. PyMOL `pair_fit`-paritet/hardening er eventuelt senere arbeid, ikke dagens operative backend.

Konfigurasjonsregel (gjeldende):
- Scripts skal ikke hardkode endringsbare tall, grenser, flags, filnavn, mappenavn, absolute paths eller verktøy-paths.
- Hvis samme config-verdi brukes flere ganger i én modul, kan den lastes én gang lokalt for lesbarhet/ytelse.
- Nye moduler skal enten bruke eksisterende config-fil eller utvide config-strukturen; ikke innfør nye spredte literals i kode.

---

## Statusoversikt

| Steg | Modul | Status |
|------|-------|--------|
| 1 | `io/discovery.py` | ✅ Ferdig. Verifisert med 9 targets, 156K samples, af3_only + latest_only filter |
| 2 | `io/mmcif_ingest.py` | ✅ Ferdig. Verifisert på ekte AF3-CIF (CEL6, NAG6, STA8) |
| 3 | `io/normalize_mmcif.py` | ✅ Ferdig. Chain mapping A→A, C→B, B→E. Verifisert på alle 3 substrattyper. Per 2026-05-19 filtreres identitetsremap bort, og `atom_map.tsv`/`rename_log.json` skrives bare ved uventet mapping/valideringsfeil for å redusere I/O. |
| 4 | `io/ccd_lookup.py` | ✅ Verifisert via test_io_contracts.sh. CCD-gate + rapportering fungerer på ekte AF3-data |
| 5–6 | `mapping/` | ✅ Verifisert via test_mapping_contracts.sh. Coverage=100%, round-trip rename validert |
| 7 | `io/protonate_export.py` | ✅ Verifisert 2026-04-23. PDBFixer/OpenMM primær-backend, ingen stille kopi-fallback. Name-based glycan-linking (C1→O4 for NAG/BGC/GLC) før protonering, og eksplisitt CONECT rewrite for komplett ligand-konnektivitet i `complex_H.pdb`. Verifisert via sbatch (572928, 572982, 572996, 573050). |
| 7b | `io/cif_to_pdb.py` | ✅ Verifisert 2026-04-23. PDBFixer-backend er primær (AF3-riktig), gemmi fallback med `backend_fallback_reason`. Verifisert via test_protonation_contracts.sh (backend=pdbfixer). |
| 8 | `qc/active_site_proximity.py` | ⚠️ Integrert og real-data-verifisert i hard-QC/analysis-core (2026-05-03), men endelig pose-manifestflate for metrikker fra droppede poser er fortsatt uavklart |
| 9–12 | `qc/` (PoseBusters, Privateer, Cu-His, QC report) | ✅ Reell 3-pose hard-QC-kjøring er verifisert 2026-05-03 via `tests/run_tests_scripts/test_hard_qc_real_cifs.sh` (jobb 612234, run `hard_qc_real_cifs_612234`). PoseBusters auto-splitter kombinerte AF3 protein+glykan-input til ligand+protein og kjører i `dock`-modus med bibliotekets standard `intermolecular_distance`-terskler (`max_distance=5.0 Å`, `search_distance=6.0 Å`). Per 2026-05-19 kan per-pose hard QC og Privateer-dispatch bruke `max_workers` fra `n_jobs`. Resultat: 1 pass, 1 soft_flag, 1 hard_fail; eneste gjenstående PB-hardfail er `minimum_distance_to_protein` (PoseBusters' renamed `no_clashes`-resultat), og samme STA4-case har fortsatt en reell Privateer anomer-feil. |
| 13 | `analysis/prolif_ifp.py` + `pose_ifp_table.tsv` | ⚠️ Standalone slice verifisert på ekte data 2026-05-04 via `test_prolif_real_cifs.sh` (jobb 617120), og koblet inn i analysis-core produksjonsstien 2026-05-04 med focused pytest. Per 2026-05-19 støtter batchen `max_workers` via `n_jobs`; ny real-data verifikasjon av den integrerte stien gjenstår |
| 13b | `analysis/residue_contact_extraction.py` + `pose_residue_contact_table.tsv` | ✅ Verifisert på ekte data 2026-05-04 via `test_residue_contact_real_cifs.sh` (jobb 619027); downstream bruk gjenstår |
| 14 | `analysis/mdanalysis_metrics.py` | ⚠️ Implementert, koblet inn i analysis-core produksjonssti, og verifisert på ekte data; RMSD-felter og videre geometry-hardening gjenstår |
| 14b | `analysis/convergence_metrics.py` | ⚠️ Standalone slice verifisert på ekte data 2026-05-04 via `test_convergence_real_cifs.sh` (jobb 619047), og koblet inn i analysis-core produksjonssti med focused pytest; ny integrert real-data verifikasjon gjenstår |
| 15 | `analysis/clustering_agglomerative.py` (`clustering_hdbscan.py` beholdt som sensitivitet) | ⚠️ Primær Stage 6-metode er låst 2026-05-20 til condition-wise agglomerative Jaccard med `linkage=average`, `distance_threshold=0.55` og `min_cluster_size=3`. HDBSCAN `min_cluster_size=3`, `min_samples=null` beholdes som sensitivitet. Focused clustering unit tests ble kjørt på nytt 2026-05-16; real-data clustering-output inspeksjon gjenstår før senere analysebygging |
| 16 | `analysis/cluster_signatures.py` | ✅ Koblet inn i analysis-core produksjonsstien og verifisert 2026-05-18 med focused pytest + real-data sbatch smoke (jobb 1105024). Produksjonsstien skriver `cluster_table.tsv`, `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, `cluster_signatures.json` og `cluster_annotation_stage_completed`. |
| 16b | `analysis/residue_importance.py` | ✅ Koblet inn i analysis-core produksjonsstien og verifisert 2026-05-18 med focused pytest + real-data sbatch smoke (jobb 1105024). 16b konsumerer Stage 16-signaturene, skriver eksplisitte nullrader ved observerte kontakter uten retained clusters, og bruker ikke loop-fraksjon. |
| 17 | `placer/run_placer.py` | ❌ FJERNET — PLACER er fjernet fra analysen (beslutning 2026-04-21) |
| 18–25 | Analyse, aktivitetsmapping, modeller, rapportering | ⚠️ Delvis påbegynt. `report/` + `cli.py` analysis-core produksjonssti skriver nå også ProLIF-, Stage 6-clustering-, Stage 7 `cluster_table`/signatur-, residue-importance-, condition/protein-summary-, crystal-anchoring- og implementerte pose/QC TSV-artefakter (`pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv`, `qc_attrition_table.tsv`, `cluster_table.tsv`, `condition_table.tsv`, `protein_summary_table.tsv`, `crystal_anchor_table.tsv`, `crystal_geometry_table.tsv`, `crystal_ifp_diagnostic_summary.tsv`). En bred integrert real-data kjøring som faktisk gir medoid-vs-crystal-sammenligninger og senere analyser gjenstår fortsatt. |

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
   Status 2026-05-19: identitetsremaps filtreres bort, CIF-remap hoppes over når ingen reell remap trengs, og rutinemessig vellykket normalisering skriver ikke lenger `atom_map.tsv`/`rename_log.json`. De filene beholdes som debug/failure-artefakter ved ufullstendig atom-mapping eller CCD-valideringsfeil.

4. ✅ **`io/ccd_lookup.py`** — CCD-cache + monosakkarid-validering.
    Status 2026-04-22: Verifisert via test_io_contracts.sh på ekte AF3-data. Hard fail ved ugyldig glykan-CCD fungerer korrekt.

5. ✅ **`mapping/cross_model_atom_mapping.py`** — 3-tier atom-matching.
   Status 2026-04-22: Verifisert via test_mapping_contracts.sh. Coverage=100%, atom_map.tsv generert korrekt.
   Status 2026-05-19: mapping-kontrakten er uendret, men per-pose `atom_map.tsv` er ikke lenger normal produksjonsoutput fra vellykket normalisering.

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
    SIF-kandidater styres fra `configs/runtime_paths.yaml`.
    Status 2026-05-03: wrapperen er verifisert både med targeted pytest og i
    full hard-QC real-run. Kombinerte AF3 protein+glykan `for_posebusters.pdb`
    auto-splittes til ligand + protein før PoseBusters, og real-run
    `hard_qc_real_cifs_612234` kjørte `dock`-modus for alle 3 caser. De tidligere
    falske real-case hard-failene `all_atoms_connected` og
    `internal_steric_clash` er borte; eneste gjenstående PB-hardfail på ekte data
    er `minimum_distance_to_protein` i STA4-caset. Dette er PoseBusters'
    renamed `no_clashes`-resultat, ikke `protein-ligand_maximum_distance`-testen.
   Test: `pytest tests/test_qc_gates.py::TestPoseBustersGate`.

10. **`qc/privateer_runner.py`** — Privateer-wrapper + 100%-recog gate.
    Status 2026-04-30: wrapperen kjøres via SIF definert i `configs/runtime_paths.yaml` med `apptainer run --cleanenv` og eksplisitte bind mounts. Avklart parse-kilde er `validation_data-privateer` fra `-mode ccp4i2`, ikke JSON stdout. Dry-run, batch-kjøring, SIF-basert versjonsdeteksjon, og filtrert artefakt-retensjon er implementert. Targeted pytest passer i `analyse_env`. Verifikasjon i full QC på ekte poser gjenstår.
    Status 2026-05-19: `run_privateer_batch` kan kjøre flere eligible poser parallelt via `max_workers`, koblet til production `n_jobs`. Oppstartskost for selve SIF-kjøringen betales fortsatt per Privateer-prosess.
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
        Status 2026-05-04 (senere): samme ProLIF-slice er nå koblet inn i
        `analysis/analysis_orchestrator.py` / `lpmo-pipeline run` for alle QC-pass/
        flagged poser. Produksjonsstien skriver `pose_ifp_table.tsv` og per-condition
        `ifp_matrix.csv`, og focused pytest dekker integrasjonen.
        Status 2026-05-19: `compute_ifp_batch` støtter `max_workers` via production
        `n_jobs`, slik at uavhengige poser kan beregnes i parallelle prosesser.
        Gjenstår:
        - kjøre den integrerte stien på ekte data og inspisere outputene
        - oppdatere downstream signatur/parsing til det nye feature-navneskjemaet
        ✅ STOPP-PUNKT verifisert 2026-05-04: `pose_ifp_table.tsv` og `ifp_matrix.csv`
        skrevet med ligand-residue-oppløste features på ekte AF3-data (jobb 617120).

13b. **`analysis/residue_contact_extraction.py`** — Residue-nivå kontakttabell (ny i v1.0).
    Kjøres parallelt med IFP-generering.
    Output: `pose_residue_contact_table.tsv` (én rad per pose × residue × interaction_type).
    Kolonner: pose_id, protein_id, condition_id, residue_chain, residue_number, residue_name,
              interaction_type, contact_present, ligand_residue_label, distance_if_available,
              is_catalytic_surface_region, is_cbm_region, is_linker_region.
    Status 2026-05-04: implementert direkte fra de ligand-residue-oppløste
    ProLIF-featurelabelene (`ligand_residue|protein_residue|interaction`) og
    koblet inn i `analysis/analysis_orchestrator.py` / `lpmo-pipeline run`.
    Focused pytest dekker extractor + integrasjon i produksjonsstien.
    Status 2026-05-04 (senere): verifisert på ekte data via
    `tests/run_tests_scripts/test_residue_contact_real_cifs.sh` (jobb 619027).
    Denne kjøringen ga 3/3 `status = ok`, 252 totale rader i
    `pose_residue_contact_table.tsv`, ingen kollapsede `UNL`-ligandetiketter,
    og eksakt samsvar mellom rekonstruerte residue-contact features og
    `pose_ifp_table.tsv` på alle tre caser.
    ✅ STOPP-PUNKT verifisert 2026-05-04: residue-mapping er konsistent med IFP-features på ekte AF3-data.

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
    - utvide den samme produksjonsstien videre inn i crystal anchoring
    Status 2026-05-04 (senere): den samme produksjonsstien kjører nå også
    protonering → ProLIF-batch → condition-wise clustering og skriver
    `pose_ifp_table.tsv`, `cluster_assignments.tsv`, `medoid_manifest.tsv` og
    `condition_cluster_summary.tsv`. Focused pytest dekker denne integrasjonen.
    ✅ STOPP-PUNKT verifisert 2026-05-03: integrert analysis-core-kjøring skrev `pose_geometry.tsv` fra analysepipen uten å eksportere `geometry_debug.pdb` i produksjon (jobb 613251).

14b. **`analysis/convergence_metrics.py`** — Konvergensmetrikker per betingelse (ny i v1.0).
    Per pose: `ligand_rmsd_to_reference`, `convergent_flag` (RMSD < 2.0 Å anbefalt).
    Per betingelse: `convergence_fraction`, `median_ligand_rmsd`, `iqr_ligand_rmsd`.
    Referansepose: bruk QC-passing seed-1 pose under pre-clustering; cluster-medoider
    brukes senere som representative strukturer for crystal anchoring, men
    konvergensmodulen er fortsatt deskriptiv og filtrerer ikke poser.
    Status 2026-05-04: modul implementert med protein-CA-alignering + ligand-heavy-atom RMSD,
    seed-1-først referansevalg, og TSV-outputene `pose_convergence.tsv` og
    `condition_convergence_summary.tsv`. Koblet inn i
    `analysis/analysis_orchestrator.py` / `lpmo-pipeline run`, og focused pytest
    dekker både modulkontrakt og produksjonssti-integrasjon.
    Status 2026-05-04 (senere): standalone real-data-verifikasjon er kjørt via
    `tests/run_tests_scripts/test_convergence_real_cifs.sh` (jobb 619047) på to
    flerpose-betingelser (`Q7SCE9_NAG4`, `Q59930_STA6`; 6 poser totalt).
    Resultat:
    - 6/6 poser `status = ok`
    - referansepose ble valgt som `seed-1_sample-0` i begge betingelser
    - `pose_convergence.tsv` skrev 6 rader og `condition_convergence_summary.tsv` skrev 2 rader
    - observerte ligand-RMSD-er var ikke trivielle (f.eks. ~6.73 Å og ~27.13 Å), så harnessen bekrefter at metrikken faktisk skiller mellom alternative positurer
    Gjenstår:
    - kjøre den integrerte stien på ekte data og inspisere konvergensoutputene
    - avklare/implementere endelig medoid-basert referanseoppdatering etter clustering for sluttrapportering

15. **`analysis/clustering_agglomerative.py`** — agglomerative Jaccard som primær Stage 6-metode.
    Kjør innen-modell clustering per protein-ligand-betingelse.
    Status 2026-05-20: clustering-piloten og parameter-sensitivitetskjøringen valgte
    agglomerative Jaccard som global primærmetode for full analyse med
    `linkage=average`, `distance_threshold=0.55` og `min_cluster_size=3`.
    Produksjonsstien i `analysis/analysis_orchestrator.py` er nå låst til denne
    metoden. `analysis/clustering_hdbscan.py` beholdes som sensitivitet med
    `min_cluster_size=3`, `min_samples=null` og `cluster_selection_method=eom`.
    Begge metodene bruker binær Jaccard og eksakte medoids med minimum summert
    within-cluster Jaccard-avstand. Nye hjelpefunksjoner bygger rå-rader for
    `cluster_assignments.tsv`, `medoid_manifest.tsv` og `condition_cluster_summary.tsv`.
    Focused pytest passer, og `tests/run_tests_scripts/run_analysis_core_real_cifs.py`
    er utvidet til å samle de nye artefaktene; selve real-data kjøringen gjenstår.
    Test: `pytest tests/test_clustering.py`.
    ⛔ STOPP: Inspiser klynge-output for ≥2 systemer.

16. **`analysis/cluster_signatures.py`** — clusterannotering med median/IQR.
    Krav: geometri brukes som annotering etter clustering, ikke som clusterinput.
    Output: `cluster_table.tsv`, `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, `cluster_signatures.json`.
    Cluster type-etiketter: C1_compatible, C4_compatible, mixed_compatible, non_plausible, uncertain.
    Terskler defineres i `configs/thresholds.yaml` og låses før full analyse.
    Status 2026-05-18: Stage 7-tabellbyggeren aggregerer eksisterende
    Stage 6-clustermedlemskap og medoids uten å endre clustering, parser
    ligand-residue-oppløste IFP-features (`ligand_residue|protein_residue|interaction`),
    beregner per-cluster IFP-/residue-contact-frekvenser og lagrer
    median/IQR-geometri i `cluster_signatures.json`. Koblet inn i
    `analysis/analysis_orchestrator.py` og CLI-manifestet som
    `cluster_annotation_stage_completed`.
    Status 2026-05-20: samme Stage 7-overflate skriver nå også `cluster_table.tsv`
    som en flat TSV-eksport av de samme retained non-noise `cluster_summaries`
    som ligger i `cluster_signatures.json`. Exporten gjenbruker eksisterende
    Stage 7-aggregater, innfører ingen ny clusterlogikk, og lager ikke
    syntetiske nullrader for no-cluster conditions.
    Status 2026-05-21: `cluster_table.tsv` er beriket bakoverkompatibelt med
    condition metadata, `cluster_size`, C1/C4
    computable/plausible/highly-plausible-fraksjoner, pose-confidence
    mean/median/medoid-felter og convergence mean/median/medoid-felter.
    Verifisering:
    `/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python -m pytest tests/test_cluster_signatures.py tests/test_residue_importance.py tests/test_analysis_orchestrator.py tests/test_cli_run.py -q`
    (10 passed), og `tests/run_tests_scripts/test_analysis_core_real_cifs.sh`
    via sbatch jobb 1105024. Real-data smoken skrev Stage 7-artefaktene og
    manifest-gaten; det korte 3-CIF utvalget hadde ingen faktiske clusters
    (1 analysert pose, ikke contact-eligible), så signaturtabellene var
    kontraktmessig tomme med headers.

16b. **`analysis/residue_importance.py`** — Residueimportansanalyse (ny i v1.0).
    Steg 1 (within-protein): occupancy-weighted residue contact score per protein × betingelse.
      formula: residue_contact_score = sum(cluster_occupancy × residue_contact_frequency_in_cluster)
      Output: `protein_condition_residue_scores.tsv`
    Steg 2 (C1 vs C4 delta): delta-kontaktscore for C1-compatible vs C4-compatible clusters.
      Output: `protein_residue_regio_delta.tsv`
    Steg 3 (valgfritt, within-family): kartlegg residuer til alignmentkolonner per familie.
      Output: `family_aligned_residue_table.tsv`, `family_residue_enrichment.tsv`
    Steg 4 (cross-dataset, patch): aromatisk/polær/ladet/hydrogenbinding-kontakttetthet og eksplisitte regionflagg der de finnes.
      Output: `condition_patch_summary.tsv`, `protein_patch_summary.tsv`
        Status 2026-05-18: Stage 16b konsumerer nå Stage 16-outputene
        `cluster_residue_signature.tsv`, `cluster_ifp_signature.tsv` og
        `cluster_signatures.json` som kontrakt, i stedet for å rekonstruere
        cluster-/residuegrunnlaget fra rå Stage 6-output. Modulen bygger alle fire
        tabellene, er koblet inn i `analysis/analysis_orchestrator.py`, og
        ekskluderer noise/outlier-clusters fordi Stage 16-signaturene bare skrives
        for retained non-noise clusters. For betingelser med observerte kontaktrester
        men ingen retained clusters skrives eksplisitte nullrader i residue-tabellene
        i stedet for kun header, slik at no-valid-cluster-tilfeller blir synlige
        downstream. Loop-fraksjon er fjernet; patch-summary bruker nå residueklasse-,
        hydrogenbinding- og eksplisitte regionflagg fra Stage 16-overflatene.
        Cluster-avhengige tester bruker den deterministiske fixture-mappen
        `tests/fixtures/clustering_stage_outputs/`, slik at ikke videre
        clusterlogikk testes på en én-pose smoke.
        Verifisering:
        `/cluster/work/projects/nn1003k/eirik/conda/analyse_env/bin/python -m pytest tests/test_cluster_signatures.py tests/test_residue_importance.py tests/test_analysis_orchestrator.py tests/test_cli_run.py -q`
        (10 passed), og `tests/run_tests_scripts/test_analysis_core_real_cifs.sh`
        via sbatch jobb 1105024. Real-data smoken skrev alle Stage 16/16b-artefakter;
        den hadde ingen faktiske clusters (`cluster_signatures.clusters=[]`), men
        skrev 6 eksplisitte nullrader i både `protein_condition_residue_scores.tsv`
        og `protein_residue_regio_delta.tsv` fra observerte kontakter.
        Status 2026-05-21 (senere): en ny retained-cluster-regresjon i
        `tests/test_residue_importance.py` verifiserer nå ikke-null,
        cluster-vektede Stage 16b-output mot den real-deriverte fixture-roten
        `tests/fixtures/clustering_stage_outputs/` for `A0A0S2GKZ1`, inkludert
        ikke-null `residue_contact_score`, ikke-null C1/C4-deltaer og eksakt
        samsvar med de innskrevne Stage 16b-tabellene. Family-grenen er nå
        implementert som en valgfri, separat postprosess gjennom
        `analysis/family_alignment.py`,
        `analysis/family_enrichment_postprocess.py` og CLI-kommandoen
        `lpmo-pipeline family-enrichment`. Den filtrerer til AA9/AA10
        `domain_only`, bruker den dedupliserte catalytic-core FASTA-en sammen
        med metadata for sekvensgruppe-resolusjon, og bygger family-alignments
        via precomputed aligned FASTA eller MAFFT uten å kobles inn i
        `run_analysis_core()`. Verifisering: `tests/test_family_alignment.py`,
        `tests/test_family_residue_enrichment.py`,
        `tests/test_family_enrichment_postprocess.py` og
        `tests/test_cli_run.py::test_cmd_family_enrichment_calls_postprocess`
        kjører grønt i `tests/run_tests_scripts/test_family_enrichment_validation.sh`,
        og sbatch jobb 1150831 prosesserte både AA9 og AA10 på de mergede
        staged pilot `domain_only`-shardene med `mafft_linsi`, som skrev 219
        rader i `family_aligned_residue_table.tsv` og 83 rader i
        `family_residue_enrichment.tsv`.

17. ~~**`placer/run_placer.py`**~~ — **FJERNET. PLACER er fjernet fra analysen helt (beslutning 2026-04-21).** Steg 17 er avviklet. Se statusoversikt.

18. **`analysis/activity_mapping.py`** — bygg `predictive_cluster_table` fra cluster-rader.
    Krav: ingen primær aggregering cluster→protein. Bruk `protein_id` som kolonnenavn.
    De reviderte aktivitets-prediktive planene
    `c1_c4_predictive_analysis_plan_simplified.yaml` og
    `substrate_activity_prediction_plan.yaml` spesifiserer nå hvilke
    condition-/cluster-/metadata-tabeller som må bygges før prediktiv modellering.
    Status 2026-05-21: modulen har nå en separat postprosess-bygger for
    `predictive_cluster_table.tsv` fra beriket `cluster_table.tsv`,
    proteinmetadata og EC/activity mapping. Koden bruker dagens Stage 7-felter
    (`Cu_C1_distance_median`, `oxyl_H_C1_distance_median`, osv.) i stedet for
    eldre `median_cu_c1`-navn, og er dekket av focused pytest.

19. **`scripts/ec_activity_mapping.R`** — metadata EC# → aktivitet/substrat/regio.
    Krav: implementer eksplisitte regler for 1.14.99.53/54/55/56 og 1.14.99.- med AA17-spesialtilfelle.
    Trenger ikke å være R, kan også være python-kode

20. **`analysis/predictive_models.py` + R scripts** — cluster-baserte modeller.
    Krav: grouped CV på enzymnivå, cluster-rader beholdes, avhengighet modelleres.
        Implementasjonen skal følge de reviderte, mindre detaljerte planene:
      - `c1_c4_predictive_analysis_plan_simplified.yaml` for C1/C4-regioaktivitet
      - `substrate_activity_prediction_plan.yaml` for substrataktivitet
    Begge planene bruker `protein_id` som CV-gruppe og eksplisitt
    exploratory-only tolkning.
    Status 2026-05-21: Python-modulen har nå et testet, foreløpig scaffold for
    leakage-sikre grouped folds på `protein_id` med 5-fold som default når det
    finnes nok proteingrupper, pluss smale baseline-runnere. Dette er ikke en
    full modellfase. De to predictive YAML-planene er nå revidert til kompakte
    aktivitetsplaner med smalere scope og færre låste prediktorer enn før.

21. **`analysis/crystal_anchoring.py`, `cbm_variant.py`, `cbm_comparison.py`**
    som sekundæranalyser.
    CBM paired analysis skal følge
    `cbm_full_length_vs_domain_only_analysis_plan.yaml`, som nå er den
    detaljerte implementasjonsplanen for full-length vs domain-only
    sammenligning, regionmerking, parvise endepunkter og outputtabeller.
    Status 2026-05-21: `analysis/cbm_comparison.py` har nå condition-level
    sideanalyse som bygger `cbm_construct_condition_summary.tsv` og
    `cbm_paired_comparison_table.tsv` fra `condition_table.tsv`,
    `cluster_table.tsv` og proteinmetadata. Pose-level dual-IFP i
    `cbm_variant.py` er fortsatt lavere prioritet.
    Crystal anchoring bruker per i dag lokal gemmi/numpy Kabsch-superposisjon
    på delte pocket C-alpha-atomer, ikke operativ PyMOL `pair_fit`.
        Status 2026-05-16: standalone crystal-reference-slice er verifisert på ekte
        data via `tests/run_tests_scripts/test_crystal_anchoring_real_cifs.sh` for
        `A0A0S2GKZ1` mot begge referansene `5ACI` og `7PXW`.
        - referanser løses fra `input_data/pdb_structure_data.csv` uten å bruke
            `Oligo_Activity` eller `Comment`
        - multikjede-krystaller subsett-es til valgt proteinkjede, tilhørende
            ligand og Cu, med fallback når foretrukket kjede A ikke er holo
        - Status 2026-05-21: crystal-subsettet bevarer nå relevant kjemimetadata
            (`_entity`, `_chem_comp`, `_chem_comp_bond`, `_struct_conn`) og
            remappes eksplisitt til `A` protein, `B/C/D` glykan og `E` site-Cu.
            Ligandbundne crystal references normaliseres før protonering.
        - ProLIF-liganden fra crystal-prep er nå glykan-only: Cu blir værende i
            `complex_H.pdb` for geometri/RMSD, men filtreres ut av
            `ligand_for_prolif.mol2` sammen med solvent/ioner. 6YDC har focused
            regresjonstest for at bare valgt site-Cu beholdes og at Cu ikke
            lekker inn i ligand-only eksporten.
        - standalone-harnessen kan auto-velge en best-rangert AF3-pose for denne
            real-data-valideringen; produksjonsstien sammenligner alle beholdte
            cluster-medoider
        - representative pose og hver crystal reference skriver nå ProLIF-artifakter
            under `crystal_anchoring_output/.../ifp/`
        - `ifp_result.json` er bevisst kort og inneholder bare pose/status,
            residue_names, totals, interaction_counts og active_feature_names; full
            feature-layout ligger i `pose_ifp_table.tsv` og `ifp_matrix.csv`
        - pocket defineres nå som proteinrester innen 5 A fra ligand eller Cu
        - holo-referanser bruker pocket fra crystal-strukturen selv; apo-referanser
            kan få pocket via sekvensprojeksjon fra representant/medoid med
            residunavn-normalisering for varianter som `HIC -> HIS`
        - focused pytest dekker nå både syntetisk ligand+Cu-pocket,
            sekvensremapping, ekte holo-case og ekte apo-case
        - samme slice er nå også koblet inn i `analysis/analysis_orchestrator.py`
            som sekundær per-condition-analyse med focused pytest; produksjonsstien
            sammenligner alle beholdte cluster-medoider, og betingelser uten
            beholdte clusters kan bruke top-level AF3 model CIF som fallback bare
            hvis fallbacken passerer hard QC
        - crystal-IFP må passere samme non-vdW contact-eligibility-regel som
            pose-clustering før `ifp_tanimoto` regnes som sammenlignbar; VdW-only,
            zero-contact og low-specific-contact crystal-IFP-er beholder RMSD og
            geometry-output, men merkes som IFP-non-comparable. Etter crystal-prep
            hardeningen er slike utfall ikke lenger kjent Cu/normaliserings-
            forurensning uten ny evidens.
        - ligandbundne crystal references får nå C1/C4-geometri i
            `crystal_geometry_table.tsv`, og nøkkelfeltene joines inn i
            `crystal_anchor_table.tsv`; `crystal_ifp_diagnostic_summary.tsv`
            oppsummerer IFP-eksklusjoner per unique crystal reference og per
            medoid/fallback-sammenligning
        - analysis-core skriver eksplisitt `crystal_anchoring_stage_completed`
            i `analysis_core_summary.json` og `run_manifest.json`, og denne stage-gaten
            er smoke-validert på ekte data
        Gjenstår:
        - kjøre den integrerte produksjonsstien på et system som faktisk gir
            ikke-tomme medoid-vs-crystal-sammenligninger
        - avklare om dagens gemmi/numpy Kabsch-backend er tilstrekkelig som endelig
            operativ løsning eller om PyMOL-paritet skal implementeres senere
        - bruke `crystal_ifp_diagnostic_summary.tsv` fra en bredere medoid-kjøring
            til å kvantifisere hvor ofte crystal-IFP-ene blir for VdW-dominerte,
            og spikre endelige soft-thresholds for plausibilitet
    ⛔ STOPP: Verifiser cluster_signatures/predictive_cluster_table mot skjema.

22. **`report/`** — build_summary_json, build_metrics_csv, build_report_html.
    Status 2026-05-03: `build_summary_json`, `build_metrics_csv` og `build_report_html` brukes nå fra analysis-core produksjonsstien. Real-data smoke run `analysis_core_real_cifs_613251` skrev `summary.json`, `metrics.csv` og `report.html`; det droppede STA4-caset ble ekskludert fra `metrics.csv`, mens `passed`/`flagged` poser ble beholdt.
    Status 2026-05-04 (senere): metrics- og summary-byggerne får nå cluster- og IFP-data fra den integrerte produksjonsstien. ⛔ STOPP: verifiser schema og HTML visuelt på nytt med de rikere feltene i en ny real-data kjøring.

23. **`cli.py`** — Integrer alle steg. Kjør `lpmo-pipeline run --config ...`.
    Status 2026-05-03: `lpmo-pipeline run` er nå koblet til analysis-core produksjonsstien (`analysis/analysis_orchestrator.py`) for discovery → normalize → hard QC → downstream geometri → rapportering. Real-data smoke run via `tests/run_tests_scripts/test_analysis_core_real_cifs.sh` (jobb 613251) staged 3 CIFs gjennom discovery og produserte 1 `pass`, 1 `soft_flag`, 1 `hard_fail`; 2 poser ble analysert videre.
    Status 2026-05-04 (senere): samme CLI-sti skriver nå også `pose_ifp_table.tsv`, `cluster_assignments.tsv`, `medoid_manifest.tsv` og `condition_cluster_summary.tsv`, og `run_manifest.json` markerer egne gates for IFP- og clustering-steg.
    Status 2026-05-16: analysis-core-stien kjører nå også crystal anchoring per condition, skriver `crystal_anchoring/<condition>/crystal_reference_screen.json`, og markerer `crystal_anchoring_stage_completed` i `run_manifest.json`. Produksjonsstien er smoke-validert uten regresjon, men den siste integrerte real-data-kjøringen som fullførte crystal stage hadde ingen tilgjengelige crystal-referanser for den analyserte betingelsen. En eksplisitt produksjonsvalidering med ikke-tomme medoid-vs-crystal-sammenligninger gjenstår derfor fortsatt.
    Status 2026-05-16 (senere): production analysis-core skriver nå også de
    implementerte tabellflatene som trengs før cluster-annotering:
    `pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv`,
    `qc_attrition_table.tsv`, `crystal_anchor_table.tsv`,
    `crystal_geometry_table.tsv` og `crystal_ifp_diagnostic_summary.tsv`. Focused tester for
    orchestrator/CLI og clustering passer, men real-data inspeksjon av
    clustering-output gjenstår før de senere prediktive lagene
    (`predictive_cluster_table.tsv` og modeller) bygges.
    Status 2026-05-21: production analysis-core skriver nå også
    `condition_table.tsv` med `qc_attrition_table.tsv` som masterradsett, join
    av cluster/convergence/patch/confidence-felter og occupancy-vektede
    cluster-aggregater, samt `protein_summary_table.tsv` som rent groupby-lag
    over `condition_table.tsv`. Focused pytest for condition summary,
    orchestrator og CLI passer.
    Status 2026-05-21 (senere): `tests/run_tests_scripts/run_summary_table_validation.py`
    kan nå regenerere `cluster_table.tsv`, `condition_table.tsv` og
    `protein_summary_table.tsv` fra eksisterende production outputs, med fallback
    til `cluster_signatures.json` når eldre shard-outputs mangler skrevet
    `cluster_table.tsv`. Harnessen er verifisert på merge av 14 komplette
    `tests/tests_results/clustering_pilot_staged/domain_only_shards/*/production_output`
    (1 ufullstendig shard ble eksplisitt rapportert og hoppet over), og
    validerte 126 betingelser med både retained-cluster-, no-valid-cluster- og
    contact-sparse/null-IFP-signaltilfeller.
    Status 2026-05-20: analysis-core skriver nå også `cluster_table.tsv` fra
    Stage 7 `cluster_summaries`, og CLI-en eksponerer artefakten sammen med de
    øvrige cluster-outputene. Focused pytest for `tests/test_cluster_signatures.py`,
    `tests/test_analysis_orchestrator.py` og `tests/test_cli_run.py` passer etter denne wiring-endringen.
    Status 2026-05-19: `lpmo-pipeline run`, real-case pilot-runneren og
    `run_clustering_pilot_full.sh` sender nå `n_jobs` inn i produksjonsstien.
    Prepare, hard QC/Privateer og ProLIF kan dermed bruke flere workers; pilotens
    prepare-oppsummeringer er samtidig komprimert for å unngå mange store,
    nesten like summary-filer. En 4-pose timingprobe viste 56.23 s med `n_jobs=1`
    og 28.94 s med `n_jobs=2`.
    Status 2026-05-19 (senere): normaliseringens CIF-remap er optimalisert i
    `io/gemmi_compat.py` ved å gruppere tag-remapping per mmCIF-loop. En 2-pose
    Slurm-sanitysjekk reduserte normalisering for de samme to domain-only-posene
    fra 37.48 s til 15.30 s totalt. Full pilot skal nå startes med
    `tests/run_tests_scripts/submit_clustering_pilot_staged.sh`, som lager
    protein-level shard-manifester og submitter domain-only/full-length som
    Slurm-arrays over flere jobber/noder. Legacy-wrapperen beholdes som fallback.
    ⛔ STOPP: Full end-to-end hovedanalyse gjenstår etter slik real-data verifikasjon av den utvidede produksjonsstien og senere analyser.

24. **`tuning/` (valgfritt etteranalyse)** — implementer sweep_runner + tune_orchestrator.
    Kjør sammenlignende tuning etter baseline leveranse, dokumenter eventuell parameteroppdatering.
