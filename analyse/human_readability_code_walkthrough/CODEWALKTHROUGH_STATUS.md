# Codewalkthrough status

Runtime update 2026-05-22: Step 6 for ProLIF-facing artifacts is now
`analysis_export`, not protonation. ProLIF uses non-protonated
`complex_for_prolif.pdb` + `ligand_only_for_prolif.pdb` and computes
`ImplicitHBAcceptor`, `ImplicitHBDonor`, and `VdWContact`; main clustering
uses only the implicit H-bond features.

## Formål

Denne walkthroughen skal dekke den nåværende produksjonsnære analysepipen under `analyse/`, fra entrypoint og konfigurasjon via discovery, per-pose forberedelse, hard QC, downstream analyse, clustering, crystal anchoring og rapportering. Målet i denne runden er å kartlegge hva koden faktisk gjør, hvilke kilder som er nyttige, hvilke pipeline-steg som ser ut til å finnes, og hvor dokumentasjon og kode ser ut til å avvike.

## Kilder brukt

| Fil | Type kilde | Relevans | Hvor pålitelig den virker | Viktige linjehenvisninger |
|---|---|---|---|---|
| `src/lpmo_pipeline/cli.py` | kode | Primær CLI-entrypoint for `run`, `discover` og `tune` | Høy | `21`, `30`, `67`, `92`, `144`, `177`, `258` |
| `src/lpmo_pipeline/analysis/analysis_orchestrator.py` | kode | Hovedorkestrator for produksjonsflyten | Høy | `383`, `389`, `400`, `474`, `825`, `956`, `981`, `988`, `995`, `1583`, `1733`, `1804`, `1838`, `1948`, `1976`, `2129`, `2319`, `2336`, `2343`, `2353`, `2367`, `2380`, `2428`, `2433`, `2438`, `2474`, `2481`, `2486`, `2491`, `2513`, `2575`, `2578`, `2581`, `2594`, `2639`, `2646`, `2652`, `2659` |
| `src/lpmo_pipeline/io/discovery.py` | kode | Viser faktisk discovery-kontrakt mot `structure_pipeline` work-root | Høy | `46`, `49`, `56`, `62`, `160`, `161`, `180`, `221`, `273` |
| `src/lpmo_pipeline/io/normalize_mmcif.py` | kode | Viser normalisering, chain remap, atom-mapping gate og glykan-CCD-gate i per-pose prepare | Høy | `109`, `127`, `133`, `157`, `179`, `212`, `244`, `331`, `389`, `437`, `602` |
| `src/lpmo_pipeline/io/cif_to_pdb.py` | kode | Viser PoseBusters-PDB-export, PDBFixer/gemmi-backendvalg og rapportkontrakt | Høy | `81`, `99`, `153`, `162`, `201`, `717` |
| `src/lpmo_pipeline/qc/privateer_runner.py` | kode | Viser at `privateer_input.cif` lages som validert case-lokal kopi fra normalized CIF | Høy | `241`, `252`, `255`, `257` |
| `src/lpmo_pipeline/qc/hard_qc_orchestrator.py` | kode | Viser faktisk hard-QC-rekkefølge og fail-closed backend-policy | Høy | `65`, `74`, `113`, `183`, `203`, `267`, `350`, `382` |
| `src/lpmo_pipeline/qc/qc_report.py` | kode | Viser hvordan proximity, PoseBusters, Privateer og geometri blir samlet til verdict/status | Høy | `51`, `71`, `88`, `109`, `150`, `179`, `224` |
| `src/lpmo_pipeline/qc/active_site_proximity.py` | kode | Viser pre-QC aktiv-sete hard gate og bevaring av avstandsmålinger | Høy | `45`, `62`, `88`, `134`, `139`, `155` |
| `src/lpmo_pipeline/qc/custom_geometry_checks.py` | kode | Viser Cu-His hard/soft gate og Cu-C1/C4-målinger | Høy | `38`, `40`, `100`, `190`, `201`, `221`, `235`, `278` |
| `src/lpmo_pipeline/qc/gates.py` | kode | Viser hvilke terskler som faktisk lastes fra YAML | Høy | `268`, `269`, `270`, `271` |
| `configs/thresholds.yaml` | config | Aktiv terskelkilde for QC/geometri/clustering | Høy | `14`, `15`, `64`, `85`, `117`, `172` |
| `configs/production.analysis_core.example.yaml` | config | Viser tenkt produksjonskonfig og input-rotvalg, men kommentarene ser delvis utdaterte ut | Middels | `3`, `10`, `22`, `27`, `32`, `35`, `36`, `40` |
| `src/lpmo_pipeline/config.py` | kode | Felles lasting av runtime paths og assets | Høy | `20`, `96`, `109`, `161`, `185` |
| `configs/runtime_paths.yaml` | config | Sentral runtime-paths-fil for verktøy og delte asset-referanser | Høy | `7`, `9`, `13`, `16`, `20`, `30`, `31`, `32`, `33`, `34`, `35`, `36` |
| `pyproject.toml` | packaging/config | Registrerer faktisk konsoll-entrypoint for `lpmo-pipeline` | Høy | `54`, `55` |
| `README.md` | README | Nyttig oversikt over nåværende status og påstått produksjonsslice | Middels | `16`, `18`, `67`, `139`, `156`, `214`, `240` |
| `DECISIONS.md` | docs | Implementasjonsforankret beslutningslogikk, spesielt Stage 6 | Middels-Høy | `6`, `101`, `110`, `201` |
| `IMPLEMENTATION_PLAYBOOK.md` | playbook | Nyttig stegoversikt og implementasjonsstatus | Middels | `6`, `8`, `36`, `50`, `142`, `183`, `194`, `317` |
| `MASTERPLAN.md` | docs/plan | Nyttig for planlagt stegmodell og forventede outputflater, men ikke fullt synkron med kode | Middels-Lav | `6`, `53`, `65`, `139`, `156`, `200`, `205` |
| `AF3_LPMO_pipeline_detailed_plan.md` | docs/plan | Primær plan-/designkilde ifølge README, men beskriver også anbefalt og fremoverskuende struktur | Middels | `3`, `4`, `60`, `233` |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | TODO/docs | Nyttig for historikk og gjenstående arbeid, men blander status, TODO og fremtidige lag | Lav-Middels | `43`, `470`, `506` |
| `OPEN_QUESTIONS.md` | docs | Nyttig for eksplisitte åpne beslutninger og usikkerheter | Middels | `66` |
| `tests/test_cli_run.py` | test | Bekrefter hvilke produksjonsoutputs og manifest-gates som forventes nå | Høy | `36`, `37`, `40`, `43`, `80`, `81`, `84`, `87`, `88`, `96`, `124`, `125` |
| `tests/run_tests_scripts/run_analysis_core_real_cifs.py` | test/harness | Viser hvordan nåværende produksjonsslice kjøres på ekte data, og hvilke artefakter som forventes | Høy | `32`, `72`, `110`, `133`, `156`, `158`, `161`, `176`, `177`, `178`, `182`, `184`, `188`, `208`, `209`, `210`, `214`, `216` |
| `tests/run_tests_scripts/run_clustering_pilot_real_case.py` | test/harness | Resumerbar pilot-entrypoint som forbereder og eventuelt kjører produksjonspipen | Høy for pilotflyt | `25`, `66`, `83`, `116`, `135` |
| `tests/run_tests_scripts/submit_clustering_pilot_staged.sh` | shell/harness | Nåværende multi-node wrapper for full clustering-pilot | Høy for pilotorkestrering | `21`, `40`, `138`, `143`, `181` |

## Viktig advarsel om støttedokumenter

README, TODO, planer, playbooks og andre støttedokumenter er brukt som orientering og hypoteser, ikke som fasit. Når slike dokumenter og koden ikke peker samme vei, skal koden under `src/`, relevante config-filer under `configs/`, og tester/harnesser som faktisk kjører produksjonsstien vektes høyest. Se særlig `DECISIONS.md:6`, som eksplisitt sier at koden og lastede config-filer er sannhetskilden, og `README.md:16`, som peker videre til planverket uten å overstyre koden.

## Foreløpig pipeline-oversikt

| Stegnummer | Stegnavn | Kort beskrivelse | Antatte input | Antatte output | Relevante filer | Status | Usikkerhet |
|---|---|---|---|---|---|---|---|
| 1 | Pipeline-start, entrypoint og konfigurasjon | `lpmo-pipeline` parser CLI, laster YAML, initialiserer manifest, resolver `ProductionRunOptions`, og starter analysis-core-orkestratoren | CLI-args, produksjons-YAML, runtime-paths-YAML/env-var | `ProductionRunOptions`, discovery-oppsett, manifest-gates og kall inn i `run_analysis_core` | `pyproject.toml:54,55`, `src/lpmo_pipeline/cli.py:21,67,144,173,177`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:131,376,405,825,909,1583,1732,2318,2636`, `src/lpmo_pipeline/config.py:18,96,101,118`, `configs/production.analysis_core.example.yaml:10,22,27,32,35,36,39,40`, `configs/runtime_paths.yaml:7,13,30,36` | Ferdig | `cli.py`-docstringen og eksempelconfigens toppkommentarer er delvis utdaterte |
| 2 | Discovery av arbeidsrot og poseutvalg | Work-root traverseres etter target/model/run/uniprot/seed-sample-moenstre; ordinare pose-inputs bygges fra sample-CIF-er, mens top-level AF3 model-CIF holdes separat som crystal-fallback per condition | `work_root`, `af3_only`, `latest_only`, `include_targets`, eventuelt `include_proteins` | `pose_inputs`, `discovery_summary`, `discovery_errors`, `crystal_top_model_fallback_candidates` | `src/lpmo_pipeline/io/discovery.py:46,56,62,138,203,238,285`, `src/lpmo_pipeline/utils/data_models.py:398,405,423`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:825,866,909,1619,1651`, `tests/test_discovery.py:264,531,558`, `tests/test_analysis_orchestrator.py:91,115,260` | Ferdig | Reell forekomst av top-level `model.cif` som ordinart pose-input er ikke bekreftet utover at koden stotter det |
| 3 | Per-pose forberedelse | Hver oppdaget CIF normaliseres, konverteres til PoseBusters-PDB, og får generert `privateer_input.cif` og in-memory struktur | Oppdagede CIF-stier | case-mapper, `normalized.cif`, PoseBusters-PDB, `privateer_input.cif`, `PreparedPose` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1031,1055,1063,1069,1074,1134`, `src/lpmo_pipeline/io/normalize_mmcif.py:109`, `src/lpmo_pipeline/io/cif_to_pdb.py:717`, `src/lpmo_pipeline/qc/privateer_runner.py:241` | Ferdig og walkthrough-dokumentert | Normaliseringens interne topologi-/remapdetaljer kan fortsatt dypdykkes senere |
| 4 | Hard QC | Hard-QC-kjeden kjører pre-QC aktiv-sete, geometri, PoseBusters og Privateer før verdict/QC-rapport bygges; aktive backendfeil failer lukket | `PreparedPose` -> `HardQCInput` | `qc_report.json`, QC-verdict per pose, case-metadata, `qc_attrition_table.tsv` | `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:113,183,203,267,350,382`, `src/lpmo_pipeline/qc/qc_report.py:71,88,109,150,179`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1246,1871,1925,2840` | Ferdig og walkthrough-dokumentert | Endelig rapporteringsflate for alle droppede-pose-metrikker kan fortsatt strammes, men main analyse stopper droppede poser riktig |
| 5 | Downstream geometri | Ikke-droppede poser analyseres videre med geometriuttrekk for `pose_geometry.tsv`; metric-feil flagges som `geometry_not_computable` og stopper ikke analysis-export/IFP | QC-pass/flagged poser + struktur | geometri-metrikker per pose, `pose_geometry.tsv`, `geometry_metrics_error` i manifest/case ved feil | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1939,1942,1962,1975,2042,2775` | Ferdig | Full kontrakt for geometri-RMSD/hardening er ikke detaljwalket i denne runden |
| 6 | Analysis export for ProLIF | Ikke-droppede poser eksporteres uten protonering til PDB/PDB-flater for senere ProLIF-batch | normaliserte CIF-er fra analyserbare poser | `analysis_export/complex_for_prolif.pdb`, `analysis_export/ligand_only_for_prolif.pdb`, `analysis_export_report.json` | `src/lpmo_pipeline/io/analysis_export.py`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py` | Ferdig | Bred real-data validering gjenstår etter migrering |
| 7 | Konvergens per betingelse | Per betingelse beregnes konvergens mot referansepose før IFP/clustering brukes videre | ikke-protonerte `complex_for_prolif.pdb` gruppert per condition | `pose_convergence.tsv`, `condition_convergence_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1948,2343,2346` | Pågår | Referanseoppdatering etter clustering er ikke fullstendig verifisert i denne runden |
| 8 | ProLIF og residue contacts | ProLIF kjøres per condition-batch; residue-level contact-tabell bygges fra IFP-resultatene | `complex_for_prolif.pdb`, `ligand_only_for_prolif.pdb`, metadata per pose | `pose_ifp_table.tsv`, condition-wise `ifp_matrix.csv`, `pose_residue_contact_table.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/residue_contact_extraction.py:110,167` | Pågår | Feature-kontrakt og contact eligibility bør walkthroughes eksplisitt senere |
| 9 | Contact eligibility og clustering | IFP-resultater klassifiseres som clusterable eller ikke; clusterable poser gate-es på `minimum_clusterable_n` før agglomerativ clustering per condition med main-feature-filtrert IFP-matrise | IFP-batch + contact eligibility | `cluster_assignments.tsv`, `medoid_manifest.tsv`, `condition_cluster_summary.tsv` med `clustering_status`, condition-lokal `main_clustering_ifp_matrix.csv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/clustering_pilot.py:258,288`, `DECISIONS.md:101,110,201` | Implementert og smoke-validert | Bred full-pilot/regression inspeksjon gjenstår, men real-data smoke dekket både `no_clusterable_poses` og contact-eligible `insufficient_clusterable_signal` 2026-05-21 |
| 10 | Cluster annotation og residue importance | Klynger annoteres med signaturer, flat cluster-tabell, residue-/patch-sammendrag og condition/protein-summary bygges oppå dette | cluster assignments, medoids, geometri, IFP, residue contacts, QC attrition | `cluster_table.tsv`, `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, `cluster_signatures.json`, `protein_condition_residue_scores.tsv`, `protein_residue_regio_delta.tsv`, `condition_patch_summary.tsv`, `protein_patch_summary.tsv`, `condition_table.tsv`, `protein_summary_table.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/cluster_signatures.py`, `src/lpmo_pipeline/analysis/condition_summary.py` | Implementert for nåværende summary-lag | Familie-alignment og bred real-data verifisering av predictive/CBM-sideanalyser gjenstår; predictive-planene er markert for revisjon av modeller/prediksjonsvariabler |
| 11 | Crystal anchoring | Alle beholdte cluster-medoider, eller hard-QC-passert AF3 top-model fallback når condition mangler clusters, screenes mot crystal references; crystal-IFP Tanimoto gates av non-vdW contact eligibility, og ligandbundne crystals får C1/C4-geometri | representative poser per condition | `crystal_reference_screen.json`, `crystal_anchor_table.tsv`, `crystal_geometry_table.tsv`, `crystal_ifp_diagnostic_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2129,2513`, `README.md:156`, `OPEN_QUESTIONS.md:66` | Implementert | Bred integrert real-data-kjøring må fortsatt brukes til å tolke hvor stor andel crystal-IFP-er som blir VdW-/low-specific-contact-ekskludert |
| 12 | Rapportering og eksportflater | Til slutt skrives de samlede pose-/QC-/summary-tabellene, `metrics.csv`, `summary.json`, `report.html` og kjøringsoppsummering | alle foregående mellomresultater | `pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv`, `qc_attrition_table.tsv`, `condition_table.tsv`, `protein_summary_table.tsv`, `metrics.csv`, `summary.json`, `report.html`, `analysis_core_summary.json` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/condition_summary.py`, `tests/run_tests_scripts/run_analysis_core_real_cifs.py` | Ferdig for implementerte eksportflater | Bred real-data verifisering av nye summary-tabeller gjenstår |
| 13 | Tuning og pilotorkestrering | Tuning og clustering-pilot lever ved siden av produksjonsstien og gjenbruker deler av samme kontrollflate | tuning-YAML eller pilot-seleksjonsmanifester | tuning-summary eller pilot-checkpoints/artefakter | `src/lpmo_pipeline/cli.py:30,92`, `src/lpmo_pipeline/tuning/tune_orchestrator.py:28`, `tests/run_tests_scripts/run_clustering_pilot_real_case.py:83,116`, `tests/run_tests_scripts/submit_clustering_pilot_staged.sh:40,181` | Trenger avklaring | Viktig, men trolig sekundært til selve produksjons-walkthroughen |

## Foreløpige entrypoints

- Produksjon: konsollskriptet `lpmo-pipeline` registreres i `pyproject.toml:54,55`, går inn i `src/lpmo_pipeline/cli.py:21`, og starter produksjonsstien via `src/lpmo_pipeline/cli.py:67,144,177`, som går inn i `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1583`.
- Discovery-only: `lpmo-pipeline discover` via `src/lpmo_pipeline/cli.py:258`, som går inn i `src/lpmo_pipeline/io/discovery.py:62`.
- Optional tuning: `lpmo-pipeline tune` via `src/lpmo_pipeline/cli.py:30,92`, som går inn i `src/lpmo_pipeline/tuning/tune_orchestrator.py:28`.
- Real-data harness for nåværende produksjonssti: `tests/run_tests_scripts/run_analysis_core_real_cifs.py:133,161`, vanligvis kjørt via `tests/run_tests_scripts/test_analysis_core_real_cifs.sh:1,57`.
- Resumerbar clustering-pilot-entrypoint: `tests/run_tests_scripts/run_clustering_pilot_real_case.py:66,83,116`.
- Full clustering-pilot wrapper, legacy én-jobb: `tests/run_tests_scripts/run_clustering_pilot_full.sh:20,30,49,119,185,210`.
- Full clustering-pilot wrapper, nåværende staged/multi-node: `tests/run_tests_scripts/submit_clustering_pilot_staged.sh:40,181`.
- Toolchain-smoke, ikke selve analysepipelinen: `tests/run_tests_scripts/run_analysis_tools_smoke.py:16,35`.

## Relevante kodeområder

- `src/lpmo_pipeline/cli.py` (`21`, `144`, `258`) for kommandoflaten.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py` (`383`, `825`, `956`, `1583`, `2319`, `2513`, `2639`) for hovedflyt og outputflater.
- `src/lpmo_pipeline/io/` med særlig discovery og inputforberedelse: `io/discovery.py:62`, `io/normalize_mmcif.py:109`, `io/cif_to_pdb.py:717`, `io/protonate_export.py`.
- `src/lpmo_pipeline/qc/` med særlig `hard_qc_orchestrator.py:96`, `active_site_proximity.py:62`, `custom_geometry_checks.py:97`, `privateer_runner.py:220`, `qc_report.py:224`.
- `src/lpmo_pipeline/analysis/` for downstream analysemoduler: `mdanalysis_metrics`, `convergence_metrics`, `prolif_ifp`, `residue_contact_extraction`, `clustering_agglomerative`, `clustering_hdbscan`, `cluster_signatures`, `residue_importance`, `crystal_anchoring`.
- `src/lpmo_pipeline/report/` via orkestratorens rapportkall på `analysis_orchestrator.py:2578,2581,2594`.
- `configs/` for runtimeparametre, terskler og feature-sett, særlig `configs/runtime_paths.yaml:7,13,30,36`, `configs/thresholds.yaml:14,15,85,117,172` og `configs/production.analysis_core.example.yaml:22,27,35,36`.
- `tests/` og `tests/run_tests_scripts/` som operative bekreftelser på hva produksjonsstien faktisk forventes å skrive nå, særlig `tests/test_cli_run.py:124,125` og `tests/run_tests_scripts/run_analysis_core_real_cifs.py:161,176,188,208,216`.

## Støttedokumenter som bør brukes med forsiktighet

- `configs/production.analysis_core.example.yaml` bør brukes som eksempelkonfig, ikke som full runtime-spesifikasjon. Den viser nå de sentrale production-nodene, men faktisk steglogikk og outputkontrakt må fortsatt leses fra orkestrator, moduler og tester.
- `MASTERPLAN.md` bør brukes med forsiktighet for runtime-nær walkthrough, fordi den blander nåværende og fremtidige output-lag og peker på summary-tabeller som ikke ser ut til å bli skrevet i dagens produksjonskode (`MASTERPLAN.md:53,156,200,205` vs. `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2428,2433,2474,2481,2486,2491,2575`).
- `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` bør brukes med forsiktighet fordi den eksplisitt peker på fremtidige summary-lag som fortsatt må implementeres, samtidig som andre deler fungerer som historisk logg (`DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:470,506`).
- `AF3_LPMO_pipeline_detailed_plan.md` er nyttig som overordnet designkilde, men bør brukes med forsiktighet når den beskriver anbefalt prosjektstruktur eller planlagte lag utover dagens runtime, siden dokumentet selv kaller deler av strukturen «recommended» og «to be adapted» (`AF3_LPMO_pipeline_detailed_plan.md:4,233`).
- `OPEN_QUESTIONS.md` er viktig som usikkerhetsregister, men er per definisjon ikke en beskrivelse av implementert sannhet; bruk den for åpne beslutninger, ikke for å rekonstruere faktisk kjøringsrekkefølge (`OPEN_QUESTIONS.md:66`).
- `plan_implementation_spec.txt` og `plan_analyse.txt` bør behandles som historiske/arkiverte kilder, ikke aktive runtime-spesifikasjoner, fordi `MASTERPLAN.md` selv merker dem som «archived legacy» (`MASTERPLAN.md:10`, `MASTERPLAN.md:11`).

## Mulige avvik mellom dokumentasjon og kode

| Dokumentasjonsfil | Hva dokumentasjonen sier | Hva koden ser ut til å gjøre | Relevante linjer | Vurdering |
|---|---|---|---|---|
| `configs/production.analysis_core.example.yaml` | Viser dagens hovednoder for production-run, inkludert discovery, IFP, clustering, cluster signatures, crystal anchoring og main clustering feature policy | Koden leser de sentrale production-nodene, men detaljert steglogikk ligger fortsatt i orkestratoren og modulene | doc: `configs/production.analysis_core.example.yaml:3`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:387,528,1948,2063,2316,2428,2513` | Oppdatert 2026-05-21 |
| `src/lpmo_pipeline/cli.py` (moduldocstring) | Brukseksempelet viser `python -m lpmo_pipeline run --mode production --config ... --output ...` | Parseren krever også `--del`, og uten dette er eksempelet ikke kjørbart | doc: `src/lpmo_pipeline/cli.py:1,6`; kode: `src/lpmo_pipeline/cli.py:67,74` | Innebygd kodekommentar/dokumentasjon er utdatert og ufullstendig |
| `README.md` | README oppsummerer Cu-His hard/soft gate og fail-closed backendfeil | Kode-defaults, README og `thresholds.yaml` er nå synket til hard gate `1.5–3.0 Å`, mens `1.8–2.6 Å` brukes som soft/preferred QC-bånd via `qc`-seksjonen; PoseBusters/Privateer runner-feil failer lukket når backendene er aktivert | doc: `README.md`; kode: `src/lpmo_pipeline/qc/gates.py`, `src/lpmo_pipeline/qc/custom_geometry_checks.py`, `src/lpmo_pipeline/qc/hard_qc_orchestrator.py` | Oppdatert 2026-05-22 |
| `README.md` | README sier nå at produksjonsstien sammenligner alle beholdte cluster-medoider og bare bruker hard-QC-passert top-level AF3 model CIF for no-cluster conditions | Dette matcher discovery- og crystal-anchoring-koden | doc: `README.md:41`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:909,953,2106,2113`; doc2: `DECISIONS.md:592,606` | Oppdatert 2026-05-20 |
| `MASTERPLAN.md` | Stage 3 sier at `pose_residue_contact_table.tsv` fortsatt er planlagt | Produksjonskoden bygger og skriver `pose_residue_contact_table.tsv`, og playbooken markerer dette som verifisert | doc: `MASTERPLAN.md:139`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2332,2336,2337`; playbook: `IMPLEMENTATION_PLAYBOOK.md:50,183,194` | Del av planverket ser eldre ut enn implementasjonen |
| `MASTERPLAN.md` | Dokumentet lister `cluster_table.tsv`, `condition_table.tsv`, `protein_summary_table.tsv`, predictive og CBM-tabeller som nødvendige analyseflater | Produksjonsstien skriver de sentrale summary-tabellene, og sideanalysemodulene bygger nå første `predictive_cluster_table.tsv`, foreløpig grouped predictive scaffold, `cbm_construct_condition_summary.tsv` og `cbm_paired_comparison_table.tsv` | doc: `MASTERPLAN.md:53,72,83,200,205`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/cluster_signatures.py`, `src/lpmo_pipeline/analysis/condition_summary.py`, `src/lpmo_pipeline/analysis/activity_mapping.py`, `src/lpmo_pipeline/analysis/predictive_models.py`, `src/lpmo_pipeline/analysis/cbm_comparison.py` | Oppdatert 2026-05-21; predictive-modeller/prediksjonsvariabler er ikke valgt og planene er markert for revisjon med 5-fold grouped CV som default |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | TODO-filen omtaler summary-tabellene som upstream outputs for senere analyser | `condition_table.tsv` og `protein_summary_table.tsv` er implementert i produksjonsstien; TODO-filen peker nå videre mot real-data verifisering, predictive modeling, family alignment og CBM paired analysis | doc: `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:470,506`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/condition_summary.py` | Oppdatert 2026-05-21 |

## Analyserte filer i denne runden

- `pyproject.toml`
- `src/lpmo_pipeline/cli.py`
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py`
- `src/lpmo_pipeline/analysis/condition_summary.py`
- `src/lpmo_pipeline/analysis/clustering_pilot.py`
- `src/lpmo_pipeline/io/discovery.py`
- `src/lpmo_pipeline/io/normalize_mmcif.py`
- `src/lpmo_pipeline/io/cif_to_pdb.py`
- `src/lpmo_pipeline/qc/privateer_runner.py`
- `src/lpmo_pipeline/utils/data_models.py`
- `src/lpmo_pipeline/config.py`
- `configs/production.analysis_core.example.yaml`
- `configs/runtime_paths.yaml`
- `README.md`
- `IMPLEMENTATION_PLAYBOOK.md`
- `DECISIONS.md`
- `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md`
- `AF3_LPMO_pipeline_detailed_plan.md`
- `tests/test_discovery.py`
- `tests/test_normalize.py`
- `tests/test_io_contracts.py`
- `tests/test_analysis_orchestrator.py`
- `tests/test_clustering_pilot.py`
- `tests/test_cli_run.py`
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py`

## Gjenstående filer og steg

- Jobb 5: analyse av hard QC og kvalitetsgrenser (`active_site_proximity`, geometri, PoseBusters, Privateer, verdict-flater) er ferdig walkthrough-dokumentert 2026-05-22.
- Jobb 6: analyse av downstream transformasjon (`pose_geometry`, analysis_export, convergence, ProLIF, residue contacts).
- Jobb 7: analyse av clustering-flyt og filtre (contact eligibility, agglomerative path, medoids, noise/outlier-behandling).
- Jobb 8: analyse av cluster annotation, residue importance og senere summary-lag som faktisk finnes versus bare er planlagt.
- Jobb 9: analyse av crystal anchoring og fallback-logikk.
- Jobb 10: analyse av output/eksport og rapportflater (`pose_manifest`, `qc_attrition`, `metrics.csv`, `summary.json`, `report.html`).
- Jobb 11: flowchart-konsolidering og avviksopprydding mellom docs og kode.

## Neste anbefalte prompt

Forslag til neste prompt:

> Analyser steg 5-8 i analysepipen: downstream geometri, analysis_export, konvergens og ProLIF/residue contacts. Fokuser på hvordan `analysis_orchestrator.py` går fra ikke-droppede QC-verdicts til `pose_geometry.tsv`, non-protonated analysis-export artefakter, `pose_convergence.tsv`, `condition_convergence_summary.tsv`, `pose_ifp_table.tsv` og `pose_residue_contact_table.tsv`. Bruk koden som fasit, oppdater `CODEWALKTHROUGH.md` og `CODEWALKTHROUGH_STATUS.md`, og stopp før ny clustering-implementasjon.
