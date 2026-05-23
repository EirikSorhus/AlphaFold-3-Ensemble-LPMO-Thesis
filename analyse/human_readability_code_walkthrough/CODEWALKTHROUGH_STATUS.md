# Codewalkthrough status

Runtime update 2026-05-22: Step 6 for ProLIF-facing artifacts is now
`analysis_export`, not protonation. ProLIF uses non-protonated
`complex_for_prolif.pdb` + `ligand_only_for_prolif.pdb` and computes
`ImplicitHBAcceptor`, `ImplicitHBDonor`, and `VdWContact`; main clustering
uses only the implicit H-bond features.

## Formål

Denne walkthroughen skal dekke den nåværende produksjonsnære analysepipen under `analyse/`, fra entrypoint og konfigurasjon via discovery, per-pose forberedelse, hard QC, downstream analyse, clustering, crystal anchoring og rapportering. Målet i denne runden var konsolidering: gjøre walkthroughen helhetlig, konsistent og klar for sluttkontroll uten å introdusere nye antagelser.

## Kilder brukt

| Fil | Type kilde | Relevans | Hvor pålitelig den virker | Viktige linjehenvisninger |
|---|---|---|---|---|
| `src/lpmo_pipeline/cli.py` | kode | Primær CLI-entrypoint for `run`, `discover` og `tune` | Høy | `21`, `30`, `67`, `92`, `144`, `177`, `258` |
| `src/lpmo_pipeline/analysis/analysis_orchestrator.py` | kode | Hovedorkestrator for produksjonsflyten | Høy | `383`, `389`, `400`, `474`, `825`, `956`, `981`, `988`, `995`, `1583`, `1733`, `1804`, `1838`, `1948`, `1976`, `2129`, `2319`, `2336`, `2343`, `2353`, `2367`, `2380`, `2428`, `2433`, `2438`, `2474`, `2481`, `2486`, `2491`, `2513`, `2575`, `2578`, `2581`, `2594`, `2639`, `2646`, `2652`, `2659` |
| `src/lpmo_pipeline/io/discovery.py` | kode | Viser faktisk discovery-kontrakt mot `structure_pipeline` work-root | Høy | `46`, `49`, `56`, `62`, `160`, `161`, `180`, `221`, `273` |
| `src/lpmo_pipeline/io/analysis_export.py` | kode | Viser non-protonated analysis export og hva som blokkerer videre ProLIF-/convergence-bruk | Høy | `57`, `68`, `77`, `132` |
| `src/lpmo_pipeline/io/normalize_mmcif.py` | kode | Viser normalisering, chain remap, atom-mapping gate og glykan-CCD-gate i per-pose prepare | Høy | `109`, `127`, `133`, `157`, `179`, `212`, `244`, `331`, `389`, `437`, `602` |
| `src/lpmo_pipeline/analysis/mdanalysis_metrics.py` | kode | Viser downstream pose-geometri, terskellasting og `geometry_not_computable`-kontrakten | Høy | `198`, `272`, `301`, `398`, `520`, `727` |
| `src/lpmo_pipeline/analysis/convergence_metrics.py` | kode | Viser fixed reference-valg, proteinalignment og condition-lokal convergence-summary | Høy | `127`, `138`, `144`, `151`, `186`, `227`, `241` |
| `src/lpmo_pipeline/analysis/prolif_ifp.py` | kode | Viser ProLIF-batch, feature-union, active interactions og contact eligibility | Høy | `133`, `162`, `487`, `524`, `552`, `566` |
| `src/lpmo_pipeline/analysis/residue_contact_extraction.py` | kode | Viser at residue-contact-tabellen bygges direkte fra ProLIF-feature-navn og bitvektor | Høy | `110`, `127`, `129`, `167` |
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
| `tests/test_mdanalysis_metrics.py` | test | Bekrefter pose-geometrifelt, terskellogikk og TSV-kontrakt | Høy | `162`, `200`, `226`, `258`, `279`, `301` |
| `tests/test_convergence_metrics.py` | test | Bekrefter seed-1-foerst referansevalg, proteinalignment og convergence-summary | Høy | `70`, `93`, `118`, `153`, `194` |
| `tests/test_prolif_ifp.py` | test | Bekrefter feature-union, ligandresidue-separasjon og contact eligibility | Høy | `31`, `57`, `170`, `218`, `225`, `245`, `266`, `287` |
| `tests/test_residue_contact_extraction.py` | test | Bekrefter at residue-contact-tabellen speiler IFP-feature-bits direkte | Høy | `15`, `81`, `97` |
| `tests/run_tests_scripts/run_analysis_core_real_cifs.py` | test/harness | Viser hvordan nåværende produksjonsslice kjøres på ekte data, og hvilke artefakter som forventes | Høy | `32`, `72`, `110`, `133`, `156`, `158`, `161`, `176`, `177`, `178`, `182`, `184`, `188`, `208`, `209`, `210`, `214`, `216` |
| `tests/run_tests_scripts/run_clustering_pilot_real_case.py` | test/harness | Resumerbar pilot-entrypoint som forbereder og eventuelt kjører produksjonspipen | Høy for pilotflyt | `25`, `66`, `83`, `116`, `135` |
| `tests/run_tests_scripts/submit_clustering_pilot_staged.sh` | shell/harness | Nåværende multi-node wrapper for full clustering-pilot | Høy for pilotorkestrering | `21`, `40`, `138`, `143`, `181` |

## Viktig advarsel om støttedokumenter

README, TODO, planer, playbooks og andre støttedokumenter er brukt som orientering og hypoteser, ikke som fasit. Når slike dokumenter og koden ikke peker samme vei, skal koden under `src/`, relevante config-filer under `configs/`, og tester/harnesser som faktisk kjører produksjonsstien vektes høyest. Se særlig `DECISIONS.md:6`, som eksplisitt sier at koden og lastede config-filer er sannhetskilden, og `README.md:16`, som peker videre til planverket uten å overstyre koden.

## Foreløpig pipeline-oversikt

| Stegnummer | Stegnavn | Kort beskrivelse | Antatte input | Antatte output | Relevante filer | Status | Usikkerhet |
|---|---|---|---|---|---|---|---|
| 1 | Pipeline-start, entrypoint og konfigurasjon | `lpmo-pipeline` parser CLI, laster YAML, initialiserer manifest, resolver `ProductionRunOptions`, og starter analysis-core-orkestratoren | CLI-args, produksjons-YAML, runtime-paths-YAML/env-var | `ProductionRunOptions`, discovery-oppsett, manifest-gates og kall inn i `run_analysis_core` | `pyproject.toml:54,55`, `src/lpmo_pipeline/cli.py:21,67,144,173,177`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:131,376,405,825,909,1583,1732,2318,2636`, `src/lpmo_pipeline/config.py:18,96,101,118`, `configs/production.analysis_core.example.yaml:10,22,27,32,35,36,39,40`, `configs/runtime_paths.yaml:7,13,30,36` | Ferdig | `cli.py`-docstringen og eksempelconfigens toppkommentarer er delvis utdaterte |
| 2 | Discovery av arbeidsrot og poseutvalg | Work-root traverseres etter target/model/run/uniprot/seed-sample-moenstre; ordinare pose-inputs bygges fra sample-CIF-er, mens top-level AF3 model-CIF holdes separat som crystal-fallback per condition | `work_root`, `af3_only`, `latest_only`, `include_targets`, eventuelt `include_proteins` | `pose_inputs`, `discovery_summary`, `discovery_errors`, `crystal_top_model_fallback_candidates` | `src/lpmo_pipeline/io/discovery.py:46,56,62,138,203,238,285`, `src/lpmo_pipeline/utils/data_models.py:405,423`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1008,1018,1092,1819,1880`, `tests/test_discovery.py:264,284,558`, `tests/test_analysis_orchestrator.py:198,222,407` | Ferdig | Reell forekomst av top-level `model.cif` som ordinart pose-input er ikke bekreftet utover at koden stotter det; `include_proteins` er kode- og testforankret, men ikke tydelig dokumentert i README |
| 3 | Per-pose forberedelse | Hver oppdaget CIF normaliseres, konverteres til PoseBusters-PDB, og får generert `privateer_input.cif` og in-memory struktur | Oppdagede CIF-stier | case-mapper, `normalized.cif`, PoseBusters-PDB, `privateer_input.cif`, `PreparedPose` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1031,1055,1063,1069,1074,1134`, `src/lpmo_pipeline/io/normalize_mmcif.py:109`, `src/lpmo_pipeline/io/cif_to_pdb.py:717`, `src/lpmo_pipeline/qc/privateer_runner.py:241` | Ferdig og walkthrough-dokumentert | Normaliseringens interne topologi-/remapdetaljer kan fortsatt dypdykkes senere |
| 4 | Hard QC | Hard-QC-kjeden kjører pre-QC aktiv-sete, geometri, PoseBusters og Privateer før verdict/QC-rapport bygges; aktive backendfeil failer lukket | `PreparedPose` -> `HardQCInput` | `qc_report.json`, QC-verdict per pose, case-metadata, `qc_attrition_table.tsv` | `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:113,183,203,267,350,382`, `src/lpmo_pipeline/qc/qc_report.py:71,88,109,150,179`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1246,1871,1925,2840` | Ferdig og walkthrough-dokumentert | Endelig rapporteringsflate for alle droppede-pose-metrikker kan fortsatt strammes, men main analyse stopper droppede poser riktig |
| 5 | Downstream geometri | Ikke-droppede poser analyseres videre med pose-lokal geometri fra in-memory struktur; metric-feil flagges som `geometry_not_computable` og stopper ikke analysis export eller IFP | QC-pass/flagged poser + `PreparedPose.structure` + geometri-terskler | geometri-metrikker per pose, `pose_geometry.tsv`, `geometry_metrics_error` i manifest/case ved feil | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2036,2079,2881`, `src/lpmo_pipeline/analysis/mdanalysis_metrics.py:198,272,301,398,727`, `tests/test_mdanalysis_metrics.py:162,226`, `tests/test_analysis_orchestrator.py:1397` | Ferdig og walkthrough-dokumentert | RMSD-feltene finnes i kontrakten, men fylles senere av convergence/crystal-steg |
| 6 | Analysis export for ProLIF | Ikke-droppede poser eksporteres uten protonering til PDB/PDB-flater; bare analysis-export-vellykkede poser blir input til konvergens og ProLIF | normaliserte CIF-er fra analyserbare poser | `analysis_export/complex_for_prolif.pdb`, `analysis_export/ligand_only_for_prolif.pdb`, `analysis_export_report.json`, analysis-export-blockers per case | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2079,2100`, `src/lpmo_pipeline/io/analysis_export.py:57,68,77,132`, `tests/test_io_contracts.py:135`, `README.md:20,60` | Ferdig og walkthrough-dokumentert | Egen focused testfil for `analysis_export.py` finnes ikke i denne workspacen; kontrakten bekreftes via modul, io-tester, orkestrator og docs |
| 7 | Konvergens per betingelse | Per condition beregnes konvergens paa analysis-export-vellykkede `complex_for_prolif.pdb`; referansen velges seed-1-foerst og deretter ranking/plDDT/sample | condition-grupperte non-protonerte complex-PDB-er + convergence-config | `pose_convergence.tsv`, `condition_convergence_summary.tsv`, `ligand_rmsd_to_reference`, `convergent_flag` per case | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2194,2628,2631`, `src/lpmo_pipeline/analysis/convergence_metrics.py:127,138,151,186,227,241`, `tests/test_convergence_metrics.py:70,118,153`, `tests/test_analysis_orchestrator.py:882,894` | Ferdig og walkthrough-dokumentert | `MASTERPLAN.md` beskriver referansevalget grovere enn dagens kode |
| 8 | ProLIF og residue contacts | ProLIF kjores per condition-batch paa analysis-export-artefaktene; feature-union alignes per condition, contact eligibility beregnes, og residue-contact-tabellen avledes direkte fra IFP-feature-navnene | `complex_for_prolif.pdb`, `ligand_only_for_prolif.pdb`, contact-eligibility-regel og posemetadata | `pose_ifp_table.tsv`, condition-wise `ifp_matrix.csv`, `pose_residue_contact_table.tsv`, contact-eligibility-felter paa case/metrics | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2222,2239,2604,2616`, `src/lpmo_pipeline/analysis/prolif_ifp.py:133,162,487,524,552,566`, `src/lpmo_pipeline/analysis/residue_contact_extraction.py:110,127,167`, `tests/test_prolif_ifp.py:31,170,225,287`, `tests/test_residue_contact_extraction.py:15,81,97` | Ferdig og walkthrough-dokumentert | Feature-kontrakt er naa dokumentert; neste uferdige tema er clustering-bruken av contact eligibility |
| 9 | Contact eligibility og clustering | IFP-resultater klassifiseres som clusterable eller ikke; clusterable poser gate-es paa `minimum_clusterable_n` og condition-lokal main-feature-matrise foer formell HDBSCAN Jaccard-clustering | IFP-batch + contact eligibility + clustering-policy | `cluster_assignments.tsv`, `medoid_manifest.tsv`, `condition_cluster_summary.tsv` med `clustering_status`, condition-lokal `main_clustering_ifp_matrix.csv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2176,2298,2309,2321,2330,2359,2565`, `src/lpmo_pipeline/analysis/clustering_hdbscan.py:197,299,387,405`, `src/lpmo_pipeline/analysis/clustering_pilot.py:14,258,288`, `tests/test_clustering.py:18,54,119`, `tests/test_analysis_orchestrator.py:919,927` | Ferdig og walkthrough-dokumentert | `DECISIONS.md` inneholder baade eldre og nyere metodeavsnitt; enkelte eldre linjer peker paa agglomerative som produksjonssti selv om aktiv kode bruker HDBSCAN |
| 10 | Cluster annotation og residue importance | Stage 16 annoterer retained non-noise clusters (type, geometri, IFP/residue-signaturer), Stage 16b bygger residue-/patch-sammendrag, og condition/protein-summary bygges som aggregasjonslag over tabellflater | cluster assignments, medoids, geometri, IFP, residue contacts, confidence, convergence, QC attrition | `cluster_table.tsv`, `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, `cluster_signatures.json`, `protein_condition_residue_scores.tsv`, `protein_residue_regio_delta.tsv`, `condition_patch_summary.tsv`, `protein_patch_summary.tsv`, `condition_table.tsv`, `protein_summary_table.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2729,2772,2973,2983`, `src/lpmo_pipeline/analysis/cluster_signatures.py:326,388,414,544,739`, `src/lpmo_pipeline/analysis/residue_importance.py:280,367,460,628`, `src/lpmo_pipeline/analysis/condition_summary.py:261,365`, `tests/test_cluster_signatures.py:37,201`, `tests/test_residue_importance.py:89`, `tests/test_condition_summary.py:18,145` | Ferdig og walkthrough-dokumentert | Familie-alignment og bred real-data verifisering av predictive/CBM-sideanalyser gjenstår; predictive-planene er markert for revisjon av modeller/prediksjonsvariabler |
| 11 | Crystal anchoring | Alle beholdte cluster-medoider, eller hard-QC-passert AF3 top-model fallback når condition mangler clusters, screenes mot crystal references; crystal-IFP Tanimoto gates av non-vdW contact eligibility, og ligandbundne crystals får C1/C4-geometri | representative poser per condition | `crystal_reference_screen.json`, `crystal_anchor_table.tsv`, `crystal_geometry_table.tsv`, `crystal_ifp_diagnostic_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1611,2378,2809`, `src/lpmo_pipeline/analysis/crystal_anchoring.py:1260,1322,1582`, `tests/test_crystal_anchoring.py:367,626`, `tests/test_analysis_orchestrator.py:1069,1170` | Ferdig og walkthrough-dokumentert | Bred integrert real-data-kjøring trengs fortsatt for andelsestimat av eligibility-eksklusjoner/fallback-bruk |
| 12 | Rapportering og eksportflater | Til slutt skrives de samlede pose-/QC-/summary-tabellene, `metrics.csv`, `summary.json`, `report.html` og kjøringsoppsummering | alle foregående mellomresultater | `pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv`, `qc_attrition_table.tsv`, `condition_table.tsv`, `protein_summary_table.tsv`, `metrics.csv`, `summary.json`, `report.html`, `analysis_core_summary.json` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2809,2942`, `src/lpmo_pipeline/report/build_metrics_csv.py:56`, `src/lpmo_pipeline/report/build_summary_json.py:30`, `src/lpmo_pipeline/report/build_report_html.py:104`, `tests/test_build_summary_json.py:6`, `tests/test_report_html.py:11` | Ferdig og walkthrough-dokumentert | Robusthet ved delvis avbrutte run (manglende deltabeller) kan fortsatt probes med ekstra integrasjonstester |
| 13 | Tuning og pilotorkestrering | Tuning og clustering-pilot lever ved siden av produksjonsstien og gjenbruker deler av samme kontrollflate | tuning-YAML eller pilot-seleksjonsmanifester | tuning-summary eller pilot-checkpoints/artefakter | `src/lpmo_pipeline/cli.py:253,289`, `src/lpmo_pipeline/tuning/tune_orchestrator.py:28`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:700,2586`, `src/lpmo_pipeline/analysis/clustering_pilot_real_case.py:659,792`, `tests/test_clustering_pilot_real_case.py:226,291`, `tests/run_tests_scripts/run_clustering_pilot_real_case.py:75,101`, `tests/run_tests_scripts/submit_clustering_pilot_staged.sh:167,217` | Ferdig og walkthrough-dokumentert | `cmd_tune`-kallsignatur mot `run_tuning` ser inkonsistent ut i kode og boer verifiseres med dedikert runtime-test |

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
| `README.md` | README dokumenterer discovery-kontrollene `work_roots`, `af3_only`, `latest_only` og `include_targets` | Koden har i tillegg et aktivt `include_proteins`-filter i discovery-orkestreringen | doc: `README.md:290,300,303`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:448,1018,1111`; test: `tests/test_analysis_orchestrator.py:198` | README er delvis komplett for discovery-steg 2 |
| `MASTERPLAN.md` | Convergence bruker topprangert QC-pass-pose som fast referansepose | Koden prioriterer seed-1-poser foerst og bruker ranking/plDDT bare som senere sortering innen kandidatsettet | doc: `MASTERPLAN.md:163`; kode: `src/lpmo_pipeline/analysis/convergence_metrics.py:127,138,144,166`; test: `tests/test_convergence_metrics.py:70,93` | Motsagt av kode |
| `MASTERPLAN.md` | `pose_residue_contact_table.tsv` er fortsatt planlagt | Produksjonskoden bygger og skriver tabellen direkte fra ProLIF-feature-navnene | doc: `MASTERPLAN.md:149`; kode: `src/lpmo_pipeline/analysis/residue_contact_extraction.py:110,167`; orkestrator: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2616,2621` | Motsagt av kode |
| `MASTERPLAN.md` | Stage 3 sier at `pose_residue_contact_table.tsv` fortsatt er planlagt | Produksjonskoden bygger og skriver `pose_residue_contact_table.tsv`, og playbooken markerer dette som verifisert | doc: `MASTERPLAN.md:139`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2332,2336,2337`; playbook: `IMPLEMENTATION_PLAYBOOK.md:50,183,194` | Del av planverket ser eldre ut enn implementasjonen |
| `MASTERPLAN.md` | Dokumentet lister `cluster_table.tsv`, `condition_table.tsv`, `protein_summary_table.tsv`, predictive og CBM-tabeller som nødvendige analyseflater | Produksjonsstien skriver de sentrale summary-tabellene, og sideanalysemodulene bygger nå første `predictive_cluster_table.tsv`, foreløpig grouped predictive scaffold, `cbm_construct_condition_summary.tsv` og `cbm_paired_comparison_table.tsv` | doc: `MASTERPLAN.md:53,72,83,200,205`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/cluster_signatures.py`, `src/lpmo_pipeline/analysis/condition_summary.py`, `src/lpmo_pipeline/analysis/activity_mapping.py`, `src/lpmo_pipeline/analysis/predictive_models.py`, `src/lpmo_pipeline/analysis/cbm_comparison.py` | Oppdatert 2026-05-21; predictive-modeller/prediksjonsvariabler er ikke valgt og planene er markert for revisjon med 5-fold grouped CV som default |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | TODO-filen omtaler summary-tabellene som upstream outputs for senere analyser | `condition_table.tsv` og `protein_summary_table.tsv` er implementert i produksjonsstien; TODO-filen peker nå videre mot real-data verifisering, predictive modeling, family alignment og CBM paired analysis | doc: `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:470,506`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, `src/lpmo_pipeline/analysis/condition_summary.py` | Oppdatert 2026-05-21 |

## Analyserte filer i denne runden

- `human_readability_code_walkthrough/CODEWALKTHROUGH.md`
- `human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md`
- `pyproject.toml`
- `src/lpmo_pipeline/cli.py`
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py`
- `src/lpmo_pipeline/analysis/mdanalysis_metrics.py`
- `src/lpmo_pipeline/analysis/convergence_metrics.py`
- `src/lpmo_pipeline/analysis/prolif_ifp.py`
- `src/lpmo_pipeline/analysis/residue_contact_extraction.py`
- `src/lpmo_pipeline/analysis/condition_summary.py`
- `src/lpmo_pipeline/analysis/cluster_signatures.py`
- `src/lpmo_pipeline/analysis/clustering_hdbscan.py`
- `src/lpmo_pipeline/analysis/clustering_pilot.py`
- `src/lpmo_pipeline/analysis/clustering_pilot_real_case.py`
- `src/lpmo_pipeline/analysis/residue_importance.py`
- `src/lpmo_pipeline/analysis/crystal_anchoring.py`
- `src/lpmo_pipeline/tuning/tune_orchestrator.py`
- `src/lpmo_pipeline/io/analysis_export.py`
- `src/lpmo_pipeline/io/discovery.py`
- `src/lpmo_pipeline/io/normalize_mmcif.py`
- `src/lpmo_pipeline/io/cif_to_pdb.py`
- `src/lpmo_pipeline/qc/privateer_runner.py`
- `src/lpmo_pipeline/report/build_metrics_csv.py`
- `src/lpmo_pipeline/report/build_summary_json.py`
- `src/lpmo_pipeline/report/build_report_html.py`
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
- `tests/test_mdanalysis_metrics.py`
- `tests/test_convergence_metrics.py`
- `tests/test_prolif_ifp.py`
- `tests/test_crystal_anchoring.py`
- `tests/test_residue_contact_extraction.py`
- `tests/test_clustering_pilot.py`
- `tests/test_clustering_pilot_real_case.py`
- `tests/test_clustering.py`
- `tests/test_cluster_signatures.py`
- `tests/test_build_summary_json.py`
- `tests/test_condition_summary.py`
- `tests/test_cli_run.py`
- `tests/test_report_html.py`
- `tests/test_residue_importance.py`
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py`
- `tests/run_tests_scripts/build_clustering_pilot_selection_manifests.py`
- `tests/run_tests_scripts/split_clustering_pilot_selection_manifest.py`
- `tests/run_tests_scripts/collect_clustering_pilot_shards.py`
- `tests/run_tests_scripts/run_clustering_pilot_full.sh`
- `tests/run_tests_scripts/run_clustering_pilot_real_case.py`
- `tests/run_tests_scripts/submit_clustering_pilot_staged.sh`

## Konsolideringsrunde 2026-05-23

Status: gjennomført.

### Hva som ble ryddet

- `CODEWALKTHROUGH.md` har nå en tydeligere introduksjon for sluttkontroll (pipeline-/fase-/stegnivå, kode som fasit, eksplisitte usikkerheter).
- Det er lagt inn et samlet Mermaid-flowchart som viser production-stien og sidegrener.
- Begrepsbruk er strammet inn med konsekvent bruk av stegnummer 1-13 og fasegruppering.
- Filplasseringstreet er korrigert for clustering-metode (HDBSCAN som aktiv produksjonssti; agglomerative som alternativ/eldre spor).
- Redundans mellom overordnet flyt og faseinnganger er redusert uten å fjerne teknisk viktige detaljer.

### Gjenværende hull og usikkerheter

- Tune-kontraktens kallesignatur er konsistent i kode, men end-to-end runtime-verifisering av datasett-materialisering til `TuningTestCase` er fortsatt ikke gjort i denne runden.
- Flere støttedokumenter har blandede historiske og nåværende avsnitt (særlig metodebeskrivelser i `DECISIONS.md`); walkthroughen markerer dette, men dokumentkildene selv er ikke ryddet her.
- Bred real-data validering av fallback-andeler (top-model fallback i crystal anchoring) er fortsatt et åpent verifiseringspunkt.

### Siste kvalitetssjekk (gjennomført 2026-05-23)

Resultat: gjennomfort med funn.

- Representativ kontroll av paastander i flere steg viser at hovedflyten fortsatt stemmer med kode (entrypoint/discovery, downstream-flyten, clustering, crystal/rapportering).
- Tune-kontraktavviket er lukket i kode: `cmd_tune()` kaller na `run_tuning(...)` med kompatible argumentnavn (`tuning_config_path`, `test_cases`, `pipeline_config`, `max_parallel`).
- Focused testdekning for tune-CLI er lagt til: `tests/test_cli_run.py` dekker na `cmd_tune()`-kallet og manifest-gaten `tuning_completed`.

Konklusjon: walkthroughen er klar for sluttkontroll som dokument. Aapent punkt er na primart funksjonelt: hvordan tuning-datasett skal materialiseres til `TuningTestCase`-liste i faktisk runtime, siden CLI per i dag sender `test_cases=[]`.

Klar for siste kvalitetssjekk: Ja, dokumentet er konsolidert, avviksmerket og konsistent nok for en siste kvalitetssjekk.

## Neste anbefalte prompt

Forslag til neste prompt:

> Forbedre tune-runtimeflyten: (1) legg inn lasting av tuning-datasett til `TuningTestCase` i `cmd_tune`/tuning-orkestrering, (2) legg til focused test for dataset-materialisering (inkludert failure-path ved manglende/ugyldig dataset), og (3) oppdater walkthrough/status med endelig tune-kontrakt med fil-/linjereferanser.

## Oppdatering 2026-05-23 (Prompt 5: hovedprosessering)

Status: gjennomfort som avgrenset leveranse.

### Hva som ble dekket na

- Walkthroughen er oppdatert med ferdige seksjoner i paakrevd mal for:
	- Steg 5: downstream geometri for beholdte poser
	- Steg 6: analysis export som gate til condition-vis prosessering
- Seksjonene er holdt til hovedprosessering (transformasjon/gating/mellomresultater), uten aa fullfore output-/eksportanalyse utover det som er nodvendig for aa forklare datflyten.

### Gjenstaaende hovedprosessering (for neste runde)

- Steg 7: konvergens per condition
- Steg 8: ProLIF og residue-contact-avledning
- Steg 9: contact eligibility + formell clustering
- Steg 10: cluster-annotasjon, residue-importance og condition/protein-aggregering

### Flowchart

- Ingen rekkefolgeendring identifisert i denne avgrensningen.
- Eksisterende flowchart beholdes uendret i denne runden.

### Filplasseringstre

- Ingen nye filer kreves for denne avgrensede dokumentrunden.
- Filplasseringstreet beholdes uendret.

### Neste anbefalte prompt

Neste i rekken basert paa brukerens promptplan: **Prompt 6 - Gjentas til all stor prosessering er dekket**.

Konkrett forslag til Prompt 6:

> Fortsett hovedprosesseringen med Steg 7-10 i samme mal: forklar input, hovedtransformasjoner, filtrering/gating, beriking/mapping, mellomresultater og hvorfor valgene gjoeres. Hold fokus paa prosessering (ikke slutt-output), og oppdater `CODEWALKTHROUGH_STATUS.md` med hva som fortsatt eventuelt gjenstaar.

## Oppdatering 2026-05-23 (Prompt 6, del 1)

Status: gjennomfort som smal avgrensning.

### Hva som ble ferdig

- Steg 7 (konvergens per condition) er dokumentert i `CODEWALKTHROUGH.md` med full mal:
	- input, hovedtransformasjoner, begrunnelser, mellomresultater, output, linjereferanser, dokument-samsvar og usikkerheter.

### Hva som gjenstar

- Steg 8: ProLIF og residue-contact-avledning
- Steg 9: contact eligibility + formell clustering
- Steg 10: cluster-annotasjon, residue-importance og condition/protein-aggregering

### Flowchart

- Ingen flytendring identifisert i denne delrunden.
- Eksisterende flowchart beholdes uendret.

### Filplasseringstre

- Ingen nye relevante filer identifisert for Steg 7.
- Filplasseringstreet beholdes uendret.

### Anbefalt neste prompt

**Prompt 6 (neste del):** Fortsett hovedprosessering med Steg 8 i samme mal (evt. Steg 8+9 hvis du onsker en liten tett gruppe), og oppdater status med nytt gjenstaende scope.

## Oppdatering 2026-05-23 (Prompt 6, del 2)

Status: gjennomfort som avgrenset hovedprosessering.

### Hva som ble ferdig

- Steg 8 (ProLIF og residue contacts) er dokumentert i `CODEWALKTHROUGH.md` med full mal.
- Steg 9 (contact eligibility og clustering) er dokumentert i `CODEWALKTHROUGH.md` med full mal.
- Steg 10 (cluster annotation, residue importance og condition/protein-aggregering) er dokumentert i `CODEWALKTHROUGH.md` med full mal.
- Ingen flytendring identifisert; eksisterende flowchart beholdes.
- Ingen nye relevante filer identifisert; filplasseringstreet beholdes.

### Hva som gjenstar

- Hovedprosessering steg 1-10 er naa dekket i walkthroughen.
- Neste naturlige scope i production-stien er steg 11-12 (crystal anchoring og sluttrapportering), eventuelt en konsistenskontroll av hele dokumentet mot aktiv kode.

### Anbefalt neste prompt med nummer

**Prompt 7:** Fortsett walkthroughen med Steg 11-12 i samme mal (crystal anchoring og rapportering/eksportflater), bekreft gating/fallback-logikk og oppdater status med eventuelle rest-usikkerheter.

## Oppdatering 2026-05-23 (Prompt 7: validering, kvalitetssjekker, feilhaandtering, logging)

Status: gjennomfort som tverrgaaende analyse.

### Hva som ble ferdig

- `CODEWALKTHROUGH.md` er oppdatert med ny seksjon:
	- `### Steg X: Validering, kvalitetssjekker og feilhaandtering`
- Seksjonen dekker naa eksplisitt:
	- tidlig input-/configvalidering i CLI og `load_production_options(...)`
	- hard/soft fail policy gjennom prepare, hard QC, downstream geometri, analysis export, clustering og crystal anchoring
	- hvordan feil materialiseres i tabeller/sammendrag (ikke bare i logger)
	- strukturert logging (`StructuredLogger`, `FailureLog`) og manifest-gates
- Overordnet flytdiagram er justert med en tverrgaaende node for validering/feilhaandtering.
- Filplasseringstreet er utvidet med manglende tverrgaaende moduler:
	- `src/lpmo_pipeline/utils/logging.py`
	- `src/lpmo_pipeline/utils/manifest.py`

### Hva som gjenstar

- Walkthroughen dekker naa hovedflyt (steg 1-13) og tverrgaaende validering/feilpolicy.
- Naturlig neste restarbeid er ikke ny stegdekning, men kvalitetssikring:
	- konsistenskontroll av linjehenvisninger i hele dokumentet
	- eventuell oppstramming av avviksseksjoner mot README/MASTERPLAN/DECISIONS etter siste kodeendringer.

### Flowchart

- Oppdatert i `CODEWALKTHROUGH.md` med egen tverrgaaende node:
	- `Steg X: Validering, QC, feilhaandtering og logging`
- Endringen tydeliggjor at validering/feilpolicy ikke er ett enkelt sekvenssteg, men styrer hele produksjonskjeden.

### Filplasseringstre

- Oppdatert i `CODEWALKTHROUGH.md` for aa synliggjore tverrgaaende observability/reproduserbarhet:
	- `src/lpmo_pipeline/utils/logging.py`
	- `src/lpmo_pipeline/utils/manifest.py`

### Anbefalt neste prompt med nummer

**Prompt 8:** Kjor en streng konsistenskontroll av hele walkthroughen mot aktiv kode: finn og rett eventuelle foreldede linjehenvisninger, fjern dupliserte eller motstridende paastander, og oppdater status med en kort "sluttkontroll bestatt/ikke bestatt"-oppsummering med konkrete avvik.
