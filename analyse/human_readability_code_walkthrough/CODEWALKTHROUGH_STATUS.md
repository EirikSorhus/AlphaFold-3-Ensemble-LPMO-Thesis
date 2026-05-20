# Codewalkthrough status

## Formål

Denne walkthroughen skal dekke den nåværende produksjonsnære analysepipen under `analyse/`, fra entrypoint og konfigurasjon via discovery, per-pose forberedelse, hard QC, downstream analyse, clustering, crystal anchoring og rapportering. Målet i denne runden er å kartlegge hva koden faktisk gjør, hvilke kilder som er nyttige, hvilke pipeline-steg som ser ut til å finnes, og hvor dokumentasjon og kode ser ut til å avvike.

## Kilder brukt

| Fil | Type kilde | Relevans | Hvor pålitelig den virker | Viktige linjehenvisninger |
|---|---|---|---|---|
| `src/lpmo_pipeline/cli.py` | kode | Primær CLI-entrypoint for `run`, `discover` og `tune` | Høy | `21`, `30`, `67`, `92`, `144`, `177`, `258` |
| `src/lpmo_pipeline/analysis/analysis_orchestrator.py` | kode | Hovedorkestrator for produksjonsflyten | Høy | `383`, `389`, `400`, `474`, `825`, `956`, `981`, `988`, `995`, `1583`, `1733`, `1804`, `1838`, `1948`, `1976`, `2129`, `2319`, `2336`, `2343`, `2353`, `2367`, `2380`, `2428`, `2433`, `2438`, `2474`, `2481`, `2486`, `2491`, `2513`, `2575`, `2578`, `2581`, `2594`, `2639`, `2646`, `2652`, `2659` |
| `src/lpmo_pipeline/io/discovery.py` | kode | Viser faktisk discovery-kontrakt mot `structure_pipeline` work-root | Høy | `46`, `49`, `56`, `62`, `160`, `161`, `180`, `221`, `273` |
| `src/lpmo_pipeline/qc/hard_qc_orchestrator.py` | kode | Viser faktisk hard-QC-rekkefølge | Høy | `96`, `166`, `189`, `344` |
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
| 3 | Per-pose forberedelse | Hver oppdaget CIF normaliseres, konverteres til PoseBusters-PDB, og får generert `privateer_input.cif` og in-memory struktur | Oppdagede CIF-stier | case-mapper, `normalized.cif`, PoseBusters-PDB, `privateer_input.cif`, `PreparedPose` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:956,981,988,995` | Ferdig | Normaliseringens interne detaljregler er ikke gått modul for modul ennå |
| 4 | Hard QC | Hard-QC-kjeden kjører pre-QC aktiv-sete, geometri, PoseBusters og Privateer før verdict/QC-rapport bygges | `PreparedPose` -> `HardQCInput` | `qc_report.json`, QC-verdict per pose, case-metadata | `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:96,166,189,344`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1733` | Ferdig | Dokumentasjon og lastet terskelverdi for Cu-His ser ikke helt synkron ut |
| 5 | Downstream geometri | Ikke-droppede poser analyseres videre med geometriuttrekk for `pose_geometry.tsv` | QC-pass/flagged poser + struktur | geometri-metrikker per pose, `pose_geometry.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1804,2575` | Ferdig | Full kontrakt for geometri-RMSD/hardening er ikke detaljwalket i denne runden |
| 6 | Protonering og eksport for analyse | Ikke-droppede poser protoneres og eksporteres til PDB/MOL2-flater for senere ProLIF-batch | normaliserte CIF-er fra analyserbare poser | `complex_H.pdb`, `ligand_for_prolif.mol2`, protoneringsrapport | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1838` | Ferdig | Fallback/eksterne backends er ikke sporet i detalj utover orkestrator-kallet |
| 7 | Konvergens per betingelse | Per betingelse beregnes konvergens mot referansepose før IFP/clustering brukes videre | protonerte komplekser gruppert per condition | `pose_convergence.tsv`, `condition_convergence_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1948,2343,2346` | Pågår | Referanseoppdatering etter clustering er ikke fullstendig verifisert i denne runden |
| 8 | ProLIF og residue contacts | ProLIF kjøres per condition-batch; residue-level contact-tabell bygges fra IFP-resultatene | `complex_H.pdb`, `ligand_for_prolif.mol2`, metadata per pose | `pose_ifp_table.tsv`, condition-wise `ifp_matrix.csv`, `pose_residue_contact_table.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:1976,2319,2332,2336,2337`, `src/lpmo_pipeline/analysis/residue_contact_extraction.py:110,167` | Pågår | Feature-kontrakt og kontakt-eligibility bør walkthroughes eksplisitt senere |
| 9 | Contact eligibility og clustering | IFP-resultater klassifiseres som clusterable eller ikke; clusterable poser sendes til agglomerativ clustering per condition | IFP-batch + contact eligibility | `cluster_assignments.tsv`, `medoid_manifest.tsv`, `condition_cluster_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2353,2367,2380`, `DECISIONS.md:101,110,201` | Pågår | Produksjonsstien bruker rå `batch.matrix`; forskjellen mellom pilotens feature-pruning og produksjon bør forklares senere |
| 10 | Cluster annotation og residue importance | Klynger annoteres med signaturer, og residue-/patch-sammendrag bygges oppå dette | cluster assignments, medoids, geometri, IFP, residue contacts | `cluster_ifp_signature.tsv`, `cluster_residue_signature.tsv`, `cluster_signatures.json`, `protein_condition_residue_scores.tsv`, `protein_residue_regio_delta.tsv`, `condition_patch_summary.tsv`, `protein_patch_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2417,2428,2433,2438,2467,2474,2481,2486,2491` | Pågår | Senere summary-lag som `cluster_table.tsv` og `condition_table.tsv` er ikke del av nåværende produksjonssti |
| 11 | Crystal anchoring | Alle beholdte cluster-medoider, eller hard-QC-passert AF3 top-model fallback når condition mangler clusters, screenes mot crystal references; crystal-IFP Tanimoto gates av non-vdW contact eligibility, og ligandbundne crystals får C1/C4-geometri | representative poser per condition | `crystal_reference_screen.json`, `crystal_anchor_table.tsv`, `crystal_geometry_table.tsv`, `crystal_ifp_diagnostic_summary.tsv` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2129,2513`, `README.md:156`, `OPEN_QUESTIONS.md:66` | Implementert | Bred integrert real-data-kjøring må fortsatt brukes til å tolke hvor stor andel crystal-IFP-er som blir VdW-/low-specific-contact-ekskludert |
| 12 | Rapportering og eksportflater | Til slutt skrives de samlede pose-/QC-tabellene, `metrics.csv`, `summary.json`, `report.html` og kjøringsoppsummering | alle foregående mellomresultater | `pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv`, `qc_attrition_table.tsv`, `metrics.csv`, `summary.json`, `report.html`, `analysis_core_summary.json` | `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2578,2581,2594,2639,2646,2652,2659`, `tests/run_tests_scripts/run_analysis_core_real_cifs.py:176,177,178,182,184,188,208,209,210,214,216` | Ferdig | Det finnes ikke ett samlet høyere lag som `condition_table.tsv` i dagens produksjonskode |
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
- `src/lpmo_pipeline/io/` med særlig discovery og inputforberedelse: `io/discovery.py:62`, `io/normalize_mmcif.py` (ikke detaljwalket ennå), `io/cif_to_pdb.py`, `io/protonate_export.py`.
- `src/lpmo_pipeline/qc/` med særlig `hard_qc_orchestrator.py:96`, `active_site_proximity.py:62`, `custom_geometry_checks.py:97`, `privateer_runner.py:220`, `qc_report.py:224`.
- `src/lpmo_pipeline/analysis/` for downstream analysemoduler: `mdanalysis_metrics`, `convergence_metrics`, `prolif_ifp`, `residue_contact_extraction`, `clustering_agglomerative`, `clustering_hdbscan`, `cluster_signatures`, `residue_importance`, `crystal_anchoring`.
- `src/lpmo_pipeline/report/` via orkestratorens rapportkall på `analysis_orchestrator.py:2578,2581,2594`.
- `configs/` for runtimeparametre, terskler og feature-sett, særlig `configs/runtime_paths.yaml:7,13,30,36`, `configs/thresholds.yaml:14,15,85,117,172` og `configs/production.analysis_core.example.yaml:22,27,35,36`.
- `tests/` og `tests/run_tests_scripts/` som operative bekreftelser på hva produksjonsstien faktisk forventes å skrive nå, særlig `tests/test_cli_run.py:124,125` og `tests/run_tests_scripts/run_analysis_core_real_cifs.py:161,176,188,208,216`.

## Støttedokumenter som bør brukes med forsiktighet

- `configs/production.analysis_core.example.yaml` bør brukes med forsiktighet fordi kommentarene på toppen fortsatt beskriver en eldre produksjonsslice uten ProLIF, clustering og crystal anchoring, mens nåværende kode og tester forventer disse stegene (`configs/production.analysis_core.example.yaml:3,10` vs. `README.md:156` og `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2319,2353,2428,2513`).
- `MASTERPLAN.md` bør brukes med forsiktighet for runtime-nær walkthrough, fordi den blander nåværende og fremtidige output-lag og peker på summary-tabeller som ikke ser ut til å bli skrevet i dagens produksjonskode (`MASTERPLAN.md:53,156,200,205` vs. `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2428,2433,2474,2481,2486,2491,2575`).
- `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` bør brukes med forsiktighet fordi den eksplisitt peker på fremtidige summary-lag som fortsatt må implementeres, samtidig som andre deler fungerer som historisk logg (`DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:470,506`).
- `AF3_LPMO_pipeline_detailed_plan.md` er nyttig som overordnet designkilde, men bør brukes med forsiktighet når den beskriver anbefalt prosjektstruktur eller planlagte lag utover dagens runtime, siden dokumentet selv kaller deler av strukturen «recommended» og «to be adapted» (`AF3_LPMO_pipeline_detailed_plan.md:4,233`).
- `OPEN_QUESTIONS.md` er viktig som usikkerhetsregister, men er per definisjon ikke en beskrivelse av implementert sannhet; bruk den for åpne beslutninger, ikke for å rekonstruere faktisk kjøringsrekkefølge (`OPEN_QUESTIONS.md:66`).
- `plan_implementation_spec.txt` og `plan_analyse.txt` bør behandles som historiske/arkiverte kilder, ikke aktive runtime-spesifikasjoner, fordi `MASTERPLAN.md` selv merker dem som «archived legacy» (`MASTERPLAN.md:10`, `MASTERPLAN.md:11`).

## Mulige avvik mellom dokumentasjon og kode

| Dokumentasjonsfil | Hva dokumentasjonen sier | Hva koden ser ut til å gjøre | Relevante linjer | Vurdering |
|---|---|---|---|---|
| `configs/production.analysis_core.example.yaml` | Tidligere sa toppkommentarene at ProLIF, clustering og crystal anchoring ikke var del av produksjonsslicen | Kommentarene er nå oppdatert til å dekke ProLIF, convergence, clustering, cluster signatures og crystal anchoring med all-medoid/top-model-fallback-logikk | doc: `configs/production.analysis_core.example.yaml:3`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2319,2353,2428,2513` | Oppdatert 2026-05-20 |
| `src/lpmo_pipeline/cli.py` (moduldocstring) | Brukseksempelet viser `python -m lpmo_pipeline run --mode production --config ... --output ...` | Parseren krever også `--del`, og uten dette er eksempelet ikke kjørbart | doc: `src/lpmo_pipeline/cli.py:1,6`; kode: `src/lpmo_pipeline/cli.py:67,74` | Innebygd kodekommentar/dokumentasjon er utdatert og ufullstendig |
| `README.md` | README oppsummerer Cu-His hard gate som `1.9–2.6 Å` | Hard gate lastes fra `hard_gates` i `thresholds.yaml` og ser ut til å være `1.5–3.0 Å`, mens `1.8–2.6 Å` brukes som soft/preferred QC-bånd via `qc`-seksjonen | doc: `README.md:240`; kode: `src/lpmo_pipeline/qc/gates.py:268,269,270,271`; config: `configs/thresholds.yaml:14,15,172` | Trolig dokumentasjon som ikke er oppdatert etter senere terskelendring/splitt mellom hard og soft range |
| `README.md` | README sier nå at produksjonsstien sammenligner alle beholdte cluster-medoider og bare bruker hard-QC-passert top-level AF3 model CIF for no-cluster conditions | Dette matcher discovery- og crystal-anchoring-koden | doc: `README.md:41`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:909,953,2106,2113`; doc2: `DECISIONS.md:592,606` | Oppdatert 2026-05-20 |
| `MASTERPLAN.md` | Stage 3 sier at `pose_residue_contact_table.tsv` fortsatt er planlagt | Produksjonskoden bygger og skriver `pose_residue_contact_table.tsv`, og playbooken markerer dette som verifisert | doc: `MASTERPLAN.md:139`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2332,2336,2337`; playbook: `IMPLEMENTATION_PLAYBOOK.md:50,183,194` | Del av planverket ser eldre ut enn implementasjonen |
| `MASTERPLAN.md` | Dokumentet lister `cluster_table.tsv`, `condition_table.tsv` og `protein_summary_table.tsv` som nødvendige/tabellag i hovedløpet | Nåværende produksjonssti skriver lavere nivå-flater som signatures, residue scores og patch summaries, men jeg har ikke funnet tilsvarende skriving av disse summary-tabellene i orkestratoren | doc: `MASTERPLAN.md:53,72,83,200,205`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2428,2433,2474,2481,2486,2491,2575` | Ser ut som planlagte senere lag, ikke nåværende produksjonsoutput |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | TODO-filen omtaler `cluster_table.tsv`, `condition_table.tsv` og `protein_summary_table.tsv` som upstream outputs som må finnes eller implementeres | Produksjonskoden ser per nå ut til å stoppe på cluster signatures, residue scores og patch summaries | doc: `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:470,506`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:2428,2433,2474,2481,2486,2491` | Støtter tolkningen av at disse summary-lagene fortsatt ligger foran, ikke bak, dagens produksjonssti |

## Analyserte filer i denne runden

- `pyproject.toml`
- `src/lpmo_pipeline/cli.py`
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py`
- `src/lpmo_pipeline/io/discovery.py`
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
- `tests/test_analysis_orchestrator.py`
- `tests/test_cli_run.py`
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py`

## Gjenstående filer og steg

- Jobb 4: analyse av input/loading og forberedelse (`normalize_mmcif`, `cif_to_pdb`, `prepare_privateer_input`, case-struktur).
- Jobb 5: analyse av hard QC og kvalitetsgrenser (`active_site_proximity`, geometri, PoseBusters, Privateer, verdict-flater).
- Jobb 6: analyse av downstream transformasjon (`pose_geometry`, protonering, convergence, ProLIF, residue contacts).
- Jobb 7: analyse av clustering-flyt og filtre (contact eligibility, agglomerative path, medoids, noise/outlier-behandling).
- Jobb 8: analyse av cluster annotation, residue importance og senere summary-lag som faktisk finnes versus bare er planlagt.
- Jobb 9: analyse av crystal anchoring og fallback-logikk.
- Jobb 10: analyse av output/eksport og rapportflater (`pose_manifest`, `qc_attrition`, `metrics.csv`, `summary.json`, `report.html`).
- Jobb 11: flowchart-konsolidering og avviksopprydding mellom docs og kode.

## Neste anbefalte prompt

Forslag til neste prompt:

> Analyser steg 3 i analysepipen: per-pose forberedelse og case-bygging. Fokuser paa `NormalizeMMCIFRunner` i `src/lpmo_pipeline/io/normalize_mmcif.py`, `convert_cif_to_pdb()` i `src/lpmo_pipeline/io/cif_to_pdb.py`, `prepare_privateer_input()` og `_prepare_pose_case()` i `src/lpmo_pipeline/analysis/analysis_orchestrator.py`, samt de naermeste testene/harnessene som viser forventet case-struktur. Bruk koden som fasit. Oppdater baade `human_readability_code_walkthrough/CODEWALKTHROUGH.md` og `human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md` med en ferdig seksjon for forberedelsessteget, presise fil-/linjehenvisninger, samsvarstabell mot stoettedokumenter og tydelige usikkerheter. Ikke ga videre til hard QC ennå.
