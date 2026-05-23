# Codewalkthrough

Runtime update 2026-05-22: ProLIF-facing downstream preparation now writes
non-protonated `analysis_export/complex_for_prolif.pdb` and
`analysis_export/ligand_only_for_prolif.pdb`. ProLIF computes implicit
H-bonds plus VdW; main clustering uses implicit H-bonds only.

## 1. Introduksjon

Denne walkthroughen dekker den produksjonsnaere analysepipen under `analyse/`, fra entrypoint og konfigurasjon via discovery og poseforberedelse videre til QC, downstream analyse, clustering, crystal anchoring og rapportering. Denne forste versjonen er bevisst en struktur- og oversiktsversjon; senere revisjoner skal fylle ut pipeline-stegene mer presist og med tettere kodeforankring (`human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:3-5`).

Gjennomgangen er kodebasert. Runtime-kode under `src/` og config-filer som faktisk lastes av koden prioriteres, mens README, TODO, dokumentasjon og planer brukes som stotte bare der de stemmer med implementasjonen (`DECISIONS.md:3-8`; `human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:31-33`).

## 2. Hvordan lese dette dokumentet

Denne walkthroughen beskriver pipen paa hoeyt nivaa, ikke funksjon for funksjon. Hensikten er aa rekonstruere faktisk flyt, ikke aa lage en API-katalog.

For hvert pipeline-steg skal teksten etter hvert beskrive:

- input
- hva som skjer
- hvorfor det gjoeres, hvis koden eller lastet config faktisk forklarer det
- output

Alle viktige paastander skal forankres med fil- og linjehenvisninger. Hvis dokumentasjon, config-kommentarer og kode peker i ulike retninger, skal avviket sies eksplisitt i stedet for aa glattes over (`DECISIONS.md:6-8`; `human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:31-33`).

Usikkerheter markeres eksplisitt. Der walkthroughen forelopig bygger paa en arbeidsantakelse, skal det staa at det er en arbeidsantakelse og hva som maa verifiseres i kode senere.

## 3. Kildegrunnlag

Kildegrunnlaget bygger forelopig paa de kildekategoriene som allerede er kartlagt i statusdokumentet, med kode som primaergrunnlag og stottedokumenter som sekundaergrunnlag (`human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:7-29`; `DECISIONS.md:6-8`).

- kode: runtime-kode under `src/`, spesielt CLI, orkestrering, I/O, QC, analyse og rapportering
- config: YAML-filer under `configs/` som faktisk lastes av koden
- README: overordnet oversikt over pipens naavaerende flate og status (`README.md:6-18`)
- TODO: dokumenter som blander status, restarbeid og manuelle kontroller
- documentation: beslutningsdokumenter, playbooks og eksplisitte open-questions-dokumenter
- planer: plan- og designfiler som forklarer intendert struktur, men som kan vaere fremoverskuende eller utdaterte

Stottedokumenter kan vaere utdaterte. Hvis de avviker fra runtime-koden, prioriteres kode under `src/`, relevante config-filer under `configs/`, og tester eller harnesser som faktisk bekrefter produksjonsstien (`DECISIONS.md:6-8`; `human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:31-33`).

## 4. Filplasseringstre

Treverket under er et forelopig arbeidsutvalg over filene og mappene som ser mest relevante ut for aa rekonstruere den produksjonsnaere analyseflyten. Det er med vilje filtrert og er ikke et fullstendig repo-tre (`human_readability_code_walkthrough/CODEWALKTHROUGH_STATUS.md:11-29`).

```text
analyse/
├── pyproject.toml                                  # pakkeoppsett og CLI-registrering
├── src/
│   └── lpmo_pipeline/
│       ├── cli.py                                  # CLI entrypoints for run, discover og tune
│       ├── config.py                               # lasting av runtime paths og felles assets
│       ├── analysis/
│       │   ├── analysis_orchestrator.py            # hovedorkestrator for produksjonsflyten
│       │   ├── convergence_metrics.py              # konvergensberegning per condition
│       │   ├── prolif_ifp.py                       # ProLIF/IFP-bygging og featureflate
│       │   ├── residue_contact_extraction.py       # residue-level contact-tabeller
│       │   ├── clustering_agglomerative.py         # naavaerende produksjonsclustering
│       │   ├── clustering_hdbscan.py               # sensitivitetsspor og delt clusteringlogikk
│       │   ├── cluster_signatures.py               # clusterannotasjon og signaturtabeller
│       │   ├── residue_importance.py               # residue- og patch-scorer
│       │   └── crystal_anchoring.py                # screening mot crystal references
│       ├── io/
│       │   ├── discovery.py                        # discovery av poser fra work-root
│       │   ├── mmcif_ingest.py                     # innlasting av mmCIF-strukturer
│       │   ├── normalize_mmcif.py                  # normalisering av mmCIF
│       │   ├── cif_to_pdb.py                       # eksport til PoseBusters-kompatibel PDB
│       │   └── protonate_export.py                 # legacy protonering; ikke aktiv ProLIF-produksjonssti
│       ├── qc/
│       │   ├── hard_qc_orchestrator.py             # hard-QC-rekkefolge og verdict
│       │   ├── active_site_proximity.py            # pre-QC sjekker rundt aktivt sete
│       │   ├── custom_geometry_checks.py           # Cu-geometri og andre strukturregler
│       │   ├── gates.py                            # lasting av terskler fra YAML
│       │   ├── posebusters_runner.py               # PoseBusters-wrapper
│       │   ├── privateer_runner.py                 # Privateer-wrapper
│       │   └── qc_report.py                        # bygging av QC-rapport
│       └── report/
│           ├── build_metrics_csv.py                # skriver metrics.csv
│           ├── build_summary_json.py               # skriver summary.json
│           └── build_report_html.py                # bygger report.html
├── configs/
│   ├── production.analysis_core.example.yaml       # eksempelkonfig for analysekjoring
│   ├── thresholds.yaml                             # terskler for QC, geometri og clustering
│   ├── runtime_paths.yaml                          # runtime path-oppslag
│   ├── prolif_features.yaml                        # featuredefinisjoner for IFP
│   ├── residue_rules.yaml                          # regler for residue-sammendrag
│   └── geometry_rules.yaml                         # geometri- og valideringsregler
├── tests/
│   ├── test_cli_run.py                             # forventede produksjonsoutputs fra CLI
│   ├── test_analysis_orchestrator.py               # tester for hovedorkestrator
│   └── run_tests_scripts/
│       ├── run_analysis_core_real_cifs.py          # real-data harness for analysekjeden
│       ├── run_clustering_pilot_real_case.py       # resumerbar pilot-entrypoint
│       └── submit_clustering_pilot_staged.sh       # staged multi-node wrapper for pilotkjoring
├── human_readability_code_walkthrough/
│   ├── CODEWALKTHROUGH_STATUS.md                   # arbeidsstatus, scope og kildeinventar
│   └── open_questions.md                           # egne sporsmal for walkthrougharbeidet
├── README.md                                       # overordnet oversikt, ikke nodvendigvis fasit
├── DECISIONS.md                                    # implementasjonsforankret beslutningslogikk
├── IMPLEMENTATION_PLAYBOOK.md                      # stegoversikt og implementasjonsstatus
├── MASTERPLAN.md                                   # planverk, delvis fremoverskuende
├── AF3_LPMO_pipeline_detailed_plan.md              # hovedplan omtalt i README
├── DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md         # TODO, historikk og manuelle oppfolginger
└── OPEN_QUESTIONS.md                               # eksplisitte apne usikkerheter
```

Dette treet vil bli strammet inn eller utvidet etter hvert som walkthroughen gaar dypere i de enkelte delene av pipen.

## 5. Pipeline-steg

### Steg 1: Pipeline-start, entrypoint og konfigurasjon

#### Status
Ferdig

#### Formaal
Dette steget etablerer en konkret produksjonskjoring av analysepipen. Det er her et skallkall eller en dokumentert kommando blir gjort om til en faktisk kjoring: CLI-entrypointen parser argumenter, hovedkonfigurasjonen lastes fra YAML, output-roten opprettes, en manifestramme initialiseres, og produksjonskonfigurasjonen blir oversatt til en intern `ProductionRunOptions` som styrer discovery og videre orkestrering. I praksis er dette broen mellom en manuell kommando og resten av analysis-core-flyten.

#### Input
- CLI-argumenter: den kanoniske produksjonsinngangen er `lpmo-pipeline run --mode production --config ... --output ... --del ... [--n-jobs N]`. Konsollskriptet `lpmo-pipeline` registreres i `pyproject.toml`, og `main()` i `src/lpmo_pipeline/cli.py` eksponerer tre subcommands: `run`, `discover` og `tune`. Bare `run` starter den integrerte produksjonsstien. `--mode` er obligatorisk, men har i praksis bare ett lovlig valg (`production`), `--config`, `--output` og `--del` er obligatoriske, og `--n-jobs` defaultes til `1`.
- config-filer: `--config` peker til hoved-YAML-en som leses med `yaml.safe_load()` i `cmd_run()`. Ved produksjonskjoring er det feltene under `production` som faktisk brukes til oppstarten: `construct_type`, `work_roots` eller `work_root`, `run_id`, `af3_only`, `latest_only`, `max_cases`, `include_targets`, `include_proteins`, `run_posebusters`, `run_privateer`, `n_jobs`, `collect_timing_events` og eventuelt `clustering_pilot.*`. `configs/production.analysis_core.example.yaml` matcher disse navna, men toppkommentarene i filen beskriver en eldre og smalere produksjonsslice enn den som kjorer i dagens kode.
- miljoevariabler: jeg fant ingen miljoevariabel i `cmd_run()` eller `run_analysis_core()` som velger produksjonssteg. Den relevante startnaere miljoevariabelen er `LPMO_PIPELINE_RUNTIME_PATHS_CONFIG`, som lar downstream-moduler bytte ut `configs/runtime_paths.yaml` nar de resolver eksterne verktoypaths og sentrale asset-filer. I tillegg ekspanderer `lpmo_pipeline.config` miljoevariabler og `~` inne i path-strenger i runtime-paths-YAML-en.
- hardkodede defaults: produksjonsopsjonene defaultes i kode dersom YAML-en ikke setter dem. Viktige defaults er `construct_type="domain_only"`, `af3_only=True`, `latest_only=True`, `run_posebusters=True`, `run_privateer=True`, `n_jobs=1` og `collect_timing_events=False`. Hvis `latest_only` er aktivt og `latest`-symlinken i en modellmappe er ubrukelig, faller discovery-koden tilbake til hoyeste numeriske run-ID. Runtime-paths loaderen faller ogsa tilbake til den bundlete `configs/runtime_paths.yaml` hvis ingen eksplisitt config-path eller miljoevariabel er satt.
- filer eller mapper: `--output` bestemmer output-roten og opprettes hvis den mangler. Data-roten velges fra `production.work_roots[construct_type]`, med fallback til `production.work_root`. Ekstra runtime-konfig ligger separat i `configs/runtime_paths.yaml`, som peker til felles assets som `configs/thresholds.yaml`, `configs/prolif_features.yaml`, `configs/defaults.yaml`, `configs/geometry_rules.yaml`, `configs/residue_rules.yaml`, `configs/cv_hierarchy.yaml` og `schemas/qc_report_schema.json`.
- dokumenterte kommandoer som stemmer med kode: README-kommanden `lpmo-pipeline run --mode production --config configs/production.analysis_core.example.yaml --output results/del_a --del del_a --n-jobs 4` matcher parseren i `cli.py`. I tillegg finnes en reell smoke-harness under `tests/run_tests_scripts/run_analysis_core_real_cifs.py` som skriver en liten produksjons-YAML, peker den mot en staged `work_root`, og kaller `cmd_run()` direkte.

#### Hva skjer
1. `pyproject.toml` registrerer `lpmo-pipeline` som konsollskript som peker til `lpmo_pipeline.cli:main`.
2. `main()` bygger en `argparse`-basert kommandoflate med `run`, `discover` og `tune`. `discover` er en separat discovery-only flate, og `tune` er en sidegren for post-analyse-tuning. Den integrerte analysepipen starter bare via `run`.
3. `cmd_run()` oppretter output-mappen, initialiserer `run_manifest.json`, laster YAML-konfigen, og registrerer eventuell `tuning_reference` dersom den finnes i configen.
4. Deretter kalles `run_analysis_core(config, output_dir, del_variant, n_jobs)`. Her blir den ra YAML-konfigurasjonen og CLI-verdiene oversatt til en `ProductionRunOptions` med resolverte stier, discoveryfiltre og worker-innstillinger. `--del del_a` eller `--del del_b` blir samtidig normalisert internt til `a` eller `b`.
5. `run_analysis_core()` starter saa discovery over valgt `work_root`, bygger ogsa en separat AF3 top-model fallback-map for senere crystal anchoring, og forbereder case-mapper per pose. Dette er fortsatt en del av pipeline-starten fordi det er her den generelle produksjonskjoringen blir konkretisert til et sett poser som resten av pipen kan arbeide videre med.
6. Videre steg velges ikke gjennom et generelt stage-register eller en eksplisitt `run_ifp/run_clustering/run_reporting`-familie av flagg. I stedet er flyten fast, mens data-tilgjengelighet og tidligere resultater gate'r hva som faktisk kjorer videre:
- bare hvis minst en pose ble forberedt (`prepared_poses`) gaar kjoringen videre til hard QC og resten av produksjonsstien
- bare poser som ikke er `dropped` av QC gaar videre til downstream geometri og non-protonated analysis export
- bare poser med vellykket analysis export samles opp som input til condition-vis konvergens og ProLIF/IFP
- bare hvis det finnes IFP-resultater skrives `pose_ifp_table.tsv` og `pose_residue_contact_table.tsv`
- bare hvis det finnes QC-pass-poser per condition skrives clustering-relaterte tabeller og cluster-annotasjon
- `clustering_pilot.enabled` legger til ekstra pilot-output, men erstatter ikke den primare HDBSCAN-produksjonsclusteringstien
7. Feil og manglende input handteres relativt sent og pragmatisk, ikke gjennom en samlet schema-validering i starten. YAML-lesefeil stoppes i `cmd_run()`. Ugyldig `work_roots`-type eller manglende `work_root` utloser unntak i `load_production_options()`, som saa fanges av `cmd_run()` som en generell produksjonsfeil. Hvis ingen poser blir forberedt, skriver `run_analysis_core()` fortsatt en oppsummering med `qc_report_error`, men returnerer `success=False`, og CLI-en avslutter med feilstatus. Discovery-only kommandoen validerer eksplisitt at `--work-root` er en katalog.

#### Hvorfor det gjoeres
- CLI-en holder et smalt offentlig grensesnitt, mens den storre styringen ligger i YAML. Det ser ut til aa vaere et bevisst skille: output-rot og DEL-branch holdes paa kommandolinjen, mens datakilde, filtre og flere runtimevalg holdes i produksjonsconfigen. README beskriver dette, og koden bekrefter det.
- Runtime-paths ligger i en separat YAML fordi containerpaths, eksterne verktoy og asset-referanser skal kunne endres uten aa hardkode dette i hver modul. Dette er eksplisitt forklart i `src/lpmo_pipeline/config.py` og i kommentarene i `configs/runtime_paths.yaml`, og koden bekrefter at dette er den faktiske mekanismen.
- Discovery er dynamisk i stedet for bygget paa hardkodede target- eller proteinlister. `io/discovery.py` sier dette eksplisitt i moduldocstringen, og filtrene (`include_targets`, `include_proteins`, `latest_only`) brukes faktisk som runtimevalg i oppstartsfasen.
- `latest_only` finnes for aa holde produksjonskjoringer fokusert paa ett valgt run per modelltre, men uten aa gjore kjoringen skjore dersom `latest`-symlinken mangler eller peker feil. Denne begrunnelsen kommer fra discovery-docstringen og bekreftes av fallback-koden til hoyeste numeriske run-ID.
- Jeg fant ingen kode som tilsier at IFP, clustering eller crystal anchoring er valgfrie hovedsteg i dagens `run`-modus. Tvert imot ser `run` ut til aa representere en fast analysis-core-sti der senere steg bare faller bort hvis upstream-data mangler. Det er en kodebasert observasjon, ikke bare en planantakelse.

#### Output
- `cmd_run()` sender videre en ra `config`-dict, en konkret `output_dir`, valgt DEL-branch og worker-antall til orkestratoren.
- `load_production_options()` produserer en `ProductionRunOptions` som samler de startrelevante valgene i en intern struktur.
- `run_analysis_core()` bygger tidlig en `analysis_summary` med blant annet `run_id`, `work_root`, discoveryfiltre, `discovery_summary`, `discovery_errors` og `crystal_top_model_fallback_candidates`.
- Pipeline-starten produserer ogsaa de sentrale mellomobjektene som senere steg bruker videre: `pose_inputs`, `cases`, `prepared_poses`, og en condition-keying som senere styrer konvergens, IFP og clustering.
- Manifestlagene oppdateres allerede fra start: `run_manifest.json` far gates som `production_pipeline_started`, `analysis_core_completed`, `hard_qc_completed`, `geometry_stage_completed`, `ifp_stage_completed`, `clustering_stage_completed`, `cluster_annotation_stage_completed` og `crystal_anchoring_stage_completed`.
- Selv nar kjoringen ikke kommer langt downstream, skrives det fortsatt en `analysis_core_summary.json` og de basisnare tabellene `pose_manifest.tsv`, `pose_confidence.tsv`, `structure_index.tsv` og `qc_attrition_table.tsv`. De rikere outputflatene blir bare satt dersom de tilsvarende stegene faktisk fikk brukbare upstream-data.

#### Relevante filer og linjer
- `pyproject.toml:L54-L55` - registrerer `lpmo-pipeline` som konsollskript mot `lpmo_pipeline.cli:main`
- `src/lpmo_pipeline/cli.py:L21-L88` - definerer CLI-flaten med `run`, `discover` og `tune`, og viser hvilke argumenter `run` faktisk krever
- `src/lpmo_pipeline/cli.py:L144-L205` - produksjonsentrypointen `cmd_run()`: oppretter output-dir, laster YAML, kaller `run_analysis_core()` og skriver manifest-gates
- `src/lpmo_pipeline/cli.py:L208-L254` - feilhandtering for tom produksjonskjoring og generelle exceptions i `cmd_run()`
- `src/lpmo_pipeline/config.py:L18-L18` - miljoevariabelnavnet `LPMO_PIPELINE_RUNTIME_PATHS_CONFIG`
- `src/lpmo_pipeline/config.py:L96-L125` - runtime-paths loaderen og prioriteringsrekkefolgen mellom eksplisitt path, miljoevariabel og bundlet standardfil
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L131-L146` - `ProductionRunOptions` og hvilke startvalg som finnes med defaults
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L376-L475` - `load_production_options()`, inkludert `construct_type`, `work_roots`/`work_root`, filtre, QC-toggle og `clustering_pilot`
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L825-L951` - discovery av ordinare pose-inputs og separat AF3 top-model fallback-input for senere crystal anchoring
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1583-L1671` - starten av `run_analysis_core()`: option-loading, discovery, summary-oppsett og tidlig pilot-metadata
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1732-L2037` - overgangen fra forberedte poser til hard QC, downstream geometri, non-protonated analysis export, konvergens, IFP og contact eligibility
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L2302-L2466` - ekstra pilot-output, betingede IFP-/clustering-flater og cluster-annotasjon nar upstream-data finnes
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L2636-L2679` - fallback nar ingen poser ble forberedt, og skriving av basis-tabeller som alltid henger paa sluttfasen
- `src/lpmo_pipeline/io/discovery.py:L62-L105` - offentlig discovery-kontrakt med `af3_only`, `latest_only` og `include_targets`
- `src/lpmo_pipeline/io/discovery.py:L143-L180` - `latest_only`-logikken og fallback til hoyeste numeriske run
- `configs/production.analysis_core.example.yaml:L1-L40` - eksempelkonfig for produksjonsstart; toppkommentarene er delvis utdaterte, men nodene og feltnavna matcher dagens loader
- `configs/runtime_paths.yaml:L7-L36` - sentrale runtime-paths for verktoy og asset-filer
- `README.md:L125-L199` - dokumentert produksjonskommando og beskrivelse av hva som styres av CLI versus YAML
- `IMPLEMENTATION_PLAYBOOK.md:L386-L401` - playbook-status for at `lpmo-pipeline run` er koblet til analysis-core-stien og senere fikk `n_jobs`
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py:L110-L166` - smoke-harness som bygger minimal produksjonsconfig og kaller `cmd_run()` direkte
- `tests/test_cli_run.py:L11-L125` - focused test som bekrefter at `cmd_run()` skriver manifest-gates for downstream-steg

#### Samsvar med stottedokumenter

| Dokument/kommentar | Hva den sier | Bekreftet i kode? | Avvik/usikkerhet | Linjereferanser |
|---|---|---|---|---|
| `README.md` | Dokumenterer `lpmo-pipeline run` som kanonisk produksjonskommando, og sier at `--output` og `--del` fortsatt styres fra CLI mens datavalg ligger i YAML | Ja | README ser ut til aa vaere oppdatert for dagens kommandoflate | `README.md:L125-L199`; `src/lpmo_pipeline/cli.py:L67-L76`; `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L376-L475` |
| `configs/production.analysis_core.example.yaml` | Viser riktige produksjonsnoder som `construct_type`, `work_roots`, `run_id`, `af3_only`, `latest_only`, `max_cases` og `include_targets` | Delvis | Selve feltnavna matcher loaderen, men toppkommentarene sier fortsatt at ProLIF, clustering og crystal anchoring ikke er med i produksjonsslicen; det stemmer ikke med dagens kode | `configs/production.analysis_core.example.yaml:L1-L40`; `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1934-L2595`; `src/lpmo_pipeline/cli.py:L183-L202` |
| `configs/runtime_paths.yaml` + `src/lpmo_pipeline/config.py` | Beskriver en separat, sentral mekanisme for eksterne verktoypaths og asset-filer | Ja | Dette er ikke hoved-YAML-en for `cmd_run()`, men en parallell runtime-konfig som downstream-moduler bruker | `configs/runtime_paths.yaml:L7-L36`; `src/lpmo_pipeline/config.py:L96-L125` |
| `IMPLEMENTATION_PLAYBOOK.md` | Beskriver at `lpmo-pipeline run` er koblet til analysis-core-stien, senere utvidet med IFP, clustering, crystal anchoring og `n_jobs` | Ja, med tidslinjeforbehold | Eldre statuslinjer beskriver en tidligere smalere slice; de senere linjene stemmer bedre med dagens kode enn de tidlige | `IMPLEMENTATION_PLAYBOOK.md:L386-L401`; `src/lpmo_pipeline/cli.py:L173-L205`; `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1732-L2595` |
| `src/lpmo_pipeline/cli.py` moduldocstring | Viser et `python -m lpmo_pipeline run ...`-eksempel for produksjon | Nei | Eksempelet mangler det obligatoriske `--del`-argumentet og er derfor ikke kjorbart slik det star | `src/lpmo_pipeline/cli.py:L1-L7`; `src/lpmo_pipeline/cli.py:L67-L76` |

#### Usikkerheter
- Jeg fant ingen egen schema-validering for hoved-YAML-en ved pipeline-start. Det betyr at beskrivelsen over bygger paa feltene som faktisk leses i `cmd_run()` og `load_production_options()`, ikke paa en separat kontraktfil.
- `--mode` er obligatorisk i CLI-en, men den styrer ikke noen reell modusgren i dagens produksjonskode fordi eneste lovlige verdi er `production`. Begrunnelsen for aa beholde flagget fremgar ikke tydelig av kommenterer eller dokumentasjon.
- Jeg har bekreftet at `configs/runtime_paths.yaml` er den delte runtime-path-mekanismen, men jeg har ikke i dette steget kartlagt alle downstream-moduler som laster den. Denne seksjonen behandler derfor runtime-paths som en startnaer konfigurasjonsflate, ikke som en fullstendig dependency-gjennomgang.
- README dokumenterer ogsaa sbatch-baserte smoke-kommandoer rundt produksjonsstien. I denne runden har jeg verifisert den Python-baserte harnessen som kaller `cmd_run()` direkte, men jeg har ikke walket alle wrapper-shellskriptene ende til ende.

### Steg 2: Discovery av arbeidsrot og poseutvalg

#### Status
Ferdig

#### Formaal
Discovery-steget gjoer produksjonskonfigens `work_root` om til et konkret sett pose-inputs som resten av analysepipen kan arbeide paa. I praksis betyr det to ting: for det forste traverseres work-rooten etter target-, modell-, run-, protein- og seed/sample-strukturen som `structure_pipeline` skriver; for det andre separeres top-level AF3 model-CIF-er ut i en egen fallback-struktur som bare brukes senere i crystal anchoring hvis en condition ikke ender opp med en brukbar cluster-medoid. Koden viser at discovery derfor ikke bare er en filskann, men et eksplisitt utvalgstrinn som bestemmer hva som teller som ordinare genererte poser og hva som holdes utenfor de vanlige denominatorene.

#### Input
- arbeidsrot/work root: `ProductionRunOptions.work_root` resolves fra produksjonsconfigen og sendes inn til `discover_work_root()` som rotnode for traverseringen. README og eksempelconfigen beskriver dette som AF3-/structure_pipeline-work-rooten valgt via `work_roots.domain_only` eller `work_roots.full_length`.
- CLI-argumenter eller config: discovery bruker ikke egne CLI-flagg i production-pathen, men arver `af3_only`, `latest_only`, `max_cases`, `include_targets`, `include_proteins` og `construct_type` fra `ProductionRunOptions`.
- mapper/filer som skannes: koden forventer et tre paa formen `TARGET/model/runs/run_id/UNIPROT_TARGET/`, med eventuell top-level `*_model.cif`, top-level JSON-filer og `seed-N_sample-M/`-mapper med CIF/JSON. Regexene i `io/discovery.py` er den faktiske kontrakten for target-navn, run-ID-er og uniprot-target-kataloger.
- pose-/modellkandidater: ordinare kandidater er CIF-filer inne i `seed-*_sample-*`-mapper. Hvis en uniprot-target-mappe ikke har noen sample-CIF-er, men har en top-level `model.cif`, kan denne top-level modellen brukes som ordinart pose-input for akkurat den katalogen.
- fallback-input: `_discover_top_model_fallback_inputs()` kjores separat og ser bare etter top-level AF3 `model.cif` per uniprot-target, uten aa blande disse inn i den ordinare poselista.
- testdata eller harness-input som viser forventet bruk: `tests/test_discovery.py` bygger minimale work-root-traer for `discover_work_root()`, `tests/test_analysis_orchestrator.py` bekrefter hvordan discovery-resultatet blir oversatt til `PoseInputRecord`, og `tests/run_tests_scripts/run_analysis_core_real_cifs.py` viser en reell harness som lager en staged work-root og setter `af3_only=True`, `latest_only=True` og `include_targets` i configen.

#### Hva skjer
Discovery begynner med at `discover_work_root()` gaar gjennom work-rooten dynamisk. Bare kataloger som matcher target-regexen tas med videre, og `include_targets` fungerer som en eksplisitt allowlist dersom den er satt. For hver target-mappe skannes bare modellnavn i `ALLOWED_MODELS`, slik at AF3 og eventuelt RF3 kan inngaa, mens Boltz holdes utenfor av modellfilteret. Dette er discovery av arbeidsrot i streng forstand: koden bestemmer hvilken del av katalogtreet som i det hele tatt er gyldig input for resten av pipen.

Deretter velges hvilke runs som skal representere hver target/model-kombinasjon. Nar `latest_only=True`, prover discovery-koden forst aa resolvere `latest`-symlinken og mappe den tilbake til en lokal katalog under `runs/`. Hvis det lykkes, skannes bare det ene runnet. Hvis symlinken mangler eller ikke kan mappes tilbake til en lokal run-katalog, velges i stedet hoyeste numeriske run-ID. Koden skanner altsaa ikke alle historiske runs nar `latest_only=True`; den komprimerer eksplisitt til ett representativt run per modelltre. Nar `latest_only=False`, traverseres alle numeriske run-kataloger.

Innenfor hvert valgt run filtreres kataloger videre paa `{UNIPROT}_{TARGET}`-moensteret. Her bygges den faktiske manifeststrukturen som analysis-orchestratoren bruker videre. Top-level `model.cif` og tilhorende JSON-filer registreres separat fra seed/sample-mappene, og hver `seed-N_sample-M`-mappe klassifiseres etter hvilke filer som finnes der. Sample-statusen blir derfor en del av discovery-resultatet, ikke noe som utledes senere.

Nar orkestratoren bygger vanlige pose-inputs i `_discover_pose_inputs()`, brukes manifestet som fasit. For hver uniprot-target-katalog filtreres eventuelt proteinet bort hvis `include_proteins` er satt og UniProt-ID-en ikke er med. Deretter bygges vanlige pose-inputs for alle sample-CIF-er som finnes i `seed-*_sample-*`-mappene. Hver pose faar med resolved CIF-path, valgt confidence-JSON, protein-ID, target/ligand-ID, modellnavn, seed/sample-indekser, sample-status og baade `source_run_id` og `discovered_run_id`. Det viktige skillet er at `source_run_id` infereres fra den resolverte filstien, slik at en symlinket staging-work-root kan beholde upstream-run-ID-en selv om discovery skjedde under en lokal staging-runmappe.

Poseutvalget bygges altsaa sample-for-sample nar sample-CIF-er finnes. Bare dersom en uniprot-target-katalog ikke har noen sample-CIF-er i det hele tatt, men faktisk har en top-level `model.cif`, legges denne modellen inn som et ordinart pose-input med `run_status="model_cif"`. Det betyr at top-level model-CIF ikke er den normale produksjonsposen nar seed/sample-poser finnes; den er bare en ordinær fallback inne i pose discovery for kataloger uten sample-CIF-er.

Etter at den ordinare poselista er bygd, sorteres den deterministisk paa target, protein, modell, source-run, seed, sample og sti. `max_cases` trunkerer deretter listen etter sorteringen. Dette er siste poseutvalgstrinnet foer neste pipeline-steg tar over.

Parallelt bygger `_discover_top_model_fallback_inputs()` en egen fallback-map for crystal anchoring. Den kaller `discover_work_root()` paa nytt, men hardkoder `af3_only=True` for denne delen og krever at `ut.model_cif_path` finnes. Resultatet indekseres per `condition_id`, ikke som en flat pose-liste. Disse recordene merkes med `run_status="af3_top_model_fallback"` og holdes eksplisitt utenfor den ordinare poselista. `run_analysis_core()` lagrer baade `discovery_summary`, kombinerte discovery-feil og `crystal_top_model_fallback_candidates` i `analysis_summary` foer den sender de vanlige pose-inputene videre til poseforberedelse.

#### Hvorfor det gjoeres
- bekreftet av kode: discovery-modulen sier eksplisitt at alt skal oppdages dynamisk uten hardkodede target- eller uniprotlister, og implementasjonen bruker regex- og katalogfiltre i stedet for forhåndsdefinerte objekter.
- bekreftet av kode: `latest_only` er laget for aa velge ett representativt run per modelltre uten aa bli skjort hvis `latest`-symlinken peker paa et annet work-tre eller mangler. Dette bekreftes baade av `_resolve_latest_run_dir()` og testene som dekker ekstern symlink og fallback til hoyeste numeriske run-ID.
- bekreftet av kode: top-level AF3 fallback holdes utenfor ordinare generated-pose/QC/IFP/clustering-denominatorer. Dette er ikke bare en implisitt effekt; `_discover_top_model_fallback_inputs()` og senere crystal-anchoring-flyt bruker en separat map, og `DECISIONS.md` dokumenterer samme skille.
- dokumentert i stottedokumenter og bekreftet i kode: README beskriver `work_roots`, `af3_only`, `latest_only` og `include_targets` som discovery-kontroller for production-run configen, og disse feltene brukes faktisk i `load_production_options()` og discovery-funksjonene.
- dokumentert i stottedokumenter og bekreftet i kode: playbooken omtaler `io/discovery.py` som ferdig og verifisert med `af3_only` og `latest_only`; testene bekrefter at akkurat disse filtrene er operative.
- usikkert eller antatt: jeg kan ikke fra koden alene fastslaa hvor ofte top-level `model.cif` faktisk brukes som ordinart pose-input i reelle produksjonskjoringer. Koden stotter det, men harnesser og de naermeste testene fokuserer hovedsakelig paa sample-CIF-er som ordinare poser og top-level model-CIF som crystal-fallback.

#### Output
- identifisert arbeidsrot: `WorkRootManifest.summary` gir en kort oppsummering med `work_root`, antall targets, modellmapper, runs, uniprot-target-mapper, totale CIF-filer og antall errors.
- liste over pose-inputs: `_discover_pose_inputs()` returnerer en sortert liste `PoseInputRecord`-objekter som representerer de ordinare posene neste steg skal forberede.
- fallback-inputs: `_discover_top_model_fallback_inputs()` returnerer en separat `dict[condition_id, PoseInputRecord]` for top-level AF3 fallback-kandidater.
- metadata/paths: hver pose inneholder resolved CIF-path, valgt confidence-JSON, modell, protein, ligand/target, sample-status, seed/sample og run-metadata (`source_run_id`, `discovered_run_id`).
- signaler til senere normalisering, QC eller analyse: `run_analysis_core()` legger discovery-resultatet inn i `analysis_summary` som `discovery_summary`, `discovery_errors` og `crystal_top_model_fallback_candidates`, og sender `pose_inputs` videre til case-preparering. Selve fallback-mappen holdes igjen til crystal-anchoring-stien senere i pipen.

#### Relevante filer og linjer
- `src/lpmo_pipeline/io/discovery.py:L1-L15` - moduldocstringen beskriver den faktiske work-root-strukturen discovery forventer.
- `src/lpmo_pipeline/io/discovery.py:L46-L56` - regex-kontraktene for target-navn og `{UNIPROT}_{TARGET}`-kataloger.
- `src/lpmo_pipeline/io/discovery.py:L62-L108` - `discover_work_root()` filtrerer work-rooten paa target-navn, `af3_only` og `include_targets`.
- `src/lpmo_pipeline/io/discovery.py:L138-L201` - `_discover_model()` implementerer `latest_only` med symlinkpreferanse og fallback til hoyeste numeriske run-ID.
- `src/lpmo_pipeline/io/discovery.py:L203-L208` - `_resolve_latest_run_dir()` forklarer hvorfor eksterne `latest`-symlinker fortsatt kan mappes til lokal `runs/`-katalog.
- `src/lpmo_pipeline/io/discovery.py:L214-L283` - run- og uniprot-target-discovery, inkludert top-level `model.cif`, confidence-JSON og seed/sample-kataloger.
- `src/lpmo_pipeline/io/discovery.py:L285-L316` - `_discover_sample()` klassifiserer sample-status fra filinnholdet.
- `src/lpmo_pipeline/utils/data_models.py:L398-L436` - `WorkRootManifest` eksponerer `iter_all_uniprot_targets()` og `summary`, som analysis-orchestratoren bygger videre paa.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L825-L907` - `_discover_pose_inputs()` bygger den ordinare poselista fra sample-CIF-er, med top-level `model.cif` bare som lokal fallback nar sample-CIF-er mangler.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L909-L953` - `_discover_top_model_fallback_inputs()` bygger separat AF3-only fallback per condition.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1619-L1658` - `run_analysis_core()` kobler discovery-resultatet inn i `analysis_summary` og lagrer fallback-kandidater separat.
- `tests/test_discovery.py:L264-L291` - tester at `latest_only` bruker `latest`-symlink eller faller tilbake til hoyeste numeriske run-ID.
- `tests/test_discovery.py:L520-L565` - bekrefter at `latest_only`, `af3_only` og `include_targets` faktisk avgrenser discovery-scope.
- `tests/test_analysis_orchestrator.py:L91-L112` - viser at `include_proteins` filtrerer bort andre UniProt-kandidater i den ordinare poselista.
- `tests/test_analysis_orchestrator.py:L115-L138` - viser at top-level AF3-modellen holdes separat fra vanlige sample-poser som crystal-fallback.
- `tests/test_analysis_orchestrator.py:L260-L299` - viser at discovery bevarer upstream `source_run_id` gjennom symlinket staging-work-root.
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py:L110-L118` - harnessen skriver produksjonsconfig med `work_root`, `af3_only`, `latest_only` og `include_targets`.
- `tests/run_tests_scripts/run_analysis_core_real_cifs.py:L155-L158` - harnessen lager en staged work-root og lar discovery jobbe mot denne strukturen.
- `README.md:L156-L189` - README beskriver discovery som del av production-run og dokumenterer `work_roots`, `af3_only`, `latest_only` og `include_targets`.
- `IMPLEMENTATION_PLAYBOOK.md:L40-L40` - playbooken oppsummerer `io/discovery.py` som ferdig og verifisert med `af3_only` og `latest_only`.
- `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:L77-L81` - dokumenterer samme `latest_only`-oppstramming som testene dekker.
- `DECISIONS.md:L592-L606` - dokumenterer at top-level AF3 fallback er separat fra generated-pose- og clustering-denominatorer.

#### Samsvar med stottedokumenter

| Stottedokument | Hva dokumentet sier | Hva koden viser | Vurdering | Relevante linjer |
|---|---|---|---|---|
| `README.md` | Production-run bruker discovery over `work_roots`, med `af3_only`, `latest_only` og `include_targets` som discovery-kontroller | `load_production_options()` leser disse feltene, og discovery-funksjonene bruker dem direkte | Bekreftet | doc: `README.md:L156-L189`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L376-L475`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L825-L953` |
| `configs/production.analysis_core.example.yaml` | Eksempelconfigen viser `work_roots`, `af3_only`, `latest_only`, `max_cases` og `include_targets` som discovery-scope | Disse nodene matcher faktisk feltene som loaderen og discovery-koden leser | Bekreftet | doc: `configs/production.analysis_core.example.yaml:L24-L40`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L376-L475` |
| `IMPLEMENTATION_PLAYBOOK.md` | `io/discovery.py` er ferdig og verifisert med `af3_only` + `latest_only` | Testene og implementasjonen bekrefter at begge filtrene er operative | Bekreftet | doc: `IMPLEMENTATION_PLAYBOOK.md:L40-L40`; kode/test: `src/lpmo_pipeline/io/discovery.py:L62-L108`, `src/lpmo_pipeline/io/discovery.py:L138-L201`, `tests/test_discovery.py:L520-L549` |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | `latest_only=True` skal respektere `latest` og ellers falle tilbake til hoyeste numeriske run-ID | Discovery-koden gjor akkurat dette, og testene dekker begge grenene | Bekreftet | doc: `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:L77-L81`; kode: `src/lpmo_pipeline/io/discovery.py:L152-L180`, `src/lpmo_pipeline/io/discovery.py:L203-L208`; test: `tests/test_discovery.py:L264-L291` |
| `DECISIONS.md` | Top-level AF3 fallback skal oppdages separat og ikke inngaa i generated-pose, QC, IFP, clustering eller medoid-denominatorer | `_discover_top_model_fallback_inputs()` bygger en separat map, og den ordinare pose discovery bruker sample-CIF-er som hovedkilde | Bekreftet | doc: `DECISIONS.md:L592-L606`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L825-L953`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1619-L1658` |
| `README.md` | Production-pathen beskriver naa at alle beholdte cluster-medoider screenes, og at no-cluster conditions bare kan bruke top-level AF3 model-CIF etter hard QC | Dette matcher discovery- og crystal-anchoring-koden, der fallbacken holdes i en separat map og ikke blandes inn i ordinare pose-/clustering-denominatorer | Bekreftet etter doc-oppdatering 2026-05-20 | doc: `README.md:L41-L59`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L909-L953`, `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L2106-L2113`; doc2: `DECISIONS.md:L592-L606` |

#### Usikkerheter
- Jeg kan ikke fra discovery-koden alene fastslaa om top-level `model.cif` som ordinart pose-input faktisk forekommer i dagens reelle produksjonsdata, bare at koden stotter det nar sample-CIF-er mangler.
- `include_proteins` filtrerer paa `ut.uniprot_id` etter at manifestet allerede er bygd. Det er klart hva koden gjor, men det fremgar ikke eksplisitt av README hvorfor proteinfilteret ligger i orkestratoren og ikke i `discover_work_root()`.
- Discovery returnerer `manifest.errors`, men de naermeste testene fokuserer mest paa vellykket traversering. Jeg har derfor ikke i denne runden kartlagt hele feilflatens praktiske variasjon utover target-mismatch, permission-feil og ugyldig work-root.

### Steg 3: Per-pose forberedelse og case-bygging

#### Status
Ferdig for dagens production-path; detaljer om normaliseringens interne kjemiregler kan fortsatt walkthroughes dypere senere.

#### Formaal
Dette steget gjoer hver oppdaget `PoseInputRecord` om til en case-mappe og et `PreparedPose`-objekt som hard QC kan bruke. Det er her den raa AF3-CIF-en blir normalisert, eksportert til en PoseBusters-spesifikk PDB, kopiert/validert som Privateer-input og lest tilbake som in-memory gemmi-struktur. Steget endrer ikke hvilke poser som er QC-pass eller analysert videre; det etablerer bare de case-lokale artefaktene og markerer tidlige prep-feil.

#### Input
- `PoseInputRecord`: kommer fra discovery og inneholder resolved CIF-path, pose-ID, protein-ID, ligand/target, modell, run-ID, seed/sample og eventuell confidence-JSON.
- output-rot: `ProductionRunOptions.output_dir` brukes til aa lage `cases/{index:04d}_{slug(pose_id)}`.
- normalisering: `NormalizeMMCIFRunner(input_cif, case_dir / "normalize")` leser raa CIF, bygger chain remap, validerer atom-mapping og glykan-CCD, og skriver `normalized.cif` ved suksess.
- PoseBusters-input: `convert_cif_to_pdb(normalized.cif, case_dir / "posebusters_input")` skriver `for_posebusters.pdb` og `cif_to_pdb_report.json`.
- Privateer-input: `prepare_privateer_input(normalized.cif, case_dir / "privateer_input.cif")` validerer CCD-koder og kopierer kontrollert til case-lokal CIF.

#### Hva skjer
1. `_prepare_pose_case()` lager case-mappen og starter en case-dict med posemetadata, `case_dir` og `status="preparing"`.
2. Normalisering kjorer forst. `NormalizeMMCIFRunner.run()` parser CIF-en med gemmi, bygger chain mapping fra `_entity`/`_struct_asym`, remapper protein/glykan/metall til canonical chain schema, sjekker atom-mapping coverage, validerer glykanresiduer mot CCD-reglene, sjekker `_struct_conn` og `_chem_comp_bond`, skriver `normalize_report.json`, og returnerer `normalized.cif` bare ved suksess.
3. Dersom atom mapping coverage er under 100 %, eller glykanvalideringen feiler, stopper normalisering med `False, None`. Missing glycan chain er en hard normalization failure, ikke en soft downstream-warning.
4. CIF-til-PDB-konverteringen kjorer etter normalisering. Den prover PDBFixer/OpenMM forst, faller tilbake til gemmi hvis PDBFixer ikke kan brukes, og skriver alltid en konverteringsrapport ved suksess. PoseBusters-eksporten er metall-strippet som default fordi Cu/metaller kan vaere inkompatible med PoseBusters' ligand/protein-splitting.
5. Privateer-input bygges fra den normaliserte CIF-en. Dagens `prepare_privateer_input()` er ikke en kjemisk transformasjon; den validerer at CIF-teksten bare bruker tillatte monosakkarid-CCD-koder og kopierer filen til `privateer_input.cif`.
6. Den normaliserte CIF-en leses inn med `gemmi.read_structure()`, og et `PreparedPose` bygges med `pose`, `case_dir`, `normalized_cif`, `posebusters_pdb`, `privateer_input_cif` og `structure`.
7. Ved suksess oppdateres case-dicten til `status="prepared"` og faar paths til `normalized_cif`, `posebusters_pdb` og `privateer_input_cif`. Ved exception settes `status="prep_error"`, med error-string og traceback, og posen blir ikke med i `prepared_poses`.
8. `_prepare_pose_cases()` kjorer enten seriell loop eller `ProcessPoolExecutor`, avhengig av `n_jobs` og antall poser. I parallellgrenen returnerer workerne bare path-payloads; parent-prosessen bygger `PreparedPose` pa nytt og leser gemmi-strukturen lokalt. Outputrekkefolgen beholdes etter original poseindeks.

#### Hvorfor det gjoeres
- Case-mappen samler alle per-pose artefakter paa ett sted, slik at senere hard QC, analysis export, IFP og feilsoking kan referere til samme `case_dir`.
- Normalisering maa skje foer hard QC fordi downstream-stegene forventer canonical chain schema: protein `A`, glykaner `B..D`, metall `E`.
- Missing glycan-chain stoppes tidlig fordi protein+metal-only CIF-er ellers ville sett ut som vellykket normalisering og feilet mye senere som "no ligand atoms".
- `for_posebusters.pdb` er et PoseBusters-spesifikt artefakt, ikke masterstrukturen. Masterformatet for analyse er fortsatt normalisert mmCIF.
- Privateer-input holdes case-lokalt selv om den for oyeblikket er en validert kopi; det gir en stabil kontrakt for hard-QC-orchestratoren og Privateer-batchen.
- Parallell prepare holder de uavhengige per-pose I/O-stegene skalerbare, men rebuild i parent-prosessen gjor at `PreparedPose.structure` ikke maa serialiseres mellom prosesser.

#### Output
- per case-mappe:
  - `normalize/normalized.cif`
  - `normalize/normalize_report.json`
  - eventuelle debug/failure-artefakter som `atom_map.tsv`, `rename_log.json` og `normalize_failures.json`
  - `posebusters_input/for_posebusters.pdb`
  - `posebusters_input/cif_to_pdb_report.json`
  - `privateer_input.cif`
- in-memory:
  - `cases`: en liste med case-dicter for alle pose-inputs, inkludert prep-feil
  - `prepared_poses`: bare de posene som faktisk fikk alle prep-artefaktene
  - `PreparedPose`: objektet som `_hard_qc_input()` senere oversetter til `HardQCInput`
- tabellflater senere i samme orkestrator:
  - `pose_manifest.tsv` bruker case-status, `case_dir` og `normalized_cif`
  - `structure_index.tsv` bruker raw path, normalized CIF, PoseBusters-PDB, Privateer-input og senere analysis-export paths

#### Relevante filer og linjer
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L201-L208` - `PreparedPose`-dataklassen og feltene som hard QC bruker videre.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1031-L1099` - `_prepare_pose_case()` bygger case-mappe, normaliserer, konverterer til PDB, forbereder Privateer-input og handterer prep-feil.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1102-L1131` - worker-payload og rebuild av `PreparedPose` etter parallell prepare.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1134-L1194` - `_prepare_pose_cases()` velger seriell eller prosessbasert prepare og bevarer inputrekkefolgen.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1197-L1203` - `_hard_qc_input()` viser hvilke prepared artefakter som gaar inn i hard QC.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1206-L1270` - `pose_manifest.tsv` bruker case-status og normalized-CIF-path fra prepare.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1324-L1362` - `structure_index.tsv` samler raw, normalized, PoseBusters- og Privateer-paths.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L109-L232` - `NormalizeMMCIFRunner.run()` og normaliseringsrekkefolgen.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L244-L319` - chain mapping fra AF3/entity-layout til canonical chain IDs.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L331-L384` - remap av gemmi-struktur og relevante CIF-loop-tags.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L389-L435` - AF3 identity atom-mapping og 100 % coverage-kontrakt.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L437-L494` - glykan-CCD-validering og hard fail ved manglende glykan.
- `src/lpmo_pipeline/io/normalize_mmcif.py:L602-L678` - debug-artefakter og `normalize_report.json`.
- `src/lpmo_pipeline/io/cif_to_pdb.py:L81-L151` - `CIFToPDBRunner.run()` og rapportskriving.
- `src/lpmo_pipeline/io/cif_to_pdb.py:L153-L239` - PDBFixer-primaerbackend og gemmi-fallback.
- `src/lpmo_pipeline/io/cif_to_pdb.py:L717-L728` - `convert_cif_to_pdb()` wrapperen som orkestratoren kaller.
- `src/lpmo_pipeline/qc/privateer_runner.py:L241-L258` - `prepare_privateer_input()` validerer og kopierer normalized CIF til case-lokal Privateer-input.
- `tests/test_analysis_orchestrator.py:L305-L394` - fake prepare/QC-harness som bekrefter orkestratorens case-byggingskontrakt.
- `tests/test_normalize.py:L278-L480` - normaliseringsregresjoner for success, remap, invalid glykan og missing glykan.
- `tests/test_io_contracts.py:L147-L288` - CIF-to-PDB-kontrakter og PDBFixer/gemmi-fallbackrelaterte forventninger.

#### Samsvar med stottedokumenter

| Stottedokument | Hva dokumentet sier | Hva koden viser | Vurdering | Relevante linjer |
|---|---|---|---|---|
| `IMPLEMENTATION_PLAYBOOK.md` | Normalisering, mapping, non-protonated analysis export og CIF-to-PDB er implementert/verifisert, og rutinemessig atom-map debug-output er redusert | Koden skriver `normalize_report.json` alltid, men skriver atom-map/rename debug bare ved mapping- eller CCD-failure | Bekreftet | doc: `IMPLEMENTATION_PLAYBOOK.md:L42-L46`; kode: `src/lpmo_pipeline/io/normalize_mmcif.py:L602-L678` |
| `MASTERPLAN.md` | Canonical chain schema er protein `A`, glycans `B..D`, metal `E`; missing glycan etter remap skal hard-faile | `_build_chain_mapping()` og `_validate_glycan_residues()` implementerer dette i prepare-steget | Bekreftet | doc: `MASTERPLAN.md:L44-L47`; kode: `src/lpmo_pipeline/io/normalize_mmcif.py:L244-L319`, `src/lpmo_pipeline/io/normalize_mmcif.py:L437-L494` |
| `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md` | Missing-glycan normalization false-success path er resolved, og no-glycan inputs skal stoppe ved normalisering | `_validate_glycan_residues()` returnerer `MISSING_GLYCAN_CHAIN` som hard fail, og `run()` returnerer `False, None` | Bekreftet | doc: `DOCUMENTATION_TODO_AND_MANUAL_CHECKS.md:L43-L81`; kode: `src/lpmo_pipeline/io/normalize_mmcif.py:L179-L207`, `src/lpmo_pipeline/io/normalize_mmcif.py:L450-L463` |
| `IMPLEMENTATION_PLAYBOOK.md` | CIF-to-PDB bruker PDBFixer som primaer backend og gemmi som fallback | `_convert_auto()` prover PDBFixer og faller tilbake til gemmi med `backend_fallback_reason` | Bekreftet | doc: `IMPLEMENTATION_PLAYBOOK.md:L46-L46`; kode: `src/lpmo_pipeline/io/cif_to_pdb.py:L153-L239` |
| `README.md` | Production-run skriver `pose_manifest.tsv`, `structure_index.tsv` og QC-/analysis-artefakter fra samme production-path | Orkestratoren skriver manifest og structure index fra case-dictene etter prepare | Bekreftet | doc: `README.md:L156-L214`; kode: `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1206-L1362` |

#### Usikkerheter
- Normaliseringens topologi- og CIF-remap-detaljer er bare walket paa kontrollflytnivaa her. En egen modulgjennomgang av `gemmi_compat.py`, atom-map-debugging og CCD-oppslag kan fortsatt vaere nyttig.
- `prepare_privateer_input()` er i dag en validert kopi, men navnet gir rom for mer transformasjon senere. Hvis Privateer-input senere divergerer fra `normalized.cif`, maa denne seksjonen oppdateres.
- `convert_cif_to_pdb()` stripper metaller for PoseBusters som default. Det er riktig for dette artefaktet, men maa ikke forveksles med analyse-/geometri-strukturen, som fortsatt bruker normalisert CIF og non-protonated `analysis_export/complex_for_prolif.pdb`.

### Steg 4: Hard QC og kvalitetsgrenser

#### Status
Ferdig

#### Formaal
Hard QC bestemmer hvilke forberedte poser som kan gaa videre til downstream geometri, analysis export, IFP og clustering. Steget er en hard gate, men det sletter ikke maalinger: pre-QC- og geometriavstander bevares i verdict-data selv naar posen droppes. Koden er fasit for rekkefolgen: aktiv-sete-proximity, Cu-His/substratgeometri, PoseBusters og deretter Privateer-batch for eligible poser.

#### Input
- `PreparedPose` fra steg 3, oversatt til `HardQCInput` med `pose_id`, `for_posebusters.pdb`, in-memory struktur og `privateer_input.cif`.
- QC-flagg fra production-config: `run_posebusters`, `run_privateer`, `n_jobs` og `collect_timing_events`.
- Terskler fra `configs/thresholds.yaml`, lastet via `load_gate_config_from_yaml()`: aktiv-sete hard cutoff, Cu-His hard/soft band, Cu-C soft flag og his-brace search radius.

#### Hva skjer
1. `run_analysis_core()` bygger en liste `HardQCInput` og kaller `run_hard_qc()` for alle prepared poser. Rapporten skrives til `qc_report.json`, valideres mot schema, og verdicts legges tilbake paa `case_by_pose_id`.
2. `run_hard_qc()` kjorer `check_active_site_proximity()` for hver pose. Manglende Cu, manglende glykanatomer eller `min_cu_ligand_distance` over hard cutoff dropper posen foer geometri, PoseBusters og Privateer.
3. Hvis pre-QC passerer, kjorer `check_geometry()`. Cu-His hard range er `1.5-3.0 A`; preferred/soft band er `1.8-2.6 A`. Hard range-feil dropper posen. Preferred-band-feil blir soft warning dersom hard range fortsatt passerer.
4. Bare poser som passerer pre-QC og geometri er `qc_eligible` for PoseBusters og Privateer. PoseBusters kjorer per pose hvis `run_posebusters=True`. Uventet PoseBusters-exception blir naa en eksplisitt hard fail med `posebusters_runner_error`.
5. Privateer kjorer som batch for eligible poser med `privateer_input.cif` hvis `run_privateer=True`. Batch-exception eller manglende resultat for en eligible pose failer lukket med `privateer_batch_runner_error` eller `privateer_missing_result`.
6. `compute_verdict()` samler proximity, PoseBusters, Privateer og geometri til `passed`, `flagged` eller `dropped`. `flagged` teller som beholdt pose videre; bare `dropped` stoppes foer downstream analyse.

#### Hvorfor det gjoeres
- Aktiv-sete-proximity ligger forst for aa unngaa dyre og misvisende kjemiverktoy paa poser der ligand/Cu allerede er fysisk for langt unna.
- Cu-His-geometri kjorer foer PoseBusters i dagens kode fordi det er en LPMO-spesifikk hard gate og avgjor `qc_eligible` for resten av hard-QC-stakken.
- Fail-closed for aktive backendfeil er viktig for hovedanalysen: en manglende PoseBusters- eller Privateer-kjoring skal ikke se ut som en kjemisk pass.
- `run_posebusters=False` og `run_privateer=False` er fortsatt bevisste deaktiveringer, ikke fail-closed-feil.

#### Output
- `qc_report.json` med `poses` og `verdicts`, validert mot `schemas/qc_report_schema.json`.
- `analysis_summary["qc_counts"]` med total/pass/flagged/dropped.
- `case["qc_verdict"]` per pose, som senere brukes av `pose_manifest.tsv`, `qc_attrition_table.tsv` og downstream-gating.
- `qc_attrition_table.tsv` aggregerer hard-fail og soft-flag reasons per condition.
- Timing-events for hard-QC-delsteg hvis `collect_timing_events=True`.

#### Relevante filer og linjer
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1246-L1252` - oversetter `PreparedPose` til `HardQCInput`.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1871-L1907` - kjorer hard QC, skriver/validerer `qc_report.json` og lagrer `qc_counts`.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1925-L1940` - legger verdict paa case og stopper bare `dropped` poser fra downstream geometri.
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1414-L1495` - bygger `qc_attrition_table.tsv` fra verdict-status, drop reasons og warnings.
- `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:L113-L152` - offentlig hard-QC-kontrakt og fail-closed-dokumentasjon.
- `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:L183-L301` - per-pose proximity, geometri og PoseBusters-rekkefolge.
- `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:L350-L415` - Privateer-batch, missing-result og batch-exception-policy.
- `src/lpmo_pipeline/qc/qc_report.py:L71-L180` - verdict-aggregasjon og endelig status.
- `src/lpmo_pipeline/qc/active_site_proximity.py:L62-L164` - pre-QC proximity-metrikker og hard cutoff.
- `src/lpmo_pipeline/qc/custom_geometry_checks.py:L100-L294` - Cu-His hard/soft geometri og Cu-C1/C4-maalinger.
- `src/lpmo_pipeline/qc/gates.py:L243-L280` - faktisk YAML-mapping for hard/soft QC-terskler.
- `tests/test_hard_qc_orchestrator.py:L194-L337` - dekker soft geometri, config-overstyring og fail-closed backendfeil.

#### Samsvar med stottedokumenter

| Dokument/kommentar | Hva den sier | Bekreftet i kode? | Avvik/usikkerhet | Linjereferanser |
|---|---|---|---|---|
| `README.md` | Hard QC har aktiv-sete pre-gate, PoseBusters dock-mode og fail-closed backendfeil | Ja | README beskriver ikke alle verdictfeltene; schema/kode er fasit | `README.md`; `src/lpmo_pipeline/qc/hard_qc_orchestrator.py:L113-L152` |
| `MASTERPLAN.md` | Pre-QC, PoseBusters, Privateer og Cu-His er hard/soft gates med metrics beholdt | Ja | Stage-navnene er grovere enn dagens kode, men policyen stemmer etter oppdatering | `MASTERPLAN.md`; `src/lpmo_pipeline/qc/qc_report.py:L71-L180` |
| `IMPLEMENTATION_PLAYBOOK.md` | Pre-QC, PoseBusters, Privateer og QC-report er implementert og real-data-verifisert | Ja | Eldre stoppunkt om pose-manifest for droppede metrics er delvis avklart av dagens `qc_report`/attrition-flater, men final rapportflate kan fortsatt diskuteres | `IMPLEMENTATION_PLAYBOOK.md`; `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1414-L1495` |
| `schemas/qc_report_schema.json` | `qc_report.json` krever per-pose PoseBusters, Privateer, Cu-geometri og overall status | Ja | Schemaet er bredt og tillater ikke alle interne verdict-detaljer; `verdicts`-delen er runtime-ekstra | `schemas/qc_report_schema.json`; `src/lpmo_pipeline/qc/qc_report.py:L224-L287` |

#### Usikkerheter
- Privateer-input er fortsatt en validert kopi av normalized CIF. Hvis dette senere blir en reell transformasjon, maa hard-QC-beskrivelsen oppdateres.
- `qc_report_schema.json` beskriver den offentlige rapportflaten, men de rikeste interne feltene ligger i `verdicts` og case-metadata. En endelig rapporteringskontrakt for alle droppede-pose-metrikker kan fortsatt strammes.
- Hard-QC-rekkefolgen er naa dokumentert som kode faktisk kjorer den. Eventuelle eldre dokumenter som sier PoseBusters foer geometri er utdaterte.
