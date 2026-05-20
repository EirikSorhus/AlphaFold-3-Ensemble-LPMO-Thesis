# Codewalkthrough

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
│       │   └── protonate_export.py                 # protonering og eksport for analyseverktoy
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
- bare poser som ikke er `dropped` av QC gaar videre til downstream geometri og protonering
- bare poser som protoneres vellykket samles opp som input til condition-vis konvergens og ProLIF/IFP
- bare hvis det finnes IFP-resultater skrives `pose_ifp_table.tsv` og `pose_residue_contact_table.tsv`
- bare hvis det finnes QC-pass-poser per condition skrives clustering-relaterte tabeller og cluster-annotasjon
- `clustering_pilot.enabled` legger til ekstra pilot-output, men erstatter ikke den primare agglomerative produksjonsclusteringstien
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
- `src/lpmo_pipeline/analysis/analysis_orchestrator.py:L1732-L2037` - overgangen fra forberedte poser til hard QC, downstream geometri, protonering, konvergens, IFP og contact eligibility
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
