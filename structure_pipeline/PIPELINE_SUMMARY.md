# Structure Pipeline Summary
Forfatter: OpenAI ChatGPT 5.4 high thinking
## Sammendrag

Structure Pipeline er en manifestdrevet orkestreringspipeline for strukturprediksjon med AlphaFold3, Boltz-2 og RoseTTAFold3. Den er bygget for HPC-kjoring med SLURM og containere, og hovedideen er at raadata ikke brukes direkte i submit-steget. I stedet oversettes FASTA, ligandkilder og MSA-kilder forst til eksplisitte manifestfiler, deretter til standardiserte `Case`-objekter, og til slutt til grupperte jobber per ligand og modell. Det gjor flyten sporbar, enklere a resumere, og lettere a inspisere nar noe gaar galt.

Det viktigste a forsta er at pipelinen har en tydelig arbeidsdeling. `cli.py` styrer hele orkestreringen. `config.py` gjor miljo og ressursvalg eksplisitte. `manifest/*` normaliserer inputdata til CSV-manifester. `cases.py` utvider dette til konkrete kjoreoppgaver. `runners/*` bygger modellspesifikke inputfiler og SLURM-skript. `executors/*` skiller submit-mekanikken fra jobbinnholdet. `resume.py` gir en enkel, filbasert mekanisme for a hoppe over ferdige grupper. Samlet er dette en kodebase som er mer opptatt av robust batch-kjoring enn av interaktiv modellbruk.

AF3, Boltz og RF3 behandles ikke likt. AF3 har et bevisst to-stegs oppsett med MSA-forbehandling paa CPU og inferens paa GPU, fordi det gir gjenbruk av MSA-arbeid og bedre ressursutnyttelse. Boltz og RF3 kjores som ett jobbtrinn per ligand, der alle proteiner for samme ligand batch-es sammen. Pipelinen legger ogsaa inn modellspesifikke detaljer som obligatorisk `CU`-ligand, oligosakkaridstotte for AF3, og forskjellig inputformat for JSON versus YAML. Disse valgene virker smaa isolert sett, men de er helt sentrale for at hele systemet skal oppfore seg konsistent paa tvers av mange jobber.

Denne oppsummeringen forklarer forst hva pipelinen gjor ende til ende, deretter de viktigste modulene og funksjonene, og til slutt hvorfor designvalgene ser ut slik de gjor i dagens kode.

## Nylig Analyse-Notat

Endringene 2026-05-19 ligger i downstream `analyse`-pipelinen, ikke i selve strukturprediksjonskjoringen. Strukturpipelinen produserer fortsatt AF3-output under `work_core` og `work_full_length`; analysepipelinen bruker disse outputene videre.

For clustering-piloten er analysewrapperen gjort Slurm-first: `submit_clustering_pilot_staged.sh` splitter domain-only og full-length selections i protein-level shards og submitter dem som Slurm-arrays, slik at pilotarbeidet kan fordeles over flere jobber/noder. Normalisering i analysepipelinen er ogsaa raskere fordi CIF chain-remap naa grupperes per mmCIF-loop i stedet for a skrive samme loop om igjen for hver tag. Rutinemessige vellykkede normaliseringer skriver bare `normalize_report.json`; `atom_map.tsv` og `rename_log.json` er debug/failure-artefakter.

## Hva Pipelinen Gjor

Pipelinen tar tre typer input og gjor dem om til kjoreplaner for strukturprediksjon:

- proteiner fra FASTA
- ligander enten fra CIF-filer eller en ren CCD-liste
- MSA-kilder for modeller som krever pre-komputert MSA

Ut av dette produserer den fire manifestfiler:

- `proteins.csv`
- `ligands.csv`
- `msa.csv`
- `cases.csv`

Disse filene fungerer som et kontraktslag mellom inputbehandling og eksekvering. Det er et viktig arkitekturvalg. Ved a materialisere mellomresultater til CSV blir det mulig a se noyaktig hvilke proteiner, ligander og MSA-koblinger systemet tror finnes, og hvilke faktiske kombinasjoner som skal kjores. Det er mye enklere a debugge en pipeline der planleggingen er eksplisitt skrevet ut enn en pipeline som genererer alle kombinasjoner i minnet ved submit-tid.

Et `Case` representerer en konkret prediksjonsoppgave: ett protein, ett ligandvalg og en modell. En stor del av systemets enkelhet kommer av at nesten all videre logikk opererer paa lister av `Case`-objekter. Sa snart input er normalisert til `Case`, blir resten av problemet hovedsakelig batching, filbygging og submit-logikk.

## Ende-til-Ende Flyt

Den praktiske flyten er styrt fra `cli.py`, og kan leses som en serie tydelige transformasjoner.

Forst brukes `init-config()` til a lage et eksempel paa konfigurasjonsfil. Dette er ikke avansert logikk, men funksjonen er viktig fordi den eksplisitt viser hvilke stier og ressurser pipelinen forventer. Konfigurasjonen er laget slik at runtime-miljoet skal beskrives i YAML i stedet for a vaere skjult i shellskript eller hardkodede antakelser.

Deretter brukes `manifest()` som pipeline-start. Den gjor fire ting i fast rekkefolge:

1. Leser FASTA og bygger proteinmanifest.
2. Leser ligander enten ved katalogskanning eller fra en CCD-liste.
3. Matcher MSA-er mot proteinidentiteter.
4. Genererer alle cases for alle aktiverte modeller.

Det sentrale poenget er at `manifest()` ikke submitter noe. Kommandoen bygger bare et eksplisitt execution plan. Denne separasjonen er klok fordi den skiller inputkvalitet fra ressursbruk. Hvis inputene er feil, finner man det i manifestlaget for man bruker GPU-tid.

`validate()` er neste stopp. Denne kommandoen sjekker at inndata, containere, checkpoints, vektfiler og andre modellspesifikke stier finnes. Den gjor i praksis en miljoverifisering. Det ser trivielt ut, men paa HPC-systemer er dette ofte forskjellen mellom rask avklaring og lange koer for jobber som feiler direkte ved oppstart.

Selve kjoringen skjer i `run()`. Her lastes manifestfilene inn igjen, cases kan filtreres paa modell eller protein, resume-filteret brukes, og de gjenværende casene grupperes med `group_cases_by_ligand_model()`. Det betyr at submit-enheten i dagens kode i hovedsak er `(ligand_ccd_code, model)`, ikke `(protein, ligand, model)`. Dette er et bevisst valg som reduserer antall jobber og gjor modellinnlasting billigere per protein.

AF3 er et spesielt tilfelle. `run()` samler forst opp hvilke proteiner som trenger AF3-MSA, submitter MSA-jobber separat, og submitter deretter inferensjobber per ligand med `afterok`-avhengigheter naar `SlurmExecutor` brukes. Boltz og RF3 submitter ett jobbskript per ligand/modell-gruppe uten denne ekstra fan-in-mekanikken.

Til slutt brukes `status()` til a oppsummere hva som er ferdig. I dagens implementasjon er dette en grov gruppestatus basert paa mappestruktur og `DONE.ok`, ikke en komplett case-for-case rapportering.

## Viktige Moduler

### CLI og Orkestrering

`cli.py` er den reelle kontrollflaten for hele systemet. Den er viktigere enn de fleste andre filer fordi den bestemmer hvilke objekter som kobles sammen og i hvilken rekkefolge.

De viktigste funksjonene er:

- `get_config()`: wrapper rundt `load_config()` med brukerrettet feilhandtering.
- `manifest()`: bygger de fire manifestfilene og viser en modellvis oppsummering av pending og skipped cases.
- `validate()`: validerer input- og modellstier via runnernes `validate_paths()`.
- `run()`: laster manifests, filtrerer, grupperer, setter opp runner- og executor-objekter og submitter jobber.
- `status()`: viser aggregert status per modell.
- `init_config()`: skriver en standard YAML-mal.

`run()` er den klart viktigste funksjonen i hele pipelinen. Den inneholder de viktigste kontrollvalgene:

- CLI-overstyringer skrives direkte tilbake til configobjektet for kjoreplanen bygges.
- resume skjer for submit, ikke etter submit.
- gruppering skjer per ligand og modell.
- AF3 har egen MSA-fase og separat inferensfase.
- oligosakkaridregister lastes bare dersom bruker eller config peker til en YAML-fil.

Et viktig, konkret detaljpunkt er at AF3-MSA-jobbene i dagens kode grupperes per protein, ikke per protein og ligand. `af3_msa_jobs` bruker `protein_id` som nokkel, og `af3_msa_ligand.setdefault(case.protein_id, ligand_ccd)` velger forste ligand som ble sett for det proteinet. Dette betyr at MSA-fasen bruker ett valgt ligandoppsett per protein, mens inferensfasen senere patcher JSON-en tilbake til riktig ligand for hver jobb. Det er effektivt, men det er ogsaa en subtil implementasjonsdetalj som er viktig a kjenne til hvis man forventer full ligandspesifikk separasjon i MSA-steget.

### Konfigurasjon

`config.py` bruker Pydantic-modeller til a beskrive hele runtime-miljoet. Dette er en sterk del av designet fordi det tvinger frem eksplisitthet. I stedet for a spre stier, partinavn og modellparametre rundt i runnerne, samles de i `PipelineConfig` med underseksjoner for `paths`, `slurm`, `af3`, `boltz`, `rf3`, `inputs` og `outputs`.

Noen viktige funksjoner og detaljer:

- `PathsConfig` samler containere, vekter, checkpoints og installasjonsrotter.
- `SlurmConfig.get_mem_per_gpu(model)` lar hver modell overstyre standard GPU-minne uten at runnerne maa vite om fallback-logikk.
- `InputsConfig.resolve_relative_paths()` og `OutputsConfig.resolve_relative_paths()` gjor relative stier absolutte etter lasting.
- `PipelineConfig.create_output_dirs()` oppretter outputmapper automatisk.
- `PipelineConfig.validate_models()` begrenser modeller til `af3`, `boltz` og `rf3`.

Hvorfor er dette gjort slik? Fordi pipelinen er avhengig av mange miljo-spesifikke stier. Hvis disse ikke normaliseres tidlig, blir runnerne fulle av lokal sti-logikk og feilsoking blir mye vanskeligere. Ved a sentralisere denne delen i Pydantic kan runnerne fokusere paa input- og kommandobygging.

Det finnes ogsaa noen nyere CLI-overstyringer i `run()` som det er verdt a merke seg: `af3_num_recycles`, `boltz_sampling_steps` og `boltz_use_potentials`. De viser at konfigurasjonen ikke bare er statisk YAML, men ogsaa en struktur som kan modifiseres ved submit-tid.

### Manifestlag

Manifestlaget bestaar av `manifest/proteins.py`, `manifest/ligands.py` og `manifest/msa.py`. Disse modulene er viktige fordi de oversetter varierende inputformater til stabile, serialiserbare records.

I `proteins.py` er `parse_fasta()` hovedfunksjonen. Den leser FASTA sekvensielt og bygger `ProteinRecord` med baade sekvens, lengde, SHA256-hash og metadata fra headeren. `_parse_header()` er mer interessant enn den ser ut: den stotter baade foretrukket pipe-separert format, eldre underscore-format, klassiske UniProt-headere og en siste fallback. Poenget er ikke bare robusthet, men at pipeline-identiteten `protein_id` maa kunne utledes paalitelig selv om inputformatet varierer.

I `ligands.py` er `scan_ligands()` den viktigste funksjonen. Den stotter flere katalogkonvensjoner, hopper over uinteressante mapper som `boltz_ccd_lib`, trekker ut CCD-kode fra CIF-innhold hvis mulig, og faller ellers tilbake til filnavn. `load_ccd_list()` gir et alternativt spor der ligander beskrives uten CIF-filer. Hvorfor er dette nyttig? Fordi noen kjoreoppsett bare trenger standard CCD-koder, og da er det unodvendig og dyrt a vedlikeholde en lokal CIF-katalog.

`msa.py` er det mest subtile manifestlaget. `match_msa()` bygger et indekskart fra potensielle UniProt-ID-er til MSA-kilder ved a bruke baade filnavn og A3M-header. Dette er et fornuftig valg, fordi MSA-filer i praksis ofte er ujevnt navngitt. Koden er tolerant: ved flere treff logger den warning og velger en deterministisk sortert foerste kandidat, med mindre `strict=True` brukes.

Det er ogsaa viktig at `match_msa()` stotter to kilder:

- vanlige `.a3m`-filer i kataloger
- `.a3m` eller `.a3m.gz` inne i squashfs-bilder via `unsquashfs -l`

Dette viser at manifestlaget er mer kapabelt enn runner-laget paa MSA-siden. Dagens kode kan indeksere squashfs-baserte MSA-kilder, men Boltz- og RF3-runnerne bruker forelopig ikke `msa_sqsh_file` aktivt i jobbskriptene. Det betyr at inputmatching er mer ferdig enn full runtime-integrasjon for squashfs-sporet.

### Cases og Gruppering

`cases.py` er stedet der abstrakt input blir til konkret arbeid. Dette er en liten fil med stor betydning.

`Case`-dataklassen definerer standardarbeidsenheten. Den inneholder identitet, modell, MSA-kilde, MSA-sti, status og eventuelt feiltekst. `output_dir_name` gir et standardnavn i formatet `{protein_id}_{ligand_ccd_code}`, som flere downstream antakelser bygger paa.

`_generate_case_id()` bruker SHA256 over `protein_id`, `ligand_id`, modell og MSA-kilde, og kutter ned til 12 heksadesimale tegn. Testene bekrefter at dette er deterministisk og stabilt. Hvorfor hash her? Fordi man vil ha en kort, reproducerbar ID uten a basere seg paa skjore posisjonsindekser eller mutable filnavn.

`generate_cases()` er sentral. Den gaar gjennom produktet av proteiner, ligander og modeller, og legger inn modellspesifikk MSA-logikk:

- AF3 faar alltid `msa_source="native"` og tom `msa_path`.
- Boltz og RF3 forventer pre-komputert MSA.
- proteiner i `missing_msa_proteins` blir generert som `SKIPPED` for ikke-AF3-modeller.

Dette er et godt eksempel paa at business-logikken ligger tidlig i flyten. I stedet for a la runnerne krasje paa manglende MSA, materialiseres dette som eksplisitt `SKIPPED` allerede i case-laget. Det er enklere a forklare, teste og oppsummere.

To grupperingsfunksjoner finnes:

- `group_cases_by_protein_model()`
- `group_cases_by_ligand_model()`

I dagens submitflyt er det `group_cases_by_ligand_model()` som er strategisk viktig. En jobb representerer ett ligand og en modell, med alle relevante proteiner inni. Testene bekrefter at `SKIPPED`-cases filtreres bort fra gruppering. Det reduserer antall submit-enheter og matcher hvordan runnerne er designet til a batch-e input.

### Runners

Runnerlaget er arkitekturens viktigste abstraksjon etter `Case`. `RunnerInterface` i `runners/base.py` beskriver kontrakten alle modeller maa oppfylle:

- `build_input_file()`
- `build_command()`
- `build_slurm_script()`
- `parse_outputs()`
- `validate_paths()`

Dette er et klassisk og fornuftig skille. CLI-laget skal vite at en runner finnes, men ikke hvordan AF3-JSON, Boltz-YAML eller RF3-JSON ser ut. Ved a tvinge alle modeller inn i samme interface kan submitlogikken holde seg forholdsvis generisk selv om detaljene under er veldig forskjellige.

`RunnerResult` viser ogsaa hva pipelineforfatteren anser som viktige utdata: suksess, outputfiler, confidence-scorer, runtime og metadata. Det er i praksis den strukturerte kontrakten mellom modelloutput og eventuell senere analyse.

Det er ogsaa verdt a se runnerlaget som mer enn bare de tre hovedklassene. `ligand_utils.py` samler modellspesifikk bygging av ligandrepresentasjoner og haandterer den delte antakelsen om obligatorisk `CU`. `oligo.py` beskriver et eget lite domenelag for oligosakkarider, med `OligoRegistry`, `OligoSpec`, `build_af3_oligo_ligand_entry()` og `build_af3_bonded_atom_pairs()`. `mounts.py` inneholder `compute_bind_mounts()`, som beregner minimale bind-prefix for dry-run-inspeksjon. Disse helperfilene er viktige fordi de trekker detaljarbeid ut av runnerklassene og hindrer at samme spesiallogikk blir kopiert flere steder.

### Executors

Executorlaget i `executors/base.py`, `executors/slurm.py` og `executors/local.py` er laget for a skille jobbinnhold fra submit-mekanikk. Det er et nyttig abstraheringsnivaa fordi runnerne skal bygge skript, ikke vite hvordan de sendes inn.

`SubmittedJob` holder paa `job_id`, `work_dir`, `script_path` og status, og har en enkel `done_marker`-egenskap som peker til `DONE.ok`. Dette passer godt med resten av systemet, som i stor grad er fil- og mappebasert snarere enn databasebasert.

`SlurmExecutor.submit()` skriver skript til `work_dir`, kaller `sbatch`, parser jobb-ID fra standardoutput og stotter `--dependency`. Det er dette som gjor AF3s to-stegs strategi mulig i praksis.

`LocalExecutor.submit()` skriver ogsaa et skript, men stripper `#SBATCH`-linjer og starter det som lokal subprocess. Dette er nyttig for testing og utvikling, men det er viktig a forsta begrensningen: lokal executor stotter ikke dependencies. For AF3 betyr det at lokal kjoremodus er mindre tro mot den reelle HPC-flyten enn SLURM-modus.

### Resume og Status

`resume.py` implementerer en enkel statusmodell med to primitive byggesteiner:

- append-only `status.jsonl`
- `DONE.ok` som grov ferdigmarkor per ligand/modell

`StatusTracker.log_event()` kan skrive detaljerte hendelser som `started`, `completed` og `failed`, og `get_completed_cases()` kan rekonstruere et sett av ferdige case-ID-er ved a lese JSONL-filen. Det er et bra grunnlag for finmasket resume.

Det viktige i dagens kode er likevel at denne infrastrukturen nesten ikke brukes utenfor `resume.py`. Et sok i kildekoden viser at `log_event()` og `mark_ligand_model_done()` ikke kalles fra `cli.py` eller runnerne. I praksis blir resume derfor dominert av `DONE.ok`-sjekken og av om en case finnes i tidligere JSONL hvis slike logger skulle bli skrevet av annen kode senere.

`get_pending_cases()` er likevel nyttig og riktig plassert. Den tar inn hele case-listen og returnerer bare det som fortsatt trenger submit. Det er bedre enn a bake resume-logikken inn i hver runner, fordi resume her er en pipeline-egenskap, ikke en modellegenskap.

`generate_summary()` er et bevisst enkelt sammendrag paa mappenivaa. Det teller fullforte grupper via `DONE.ok` og pending grupper via eksistensen av `status.jsonl`. Feltet `failed` finnes i outputstrukturen, men inkrementeres ikke i dagens implementasjon. Statuskommandoen er derfor mer en oversikt over fremdrift enn en full feilrapport.

## Modellspesifikke Forskjeller

### AlphaFold3

AF3-runneren i `runners/af3.py` er den mest komplekse modellen i pipelinen, og det er der mye av den mest interessante designlogikken ligger.

Det sentrale er at AF3 er splittet i to trinn:

1. `build_msa_slurm_script()` lager CPU-jobber som bare genererer MSA-relaterte `_data.json`-filer.
2. `build_slurm_script()` lager GPU-jobber som bruker disse filene som grunnlag for inferens.

Dette er et sterkt designvalg av to grunner. For det forste er MSA-generering og inferens veldig ulike ressursproblemer. CPU-jobber og GPU-jobber bor ikke blandes hvis man vil bruke clusterressurser effektivt. For det andre kan MSA-arbeid gjenbrukes mellom flere ligandgrupperinger dersom proteinet er det samme. Derfor er det smart a skille dem.

`build_input_file()` og `build_msa_input_file()` bygger AF3 JSON i versjon 4-format. Begge bruker `build_af3_ligand_entries()`, som alltid legger inn `CU` naar hovedliganden ikke allerede er `CU`. Dette er verifisert i testene og er en gjennomgaaende kontrakt paa tvers av modellene. For vanlige ligander kan `userCCDPath` settes til en custom CIF. For oligosakkarider brukes ikke custom CIF; i stedet bygges en liste med monomer-CCD-er og `bondedAtomPairs`.

En viktig detalj i inferensfasen er `_write_patch_script()`. I stedet for a generere en helt ny inferens-JSON fra bunnen av, skriver runneren et separat Python-skript som patcher AF3s produserte `_data.json`. Dette er et veldig godt teknisk valg. Hvis inferensfasen laget helt nye JSON-filer, ville man risikere a miste MSA-innholdet som CPU-trinnet allerede har beregnet. Patching er derfor ikke bare en implementasjonsdetalj, men selve mekanismen som gjor to-stegsarkitekturen korrekt.

MSA-jobben rsync-er input til scratch, kjoerer AF3 i container med `--run_inference=false`, og kopierer tilbake den produserte `_data.json`-filen som `produced_msa.json` til `work_dir/input/`. Inferensjobben leser denne filen, patcher navn og ligandfelt, kjoerer AF3 med GPU, og skriver output til `runs/{SLURM_JOB_ID}` med `latest`-symlink og `DONE.ok`. Denne per-run-isolasjonen er viktig fordi den unngaar at flere kjoringer trakker i samme outputmappe.

`parse_outputs()` for AF3 leter etter `*_summary_confidences_*.json` og `*_model.cif`. Den henter ut `ptm`, `iptm` og `ranking_score`. Det viser hva pipelineforfatteren oppfatter som de viktigste AF3-metrikkene for videre bruk.

AF3 er ogsaa den eneste modellen som faktisk bruker oligosakkaridregisteret i dagens kode. Det er en viktig domenedetalj: oligo-YAML er ikke et generelt ligandlag for alle modeller, men et AF3-spesifikt hjelpemiddel for aa uttrykke multimonomere ligander paa den maten AF3 forventer.

### Boltz-2

`runners/boltz.py` er enklere enn AF3, men designet er tydelig og praktisk. Hovedideen er ett jobbskript per ligand, og inni dette skriptet ett YAML-input per protein. `build_input_file()` lager derfor en katalog med YAML-filer, ikke en enkelt inputfil.

For hvert protein bygges `sequences` med en proteindel og en eller flere liganddeler. Igjen legges `CU` automatisk inn naar hovedliganden ikke allerede er `CU`, og testene bekrefter dette. Hvis `use_affinity` er aktiv, legges affinity-prediksjon til via YAML-feltet `properties`, der binder peker til hovedligandens ID. Det er et viktig detaljpoeng at affinity her uttrykkes i YAML, ikke som et CLI-flagg. Kommentarene i filen viser at dette er en bevisst korreksjon fra en tidligere antakelse.

`build_command()` og `build_slurm_script()` viser ogsaa den operative modellen for Boltz:

- input er katalogen `/work/input`
- output skrives til `/work/result`
- vekter bind-mountes under `/weights`
- HuggingFace-cache peker til `/weights/huggingface`

Runneren lager ett relativt enkelt GPU-skript som kaller `boltz predict /work/input`. Det passer godt med batch-ideen: modellen lastes en gang, flere YAML-er prosesseres i samme jobb.

`parse_outputs()` forventer i hovedsak en katalogstruktur `predictions/{protein}_{ligand}/` og leter etter `*_model_*.cif`, `confidence_*.json` og eventuelt `affinity_*.json`. Det hentes ut `ptm`, `iptm`, `ligand_iptm`, `confidence_score` og mulig affinity-verdi.

Det finnes likevel en viktig praktisk detalj her: runneren tar imot `msa_sqsh_file`, men bruker det ikke aktivt i skriptet. Kommentarene sier eksplisitt at squashfs ikke brukes i denne versjonen. Det betyr at Boltz-delen i praksis er sterkest for vanlige MSA-filer, selv om manifestlaget kan indeksere squashfs-kilder.

### RoseTTAFold3

`runners/rf3.py` bruker den samme batchideen som Boltz, men med et annet inputformat og et tyngre runtime-miljo. `build_input_file()` lager ett JSON-array med ett eksempel per protein. Hvert eksempel har et navn og en liste `components`, der proteinet alltid er komponent `A`, mens ligandene uttrykkes som CCD-koder eller CIF-bane.

Ogsaa her legger helperlaget automatisk til `CU` naar hovedliganden ikke allerede er `CU`, og testene bekrefter dette. Hvis en ligand-CIF finnes, foretrekkes denne fremfor bare CCD-kode. Det viser at RF3-sporet er mer fleksibelt enn Boltz naar det gjelder custom ligandrepresentasjon.

`build_slurm_script()` er mer kompleks enn Boltz sin fordi RF3 krever Foundry-miljo, checkpoint, Hydra-lignende argumenter og ved behov en sideinstallasjon av `triton==3.5.0` i et temporart site directory. Dette ser litt tungt ut, men det gjenspeiler den faktiske kostnaden ved a fa RF3 til aa kjoere stabilt i container paa cluster. Det er ikke pyntelogikk; det er miljofiksing som er nodvendig for reproducerbar kjoreevne.

Runneren setter flere `APPTAINERENV_*`-variabler, binder Foundry-root og temporarbasen, og kjoerer deretter `rf3/cli.py fold` med eksplisitte parametre for recycles, diffusion batch size, steps, seed og early stopping-threshold. Dette er et godt eksempel paa at pipelinen ikke bare orkestrerer jobber, men ogsaa innkapsler ekspertkunnskap om hvordan det aktuelle modellmiljoet maa settes opp.

`parse_outputs()` leter etter CIF-filer, metrics-CSV og eventuelle `.score`-filer som indikerer early stopping. Den henter primart ut `pTM` og `pLDDT`. Sammenlignet med AF3 og Boltz er denne outputtolkningen enklere, men ogsaa mer avhengig av at filnavngivningen holder seg stabil.

Som med Boltz aksepterer runneren `msa_sqsh_file`, men bruker det forelopig ikke aktivt i jobbskriptet. Det er derfor riktig a se squashfs-stotte i RF3-sporet som delvis planlagt, men ikke fullt realisert i runtime-delen.

## Viktige Designvalg Og Hvorfor De Er Gjort

Flere designvalg gaar igjen gjennom hele pipelinen.

For det forste er systemet manifestdrevet. Det er gjort for sporbarhet og debugbarhet. En CSV-basert kjoreplan er ikke like elegant som en ren in-memory pipeline, men den er mye lettere a verifisere og gjenbruke i praksis.

For det andre er submit-enheten valgt som `(ligand, modell)` heller enn individuelle cases. Det reduserer jobbtall, minimerer gjentatt modellinnlasting og matcher hvordan runnerne bygger batch-input. Det er en ytelsesbeslutning som ogsaa forenkler ressursstyring.

For det tredje er AF3 splittet i CPU-MSA og GPU-inferens. Dette er gjort fordi MSA og inferens har ulike ressursprofiler, og fordi MSA-resultatet kan gjenbrukes. Det er trolig den viktigste optimaliseringen i hele kodebasen.

For det fjerde skilles runner og executor. Det er gjort for at modellspesifikk inputbygging ikke skal blandes med SLURM- eller lokal kjoremekanikk. Denne separasjonen er enkel, men veldig verdifull ved testing og videre utvidelser.

For det femte brukes `DONE.ok`, `latest`-symlink og `runs/{jobid}` som filbasert kontrakt rundt output. Det er en pragmatisk losning som passer godt til batch-kjoring paa cluster, der robust filsystemtilstand ofte er viktigere enn raffinerte metadataoppsett.

For det sjette finnes det et eget helperlag for ligander. `ligand_utils.py` sikrer at alle modellene anvender samme domeneantakelse om obligatorisk `CU`, og at modellspesifikke ligandrepresentasjoner bygges ett sted. Dette er en liten, men god abstrahering som hindrer at samme kjemi-antakelse dupliseres i tre runners.

## CLI-Bruk Og Praktiske Konsekvenser

Fra brukerperspektiv er det verdt a tenke paa pipelinen som tre forskjellige aktiviteter i stedet for en enkelt kommando:

- planlegging med `manifest`
- miljoverifisering med `validate`
- kjoring med `run`

Det er et godt arbeidsmonster. Hvis `manifest` ser riktig ut og `validate` passerer, er det mye mindre sannsynlig at `run` feiler av trivielle grunner.

`run --dry-run` er spesielt nyttig fordi det viser mount-planen fra `compute_bind_mounts()`. Det gir et raskt bilde av hvilke stier som maa vaere tilgjengelige fra containerne. Det er likevel viktig a vite at denne mount-beregningen i dag hovedsakelig brukes til dry-run-inspeksjon. De faktiske Boltz- og RF3-jobbskriptene har eksplisitte bind-mounts i skriptteksten og bruker ikke `compute_bind_mounts()` direkte under vanlig submit.

CLI-overstyringer er ogsaa viktige i praksis. De lar brukeren tune modellparametre uten a endre YAML-filen permanent. Det passer godt til eksperimentell HPC-bruk, der man ofte vil prove en annen seed, flere diffusion samples eller justert recycle-antall uten aa lage nye konfigurasjonsfiler for hver variant.

Det er ogsaa verdt aa merke seg at `ligand_ccd_list`-modus endrer hva runnerne har tilgang til. I denne modusen finnes det ikke lokale CIF-baner i ligandmanifestet, sa modellene maa basere seg paa standard CCD-koder. Det er praktisk og lettvekts, men mindre fleksibelt enn CIF-modus for modeller eller ligander som trenger custom definisjoner.

## Begrensninger Og Ting A Vaere Oppmerksom Pa

Det er flere viktige begrensninger og skarpe kanter i dagens kode som er verdt a kjenne til.

- Pipelinen er bygget rundt single-chain proteiner. Det er ikke generell stotte for paired MSA eller mer komplekse multikjedeoppsett.
- Oligosakkaridregisteret brukes bare av AF3. Boltz og RF3 har ikke tilsvarende oligo-ekspansjon i dagens implementasjon.
- `msa_sqsh_file` matches i manifestlaget, men runtime-stotte i Boltz og RF3 er forelopig ikke fullfort.
- `status()` viser en `failed`-kolonne, men `generate_summary()` fyller den i praksis ikke ut. Statusvisningen er derfor enklere enn den ser ut.
- `StatusTracker.log_event()` finnes, men brukes ikke i dagens submitflyt. Resume-logikken hviler derfor mest paa `DONE.ok` og ikke paa rike hendelseslogger.
- `LocalExecutor` er nyttig for testing, men mangler dependency-stotte og er derfor spesielt mindre troverdig for AF3s to-stegsmodell.

Til tross for disse begrensningene er hovedarkitekturen sterk. Det som er best ved denne pipelinen er ikke at hver enkelt runner er vakker isolert sett, men at systemet som helhet er tydelig strukturert rundt sporbarhet, batching, ressursdeling og praktisk clusterdrift. Det gjor den til en god base for videre arbeid, sa lenge man er bevisst paa hvor implementasjonen allerede er moden, og hvor den fortsatt er litt halvferdig eller sterkt miljoavhengig.
