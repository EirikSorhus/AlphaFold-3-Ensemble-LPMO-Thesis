---
applyTo: "Masteroppgave/pipe_test/scripts/module_2_new/**"
---

# module_2_new – Repo-instruksjoner for Copilot og utviklere

> **FROSSEN MODUL** – Versjon 1.2 (2026-02-20). Kun kritiske bugfix aksepteres.
> Ny funksjonalitet og refaktorering skal IKKE gjøres uten eksplisitt godkjenning.

---

## Formål

Module 2 er en automasjonsmodul for metadata-berikelse av LPMO-proteiner (Lytic Polysaccharide Monooxygenases). Den henter og annoterer proteindata fra flere offentlige databaser:

- **UniProt** – organisme, proteinnavn, EC-nummer, sekvens
- **InterPro** – domeneannotasjoner (Pfam, CDD, SMART, PROSITE)
- **CAZy** – familietilhørighet (AA9, AA10, AA11, AA13, AA14, AA15, AA16, AA17 m.fl.)
- **NCBI BLAST** – sekvenssøk for ukjente FASTA-headere (valgfritt, kun for små batch)

Modulen støtter fire inputmodi og produserer en strukturert metadata-tabell (TSV) + FASTA-sekvenser for videre analyse.

---

## Nøkkelstier og struktur

```
Masteroppgave/pipe_test/
├── scripts/module_2_new/
│   ├── main_driver.py          # Entrypoint og orkestrering (4 modi)
│   ├── input_handler.py        # Parsing: FASTA, CAZy TSV, characterized CSV, ID-liste
│   ├── uniprot_client.py       # UniProt REST API-klient + JSON-cache
│   ├── interpro_client.py      # InterPro REST API-klient + JSON-cache
│   ├── feature_parser.py       # Domeneutvinning, LPMO-kjerne, H1-verifikasjon
│   ├── blast_client.py         # NCBI BLAST API (polling, maks 60 forsøk = 10 min)
│   ├── discovery_cache.json    # UniProt-cache (ligge i skriptmappen, absolutt stioppslag)
│   └── interpro_cache.json     # InterPro-cache (ligge i skriptmappen, absolutt stioppslag)
├── run_module_2_new_cazy.sh         # SBATCH-runner: cazy-modus
├── run_module_2_new_fasta.sh        # SBATCH-runner: fasta-modus
├── run_module_2_new_list.sh         # SBATCH-runner: list-modus
├── run_module_2_new_characterized.sh # SBATCH-runner: characterized-modus
└── data/
    ├── metadata/               # Output: metadata_expanded_{timestamp}.tsv
    ├── sequences/              # Output: all_sequences_{timestamp}.fasta
    ├── run/                    # Output: failed_ids_{timestamp}.txt, run_metadata_{timestamp}.json
    ├── cazy_raw/               # Nedlastede CAZy-råfiler (kun cazy-modus)
    └── test_data/              # Testinputfiler (AA15.txt, characterized_AA17.csv, m.fl.)
```

---

## Miljø og avhengigheter

**Miljø:** `conda/lpmo_pipe_env/` (aktiveres via `export PATH=.../lpmo_pipe_env/bin:$PATH`)

**Python:** 3.8 eller nyere

**Tredjepartsbiblioteker:**
| Pakke | Minimumsversjon |
|-------|----------------|
| `requests` | 2.31.0 |
| `pandas` | 2.0.0 |
| `tqdm` | 4.65.0 |

**Standardbibliotek brukt:** `argparse`, `logging`, `os`, `sys`, `re`, `json`, `hashlib`, `datetime`, `xml.etree.ElementTree`, `tempfile`

**SLURM-ressurser (typisk):**
```
#SBATCH --account=nn1003k
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
```

---

## Kjøring og bruk

### Grunnleggende syntaks

```bash
cd Masteroppgave/pipe_test
python scripts/module_2_new/main_driver.py --mode <MODE> --input <INPUT_FILE> [VALG]
```

### Alle parametre

| Parameter | Påkrevd | Standard | Beskrivelse |
|-----------|---------|---------|-------------|
| `--mode` | Ja | – | `fasta`, `cazy`, `list`, `characterized` |
| `--input` | Ja | – | Inputfil. I `cazy`-modus: kommaseparerte familienavn (`AA9,AA11`) eller filsti |
| `--output_dir` | Nei | `data` | Rotmappe for output |
| `--cazy-family` | Betinget | – | Påkrevd for `fasta`, `list`, `characterized`. Valgfri for `cazy` (auto-detekteres fra filnavn) |
| `--allow-ncbi-fallback` | Nei | Av | Aktiverer NCBI fallback for GenBank-IDer (kun `characterized`-modus) |
| `--allow-sequence-search` | Nei | Av | Aktiverer BLAST-søk for ukjente FASTA-headere (`fasta`-modus) |
| `--max-sequence-searches` | Nei | 20 | Maks antall BLAST-søk. Hold lavt (<20). |

### Modi

**cazy** – henter direkte fra CAZy.org eller fra lokal CAZy-eksportfil:
```bash
python scripts/module_2_new/main_driver.py \
  --mode cazy \
  --input AA15 \            # Eller: data/test_data/AA15.txt
  --output_dir data
```

**fasta** – metadata for proteiner i FASTA med UniProt-headere:
```bash
python scripts/module_2_new/main_driver.py \
  --mode fasta \
  --input sequences.fasta \
  --cazy-family AA9 \
  --output_dir data
```

**list** – metadata for en liste UniProt-accessions (én per linje):
```bash
python scripts/module_2_new/main_driver.py \
  --mode list \
  --input protein_ids.txt \
  --cazy-family AA9 \
  --output_dir data
```

**characterized** – semikolonseparert CSV (CAZy characterized-format):
```bash
python scripts/module_2_new/main_driver.py \
  --mode characterized \
  --input data/test_data/characterized_AA17.csv \
  --cazy-family AA17 \
  --output_dir data
```

### Via SBATCH

Bruk runner-skriptene i `Masteroppgave/pipe_test/`:
```bash
sbatch run_module_2_new_cazy.sh
sbatch run_module_2_new_characterized.sh
sbatch run_module_2_new_fasta.sh
sbatch run_module_2_new_list.sh
```

Skriptene setter arbeidsmappe til `pipe_test/` og aktiverer `lpmo_pipe_env` før kjøring.

---

## Input/Output og lagring

### Inputformater

| Modus | Filformat | Eksempel på ID-format |
|-------|-----------|----------------------|
| `fasta` | FASTA med UniProt-headere (`>sp|P12345|...` eller `>tr|A0A123|...`) | UniProt accession i header |
| `cazy` | CAZy TSV (tab-separert: Family, Kingdom, Organism, Accession, Source) eller familienavn | `AA15`, `XP_123456`, JGI-IDer |
| `list` | Ren tekst, én UniProt-accession per linje | `P12345`, `A0A1B2C3D4` |
| `characterized` | Semikolonseparert CSV med kolonner: `Protein Name;EC#;Reference;Organism;GenBank;UniProt;PDB/3D` | `P12345`, `WP_123456` |

**Merk:** GenBank → UniProt crossref-mapping har kjent lav suksessrate. `--allow-ncbi-fallback` anbefales for `characterized`-modus når GenBank-IDer er eneste referanse.

### Output

Alle outputfiler er tidsstemplet (`YYYYMMDD_HHMMSS`). Mapper opprettes automatisk.

| Fil | Sti | Innhold |
|-----|-----|---------|
| `metadata_expanded_{ts}.tsv` | `{output_dir}/metadata/` | 17 kolonner (se under) |
| `all_sequences_{ts}.fasta` | `{output_dir}/sequences/` | Normaliserte headere: `>UniProtID|{ID}|{Org}|{Name}` |
| `failed_ids_{ts}.txt` | `{output_dir}/run/` | IDer som ikke ble løst |
| `run_metadata_{ts}.json` | `{output_dir}/run/` | Kjørestatistikk, filbaner, duration |

### Metadata TSV – kolonner

| Kolonne | Type | Beskrivelse |
|---------|------|-------------|
| `UniProt_ID` | str | Primær UniProt-accession |
| `CAZy_family` | str | Familienavn f.eks. `AA9` |
| `InterPro_IDs` | str | Semikolonseparerte InterPro-IDer |
| `Protein_Name` | str | Fullt proteinnavn |
| `Organism` | str | Vitenskapelig artsnavn |
| `EC_Number` | str | EC-nummer |
| `Signal_End` | int | Signalpeptid-kløyvingsposisjon (0 = ingen) |
| `Transmembrane_Regions` | str | TM-regioner som `start-end; start-end` eller `None` |
| `LPMO_Core_Start` | int | LPMO-domene start (1-indeksert) |
| `LPMO_Core_End` | int | LPMO-domene slutt (1-indeksert, inklusiv) |
| `Domain_Provenance` | str | Kilde: `DATABASE:MODEL` eller `DATABASE:MODEL\|original:start-end` |
| `Binding_Modules` | str | CBM-info: `Name [Source:Model] (start-end)` eller `None` |
| `H1_Verified` | bool | True hvis H1-histidin verifisert |
| `H1_AminoAcid` | str | Faktisk aminosyre ved H1-posisjon |
| `Match_Status` | str | `Success_CAZy_NCBI`, `Success_CAZy_JGI`, `UniProt`, `BLAST_exact` m.fl. |
| `Sequence_Group` | str | Primær-ID for sekvensgruppe (`characterized`-modus) |
| `CoAccessions` | str | Semikolonseparerte ko-accessioner (`characterized`-modus) |

### Cache-håndtering

`discovery_cache.json` (UniProt) og `interpro_cache.json` (InterPro) ligger i `scripts/module_2_new/` og brukes uavhengig av arbeidsmappe (absolutt stioppslag). Skriving er atomisk via `tempfile.NamedTemporaryFile` + `os.replace()`.

---

## Konvensjoner og kvalitetskrav

**Logging:**
- Nivå INFO til stdout (tidsmerket, format: `%(asctime)s - %(levelname)s - %(message)s`)
- Kun oppsummerende INFO-meldinger; detaljert deduplisering på DEBUG-nivå
- Feilede IDer logges og samles i `failed_ids_{ts}.txt`
- Kritiske feil logger full traceback og avslutter med exit-kode > 0

**Feilhåndtering:**
- Global try-except rundt `main()` sikrer at delvise resultater skrives til disk ved krasj
- BLAST polling: maks 60 forsøk (10 min), feiler deretter med timeout-feil
- Cache-skriving er atomisk (ingen korrupte cacher ved krasj)
- InterPro signalend-validering: None-sjekk, typevalidering og grenseverdikontroll

**Konfig:**
- Ingen konfigurasjonsfil – alle parametere via CLI og hardkodede konstanter i `main_driver.py`:
  - `CAZY_BASE_URL = "https://www.cazy.org"`
  - `USER_AGENT = "LPMO-Pipeline/1.0 (eirik.sorhus@nmbu.no)"`
  - `REQUEST_TIMEOUT = 30` (sekunder)

**Navngiving:**
- Output-filer er alltid tidsstemplet for å unngå overskrivning
- Match_Status-verdier dokumenterer proveniensen til hvert treff
- `Domain_Provenance` skiller mellom justerte og ujusterte domener

---

## Do/Don't for endringer

### Do

- **Kun kritiske bugfix:** Rett opp feil som fører til krasj, datakorrupsjon, eller åpenbart feil output
- Bevar atomisk cache-skriving (tempfile + os.replace) ved endringer i `uniprot_client.py` / `interpro_client.py`
- Behold konsistens mellom `H1_Verified` og `H1_AminoAcid` (alltid synkron boolean-tildeling)
- Test mot `pipe_test/data/test_data/` etter bugfix
- Dokumenter fix i `CHANGES.md` med dato, versjonsnummer, symptom, fix og effekt

### Don't

- **Ikke legg til ny funksjonalitet** – modulen er frossen (v1.2)
- Ikke endre output-kolonner, kolonnenavn eller rekkefølge – dette bryter nedstrømsbruk i `analyse`-modulen
- Ikke endre cachefilens plassering eller format uten å migrere eksisterende cache
- Ikke øk `--max-sequence-searches` over 20 i produksjonskjøringer (BLAST-API-belastning)
- Ikke aktiver `--allow-ncbi-fallback` uten å forstå at det gjør HTTP-kall til NCBI Entrez
- Ikke refaktorer noe kun for stilens skyld – risikerer regresjoner i en ellers stabil modul

---

## Referanser i repoet

Følgende stier ble faktisk brukt som kilde for denne instruksjonsfilen:

- [Masteroppgave/pipe_test/scripts/module_2_new/README.md](../pipe_test/scripts/module_2_new/README.md) – fullstendig bruksdokumentasjon, modi, kolonner, kjente feil
- [Masteroppgave/pipe_test/scripts/module_2_new/CHANGES.md](../pipe_test/scripts/module_2_new/CHANGES.md) – v1.1 og v1.2 bugfix-historikk
- [Masteroppgave/pipe_test/scripts/module_2_new/CODE_WALKTHROUGH.md](../pipe_test/scripts/module_2_new/CODE_WALKTHROUGH.md) – funksjonsbeskrivelser, risikopunkter, arkitektur
- [Masteroppgave/pipe_test/scripts/module_2_new/main_driver.py](../pipe_test/scripts/module_2_new/main_driver.py) – entrypoint, CLI-argumenter, output-struktur, konstanter
- [Masteroppgave/pipe_test/run_module_2_new_cazy.sh](../pipe_test/run_module_2_new_cazy.sh) – SBATCH-innstillinger, miljøoppsett, arbeidsmappe
- [Masteroppgave/pipe_test/run_module_2_new_characterized.sh](../pipe_test/run_module_2_new_characterized.sh) – characterized-modus kjøringseksempel
