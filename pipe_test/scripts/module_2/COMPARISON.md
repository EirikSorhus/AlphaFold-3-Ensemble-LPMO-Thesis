# Module 2: Sammenligning Originale vs. Forbedrede Filer

## Filstruktur

### Nye filer opprettet:
```
scripts/
├── __init__.py                              [NY]
└── module_2/
    ├── __init__.py                          [NY]
    ├── config_improved.py                   [NY - forbedret versjon]
    ├── uniprot_fetch_improved.py            [NY - forbedret versjon]
    ├── cazy_fetch_improved.py               [NY - forbedret versjon]
    ├── ncbi_fetch_improved.py               [NY - forbedret versjon]
    ├── merge_sequences_improved.py          [NY - forbedret versjon]
    ├── run_module2_improved.py              [NY - forbedret versjon]
    ├── test_improved.py                     [NY - verifikasjonsscript]
    ├── README_IMPROVED.md                   [NY - dokumentasjon]
    └── COMPARISON.md                        [DENNE FILEN]
```

### Originale filer (bevart uendret):
- `config.py`
- `uniprot_fetch.py`
- `cazy_fetch.py`
- `ncbi_fetch.py`
- `merge_sequences.py`
- `run_module2.py`

## Hovedforskjeller

### 1. config.py → config_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Filstier | String-basert (`"data_raw/..."`) | pathlib.Path (`Path(__file__).parent...`) |
| NCBI email | Hardkodet konstant | Funksjon som tar parameter/miljøvariabel |
| Dependency check | Ingen | `check_dependencies()` funksjon |
| Validering | Ingen | `validate_config()` funksjon |
| PROJECT_ROOT | Ikke definert | Automatisk beregnet fra filsti |

**Nye funksjoner:**
- `check_dependencies()` - Sjekker at alle pakker er installert
- `get_ncbi_email(email)` - Henter email fra param/env/raise error
- `validate_config()` - Validerer at config er korrekt
- `ensure_data_directory()` - Sikrer at data_raw/ finnes

### 2. uniprot_fetch.py → uniprot_fetch_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Imports | Relative (`.config`) | Absolute (`scripts.module_2.config_improved`) |
| Feilhåndtering | Ingen try/except | Omfattende try/except blokker |
| Return value | Ingen (bare print) | Returnerer antall sekvenser |
| Parameters | Ingen | `query` og `output_file` som optional params |
| Validering | Ingen | Sjekker at response ikke er tom/None |
| Sekvenstelling | Ingen | Teller og rapporterer antall sekvenser |

**Nye features:**
- Funksjonen kan ta custom query og output file
- Returnerer 0 ved feil, antall sekvenser ved suksess
- Bedre feilmeldinger med context

### 3. cazy_fetch.py → cazy_fetch_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Imports | Relative | Absolute |
| Timeout | Ingen | 30 sekunder (konfigurerbar) |
| Feilhåndtering | `raise_for_status()` bare | Try/except for timeout, HTTP errors, parsing |
| Rate limiting | Ingen | `time.sleep()` mellom requests |
| Duplikathåndtering | I `fetch_all...` bare | Også i `fetch_...for_family` |
| Return ved feil | Exception eller crash | Tom liste `[]` |

**Nye features:**
- `delay_between_requests` parameter (default 1.0 sek)
- `timeout` parameter for requests
- Rapporterer antall duplikater fjernet
- Graceful handling av feil per familie

### 4. ncbi_fetch.py → ncbi_fetch_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Imports | Relative | Absolute |
| Batch fetching | Nei (1 per request) | Ja (100 per request) |
| Email validering | Ingen | Kaller `get_ncbi_email()` med validering |
| Return value | Ingen | (success_count, fail_count) tuple |
| Feilrapportering | Print per ID | Aggregert per batch med totaler |
| Rate limiting | `sleep(0.5)` hvert 10. | `sleep(0.4)` per batch (3 req/sek) |
| Progress tracking | Minimal | Viser batch X/Y med antall |

**Nye features:**
- `BATCH_SIZE` konfigurasjon (100 sekvenser per request)
- Teller vellykkede vs. feilede hentinger
- Bedre NCBI rate limiting compliance
- Email som parameter (ikke bare global config)

### 5. merge_sequences.py → merge_sequences_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Imports | Relative | Absolute |
| Deduplisering | Hardkodet (seq, id) | Konfigurerbar strategi |
| Return value | Ingen | (total, unique) tuple |
| File validation | Sjekker kun exists | Try/except ved parsing også |
| Rapportering | Basic | Detaljert med per-fil og totaler |

**Nye features:**
- `deduplicate_by` parameter: "sequence", "id", "both", "none"
- Teller hvor mange duplikater som fjernes
- Returnerer tuple med før/etter tall
- Rapporterer antall sekvenser per input-fil

### 6. run_module2.py → run_module2_improved.py

| Feature | Original | Improved |
|---------|----------|----------|
| Imports | Relative | Absolute |
| CLI support | Ingen | Full argparse-basert CLI |
| Error handling | Ingen | Try/except per steg |
| Validation | Ingen | Sjekker success mellom steg |
| Exit code | Ingen | 0 for suksess, 1 for feil |
| Skipping steps | Nei | `--skip-uniprot`, `--skip-cazy`, `--skip-merge` |
| Pre-flight checks | Ingen | Dependencies + config validering |

**Nye features:**
- `--email` argument for NCBI email
- `--skip-*` flagg for å hoppe over steg
- `--deduplicate-by` for merge-strategi
- Dependency check før pipeline starter
- Detaljert summary etter kjøring
- Exit code for SLURM/batch jobs

## Hvordan bruke

### Original versjon:
```bash
# MÅ editere config.py først for å sette NCBI_EMAIL
# Kan IKKE kjøre: python scripts/module2/run_module2.py (import errors)
# Krever riktig PYTHONPATH eller annen workaround
```

### Forbedret versjon:
```bash
# Metode 1: Miljøvariabel
export NCBI_EMAIL="din.epost@institusjon.no"
python -m scripts.module_2.run_module2_improved

# Metode 2: CLI argument
python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no

# Metode 3: Med options
python -m scripts.module_2.run_module2_improved \
    --email din.epost@institusjon.no \
    --skip-uniprot \
    --deduplicate-by sequence
```

## Testing

### Test forbedret versjon:
```bash
# Kjør verifikasjonsscript
python -m scripts.module_2.test_improved

# Test config
python -m scripts.module_2.config_improved

# Test individuelle moduler
python -m scripts.module_2.uniprot_fetch_improved
python -m scripts.module_2.cazy_fetch_improved
python -m scripts.module_2.ncbi_fetch_improved
python -m scripts.module_2.merge_sequences_improved
```

## Kritiske forbedringer oppsummert

### ✅ Fikset import-problemet
- Original: Relative imports (`.config`) som krever spesiell setup
- Forbedret: Absolute imports som fungerer med `python -m`

### ✅ NCBI email ikke hardkodet
- Original: Må editere config.py
- Forbedret: Miljøvariabel eller CLI-argument

### ✅ Robust feilhåndtering
- Original: Crashes ved feil
- Forbedret: Try/except overalt, graceful degradation

### ✅ Return values for validering
- Original: Bare print statements
- Forbedret: Returnerer status/count for å sjekke suksess

### ✅ Path-håndtering
- Original: String-paths som kan feile
- Forbedret: pathlib.Path for robusthet

### ✅ Batch fetching (NCBI)
- Original: 1 sekvens per request (tregt!)
- Forbedret: 100 sekvenser per request (100x raskere!)

### ✅ Exit codes
- Original: Ingen
- Forbedret: 0 = success, 1 = feil (viktig for SLURM)

### ✅ Pre-flight checks
- Original: Ingen
- Forbedret: Sjekker dependencies og config før start

## Anbefaling

**Bruk den forbedrede versjonen for:**
- Produksjonskjøringer på SLURM
- Når du trenger robusthet og feilhåndtering
- Når du vil ha fleksibilitet (skip steps, custom params)
- Når du vil ha tydelig status og logging

**Bruk original versjon for:**
- Rask prototyping hvis du allerede har alt satt opp
- Hvis du ikke trenger feilhåndtering
- Hvis du ikke vil endre noe i eksisterende oppsett

## Migreringsguide

1. **Installer pakker** (hvis ikke allerede gjort):
   ```bash
   pip install bioservices biopython beautifulsoup4 requests
   ```

2. **Sett NCBI email**:
   ```bash
   export NCBI_EMAIL="din.epost@institusjon.no"
   # Eller legg til i ~/.bashrc for permanent
   ```

3. **Test at alt fungerer**:
   ```bash
   cd /path/to/project_root
   python -m scripts.module_2.test_improved
   ```

4. **Kjør pipeline**:
   ```bash
   python -m scripts.module_2.run_module2_improved
   ```

5. **For SLURM**, bruk dette template:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=lpmo_fetch
   #SBATCH --time=02:00:00
   #SBATCH --mem=4G
   
   module load Python/3.12.3
   export NCBI_EMAIL="din.epost@institusjon.no"
   
   cd /path/to/project_root
   python -m scripts.module_2.run_module2_improved
   
   # Sjekk exit code
   if [ $? -eq 0 ]; then
       echo "Success"
   else
       echo "Failed"
       exit 1
   fi
   ```
