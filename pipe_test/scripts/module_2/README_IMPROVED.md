# Module 2: LPMO Sequence Fetching Pipeline - Improved Version

## Oversikt

Dette er den forbedrede versjonen av Module 2 som henter LPMO (Lytic Polysaccharide Monooxygenase) proteinsekvenser fra tre databaser:
- **UniProt**: Søk etter proteiner annotert med CAZy-familier
- **CAZy**: Web-scraping for å hente GenBank-IDer
- **NCBI**: Henting av proteinsekvenser fra GenBank-IDer

## Forbedringer i denne versjonen

### Kritiske forbedringer
1. ✅ **Fikset import-system**: Bruker absolute imports som fungerer med `python -m`
2. ✅ **pathlib for filhåndtering**: Alle stier er nå Path-objekter for robusthet
3. ✅ **Feilhåndtering**: Omfattende try/except blokker med informative feilmeldinger
4. ✅ **NCBI email som parameter**: Ikke lenger hardkodet, kan settes via miljøvariabel eller CLI
5. ✅ **Return values**: Alle funksjoner returnerer status/antall for validering
6. ✅ **Batch-henting fra NCBI**: Raskere og mer effektivt (100 sekvenser per request)
7. ✅ **Dependency check**: Sjekker at alle pakker er installert før kjøring

### Ekstra forbedringer
8. ✅ **Command-line interface**: Kan kjøres med argumenter (`--email`, `--skip-uniprot`, etc.)
9. ✅ **Type hints**: Alle funksjoner har type annotations
10. ✅ **Docstrings**: Detaljert dokumentasjon av alle funksjoner
11. ✅ **Bedre logging**: Tydelig output med status-indikatorer (✓, ✗, ⚠)
12. ✅ **Fleksibel deduplisering**: Kan velge strategi (sequence, id, both, none)
13. ✅ **Rate limiting**: Respekterer NCBI sine grenser (3 req/sek)
14. ✅ **Exit codes**: Returnerer 0 ved suksess, 1 ved feil (viktig for SLURM)

## Installasjon

### 1. Pakkestruktur
Alle nødvendige `__init__.py` filer er nå på plass:
```
scripts/
├── __init__.py
└── module_2/
    ├── __init__.py
    ├── config_improved.py
    ├── uniprot_fetch_improved.py
    ├── cazy_fetch_improved.py
    ├── ncbi_fetch_improved.py
    ├── merge_sequences_improved.py
    └── run_module2_improved.py
```

### 2. Python-pakker
Installer nødvendige pakker:
```bash
pip install bioservices biopython beautifulsoup4 requests
```

Eller med conda:
```bash
conda install -c bioconda bioservices biopython beautifulsoup4 requests
```

### 3. NCBI Email-konfigurasjon
NCBI krever en gyldig e-postadresse. Velg én av disse metodene:

**Metode 1: Miljøvariabel (anbefalt for SLURM)**
```bash
export NCBI_EMAIL="din.epost@institusjon.no"
```

**Metode 2: Command-line argument**
```bash
python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no
```

## Bruk

### Grunnleggende bruk
Fra project_root:
```bash
# Sett email først
export NCBI_EMAIL="din.epost@institusjon.no"

# Kjør full pipeline
python -m scripts.module_2.run_module2_improved
```

### Med command-line argumenter
```bash
# Spesifiser email direkte
python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no

# Hopp over UniProt (hvis du allerede har filen)
python -m scripts.module_2.run_module2_improved --skip-uniprot --email din.epost@institusjon.no

# Bruk annen dedupliseringsstrategi
python -m scripts.module_2.run_module2_improved --deduplicate-by both --email din.epost@institusjon.no
```

### Se alle alternativer
```bash
python -m scripts.module_2.run_module2_improved --help
```

### For SLURM job
I din SLURM-script:
```bash
#!/bin/bash
#SBATCH --job-name=lpmo_fetch
#SBATCH --time=02:00:00
#SBATCH --mem=4G

# Last Python-miljø
module load Python/3.12.3  # eller din versjon

# Sett NCBI email
export NCBI_EMAIL="din.epost@institusjon.no"

# Kjør pipeline
cd /path/to/project_root
python -m scripts.module_2.run_module2_improved

# Sjekk exit code
if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully"
else
    echo "Pipeline failed"
    exit 1
fi
```

## Testing

### Test konfigurasjon
```bash
python -m scripts.module_2.config_improved
```
Dette sjekker:
- At alle pakker er installert
- At konfigurasjonen er gyldig
- At stier er riktige

### Test individuelle moduler

**Test UniProt:**
```bash
python -m scripts.module_2.uniprot_fetch_improved
```

**Test CAZy scraping:**
```bash
python -c "from scripts.module_2.cazy_fetch_improved import fetch_all_cazy_genbank_ids; print(f'Found {len(fetch_all_cazy_genbank_ids())} IDs')"
```

**Test NCBI fetching:**
```bash
export NCBI_EMAIL="din.epost@institusjon.no"
python -m scripts.module_2.ncbi_fetch_improved
```

**Test merge:**
```bash
python -m scripts.module_2.merge_sequences_improved
```

## Output

Pipeline lager følgende filer i `data_raw/`:
- `uniprot_LPMO_raw.fasta` - Sekvenser fra UniProt
- `cazy_LPMO_raw.fasta` - Sekvenser fra CAZy/NCBI
- `LPMO_all_raw.fasta` - Sammenslått og deduplisert fil (**sluttprodukt**)

## Feilsøking

### Problem: Import errors
**Feil:** `ModuleNotFoundError: No module named 'scripts'`

**Løsning:** Kjør alltid fra project_root med `python -m scripts.module_2.run_module2_improved`

### Problem: NCBI email ikke satt
**Feil:** `ValueError: NCBI email not configured`

**Løsning:** 
```bash
export NCBI_EMAIL="din.epost@institusjon.no"
```
eller bruk `--email` argument

### Problem: Missing dependencies
**Feil:** Import errors for bioservices, Bio, bs4, eller requests

**Løsning:**
```bash
pip install bioservices biopython beautifulsoup4 requests
```

### Problem: Tom output
Hvis pipeline fullføres men filene er tomme eller har 0 sekvenser:
1. Sjekk nettverkstilkobling
2. Sjekk at UniProt/CAZy/NCBI er tilgjengelige
3. Sjekk query-parametere i `config_improved.py`
4. Se etter feilmeldinger i output

## Deduplication Strategies

- `sequence` (standard): Fjerner sekvenser med identisk aminosyresekvens (anbefalt)
- `id`: Fjerner entries med samme ID
- `both`: Fjerner kun entries der både ID og sekvens er identiske
- `none`: Ingen deduplisering

## Viktige forskjeller fra original versjon

| Feature | Original | Improved |
|---------|----------|----------|
| Import-system | Relative imports (`.config`) | Absolute imports (`scripts.module_2.config_improved`) |
| Filstier | String-based | pathlib.Path |
| NCBI email | Hardkodet i config | Parameter/miljøvariabel |
| Feilhåndtering | Minimal | Omfattende med try/except |
| Return values | Ofte None | Status/count returverdier |
| CLI support | Nei | Ja, med argparse |
| Batch fetching | Nei (1 per request) | Ja (100 per request) |
| Dependency check | Nei | Ja, før pipeline starter |

## Kontaktinformasjon

Ved spørsmål eller problemer, sjekk output-meldingene nøye - de inneholder vanligvis informasjon om hva som gikk galt.
