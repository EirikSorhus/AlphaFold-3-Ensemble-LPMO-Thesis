# QUICK START GUIDE - Module 2 Forbedret

## TL;DR - Kom i gang på 2 minutter

### 1. Installer pakker (hvis nødvendig)
```bash
pip install bioservices biopython beautifulsoup4 requests
```

### 2. Sett NCBI email
```bash
export NCBI_EMAIL="din.epost@institusjon.no"
```

### 3. Kjør pipeline
```bash
cd /path/to/project_root
python -m scripts.module_2.run_module2_improved
```

---

## Verifisering - Sjekkliste før første kjøring

Kjør disse kommandoene fra **project_root** (`c:\Users\eirik\Masteroppgave\VS_code_test`):

### ✓ Sjekk 1: Python-versjon
```bash
python --version
# Forventet: Python 3.12.3 (eller 3.8+)
```

### ✓ Sjekk 2: Filstruktur er korrekt
```bash
# Windows PowerShell:
Test-Path scripts/__init__.py
Test-Path scripts/module_2/__init__.py
Test-Path scripts/module_2/config_improved.py
Test-Path scripts/module_2/run_module2_improved.py

# Alle skal returnere: True
```

### ✓ Sjekk 3: Dependencies
```bash
python -m scripts.module_2.config_improved

# Forventet output:
# [CONFIG] All required packages are installed
# [CONFIG] Configuration is valid
# ✓ Configuration is ready
```

### ✓ Sjekk 4: NCBI email er satt
```bash
# Windows PowerShell:
echo $env:NCBI_EMAIL

# Skal vise din e-post, ikke "din.epost@institusjon.no"
```

### ✓ Sjekk 5: Test imports
```bash
python -m scripts.module_2.test_improved

# Forventet output:
# ✓ ALL TESTS PASSED
```

---

## Hvis noe feiler

### Problem: "ModuleNotFoundError: No module named 'scripts'"
**Løsning:** Du kjører fra feil mappe. Gå til project_root:
```bash
cd c:\Users\eirik\Masteroppgave\VS_code_test
```

### Problem: "Import 'bioservices' could not be resolved"
**Løsning:** Installer pakker:
```bash
pip install bioservices biopython beautifulsoup4 requests
```

### Problem: "ValueError: NCBI email not configured"
**Løsning:** Sett miljøvariabel:
```bash
# Windows PowerShell:
$env:NCBI_EMAIL = "din.epost@institusjon.no"

# Eller legg til i kommando:
python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no
```

### Problem: Test-scriptet feiler
**Løsning:** Sjekk spesifikk feilmelding og se README_IMPROVED.md

---

## SLURM Job Template

Kopier denne malen til en fil (f.eks. `run_lpmo_fetch.sh`):

```bash
#!/bin/bash
#SBATCH --job-name=lpmo_fetch
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=lpmo_fetch_%j.log
#SBATCH --error=lpmo_fetch_%j.err

# Last Python-miljø
module load Python/3.12.3  # Juster til din versjon

# Sett NCBI email (VIKTIG!)
export NCBI_EMAIL="din.epost@institusjon.no"

# Gå til project root
cd /cluster/home/dittbrukernavn/Masteroppgave/VS_code_test

# Kjør pipeline
echo "Starting LPMO fetch pipeline..."
python -m scripts.module_2.run_module2_improved

# Sjekk exit code
if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully!"
    exit 0
else
    echo "Pipeline failed. Check logs for details."
    exit 1
fi
```

Submit med:
```bash
sbatch run_lpmo_fetch.sh
```

---

## Forventet Output

### Vellykket kjøring ser slik ut:

```
======================================================================
MODULE 2: LPMO SEQUENCE FETCHING PIPELINE
======================================================================

Pre-flight checks:
----------------------------------------------------------------------
[CONFIG] All required packages are installed
[CONFIG] Configuration is valid
[CONFIG] Data directory ready: /path/to/data_raw

STEP 1: UniProt Sequence Fetching
----------------------------------------------------------------------
[UniProt] Starting fetch with query: (family:"AA9" OR ...)
[UniProt] SUCCESS: Saved 1234 sequences to /path/to/uniprot_LPMO_raw.fasta
[PIPELINE] ✓ UniProt: 1234 sequences fetched

STEP 2: CAZy + NCBI Sequence Fetching
----------------------------------------------------------------------
[CAZy] Fetching GenBank IDs from 4 families: ['AA9', 'AA10', 'AA11', 'AA13']
[CAZy] AA9: Found 567 GenBank IDs
[CAZy] AA10: Found 234 GenBank IDs
[CAZy] AA11: Found 123 GenBank IDs
[CAZy] AA13: Found 89 GenBank IDs
[CAZy] Total unique GenBank IDs: 987
[NCBI] Batch 1/10: Fetching 100 sequences...
[NCBI] Batch 2/10: Fetching 100 sequences...
...
[NCBI] SUCCESS: Saved 987 sequences to /path/to/cazy_LPMO_raw.fasta
[PIPELINE] ✓ CAZy/NCBI: 987 sequences fetched

STEP 3: Merge and Deduplicate Sequences
----------------------------------------------------------------------
[MERGE] Merging files: UniProt, CAZy
[MERGE] Read 1234 sequences from uniprot_LPMO_raw.fasta
[MERGE] Read 987 sequences from cazy_LPMO_raw.fasta
[MERGE] Total sequences before deduplication: 2221
[MERGE] Sequences after deduplication: 1850
[MERGE] Duplicates removed: 371
[MERGE] SUCCESS: Wrote merged FASTA to /path/to/LPMO_all_raw.fasta
[PIPELINE] ✓ Merge: 2221 → 1850 unique sequences

======================================================================
PIPELINE SUMMARY
======================================================================
✓ UniProt file: /path/to/uniprot_LPMO_raw.fasta
✓ CAZy file: /path/to/cazy_LPMO_raw.fasta
✓ Merged file: /path/to/LPMO_all_raw.fasta

✓ PIPELINE COMPLETED SUCCESSFULLY
```

### Output-filer (i `data_raw/`):
- `uniprot_LPMO_raw.fasta` - Sekvenser fra UniProt
- `cazy_LPMO_raw.fasta` - Sekvenser fra CAZy/NCBI
- `LPMO_all_raw.fasta` - **Sammenslått og deduplisert (SLUTTPRODUKT)**

---

## Hurtigreferanse - Nyttige kommandoer

### Se alle alternativer:
```bash
python -m scripts.module_2.run_module2_improved --help
```

### Kjør kun én del av pipeline:
```bash
# Kun UniProt:
python -m scripts.module_2.uniprot_fetch_improved

# Kun CAZy/NCBI:
export NCBI_EMAIL="your@email.com"
python -m scripts.module_2.ncbi_fetch_improved

# Kun merge (hvis du allerede har de to andre filene):
python -m scripts.module_2.merge_sequences_improved
```

### Hopp over steg:
```bash
# Hopp over UniProt (hvis du allerede har filen):
python -m scripts.module_2.run_module2_improved --skip-uniprot --email your@email.com

# Hopp over CAZy (hvis du bare vil ha UniProt-data):
python -m scripts.module_2.run_module2_improved --skip-cazy
```

### Bruk annen dedupliseringsstrategi:
```bash
# Dedupliser kun på sekvens (standard):
python -m scripts.module_2.run_module2_improved --deduplicate-by sequence

# Dedupliser kun på ID:
python -m scripts.module_2.run_module2_improved --deduplicate-by id

# Dedupliser på både ID og sekvens:
python -m scripts.module_2.run_module2_improved --deduplicate-by both

# Ingen deduplisering:
python -m scripts.module_2.run_module2_improved --deduplicate-by none
```

---

## Filreferanse

| Fil | Formål |
|-----|--------|
| `config_improved.py` | Konfigurasjon og konstanter |
| `uniprot_fetch_improved.py` | Henter fra UniProt API |
| `cazy_fetch_improved.py` | Scraper CAZy-websider |
| `ncbi_fetch_improved.py` | Henter fra NCBI Protein |
| `merge_sequences_improved.py` | Slår sammen FASTA-filer |
| `run_module2_improved.py` | **Hovedscript - kjør dette** |
| `test_improved.py` | Verifikasjonsscript |
| `README_IMPROVED.md` | Detaljert dokumentasjon |
| `COMPARISON.md` | Sammenligning original vs improved |
| `QUICK_START.md` | Denne filen |

---

## Hjelp og support

1. **Les output-meldingene nøye** - de er designet for å være informative
2. Sjekk `README_IMPROVED.md` for detaljert dokumentasjon
3. Sjekk `COMPARISON.md` for forskjeller fra original versjon
4. Kjør `test_improved.py` for å diagnostisere problemer

---

**VIKTIG:** Alle de originale filene (`config.py`, `uniprot_fetch.py`, etc.) er bevart uendret. De nye forbedrede filene har suffikset `_improved` i navnet.
