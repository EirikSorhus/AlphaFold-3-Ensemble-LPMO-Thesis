# KONFIGURERINGSGUIDE - Module 2 Forbedret

## Hva er hardkodet og hva kan endres?

### ✅ KAN ENDRES VIA KOMMANDOLINJE

Følgende kan styres når du kjører pipelinen uten å endre kode:

#### 1. **Proteinfamilier** (NYTT!)
```bash
# Standard (AA9, AA10, AA11, AA13):
python -m scripts.module_2.run_module2_improved

# Velg én eller flere familier:
python -m scripts.module_2.run_module2_improved --families AA9 AA10

# Enkelt familie:
python -m scripts.module_2.run_module2_improved --families AA13
```

**Hvorfor:** Lar deg teste forskjellige familier uten å editere kode.  
**MÅ velges:** NEI (bruker standard hvis utelatt)

**Output:** Pipeline lager separate FASTA-filer per familie:
- `AA9_20241202_143022_uniprot.fasta`
- `AA9_20241202_143022_cazy.fasta`
- `AA10_20241202_143022_uniprot.fasta`
- ... osv
- `20241202_143022_merged.fasta` (alle familier samlet)
- `20241202_143022_metadata.json` (info om kjøringen)

---

#### 2. **NCBI Email** (PÅKREVD)
```bash
# Via miljøvariabel:
export NCBI_EMAIL="din.epost@institusjon.no"
python -m scripts.module_2.run_module2_improved

# Via CLI-argument:
python -m scripts.module_2.run_module2_improved --email din.epost@institusjon.no
```

**Hvorfor:** NCBI krever en gyldig e-postadresse for API-tilgang.  
**MÅ velges:** JA (pipeline vil feile uten)

---

#### 3. **Hvilke steg som kjøres**
```bash
# Hopp over UniProt:
python -m scripts.module_2.run_module2_improved --skip-uniprot

# Hopp over CAZy/NCBI:
python -m scripts.module_2.run_module2_improved --skip-cazy

# Hopp over merge:
python -m scripts.module_2.run_module2_improved --skip-merge

# Kombiner (kun kjør merge av eksisterende filer):
python -m scripts.module_2.run_module2_improved --skip-uniprot --skip-cazy
```

**Hvorfor:** Hvis du allerede har noen filer eller bare vil teste deler av pipelinen.  
**MÅ velges:** NEI (standard kjører alle steg)

---

#### 4. **Dedupliseringsstrategi**
```bash
# Kun sekvens (STANDARD - fjerner identiske sekvenser):
python -m scripts.module_2.run_module2_improved --deduplicate-by sequence

# Kun ID (fjerner entries med samme ID):
python -m scripts.module_2.run_module2_improved --deduplicate-by id

# Både ID og sekvens må matche:
python -m scripts.module_2.run_module2_improved --deduplicate-by both

# Ingen deduplisering:
python -m scripts.module_2.run_module2_improved --deduplicate-by none
```

**Hvorfor:** Forskjellige strategier gir forskjellig antall sekvenser i sluttfilen.  
**MÅ velges:** NEI (standard er `sequence`)

**Anbefaling:** Bruk `sequence` (standard) for LPMO-er siden du vil ha unike sekvenser uavhengig av hvor de kommer fra.

---

### ✅ OPPDATERT: Proteinfamilier er nå valgbare!

#### 5. **Standard proteinfamilier (kan overstyres via CLI)**
```bash
# Standard (bruker alle LPMO-familier: AA9, AA10, AA11, AA13):
python -m scripts.module_2.run_module2_improved --email din@email.com

# Velg spesifikke familier:
python -m scripts.module_2.run_module2_improved --families AA9 AA10 --email din@email.com

# Kun én familie:
python -m scripts.module_2.run_module2_improved --families AA13 --email din@email.com

# Andre CAZy-familier (f.eks. GH-familier):
python -m scripts.module_2.run_module2_improved --families GH5 GH6 GH7 --email din@email.com
```

**Standard familier:** AA9, AA10, AA11, AA13 (brukes hvis --families ikke spesifiseres)

**MÅ velges:** NEI (bruker standard hvis ikke spesifisert)

---

### ⚠️ HARDKODET (må endre fil for å endre)

---

#### 6. **Output-mappe (HARDKODET)**
**Fil:** `scripts/module_2/config_improved.py`  
**Linjer:** 59-64

```python
DATA_RAW_DIR = PROJECT_ROOT / "data_raw"

UNIPROT_FASTA = DATA_RAW_DIR / "uniprot_LPMO_raw.fasta"
CAZY_FASTA = DATA_RAW_DIR / "cazy_LPMO_raw.fasta"
MERGED_FASTA = DATA_RAW_DIR / "LPMO_all_raw.fasta"
```

**For å endre:**
Endre `DATA_RAW_DIR` i config_improved.py

**Hvorfor hardkodet:** Standard plassering for rådata.

**Output-filer (genereres automatisk per kjøring):**
- Per familie: `{familie}_{run_id}_{kilde}.fasta` (f.eks. `AA9_20241202_143022_uniprot.fasta`)
- Merged: `{run_id}_merged.fasta` (f.eks. `20241202_143022_merged.fasta`)
- Metadata: `{run_id}_metadata.json` (f.eks. `20241202_143022_metadata.json`)

---

#### 7. **API-parametere (HARDKODET)**
**Fil:** `scripts/module_2/config_improved.py`  
**Linjer:** 113-122

```python
# HTTP request timeout (seconds)
REQUEST_TIMEOUT = 30

# NCBI rate limiting (requests per second)
NCBI_REQUESTS_PER_SECOND = 3
NCBI_BATCH_SIZE = 100  # Number of sequences to fetch per request
```

**For å endre:**
Rediger verdiene i config_improved.py

**Hvorfor hardkodet:** Følger best practices for API-bruk og NCBI retningslinjer.

---

## 📋 Oppsummering: Hva kan endres hvor?

| Parameter | Hvor endres | Hvordan | Må endre fil? |
|-----------|-------------|---------|---------------|
| **Proteinfamilier** | CLI | `--families AA9 AA10` | NEI |
| **NCBI Email** | CLI / miljøvariabel | `--email` eller `export NCBI_EMAIL=` | NEI |
| **Hoppe over kilder** | CLI | `--skip-uniprot` / `--skip-ncbi` | NEI |
| **Deduplisering** | CLI | `--deduplicate-by sequence` | NEI |
| **Output-mappe** | `config_improved.py` | Endre `DATA_RAW_DIR` | JA |
| **API-parametere** | `config_improved.py` | Endre `REQUEST_TIMEOUT`, etc. | JA |

---

## 🎯 Vanlige scenarier

### Scenario 1: Kun hente AA9 og AA10 familier
```bash
python -m scripts.module_2.run_module2_improved \
    --families AA9 AA10 \
    --email din.epost@institusjon.no
```

### Scenario 2: Teste kun én familie (f.eks. AA13)
```bash
python -m scripts.module_2.run_module2_improved \
    --families AA13 \
    --email din.epost@institusjon.no
```

### Scenario 3: Hoppe over UniProt (bruk bare CAZy/NCBI)
```bash
python -m scripts.module_2.run_module2_improved \
    --skip-uniprot \
    --email din.epost@institusjon.no
```

### Scenario 4: Hente kun fra UniProt (ingen CAZy/NCBI)
```bash
python -m scripts.module_2.run_module2_improved --skip-ncbi
```

### Scenario 5: Kun merge eksisterende filer
```bash
python -m scripts.module_2.run_module2_improved \
    --skip-uniprot \
    --skip-ncbi
```

---

## Se alle alternativer
```bash
python -m scripts.module_2.run_module2_improved --help
```

Dette viser komplett liste over alle CLI-argumenter.

---

## Spørsmål eller ønsker?

Hvis du vil ha:
1. **CLI-støtte for å velge proteinfamilier** - Jeg kan legge til `--families` argument
2. **CLI-støtte for å velge output-filnavn** - Jeg kan legge til `--output` argument
3. **Andre konfigurerbare parametere** - Bare spør!

La meg vite om du vil ha noen av disse tilleggene.
