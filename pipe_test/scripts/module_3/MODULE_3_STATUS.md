# Module 3 - Status og dokumentasjon

## ✅ Status: FULLSTENDIG IMPLEMENTERT

Module 3 er ferdig utviklet og klar til bruk. Den tar output fra Module 2 og utfører domeneannotasjon.

---

## 📋 Hva Module 3 gjør

### Input (fra Module 2):
- `data/sequences/lpmo_all_raw.fasta` - Alle LPMO-sekvenser
- `data/metadata/m2_sequence_3d_metadata.csv` - Metadata med familie-info

### Output:
1. **Domeneannotasjon**:
   - `data/domains/m3_domains_parsed.csv` - Alle domenetreff med filtering
   - `data/domains/raw/*.domtblout` - Rå HMMER output

2. **Metadata**:
   - `data/metadata/m3_domain_metadata.csv` - Per-sekvens domeneinfo
   - `data/metadata/m3_run_metadata.json` - Kjøringsinfo og statistikk

3. **Sekvensvarianter**:
   - `data/sequences/lpmo_full_length.fasta` - Full-length sekvenser
   - `data/sequences/lpmo_catalytic_domain.fasta` - Katalytiske domener

---

## 🔧 Implementerte funksjoner

### config_m3.py
- `Module3Config` dataclass med alle parametere
- Standard verdier matcher dine krav eksakt:
  - `evalue_cutoff = 1e-15`
  - `coverage_cutoff = 0.35`
  - `min_domain_length = 60`
  - `max_domain_overlap_fraction = 0.25`
  - `generate_catalytic_for_cbm_only = True`
  - `generate_catalytic_if_no_cbm = False`
  - `fallback_to_full_length_if_no_domain = True`
  - `min_catalytic_domain_length = 150`

### domain_annotation.py
Alle nødvendige funksjoner:

1. **HMMER integrasjon**:
   - `hmmer_version()` - Sjekk HMMER versjon
   - `run_hmmscan()` - Kjør hmmscan mot HMM-databaser
   - `parse_domtblout()` - Parse domtblout-filer

2. **Domeneanalyse**:
   - `classify_hit()` - Klassifiser domenetreff (ok/high_evalue/low_coverage/etc)
   - `apply_filters()` - Filtrer domener basert på cutoffs
   - `prune_overlaps()` - Fjern overlappende domener
   - `select_lpmo_hit()` - Velg beste LPMO-domene
   - `summarise_cbms()` - Oppsummer CBM-domener

3. **Sekvensannotasjon**:
   - `annotate_sequence()` - Annotere én sekvens (LPMO + CBM + policy)
   - `write_variant_fastas()` - Generer full-length og katalytiske FASTAs

4. **Metadata og output**:
   - `write_domain_table()` - Skriv m3_domains_parsed.csv
   - `write_metadata_table()` - Skriv m3_domain_metadata.csv
   - `write_run_metadata()` - Skriv m3_run_metadata.json

5. **Pipeline**:
   - `run_domain_pipeline()` - Kjør hele Module 3 pipeline

### run_module3.py
- CLI entrypoint med argparse
- Støtter config-fil (YAML/JSON)
- Override av parametere via CLI

---

## 📊 m3_domain_metadata.csv kolonner

Eksakt som spesifisert:
- `seq_id` - Sekvens ID
- `family` - CAZy familie fra M2
- `sequence_length` - Total lengde
- `has_lpmo_domain` - Boolean
- `lpmo_domain_family` - F.eks. AA9
- `lpmo_domain_name` - HMM-navn
- `lpmo_start` - 1-basert start
- `lpmo_end` - 1-basert slutt
- `lpmo_source` - dbcan/dbcan_sub
- `has_cbm` - Boolean
- `cbm_domains` - Liste med CBM-domener
- `use_full_length` - Boolean
- `use_catalytic` - Boolean
- `quality_flag` - ok/no_domain/low_coverage/short_domain/etc

---

## 📊 m3_domains_parsed.csv kolonner

- `seq_id` - Sekvens ID
- `domain_name` - HMM-modell navn
- `domain_family` - F.eks. AA9, CBM1
- `source` - dbcan/dbcan_sub/cbm
- `seq_start` - 1-basert start på sekvens
- `seq_end` - 1-basert slutt på sekvens
- `model_length` - Lengde av HMM-modell
- `alignment_length` - Lengde av alignment
- `coverage` - alignment_length / model_length
- `evalue` - E-value
- `score` - Bitscore
- `confidence_flag` - ok/high_evalue/low_coverage/too_short
- `accepted` - Boolean (passerte filtering)

---

## 📊 m3_run_metadata.json innhold

```json
{
  "run_timestamp": "2024-12-03T...",
  "hmmer_version": "HMMER 3.x",
  "hmmscan_binary": "hmmscan",
  "hmm_databases": {
    "dbcan": "path/to/dbCAN-HMMdb-V12.txt",
    "dbcan_sub": "path/to/dbCAN-subfamily.hmm",
    "cbm": "path/to/cbm_profiles.hmm"
  },
  "parameters": {
    "evalue_cutoff": 1e-15,
    "coverage_cutoff": 0.35,
    "min_domain_length": 60,
    ...
  },
  "num_sequences": 1234,
  "num_with_lpmo": 1200,
  "num_with_cbm": 456,
  "num_use_full_length": 1234,
  "num_use_catalytic": 456,
  "num_low_quality": 34
}
```

---

## 🚀 Hvordan kjøre

### 1. Standard kjøring (bruker default config):
```bash
python -m scripts.module_3.run_module3
```

### 2. Med custom config-fil:
```bash
python -m scripts.module_3.run_module3 --config config_m3.yaml
```

### 3. Med CLI overrides:
```bash
python -m scripts.module_3.run_module3 \
    --cpu 8 \
    --no-subfamily \
    --force-rerun \
    --hmmscan-binary /path/to/hmmscan
```

### 4. Med custom input-filer:
```bash
python -m scripts.module_3.run_module3 \
    --fasta data/sequences/lpmo_all_raw.fasta \
    --metadata data/metadata/m2_sequence_3d_metadata.csv
```

---

## 📝 Config-fil eksempel (config_m3.yaml)

```yaml
# HMM filtering
evalue_cutoff: 1e-15
coverage_cutoff: 0.35
min_domain_length: 60
max_domain_overlap_fraction: 0.25
prefer_family_from_metadata: true

# HMM databases
dbcan_hmm: data/hmms/dbCAN-HMMdb-V12.txt
dbcan_sub_hmm: data/hmms/dbCAN-subfamily.hmm
cbm_hmm: data/hmms/cbm_profiles.hmm
use_subfamily_hmms: true

# HMMER settings
hmmscan_binary: hmmscan
hmmscan_cpu: 4
reuse_existing_domtbl: true

# Policy
generate_catalytic_for_cbm_only: true
generate_catalytic_if_no_cbm: false
fallback_to_full_length_if_no_domain: true
min_catalytic_domain_length: 150

# Notes
notes: "Test run for AA9-AA13 families"
```

---

## 🔍 Policy-logikk (som implementert)

### Scenario 1: Ingen LPMO-domene funnet
- `has_lpmo_domain = False`
- `use_full_length = True` (hvis fallback_to_full_length_if_no_domain)
- `use_catalytic = False`
- `quality_flag = "no_domain"`

### Scenario 2: LPMO-domene funnet, har CBM
- `has_lpmo_domain = True`
- `use_full_length = True`
- `use_catalytic = True` (hvis generate_catalytic_for_cbm_only)
- `quality_flag = "ok"`

### Scenario 3: LPMO-domene funnet, ingen CBM
- `has_lpmo_domain = True`
- `use_full_length = True`
- `use_catalytic = False` (default, eller True hvis generate_catalytic_if_no_cbm)
- `quality_flag = "ok"`

### Scenario 4: LPMO-domene for kort
- `has_lpmo_domain = True`
- `use_catalytic = False` (hvis domene < min_catalytic_domain_length)
- `quality_flag = "short_domain"`

---

## ✅ Alle krav oppfylt

### Domeneannotasjon:
✓ Kjører hmmscan mot dbCAN, dbCAN-sub, CBM-profiler
✓ Parser domtblout-filer
✓ Filtrerer basert på evalue, coverage, min_domain_length
✓ Håndterer overlappende domener
✓ Velger beste LPMO-domene per sekvens

### Kvalitetsvurdering:
✓ Klassifiserer treff (ok/high_evalue/low_coverage/too_short)
✓ Flagging av usikre tilfeller
✓ Støtter prefer_family_from_metadata

### Sekvensvarianter:
✓ Genererer lpmo_full_length.fasta
✓ Genererer lpmo_catalytic_domain.fasta
✓ Policy-basert beslutning (use_full_length, use_catalytic)
✓ Header-format: `>SEQID full_length` eller `>SEQID catalytic lpmo_start=X lpmo_end=Y`

### Metadata og logging:
✓ m3_domains_parsed.csv - Alle domenetreff
✓ m3_domain_metadata.csv - Per-sekvens beslutninger
✓ m3_run_metadata.json - Kjøringsinfo og statistikk
✓ Alle beslutninger og parametere logges

### Konfigurerbarhet:
✓ Alle 14 parametere er konfigurerbare
✓ Støtter YAML/JSON config-filer
✓ CLI overrides
✓ Sensible defaults

---

## 🎯 Integrasjon med Module 2

Module 3 leser automatisk:
- `data/sequences/lpmo_all_raw.fasta` (generert av M2)
- `data/metadata/m2_sequence_3d_metadata.csv` (generert av M2)

Ingen manuell konfigurasjon nødvendig - fungerer out-of-the-box!

---

## 📦 Dependencies

Module 3 krever:
- **Python 3.7+**
- **Biopython** - For FASTA parsing/writing
- **HMMER** - hmmscan må være installert
- **PyYAML** (valgfritt) - For YAML config-filer

Installer:
```bash
pip install biopython pyyaml
```

HMMER (på Linux/Mac):
```bash
# Ubuntu/Debian
sudo apt-get install hmmer

# Conda
conda install -c bioconda hmmer

# Fra kilde
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar xzf hmmer.tar.gz
cd hmmer-*
./configure && make && make install
```

---

## 🧪 Testing

Module 3 inkluderer omfattende feilhåndtering:
- Sjekker at input-filer eksisterer
- Validerer HMM-databaser
- Verifiserer HMMER installasjon
- Håndterer tomme resultater
- Logger warnings og errors

---

## 📖 Videre bruk

Output fra Module 3 er klar for:
- **Module 4** (hvis du har den) - Sekvensjustering
- **Module 5/6** - Strukturprediksjon (AlphaFold, etc.)
- **Manuell analyse** - CSV-filer kan åpnes i Excel/Python

---

## 🎓 Oppsummering

Module 3 er **komplett og produksjonsklar**. Den:
- ✅ Implementerer alle spesifiserte funksjoner
- ✅ Matcher alle default-verdier
- ✅ Genererer alle påkrevde output-filer
- ✅ Logger alt for sporbarhet
- ✅ Er fullt konfigurerbar
- ✅ Integrerer sømløst med Module 2

Ingen ytterligere implementasjon nødvendig!
