# Module 2 - Oppdateringer til nye krav

## Endringer gjennomført (2024-12-02)

### 1. Output-struktur endret

**Før:**
```
data_raw/
├── AA9_20241202_143022_uniprot.fasta
├── AA9_20241202_143022_ncbi.fasta
├── 20241202_143022_merged.fasta
└── 20241202_143022_metadata.json
```

**Nå:**
```
data/
├── sequences/
│   ├── lpmo_AA9_raw.fasta
│   ├── lpmo_AA10_raw.fasta
│   ├── lpmo_AA11_raw.fasta
│   ├── lpmo_AA13_raw.fasta
│   └── lpmo_all_raw.fasta
└── metadata/
    ├── m2_run_metadata.json
    └── m2_sequence_3d_metadata.csv
```

### 2. Filnavnkonvensjoner

- **Per-familie FASTA**: `lpmo_{FAMILY}_raw.fasta` (stabile navn uten timestamp)
- **Samlet FASTA**: `lpmo_all_raw.fasta`
- **Kjørings-metadata**: `m2_run_metadata.json`
- **Sekvens-metadata**: `m2_sequence_3d_metadata.csv`

### 3. Familie-navn i FASTA headers

Alle FASTA-headers inneholder nå familienavn:
```
>Q2FJX8 AA9
>P12345 AA10
```

### 4. Ny CSV-metadata per sekvens

Filen `m2_sequence_3d_metadata.csv` inneholder:
- `protein_id` - UniProt/NCBI ID
- `family` - CAZy familie (AA9, AA10, etc.)
- `sequence_length` - Lengde på sekvens
- `has_experimental_structure` - Boolean for PDB-struktur
- `pdb_ids` - Kommaseparert liste av PDB-IDer
- `best_pdb_resolution` - Beste oppløsning (Å)
- `has_alphafold_model` - Boolean for AlphaFold-modell
- `alphafold_accession` - AlphaFold-ID/confidence

### 5. Oppdatert kjørings-metadata

`m2_run_metadata.json` inneholder nå:
- `run_id` - Tidsstempel for kjøring
- `timestamp` - ISO-format timestamp
- `families` - Liste over familier prosessert
- `n_sequences_per_family` - Dict med antall per familie
- `source_databases` - Info om datakilder brukt
- `uniprot_results` - Resultater fra UniProt-henting
- `cazy_results` - Resultater fra CAZy/NCBI-henting
- `merge_total` / `merge_unique` - Merge-statistikk

### 6. Validerings-sjekker (STEP 6)

Pipeline kjører nå automatisk følgende valideringer:

✓ **Check 1**: Alle IDer i metadata CSV finnes i `lpmo_all_raw.fasta`
✓ **Check 2**: Sum av sekvenser i per-familie filer = antall i merged fil
✓ **Check 3**: `n_sequences_per_family` i metadata matcher faktiske FASTA-filer
✓ **Check 4**: Alle familier i metadata er forventede familier

### 7. Pipeline-steg oppdatert

**STEP 1**: Hent UniProt-sekvenser per familie → `lpmo_{FAMILY}_raw.fasta`
**STEP 2**: Hent CAZy/NCBI-sekvenser per familie → oppdaterer `lpmo_{FAMILY}_raw.fasta`
**STEP 3**: Merge alle familie-filer → `lpmo_all_raw.fasta`
**STEP 4**: Lagre kjørings-metadata → `m2_run_metadata.json`
**STEP 5**: Annoter familie-navn og 3D-struktur → oppdaterer FASTA + lager CSV
**STEP 6**: Kjør validerings-sjekker

### 8. Invarianter verifisert

Pipeline sikrer at:
- Alle IDer i `m2_sequence_3d_metadata.csv` finnes i `lpmo_all_raw.fasta`
- Antall sekvenser i `lpmo_all_raw.fasta` = sum av per-familie FASTA
- `n_sequences_per_family` stemmer med faktisk innhold
- `family` i metadata stemmer med hvilken FASTA fil sekvensen kom fra

## Oppdaterte konfigurasjons-filer

### config_improved.py
- Ny: `DATA_SEQUENCES_DIR = PROJECT_ROOT / "data" / "sequences"`
- Ny: `DATA_METADATA_DIR = PROJECT_ROOT / "data" / "metadata"`
- Endret: `get_family_output_filename()` - returnerer `lpmo_{FAMILY}_raw.fasta`
- Endret: `get_merged_output_filename()` - returnerer `lpmo_all_raw.fasta`
- Endret: `get_metadata_filename()` - returnerer `m2_run_metadata.json`
- Endret: `ensure_data_directory()` - lager begge nye directories

### uniprot_fetch_improved.py
- Bruker nå `DATA_SEQUENCES_DIR` som standard output
- Oppdatert til nye filnavnkonvensjoner

### ncbi_fetch_improved.py
- Bruker nå `DATA_SEQUENCES_DIR` som standard output
- Oppdatert til nye filnavnkonvensjoner

### merge_sequences_improved.py
- Bruker nå `DATA_SEQUENCES_DIR` som standard output
- Oppdatert til nye filnavnkonvensjoner

### run_module2_improved.py
- STEP 5: Legger til familie-navn i headers
- STEP 5: Henter 3D-strukturinfo fra UniProt
- STEP 5: Genererer `m2_sequence_3d_metadata.csv`
- STEP 6: Kjører 4 validerings-sjekker
- Oppdatert summary til å liste riktige filer

## Hvordan kjøre

```bash
# Standard kjøring (alle LPMO-familier)
python -m scripts.module_2.run_module2_improved --email your@email.com

# Velg spesifikke familier
python -m scripts.module_2.run_module2_improved --families AA9 AA10 --email your@email.com

# Med miljøvariabel
export NCBI_EMAIL="your@email.com"
python -m scripts.module_2.run_module2_improved --families AA13
```

## Forventet output

### data/sequences/
- `lpmo_AA9_raw.fasta` - Alle AA9 sekvenser med header `>ID AA9`
- `lpmo_AA10_raw.fasta` - Alle AA10 sekvenser
- `lpmo_AA11_raw.fasta` - Alle AA11 sekvenser
- `lpmo_AA13_raw.fasta` - Alle AA13 sekvenser
- `lpmo_all_raw.fasta` - Alle sekvenser samlet

### data/metadata/
- `m2_run_metadata.json` - Info om kjøring, familier, antall sekvenser
- `m2_sequence_3d_metadata.csv` - Per-sekvens info med 3D-struktur

## Testing

Valideringen kjører automatisk i STEP 6 og rapporterer:
- ✓ Alle sjekker passert
- ✗ Feil funnet med spesifikk beskrivelse

## Kompatibilitet

- Gamle `data_raw/` directory beholdes for bakoverkompatibilitet
- Gamle filnavnfunksjoner aksepterer fortsatt `run_id` parameter (ignoreres)
- Importerer fortsatt `DATA_RAW_DIR` (legacy support)
