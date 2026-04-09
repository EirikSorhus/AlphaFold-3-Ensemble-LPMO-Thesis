---
applyTo: "Masteroppgave/structure_pipeline/**"
---

# structure_pipeline – Repo-instruksjoner for Copilot og utviklere

> Ingen endringer, slettinger eller laging av nye filer skal skje ved kommandolinje operasjoner.

## Formål

`structure_pipeline` er en produksjonsklar, manifest-drevet arbeidsflyt for storskala
prediksjon av proteinstruktur på HPC-klynger med SLURM-scheduling. Pipelinen støtter tre
prediksjonsmodeller: **AlphaFold3 (AF3)**, **Boltz-2** og **RoseTTAFold3 (RF3)**. Kjøring
skjer gjennom Apptainer/Singularity-konteinere for reproduserbarhet.

Pipelinen anses som **stabil**. Endringer skal kun gjøres ved feil eller konfigurasjonsbehov –
se "Do/Don't for endringer" nedenfor.

---

## Nøkkelstier og struktur

```
structure_pipeline/
├── src/structure_pipeline/        # Kildekode
│   ├── cli.py                     # CLI-entrypoint (typer-app)
│   ├── config.py                  # Pydantic-konfig (PipelineConfig, SlurmConfig)
│   ├── cases.py                   # Kasus-generering (Cartesian product)
│   ├── resume.py                  # Resume/statussporing (DONE.ok + status.jsonl)
│   ├── manifest/                  # CSV-manifest-parsing
│   ├── runners/                   # Modellspesifikk kjørelogikk
│   │   ├── af3.py                 # AF3Runner (tostegsfan-in)
│   │   ├── boltz.py               # BoltzRunner
│   │   ├── rf3.py                 # RF3Runner
│   │   ├── base.py                # Abstrakt BaseRunner
│   │   ├── ligand_utils.py        # Ligand-hjelper
│   │   ├── mounts.py              # Bind-mount-optimering
│   │   └── oligo.py               # Oligosakkarid-håndtering
│   └── executors/                 # Abstrakt jobbutlegging (SLURM / local)
├── bin/
│   ├── structure-pipeline         # Shell-wrapper (bruk denne i stedet for pip-entry)
│   └── run_pipeline.sh            # Eksempelkjøreskript
├── config/
│   ├── pipeline.yaml              # Aktiv konfigurasjonsfil (redigerbar)
│   └── oligo_definitions.yaml     # Oligosakkarid-definisjoner for AF3
├── input/
│   ├── fasta/                     # FASTA-filer med proteinsekvenser
│   ├── ligand/                    # CCD-format CIF-filer organisert per kategori
│   ├── msa/                       # Forhåndsberegnede MSA-filer (.a3m, for Boltz/RF3)
│   └── template/                  # (valgfritt) malstrukturer
├── manifest/                      # Genererte CSV-filer (cases, proteins, ligands, msa)
├── work/                          # Kjøretidsdata per (ligand, modell)
├── results/                       # Prediksjonsresultater
├── logs/                          # SLURM-jobberlogger
├── tests/                         # Pytest-tester
├── pyproject.toml                 # Pakkemetadata og avhengigheter
└── README.md                      # Fullstendig dokumentasjon (autoritativ kilde)
```

---

## Miljø og avhengigheter

### Python-pakken

- Python ≥ 3.10
- Avhengigheter (fra `pyproject.toml`): `typer[all]>=0.9.0`, `pydantic>=2.0`, `pyyaml>=6.0`, `rich>=13.0`
- Installert i `structure_pipeline_env` på Olivia (se `conda/structure_pipeline_env/`)

### PATH-oppsett (Olivia HPC)

```bash
export PATH="/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env/bin:$PATH"
export PATH="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/bin:$PATH"
```

> **OBS:** Bruk `bin/structure-pipeline`-wrapperen (eller `python -m structure_pipeline.cli`)
> fremfor pip-installert entrypoint for å unngå path-problemer inne i konteinere.

### HPC-plattform

- Klynge: **Olivia / Betzy** med SLURM
- Konteinerisering: **Apptainer/Singularity**
- SLURM-konto: `nn1003k`
- Partisjoner: `normal` (CPU), `accel` (GPU)
- Modellvekter og konteinerimager: under `/cluster/projects/nn1003k/prog/` (se `config/pipeline.yaml` for faktiske stier)

---

## Kjøring og bruk

### Typisk arbeidsflyt (i rekkefølge)

```bash
# 1. Initialiser konfig (bare én gang per prosjekt)
structure-pipeline init-config --output config/pipeline.yaml
# → Rediger config/pipeline.yaml med faktiske stier

# 2. Valider konfig og tilgjengelighet
structure-pipeline validate -c config/pipeline.yaml

# 3. Generer manifester (parsing av input)
structure-pipeline manifest -c config/pipeline.yaml

# 4. Forhåndsvis jobber (uten å sende inn)
structure-pipeline run -c config/pipeline.yaml --dry-run

# 5. Kjør alle ventende jobber
structure-pipeline run -c config/pipeline.yaml

# 6. Sjekk fremdrift
structure-pipeline status -c config/pipeline.yaml
```

### Filtrering og overrides (vanlige eksempler)

```bash
# Kun én modell
structure-pipeline run -c config/pipeline.yaml --model af3
structure-pipeline run -c config/pipeline.yaml --model boltz
structure-pipeline run -c config/pipeline.yaml --model rf3

# Kun ett protein
structure-pipeline run -c config/pipeline.yaml --protein B6EQJ6

# Override modellparametre
structure-pipeline run -c config/pipeline.yaml \
  --af3-seeds 5 --af3-diffusion-samples 3

structure-pipeline run -c config/pipeline.yaml \
  --boltz-recycling-steps 15

# Tvingen kjøring (ignorer DONE.ok)
structure-pipeline run -c config/pipeline.yaml --no-resume

# Lokal kjøring uten SLURM (for testing)
structure-pipeline run -c config/pipeline.yaml --local --model rf3
```

### AF3 tostegsfan-in (automatisk)

AF3 sender inn MSA-jobber (CPU) og inferensjobber (GPU) separat. Pipelinen håndterer
`--dependency=afterok`-kjeden automatisk. Ingen manuell interaksjon nødvendig.

---

## Input/Output og lagring

### Input

| Type | Plassering | Format |
|------|-----------|--------|
| Proteinsekvenser | `input/fasta/proteins.fasta` | Standard FASTA, header: `>UniProtIDs_{ID}_...` |
| Ligander | `input/ligand/{kategori}/{KODE}.cif` | CCD-format mmCIF; filnavn = CCD-kode |
| MSA-filer | `input/msa/{ID}.a3m` | `.a3m`; påkrevd for Boltz/RF3, ikke for AF3 |
| Oligosakkarid-def | `config/oligo_definitions.yaml` | YAML; kun for AF3 multi-monomer |

**FASTA header-format:**
```
>UniProtIDs_B6EQJ6_Aliivibrio_salmonicida_AA9_LPMO
```
Første UniProt-ID (før `;`) brukes som `protein_id`.

**CCD-kode:** Utledes av CIF-filnavnet: `STA6.cif` → kode `STA6`.

### Output

| Artefakt | Plassering | Beskrivelse |
|---------|-----------|-------------|
| Kjøringsdata | `work/{ligand_ccd}/{modell}/` | Input-JSON/YAML, SLURM-skript, status.jsonl |
| Ferdigmarkør | `work/{ligand_ccd}/{modell}/DONE.ok` | Opprettet når jobben er vellykket; brukes av resume-logikk |
| Siste kjøring | `work/{ligand_ccd}/{modell}/latest` | Symlenke til nyeste `runs/{SLURM_JOB_ID}/` |
| Jobblogger | `work/{ligand_ccd}/{modell}/slurm_*.err` / `.out` | Feil- og utdatalogg per SLURM-jobb |
| AF3 MSA-data | `work/af3_msa/{protein}_{ligand}/` | Pre-beregnede MSA-data for AF3 inference |
| Resultater | `results/` | Endelige prediksjonsresultater (antatt – strukturen avhenger av modell) |
| Manifester | `manifest/proteins.csv`, `ligands.csv`, `msa.csv`, `cases.csv` | Revisjonslogg over kjøringsplanen |

### Manifest-drevet arbeidsflyt

`cases.csv` er det autoritative kjøringsdokumentet: et kartesisk produkt av
(protein × ligand × modell). Forespørte jobber matches mot `DONE.ok`-markører for idempotent
kjøring.

### Integrasjon fra `analyse`-modulen

For å kalle pipelinen programmatisk anbefales det å:
1. Generere/oppdatere `config/pipeline.yaml` med ønskede parametre
2. Kalle `structure-pipeline run -c config/pipeline.yaml` som subprocess
3. Polle `structure-pipeline status -c config/pipeline.yaml` for fremdrift
4. Lese `manifest/cases.csv` og `work/{ligand}/{modell}/status.jsonl` for detaljert status

**Fremtidig rsync-sync til NIRD:** Anbefalt trigger er etter at `DONE.ok` er opprettet for
alle cases i et batch. Konfigurer mål-sti via miljøvariabel (f.eks. `NIRD_RESULTS_PATH`) –
ikke hardkode klyngestier. (Se `analyse`-instruksjonsfil for plan.)

---

## Konvensjoner og kvalitetskrav

**Logging:**
- Kjøringsstatus logges append-only til `work/{ligand}/{modell}/status.jsonl`
- SLURM-jobbutdata i `work/{ligand}/{modell}/slurm_*.{out,err}`
- Bruk `tail -f work/{ligand_ccd}/{modell}/slurm_*.err` for live overvåking

**Feilhåndtering:**
- `DONE.ok`-markøren opprettes kun ved vellykket kjøring
- Feila AF3-MSA-jobber blokkerer inferensjobben (manuell intervensjon nødvendig)
- Resume-logikk (`--resume`, standard) hopper over alt med `DONE.ok`
- Bruk `--no-resume` kun for tvungen rekjøring (f.eks. ved modelloppgradering)

**Konfig:**
- Alle stier i `paths:`-seksjonen MÅ være absolutte (compute-noder ser ikke relative stier)
- Støtter miljøvariabel-ekspansjon: `${ALPHAFOLD_WEIGHTS}/af3`
- `inputs:` og `outputs:` støtter relative stier (løses relativt til konfig-filens mappe)
- Modell-spesifikke parametre foretrukket via konfig eller CLI-overrides; ikke hardkod i kode

**Navngiving:**
- `protein_id`: første UniProt-ID fra FASTA-header (f.eks. `B6EQJ6`)
- `ligand_ccd`: CDD-kode fra CIF-filnavn (f.eks. `STA6`)
- SLURM-jobbnavn følger mønster: `af3_msa_{protein}_{ligand}`, `af3_{ligand}_inf`, `boltz_{ligand}`
- Kjøringskatalog: `work/{ligand_ccd}/{modell}/runs/{SLURM_JOB_ID}/`

---

## Do/Don't for endringer

### Do

- **Juster parametre via `config/pipeline.yaml`** (seeds, recycling steps, minne, tidsgrenser)
- **Bruk CLI-overrides** (`--af3-seeds`, `--boltz-recycling-steps`, osv.) for engangsbrukt
- **Legg til nye ligander** ved å plassere CIF-filer i `input/ligand/{kategori}/` og kjøre `manifest --force`
- **Legg til nye proteiner** i `input/fasta/proteins.fasta` og kjøre `manifest --force`
- **Skriv tester** i `tests/` for eventuell ny funksjonalitet
- **Oppdater `oligo_definitions.yaml`** for nye oligosakkaridtyper

### Don't

- **IKKE endre runner-logikk (`runners/*.py`)** uten å forstå AF3 tostegsfan-in og DONE.ok-semantikken
- **IKKE slett `DONE.ok`-filer** med mindre du bevisst ønsker rekjøring
- **IKKE hardkod absolutte klyngestier** i kildekode – bruk konfig
- **IKKE legg til ny funksjonalitet** uten at det er en faktisk feil eller et kritisk behov
- **IKKE endre FASTA-header-formatet** – `protein_id`-parsing er tett koblet mot manifest-generering
- **IKKE endre MSA-filnavnkonvensjonen** (`.a3m`, navngitt etter `protein_id`) uten å oppdatere `manifest/`-koden
- **IKKE kjør uten `--dry-run`** første gang etter konfigurasjonsendrigner

---

## Referanser i repoet

- [structure_pipeline/README.md](Masteroppgave/structure_pipeline/README.md) – fullstendig dokumentasjon (autoritativ)
- [structure_pipeline/pyproject.toml](Masteroppgave/structure_pipeline/pyproject.toml) – pakkemetadata og avhengigheter
- [structure_pipeline/config/pipeline.yaml](Masteroppgave/structure_pipeline/config/pipeline.yaml) – aktiv konfig med faktiske klyngestier
- [structure_pipeline/config/oligo_definitions.yaml](Masteroppgave/structure_pipeline/config/oligo_definitions.yaml) – oligosakkarid-definisjoner
- [structure_pipeline/src/structure_pipeline/cli.py](Masteroppgave/structure_pipeline/src/structure_pipeline/cli.py) – CLI-entrypoint
- [structure_pipeline/src/structure_pipeline/config.py](Masteroppgave/structure_pipeline/src/structure_pipeline/config.py) – Pydantic-konfig-modell
- [structure_pipeline/src/structure_pipeline/runners/af3.py](Masteroppgave/structure_pipeline/src/structure_pipeline/runners/af3.py) – AF3-spesifikk logikk
- [structure_pipeline/src/structure_pipeline/runners/boltz.py](Masteroppgave/structure_pipeline/src/structure_pipeline/runners/boltz.py) – Boltz-2-spesifikk logikk
- [structure_pipeline/src/structure_pipeline/runners/rf3.py](Masteroppgave/structure_pipeline/src/structure_pipeline/runners/rf3.py) – RF3-spesifikk logikk
- [structure_pipeline/src/structure_pipeline/resume.py](Masteroppgave/structure_pipeline/src/structure_pipeline/resume.py) – resume/status-sporing
- [structure_pipeline/bin/structure-pipeline](Masteroppgave/structure_pipeline/bin/structure-pipeline) – anbefalt kjøre-wrapper
