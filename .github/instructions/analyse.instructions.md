---
applyTo: "Masteroppgave/analyse/**"
---

# analyse – Repo-instruksjoner for Copilot og utviklere

> **Status: Tidlig implementasjonsstadium.** Mye av koden er pseudo-kode eller tomme stubs.
> Se `IMPLEMENTATION_PLAYBOOK.md` for prioritert 20-stegs implementasjonsrekkefølge med stopp-punkter.
> Ingen endringer, slettinger eller laging av nye filer skal skje ved kommandolinje operasjoner.

---

## Formål

`analyse/` er **orkestreringslaget** i LPMO-pipen. Den:

1. Tar predikerte strukturer (mmCIF) fra `structure_pipeline/` som input.
2. Normaliserer, validerer (QC), raffinerer (PLACER), analyserer og lager rapport.
3. Produserer reproduserbare resultater via manifest (`run_manifest.json`) med git-commit, config-hash og tool-versjoner.
4. Skal (fremtidig) auto-synkronisere resultater til NIRD via rsync etter fullfort kjøring.

**Hard regel:** Koden her kaller prediksjonskjøring via `structure_pipeline/`-adapterne.
Ingen AF3/RF3/Boltz-2-logikk skal dupliseres eller reimplementeres i `analyse/`.

---

## Nøkkelstier og struktur

```
analyse/
  src/lpmo_pipeline/
    cli.py                    ← Entrypoint: `lpmo-pipeline tune|run`
    io/                       ← mmCIF-innlesing, normalisering, protonering, CCD-oppslag
      mmcif_ingest.py
      normalize_mmcif.py
      ccd_lookup.py
      protonate_export.py
    mapping/                  ← Kryssmodell atom-mapping og omdøping
    qc/                       ← PoseBusters, Privateer, Cu-geom, QC-rapport
      gates.py                ← Pass/fail-porter
      posebusters_runner.py
      privateer_runner.py
      custom_geometry_checks.py
      qc_report.py
    placer/                   ← PLACER-wrapper (obligatorisk GPU-steg)
    analysis/                 ← MDAnalysis, ProLIF IFP, HDBSCAN, CBM, aktivitetskartlegging
    tuning/                   ← Parameter-sweep, optimalisering av HDBSCAN
    report/                   ← Bygg metrics.csv, summary.json, report.html
    utils/                    ← Manifest, logging, hashing, stier, datamodeller
  configs/
    defaults.yaml             ← Global konfig (kjede-skjema, PLACER, CCD-whitelist)
    thresholds.yaml           ← Harde porter og myke terskler
    tuning_af3.yaml           ← Tuning-grid for AF3
    tuning_boltz2.yaml        ← Tuning-grid for Boltz-2
    tuning_rf3.yaml           ← Tuning-grid for RF3
    cv_hierarchy.yaml         ← Kryss-validerings-hierarki
  schemas/                    ← JSON-skjema for I/O-validering (run_manifest, qc_report, etc.)
  tests/
    test_atom_mapping.py
    test_clustering.py
    test_io_contracts.py
    test_qc_gates.py
    fixtures/
  MASTERPLAN.md               ← Full spec (v2.1): I/O-kontrakter, 12 pipeline-steg
  IMPLEMENTATION_PLAYBOOK.md  ← Stegvis implementasjonsplan med stopp-punkter
  OPEN_QUESTIONS.md           ← Uavklarte spørsmål med foreslåtte defaults
```

---

## Miljø og avhengigheter

**Installasjon:**
```bash
pip install -e ".[dev]"   # fra analyse/
```

**Python-avhengigheter (se `pyproject.toml`):**
- `gemmi>=0.6.4` — CIF-innlesing og normalisering
- `posebusters>=0.3.0` — QC (complex mode)
- `MDAnalysis>=2.6.0` — 3D-metrikker (Cu–C1, Cu–C4, His-brace-vinkel)
- `prolif>=2.0.0` — Interaksjons-fingeravtrykk (IFP)
- `hdbscan>=0.8.33` — Klyngeanalyse med Jaccard-metrikk
- `scikit-learn>=1.3.0`, `pandas>=2.0.0`, `numpy>=1.24.0`
- `click>=8.1.0` — CLI
- `pyyaml>=6.0`, `jsonschema>=4.20.0`, `jinja2>=3.1.2`

**System-avhengigheter (må installeres separat, ikke pip):**
- `privateer` — sukkervalideringsprogram; versjon autodetektes fra `privateer --version`
- `reduce` (AmberTools) — protonering
- `PLACER` — GPU-basert raffinering; laster via `module load PLACER` (antatt, se OPEN_QUESTIONS.md pkt. 1)

---

## Kjøring og bruk

### Tuning-modus (finn optimale parametre per modell)
```bash
lpmo-pipeline tune \
  --model AF3 \
  --config analyse/configs/tuning_af3.yaml \
  --output results/tuning_af3/
```

### Produksjon (full pipeline-kjøring)
```bash
lpmo-pipeline run \
  --mode production \
  --config analyse/configs/defaults.yaml \
  --output results/del_a/ \
  --del del_a \
  [--n-jobs N]
```

Valg av `--del`:
- `del_a` — alle proteiner
- `del_b` — kun CBM-bærende LPMOer (antatt, se MASTERPLAN kap. om CBM-varianter)

### Pipen i rekkefølge (MASTERPLAN steg 1–12)
1. Innlesing og grov QC (`io/mmcif_ingest.py`)
2. Normalisering av mmCIF — kjede-omdøping A/B-D/E (`io/normalize_mmcif.py`)
3. CCD-oppslag og monosakkarid-validering (`io/ccd_lookup.py`)
4. Protonering og eksport (`io/protonate_export.py`)
5. **PLACER-raffinering (obligatorisk, GPU)** (`placer/`)
6. Full QC: PoseBusters + Privateer + Cu-geometri (`qc/`)
7. Analyse per pose: MDAnalysis + ProLIF IFP (`analysis/`)
8. Klyngeanalyse: HDBSCAN (within-run → cross-run) (`analysis/clustering_hdbscan.py`)
9. 3D-geometri per klynge
10. Krystallforankring og CBM-analyse (`analysis/crystal_anchoring.py`, `analysis/cbm_variant.py`)
11. Rapport: `metrics.csv`, `summary.json`, `report.html` (`report/`)
12. Manifest skrives (`utils/manifest.py`)

---

## Input/Output og lagring

### I/O-kontrakter (fra MASTERPLAN tabell)

| Artefakt | Format | Ansvarlig modul |
|---|---|---|
| `normalized.cif` | mmCIF | `io/normalize_mmcif.py` |
| `atom_map.tsv` | TSV | `mapping/` |
| `qc_report.json` | JSON | `qc/qc_report.py` |
| `placer_scores.json` | JSON | `placer/` |
| `ifp_matrix.csv` | CSV | `analysis/prolif_ifp.py` |
| `clusters.json` | JSON | `analysis/clustering_hdbscan.py` |
| `metrics.csv` | CSV flat | `report/build_metrics_csv.py` |
| `summary.json` | JSON | `report/build_summary_json.py` |
| `report.html` | HTML | `report/build_report_html.py` |
| `run_manifest.json` | JSON | `utils/manifest.py` |

Validering mot JSON-skjema i `schemas/` er obligatorisk for alle nøkkelartefakter.

### Stikatalog-konvensjon (fra `utils/paths.py`)
```
results/
  {del_a|del_b}/
    {model}/{protein_id}/{ligand_id}/seed_{seed}/
      normalized.cif
      atom_map.tsv
      qc_report.json
      placer_ensemble/
      ...
```

### Auto-sync til NIRD (fremtidig — ikke implementert ennå)

> **TODO — prioritert neste steg etter end-to-end-test på 1 system.**

Planen er rsync etter fullfort produksjonskjøring:
- Konfigurer kilde- og målsti via konfig (f.eks. `sync.source` og `sync.destination` i YAML).
- **Ikke hardkod cluster- eller NIRD-stier** i kode. All stilogikk via `utils/paths.py` eller konfig.
- Bruk `--checksum`-flagg i rsync (ikke bare dato/størrelse).
- Logg sync-resultat i `run_manifest.json`.
- Implementer safeguard: synkroniser kun hvis alle porter er bestått (`run_manifest.gates_passed = true`).
- Spørsmål som må avklares: NIRD-autentisering (SSH-nøkkel vs. kerberos), målstistruktur, retjobb vs. inline.

---

## Konvensjoner og kvalitetskrav

**Logging:**
- Bruk `utils/logging.py`-verktøyet (ikke bare `print`).
- Hvert pipeline-steg skal logge start, slutt og gate-resultat.
- Feil som medfører "skip"-logikk skal skrives til step-spesifikk logfil i output-mappen.

**Feilhåndtering:**
- Bruk unntaksklassene definert i `utils/exceptions.py`.
- Gate-feil (hard-fail) → skip pose/system + logg, fortsett pipeline.
- Kritisk systemfeil → re-raise, avbryt runtime med exit-kode ≠ 0.

**Konfig:**
- All konfigurasjon via YAML-filer i `configs/`. Ingen magic-verdier i kode.
- `defaults.yaml` = globale defaults; `thresholds.yaml` = harde/myke porter.
- Tuning-konfiger (`tuning_*.yaml`) er separate og skal ikke mikses med produksjonsskjøringer.
- HDBSCAN-parametre i `thresholds.yaml` merkes `locked: true` etter tuning — aldri endre i produksjon.

**Navngiving:**
- Artefaktfiler: `snake_case.{tsv,json,csv,cif,pdb,mol2,html}`.
- Python-moduler: `snake_case.py`.
- Klasser: `PascalCase`; funksjoner/variabler: `snake_case`.
- Run-ID: `{model}_{protein_id}_{ligand_id}_seed{seed}` (antatt, basert på `utils/paths.py`).

**Kjede-skjema (hard regel):**
- Protein = kjedeID `A`; glykaner = `B`, `C`, `D`; metall/Cu = `E`.
- Alle strukturer må normaliseres til dette skjemaet i steg 2 (`normalize_mmcif.py`).

**Kritiske hard-gates (fra `thresholds.yaml`):**
- `atom_mapping_coverage = 100%` — ingen unntak
- `privateer_recognized_sugars = 100%` — ingen egendefinerte oligomer-IDer
- `no_critical_posebusters_errors = true`
- `cu_his_distance`: 1.9–2.6 Å
- PLACER = obligatorisk — ingen scoring på rå prediksjoner

---

## Do/Don't for endringer

### Do
- Implementer én modul av gangen (følg IMPLEMENTATION_PLAYBOOK.md rekkefølge).
- Stopp ved hvert ⛔ STOPP-punkt og verifiser output manuelt eller med test.
- Legg all stilogikk i `utils/paths.py`.
- Valider artefakter mot skjema i `schemas/` etter generering.
- Bruk `ManifestBuilder` (`utils/manifest.py`) i alle pipeline-kjøringer.
- Skriv tester under `tests/` for nye moduler; bruk fixtures i `tests/fixtures/`.
- Nye behov fra `analyse/` som krever endringer i `structure_pipeline/` → implementer i `structure_pipeline/`-adapterne (ikke her).
- Merk uavklarte design-valg som `# ANTATT:` i kode og legg til i `OPEN_QUESTIONS.md`.

### Don't
- **Ikke reimplementer AF3/RF3/Boltz-2-kall** inne i `analyse/`. Bruk `structure_pipeline/`-adapterne.
- Ikke hardkode stier til cluster-filsystem, NIRD eller conda-miljøer.
- Ikke endre HDBSCAN-parametre etter at `locked: true` er satt i `thresholds.yaml`.
- Ikke re-tune i produksjonskjøringer.
- Ikke bruk `pipe_test/` som kilde til konfig eller importlogikk.
- Ikke hopp over PLACER-steget (selv under testing med enkelt-pose).
- Ikke legg egendefinerte oligomer-IDer (f.eks. `CELLO4`) inn i Privateer-input — kun CCD-monosakkarider.

---

## Prioriterte TODO-er (klar til implementasjon)

Basert på IMPLEMENTATION_PLAYBOOK.md (steg som er klare neste):

1. **`io/mmcif_ingest.py`** — Ekte gemmi-parsing; stopp og verifiser på én AF3-CIF.
2. **`io/normalize_mmcif.py`** — Kjede-omdøping + confidence-ekstraksjon; validér med `gemmi validate`.
3. **`io/ccd_lookup.py`** — CCD-cache + monosakkarid-whitelist; test med `test_io_contracts.py::TestCCDLookup`.
4. **`mapping/`** — Kryssmodell atom-mapping (3-tier: element → CCD → bond_graph → 3D).
5. **Rsync/NIRD-lagring** — Design konfig-styrt sync-modul etter første end-to-end-test.

---

## Referanser i repoet

Faktiske filer brukt som kilde for disse instruksjonene:

- [analyse/MASTERPLAN.md](../../analyse/MASTERPLAN.md) — Full spec: I/O-kontrakter, 12 pipeline-steg, hard-rules
- [analyse/IMPLEMENTATION_PLAYBOOK.md](../../analyse/IMPLEMENTATION_PLAYBOOK.md) — 20-stegs sekvensiert implementasjonsplan
- [analyse/OPEN_QUESTIONS.md](../../analyse/OPEN_QUESTIONS.md) — Uavklarte spørsmål med foreslåtte defaults
- [analyse/README.md](../../analyse/README.md) — Quick start og prosjektstruktur
- [analyse/pyproject.toml](../../analyse/pyproject.toml) — Avhengigheter og CLI-entrypoint (`lpmo-pipeline`)
- [analyse/configs/defaults.yaml](../../analyse/configs/defaults.yaml) — Global konfig (kjede-skjema, PLACER, CCD-whitelist)
- [analyse/configs/thresholds.yaml](../../analyse/configs/thresholds.yaml) — Harde porter og HDBSCAN-lås
- [analyse/src/lpmo_pipeline/cli.py](../../analyse/src/lpmo_pipeline/cli.py) — CLI-entrypoint med `tune` og `run` sub-kommandoer
- [analyse/src/lpmo_pipeline/utils/manifest.py](../../analyse/src/lpmo_pipeline/utils/manifest.py) — ManifestBuilder
- [analyse/src/lpmo_pipeline/utils/paths.py](../../analyse/src/lpmo_pipeline/utils/paths.py) — Kanonisk stilogikk
- [analyse/copilot.md](../../analyse/copilot.md) — Copilot-guide for mappestruktur og eierskap
- [analyse/schemas/](../../analyse/schemas/) — JSON-skjema for alle I/O-artefakter
- [analyse/tests/](../../analyse/tests/) — Test-suite med fixtures
