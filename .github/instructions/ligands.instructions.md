---
applyTo: "Masteroppgave/ligands/**"
---

# ligands – Repo-instruksjoner for Copilot og utviklere

> Ingen endringer, slettinger eller laging av nye filer skal skje ved kommandolinje operasjoner.

## Formål

Mappen inneholder ligand-molekyler (polysakkarider) i ulike filformater, brukt som input til strukturprediksjonspipelines (AlphaFold 3, Boltz 2, RoseTTAFold 3). I tillegg ligger det konverteringsskript og et hjelpemakrobibliotek (`boltz_ccd_lib`) for å produsere Boltz-kompatible pickle-filer fra CIF-input.

Ligandfamilier som finnes i repoet:
- **NAG/kitin** – N-acetylglukosamin-oligos (dimerer–oktamerer)
- **Amylose** – glukose α-1,4-kjedede oligos
- **Cellulose** – glukose β-1,4-kjedede oligos

---

## Nøkkelstier og struktur

```
ligands/
├── convert_mol_to_cif.py          # MOL → CIF-konverterer (AF3 custom CCD-format)
├── check_file_convertion.py       # Verifiserer at MOL, CIF og PDB samsvarer
├── {familie}/                     # En undermappe per ligandfamilie
│   ├── mol/                       # Råfiler: MOL V2000/V3000 (kilde)
│   ├── smiles/                    # Råfiler: SMILES-strenger
│   ├── pdb/                       # Råfiler: PDB fra 3D-minimering
│   ├── cif/                       # Avledede filer: CIF (AF3 custom CCD-format)
│   ├── af3_eksempel_*.json        # Eksempel-input til AlphaFold 3
│   └── boltz_eksempel_*.yaml      # Eksempel-input til Boltz 2
└── boltz_ccd_lib/
    ├── cif_to_pkl.py              # CIF → RDKit Mol pickle (for Boltz 2)
    ├── validate_pkl.py            # Validerer output-pkl
    ├── compare_pkl_props.py       # Sammenligner to pkl-filer
    └── test/                      # Testsuite (SLURM)
```

Navnemønster for ligandfiler: `{familie}_{n}.{ext}` eller `{FAMILIE}{n}.cif`, der `n` er oligomerlengde (2–8).

---

## Miljø og avhengigheter

### convert_mol_to_cif.py
- Kun standard-bibliotek + `re`, `pathlib`, `argparse` — **ingen spesialcontainer nødvendig**.

### check_file_convertion.py
- Kun standard-bibliotek + `math`, `re`, `dataclasses`, `pathlib`.

### boltz_ccd_lib
- Krever Python 3.11, `gemmi` og `rdkit`.
- Container: `/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/`
- Aktivering:
  ```bash
  export PATH="/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/bin:$PATH"
  ```

---

## Kjøring og bruk

### MOL → CIF (for AlphaFold 3 custom CCD)
```bash
python ligands/convert_mol_to_cif.py <mappe_med_mol_filer> --output-dir <output_cif_mappe>
```
Eksempel for NAG/kitin:
```bash
python ligands/convert_mol_to_cif.py ligands/NAG_kitin/mol --output-dir ligands/NAG_kitin/cif
```
Skriptet håndterer både V2000 og V3000 MOL-format, inkludert kjente formateringsfeil (sammenkjedede atom-/bindingstall).

### CIF → PKL (for Boltz 2)
```bash
# Krever boltz_ccd_env
export PATH="/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/bin:$PATH"
python ligands/boltz_ccd_lib/cif_to_pkl.py <input.cif> [output_dir]
```

### Verifisering av konvertering
```bash
python ligands/check_file_convertion.py <base_familje_mappe> --tol 1e-5
```
Sjekker at MOL, CIF og PDB har konsistente atomtall og koordinater for samme oligomerlengde.

### Validering av PKL
```bash
python ligands/boltz_ccd_lib/validate_pkl.py <ligand.pkl>
```

### Test av boltz_ccd_lib
```bash
cd ligands/boltz_ccd_lib/test
sbatch run_test.sh
```

---

## Input/Output og lagring

| Filtype | Rolle | Kategorisering |
|---------|-------|----------------|
| `mol/` | Rå molekylstruktur (V2000/V3000) | **Råfil – ikke overskriv** |
| `smiles/` | SMILES-representasjon | **Råfil – ikke overskriv** |
| `pdb/` | 3D-minimerte strukturer | **Råfil – ikke overskriv** |
| `cif/` | Konvertert CIF (AF3 CCD-format) | **Avledet** – kan regenereres med `convert_mol_to_cif.py` |
| `*.pkl` | RDKit Mol pickle (Boltz 2) | **Avledet** – kan regenereres med `cif_to_pkl.py` |

**Policy: Overskriv aldri råfiler (`mol/`, `smiles/`, `pdb/`).** CIF- og PKL-filer er avledede produkter og kan trygt regenereres.

AF3-eksempel-JSON og Boltz-eksempel-YAML er maler – de skal ikke overskrives; kopier dem ved behov.

---

## Konvensjoner og kvalitetskrav

**Logging:**
- `convert_mol_to_cif.py` printer `"Wrote <sti>"` til stdout for hver produsert fil.
- `cif_to_pkl.py` printer status og eventuelle advarsler til stdout.
- Feil kastes som `ValueError` med beskrivende melding.

**Feilhåndtering:**
- MOL-parsing håndterer kjente formateringsfeil (V2000 sammenkjedede teller, V3000 nøkkelverdi-props).
- CIF-parsing med gemmi validerer obligatoriske kolonner (`atom_id`, `type_symbol`, koordinater) og kaster feil ved manglende data.
- `check_file_convertion.py` matcher filer på oligomerlengde-suffiks (`trailing_int`), uavhengig av familienavn.

**Konfig:**
- Ingen konfigurasjonsfiler — alle parametere sendes som CLI-argumenter.
- `--tol` i `check_file_convertion.py` styrer koordinat-toleransen (standard 1e-5 Å).

**Navngiving:**
- Ligand-CIF-komponenter: `{FAMILIE}{n}` (f.eks. `NAGCHITIN3`, `AMYLASE4`, `CEL6`) — all caps, uten skilletegn.
- MOL/SMILES-filer: `{familie}_{n}.{ext}` (f.eks. `nag_chitin_3.mol`) — lowercase, underscore.
- PDB-filer: full kjemisk IUPAC-notasjon (autogenerert fra ekstern kilde — behold som de er).

---

## Do/Don't for endringer

### Do
- Legg til nye ligandfamilier i en ny undermappe med samme struktur: `mol/`, `smiles/`, `pdb/`, `cif/`.
- Regenerer CIF-filer med `convert_mol_to_cif.py` hvis MOL-filer oppdateres.
- Regenerer PKL-filer med `cif_to_pkl.py` etter CIF-oppdateringer.
- Kjør `check_file_convertion.py` og `validate_pkl.py` etter alle konverteringer for å verifisere konsistens.
- Oppdater AF3-JSON og Boltz-YAML-eksempler ved nye ligandfamilier (kopier eksisterende mal).

### Don't
- Ikke overskriv eller slett råfiler i `mol/`, `smiles/`, `pdb/`.
- Ikke endre CIF-komponent-ID-er manuelt etter at de er brukt i pipeline-kjøringer (bryter sporbarhet).
- Ikke bruk `boltz_ccd_lib`-skript utenfor `boltz_ccd_env`-containeren — `gemmi` og `rdkit` er ikke tilgjengelig i andre miljøer.
- Ikke legg grafisk visualisering eller tunge avhengigheter inn i konverteringsskriptene — de skal forbli lette og kjørbare uten container.
- Ikke hardkod clusterstier i koden; bruk relative stier eller CLI-argumenter.

---

## Referanser i repoet

- [ligands/convert_mol_to_cif.py](../../ligands/convert_mol_to_cif.py) — MOL-til-CIF-konverterer (V2000/V3000)
- [ligands/check_file_convertion.py](../../ligands/check_file_convertion.py) — Konverteringsverifisering
- [ligands/boltz_ccd_lib/README.md](../../ligands/boltz_ccd_lib/README.md) — Bruksdokumentasjon for CIF→PKL
- [ligands/boltz_ccd_lib/CODE_WALKTHROUGH.md](../../ligands/boltz_ccd_lib/CODE_WALKTHROUGH.md) — Detaljert kode-walkthrough
- [ligands/boltz_ccd_lib/cif_to_pkl.py](../../ligands/boltz_ccd_lib/cif_to_pkl.py) — CIF-til-pickle-konverterer
- [ligands/NAG_kitin/af3_eksempel_nag_chitin.json](../../ligands/NAG_kitin/af3_eksempel_nag_chitin.json) — AF3 input-eksempel
- [ligands/NAG_kitin/boltz_eksempel_nag_chitin.yaml](../../ligands/NAG_kitin/boltz_eksempel_nag_chitin.yaml) — Boltz 2 input-eksempel
