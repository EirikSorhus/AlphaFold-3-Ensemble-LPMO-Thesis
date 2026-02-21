# boltz_ccd_lib

Konverterer CCD-lignende mmCIF-filer (en komponent) til Boltz-kompatible RDKit Mol pickle-filer.

## Forutsetninger

Krever container med Python 3.11, gemmi, og rdkit:
- Container: `/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/`
- Build-script: `/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/build_container.sh`

## Bruk

### Konvertere CIF til PKL
```bash
export PATH="/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env/bin:$PATH"
python cif_to_pkl.py <input.cif> [output_dir]
```

Output: `<COMP_ID>.pkl` fil med RDKit Mol objekt

### Validere PKL-fil
```bash
python validate_pkl.py <ligand.pkl>
```

Sjekker:
- Antall atomer, bindinger, conformers
- Atom "name" properties
- SMILES-representasjon

### Sammenligne PKL-filer (valgfritt)
```bash
python compare_pkl_props.py <reference.pkl> <new.pkl>
```

Sammenligner mol-properties mellom to pickle-filer.

## Testing

Kjør test-suite med SLURM:
```bash
cd test
sbatch run_test.sh
```

Test-filer:
- `test/cif/` - Input CIF-filer (CEL6, NAG6, STA6)
- `test/pkl/` - Output PKL-filer
- `test/logs/` - SLURM logger

## CIF-format krav

Input CIF må inneholde:
- `_chem_comp.id` - Komponent-ID
- `loop_ _chem_comp_atom.*` med kolonner:
  - `atom_id` - Atom-navn
  - `type_symbol` - Elementtype
  - `charge` - Formal ladning
  - `pdbx_model_Cartn_x_ideal` - X-koordinat
  - `pdbx_model_Cartn_y_ideal` - Y-koordinat  
  - `pdbx_model_Cartn_z_ideal` - Z-koordinat
- `loop_ _chem_comp_bond.*` med kolonner:
  - `atom_id_1`, `atom_id_2` - Atom-par
  - `value_order` - Bindingstype (SING/DOUB/TRIP/AROM)
  - `pdbx_aromatic_flag` - Aromatisk flagg (Y/N)

## Output PKL-format

RDKit Mol objekt med:
- Atom property "name" = atom_id fra CIF
- 3D conformer med ideal-koordinater
- Mol property "MOL_NAME" = COMP_ID
- Full sanitization

