# Pseudocode - `ligands`

```text
START ligand preparation

DEFINE ligand families:
    amylose
    cellulose
    NAG_kitin

FOR each ligand family:
    COLLECT source structures from pdb/
    COLLECT molecule files from mol/
    COLLECT ligand definitions from smiles/
    COLLECT custom CCD files from cif/

    FOR each oligomer length available:
        IDENTIFY matching PDB, MOL, SMILES, and CIF files
        USE trailing number in filename as oligomer length
        GROUP files that describe the same ligand length

IF MOL files need CIF files:
    RUN convert_mol_to_cif.py
    FOR each input MOL:
        READ MOL content
        DETECT V2000 or V3000 format
        PARSE atoms, coordinates, charges, and labels
        PARSE bonds and bond orders
        CLEAN component id and atom ids
        MAP atom and bond data to CCD-style CIF fields
        WRITE one custom ligand CIF file

IF converted files need validation:
    RUN check_file_convertion.py on ligand family folder
    FOR each matched PDB/MOL/CIF triplet:
        PARSE atoms, elements, coordinates, and bonds
        COMPARE atom counts and element ordering
        COMPARE bond topology where available
        COMPARE coordinates within tolerance
        REPORT mismatches or successful conversion

IF Boltz needs ligand pickle files:
    ENTER boltz_ccd_lib/
    FOR each custom CCD CIF:
        READ CIF with gemmi
        EXTRACT component id
        EXTRACT atom loop and bond loop
        BUILD RDKit molecule
        ADD atom names, charges, bonds, and 3D conformer
        SANITIZE molecule
        SAVE as <component_id>.pkl
        VALIDATE pickle with validate_pkl.py

USE resulting ligand files in AF3 or Boltz input examples:
    af3_eksempel*.json
    boltz_eksempel*.yaml

END ligand preparation
```
