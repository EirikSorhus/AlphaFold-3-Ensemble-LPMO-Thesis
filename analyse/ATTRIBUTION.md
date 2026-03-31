# Third-Party Attribution

This project adapts patterns and approaches from the following open-source projects.
No code is copied verbatim — implementations are rewritten for the LPMO analysis pipeline.
Where adapted, source files contain inline attribution comments.

---

## PoseBench

- **Repository**: https://github.com/BioinfoMachineLearning/PoseBench
- **License**: MIT (see THIRD_PARTY_LICENSES/PoseBench-MIT.txt)
- **Copyright**: Copyright (c) PoseBench Authors
- **What we adapted**:
  - CIF→PDB conversion approach using PDBFixer (OpenMM) from `af3_output_extraction.py`.
    Biopython does not correctly handle mmCIF files produced by AlphaFold 3; PDBFixer is required.
  - PoseBusters Python API usage pattern (`PoseBusters(config="redock")`, `bust_table(mol_table)`)
    from `inference_analysis.py`.
  - mol_table DataFrame structure (mol_cond, mol_true, mol_pred columns) for batch scoring.
- **Files affected in this project**:
  - `src/lpmo_pipeline/io/cif_to_pdb.py`
  - `src/lpmo_pipeline/qc/posebusters_runner.py`

---

## benchmarking-af3

- **Repository**: https://github.com/lyulab/benchmarking-af3
- **License**: MIT (see THIRD_PARTY_LICENSES/benchmarking-af3-MIT.txt)
- **Copyright**: Copyright (c) 2025 Aakash Davasam
- **What we adapted**:
  - Pocket residue identification via spatial proximity cutoff, inspired by `out_of_sample/3_find_pocket_residues.py`.
    Our implementation uses gemmi/PyMOL instead of Biopython NeighborSearch.
  - Concept of residue mapping between reference and predicted structures.
  - Metrics aggregation pattern from `out_of_sample/9_save_metrics.py`.
- **What we do NOT use**:
  - APoc alignment tool (we use PyMOL `pair_fit` instead).
  - DockRMSD tool.
  - Biopython-based structure parsing (we use gemmi).
- **Files affected in this project**:
  - `src/lpmo_pipeline/analysis/crystal_anchoring.py`
