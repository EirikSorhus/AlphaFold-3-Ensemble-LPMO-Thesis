# Masteroppgave

Master's thesis in structural biology and bioinformatics at the Norwegian University of Life Sciences (NMBU). Created by Eirik Sørhus with Åsmund Røhr Kjendseth as supervisor.

**Thesis title:** *Ensemble-Based Analysis of AlphaFold 3-Predicted LPMO–Oligosaccharide Complexes for Substrate Recognition*
**Thesis handle:** https://hdl.handle.net/11250/5554967

---

## Repository Structure

| Directory | Description |
|-----------|-------------|
| [`structure_pipeline/`](structure_pipeline/) | Production manifest-driven pipeline for running protein structure predictions (AlphaFold 3, Boltz-2, RoseTTAFold 3) at scale on HPC clusters via SLURM. Handles input preparation, job submission, batching, and result collection. |
| [`analyse/`](analyse/) | Analysis pipeline for AF3-predicted LPMO–oligosaccharide complexes. Covers structure normalisation, quality control (PoseBusters, Privateer, Cu-geometry), interaction fingerprinting (ProLIF), binding-mode clustering, and summary reporting. |
| [`ligands/`](ligands/) | Ligand files and format-conversion helpers for oligosaccharide substrates (chitin, cellulose, amylose). Includes CCD-format CIF libraries used by the prediction pipeline. |
| [`pipe_test/`](pipe_test/) | Utility scripts for sequence retrieval, signal-peptide processing, and earlier experimental pipeline modules. The stable working area is `pipe_test/scripts/module_2_new/`. |

---

## Getting Started

Each subdirectory contains its own README with detailed usage instructions. Start with the component relevant to your task:

- **Running predictions →** see [`structure_pipeline/README.md`](structure_pipeline/README.md)
- **Analysing results →** see [`analyse/README.md`](analyse/README.md)
- **Ligand preparation →** see [`ligands/`](ligands/)
- **Sequence retrieval & helpers →** see [`pipe_test/scripts/module_2_new/README.md`](pipe_test/scripts/module_2_new/README.md)

## Environment

All heavy computation is designed for HPC execution with SLURM scheduling (tested on Olivia). Conda environments are managed under `/cluster/work/projects/nn1003k/eirik/conda/` using containerised wrappers (`hpc-container-wrapper`).
