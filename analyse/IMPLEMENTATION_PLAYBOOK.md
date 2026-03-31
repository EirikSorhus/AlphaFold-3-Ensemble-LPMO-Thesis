# IMPLEMENTATION PLAYBOOK

Prioritert rekkefølge for implementasjon. Hvert steg kan testes isolert
før neste startes. **Stopp-punkt** = manuell verifikasjon før du går videre.

**⚠️ VIKTIG: Rekkefølgen på stegene og hvilke analyser som skal kjøres skal ikke endres av AI uten godkjenning fra deg først.**

**Status 2026-03-26: Overgang fra pseudokode til implementasjon.**
Hvert steg implementeres, testes og verifiseres steg for steg.

Forutsetninger (18.03.2026, oppdatert 26.03):
- Strukturprediksjon (kun AF3) kjores separat fra analysekoden og produseres kun en gang per valgt parameteroppsett. (**RF3/RosettaFold 3 og Boltz-2 er ekskludert fra analyse.**)
- AF3 faste parametre (hovedanalyse): `num_recycles=10`, `num_seeds=15`, `num_diffusion_samples=5`.
- Denne playbooken dekker analysepipen som konsumerer ferdige prediksjonsartifakter.
- Hovedanalyse kjores forst. Tuning er valgfri etteranalyse/sensitivitet.
- Primarenhet etter posegenerering er cluster (ikke tidlig aggregering til enzym).
- Beregnede numeriske metrikker beholdes i tabeller/logg og skal ikke slettes ved gating.
- Bruk R der det er praktisk for beskrivende/prediktiv statistikk.
- PoseBusters kjøres via SIF-container: `/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif` (primær), `/cluster/projects/nn1003k/prog/posebusters/posebusters.sif` (fallback).
- Conda env: `/cluster/work/projects/nn1003k/eirik/conda/analyse_env/`.
- CIF→PDB konvertering bruker PDBFixer (OpenMM), adaptert fra PoseBench (MIT).
- Crystal anchoring alignment bruker PyMOL `pair_fit` på histidine-brace, Cu-nære residuer og substrat-recognition surface residues.

---

1. **Installer avhengigheter** — `pip install -e ".[dev]"` fra pyproject.toml.
   Verifiser at gemmi, posebusters, hdbscan, prolif, mdanalysis importeres.

2. **`io/mmcif_ingest.py`** — Implementer `ingest()` med ekte gemmi-parsing.
   ⛔ STOPP: Kjør på én AF3-CIF og inspiser resultatet manuelt.

3. **`io/normalize_mmcif.py`** — Chain-rename (A/B-D/E), confidence-ekstraksjon.
   ⛔ STOPP: Verifiser `normalized.cif` med `gemmi validate`.

4. **`io/ccd_lookup.py`** — CCD-cache + monosakkarid-validering.
   Test: `pytest tests/test_io_contracts.py::TestCCDLookup -v`.

5. **`mapping/cross_model_atom_mapping.py`** — 3-tier atom-matching.
   ⛔ STOPP: Inspiser `atom_map.tsv` for én kjent NAG-residue. Coverage=100%?

6. **`mapping/rename_atoms.py`** — Last og anvend atom_map.
   Test: round-trip test med atom_map.tsv.

7. **`io/protonate_export.py`** — Reduce + OpenBabel for protonering.
   Forutsetning: `io/cif_to_pdb.py` omgjører mmCIF→PDB via PDBFixer først.
   ⛔ STOPP: Inspiser for_posebusters.pdb (CONECT) og .mol2 (bond orders).

7b. **`io/cif_to_pdb.py`** — CIF→PDB konvertering via PDBFixer (adaptert fra PoseBench MIT).
    Bruk PDBFixer for å håndtere AF3 mmCIF-format korrekt.
    ⛔ STOPP: Verifiser at output PDB har CONECT records og korrekt bondgraf.

8. **`qc/active_site_proximity.py`** — pre-QC gate (ligand nar aktivt sete).
   Krav: beregn og logg `min_cu_ligand_distance`, `min_cu_c1`, `min_cu_c4` for alle poser.
   ⛔ STOPP: Verifiser at poser med stor avstand droppes for PB/Privateer, men at metrikker beholdes i `pose_table`.

9. **`qc/posebusters_runner.py`** — PoseBusters-wrapper.
   Bruker PoseBusters Python API (`PoseBusters(config="redock")`, `bust_table(mol_table)`)
   for batch-kjøring. Alternativt CLI via SIF-container som fallback.
   SIF: `/cluster/projects/nn1003k/prog/posebusters/build/posebusters.sif`.
   Test: `pytest tests/test_qc_gates.py::TestPoseBustersGate`.

10. **`qc/privateer_runner.py`** — Privateer-wrapper + 100%-recog gate.
   Test: `pytest tests/test_qc_gates.py::TestPrivateerGate`.

11. **`qc/custom_geometry_checks.py`** — Cu-His 1.9–2.6 Å gate.
    Test: `pytest tests/test_qc_gates.py::TestCuHisGate`.

12. **`qc/qc_report.py`** — Samle pre-QC + PB + Privateer + geometri → verdict.
    ⛔ STOPP: Kjør full QC på 3 poser, verifiser JSON-rapport mot skjema.

13. **`placer/run_placer.py`** — PLACER-wrapper (GPU).
    ⛔ STOPP: Verifiser ≥1 pose returnert, inspiser score-fil.

14. **`placer/rank_ensemble.py`** — Kompositt-scoring av PLACER-ensemble.
    Verifiser at ranked_poses.csv inneholder alle poser, korrekt sortert.

15. **`analysis/prolif_ifp.py`** — IFP-beregning med ProLIF.
    ⛔ STOPP: Visualiser IFP-matrise for 5 poser.

16. **`analysis/mdanalysis_metrics.py`** — 3D-metrikker per pose.
    ⛔ STOPP: Sjekk Cu-C1, Cu-C4, RMSD mot forventa størrelsesorden.

17. **`analysis/clustering_hdbscan.py`** — HDBSCAN med Jaccard-metrikk pa IFP-only.
    Kjor innen-modell og kryss-modell clustering separat.
    Test: `pytest tests/test_clustering.py`.
    ⛔ STOPP: Inspiser klynge-output for ≥2 systemer.

18. **`analysis/cluster_signatures.py`** — clusterannotering med median/IQR.
    Krav: geometri brukes som annotering etter clustering, ikke som clusterinput.

19. **`analysis/activity_mapping.py`** — bygg `predictive_cluster_table` fra cluster-rader.
    Krav: ingen primar aggregering cluster->enzym.

20. **`scripts/ec_activity_mapping.R`** — metadata EC# -> aktivitet/substrat/regio.
    Krav: implementer eksplisitte regler for 1.14.99.53/54/55/56 og 1.14.99.- med AA17-spesialtilfelle.

21. **`analysis/predictive_models.py` + R scripts** — cluster-baserte modeller.
    Krav: grouped CV pa enzymniva, cluster-rader beholdes, avhengighet modelleres.

22. **`analysis/crystal_anchoring.py`, `cbm_variant.py`, `cbm_comparison.py`**
    som sekundaranalyser.
    Crystal anchoring bruker PyMOL `pair_fit` for optimal lokal superposisjon:
    - Align på histidine-brace (Cα/Nε2), Cu-nære residuer, og substrat-recognition residues.
    - Substrat-recognition residues: litteraturbasert (foretrukket) eller proximity-basert (alle residuer innen cutoff fra ligand).
    - Mål pocket RMSD etter optimal lokal alignment. Dokumenter tydelig at det er optimalisert alignment.
    - Inspirert av benchmarking-af3 (MIT) sitt pocket-residue-mønster, men bruker PyMOL istedenfor APoc/Biopython.
    ⛔ STOPP: Verifiser cluster_signatures/predictive_cluster_table mot skjema.

23. **`report/`** — build_summary_json, build_metrics_csv, build_report_html.
    Verifiser summary.json mot schema, metrics.csv kolonner, HTML visuelt.

24. **`cli.py`** — Integrer alle steg. Kjor `lpmo-pipeline run --config ...`.
    ⛔ STOPP: End-to-end test på 1 system (AF3). Sammenlign med forventet output.

25. **`tuning/` (valgfritt etteranalyse)** — implementer sweep_runner + tune_orchestrator.
    Kjor sammenlignende tuning etter baseline leveranse, dokumenter eventuell parameteroppdatering.
