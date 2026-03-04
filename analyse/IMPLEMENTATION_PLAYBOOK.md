# IMPLEMENTATION PLAYBOOK

Prioritert rekkefølge for implementasjon. Hvert steg kan testes isolert
før neste startes. **Stopp-punkt** = manuell verifikasjon før du går videre.

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
   ⛔ STOPP: Inspiser for_posebusters.pdb (CONECT) og .mol2 (bond orders).

8. **`qc/posebusters_runner.py`** — PoseBusters-wrapper.
   Test: `pytest tests/test_qc_gates.py::TestPoseBustersGate`.

9. **`qc/privateer_runner.py`** — Privateer-wrapper + 100%-recog gate.
   Test: `pytest tests/test_qc_gates.py::TestPrivateerGate`.

10. **`qc/custom_geometry_checks.py`** — Cu-His 1.9–2.6 Å gate.
    Test: `pytest tests/test_qc_gates.py::TestCuHisGate`.

11. **`qc/qc_report.py`** — Samle alle QC-resultater → verdict.
    ⛔ STOPP: Kjør full QC på 3 poser, verifiser JSON-rapport mot skjema.

12. **`placer/run_placer.py`** — PLACER-wrapper (GPU).
    ⛔ STOPP: Verifiser ≥1 pose returnert, inspiser score-fil.

13. **`placer/rank_ensemble.py`** — Kompositt-scoring av PLACER-ensemble.
    Verifiser at ranked_poses.csv inneholder alle poser, korrekt sortert.

14. **`analysis/mdanalysis_metrics.py`** — 3D-metrikker per pose.
    ⛔ STOPP: Sjekk Cu-C1, Cu-C4, RMSD mot forventa størrelsesorden.

15. **`analysis/prolif_ifp.py`** — IFP-beregning med ProLIF.
    ⛔ STOPP: Visualiser IFP-matrise for 5 poser.

16. **`analysis/clustering_hdbscan.py`** — HDBSCAN med Jaccard-metrikk.
    Test: `pytest tests/test_clustering.py`.
    ⛔ STOPP: Inspiser klynge-output for ≥2 systemer.

17. **`tuning/`** — Implementer sweep_runner + tune_orchestrator.
    ⛔ STOPP: Kjør tuning på AF3 med 4 test-cases. Verifiser locked config.

18. **`analysis/` (remaining)** — cluster_signatures, crystal_anchoring,
    cbm_variant, cbm_comparison, activity_mapping, predictive_models.
    ⛔ STOPP: Verifiser cluster_signatures.json mot skjema.

19. **`report/`** — build_summary_json, build_metrics_csv, build_report_html.
    Verifiser summary.json mot schema, metrics.csv kolonner, HTML visuelt.

20. **`cli.py`** — Integrer alle steg. Kjør `lpmo-pipeline run --config ...`.
    ⛔ STOPP: End-to-end test på 1 system (AF3). Sammenlign med forventet output.
