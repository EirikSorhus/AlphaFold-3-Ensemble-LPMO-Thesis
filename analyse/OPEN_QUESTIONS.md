# OPEN QUESTIONS

Spørsmål som må avklares under implementasjon. Hvert punkt har et
foreslått default-valg slik at arbeidet kan fortsette uten blokkering.

---

1. **PLACER GPU-krav** — Trenger vi spesifikk CUDA-versjon eller holder
   standard IDUN-moduler? *Default: bruk `module load PLACER` og test.*

2. **Privateer CLI-versjon** — Privateer v1 vs v2 har ulik JSON-output.
   Hvilken versjon er installert? *Default: autodetekt fra `privateer --version`.*

3. **Reduce-versjon** — MolProbity Reduce vs AmberTools reduce?
   *Default: AmberTools `reduce` (mest tilgjengelig via conda).*

4. **Kristallstrukturer for anchoring** — Hvilke PDB-koder skal brukes
   som referanse? Trenger vi ligand-bundet + apo for alle LPMO-familier?
   *Default: bruk AA9-referanser fra litteraturen (4EIS, 5ACF, etc.).*

5. **CBM-varianter (DEL A / DEL B)** — Skal begge CBM-deletions kjøres
   for alle systemer, eller bare for CBM-bærende LPMOer?
   *Default: kun for systemer der full-length har annotert CBM.*

6. **DP-range for glykaner** — Plan sier DP2–DP8. Skal alle DP-lengder
   kjøres for hvert system, eller bare de biologisk relevante?
   *Default: DP4 + DP6 som primære, DP2/DP3/DP5/DP7/DP8 som utvidelse.*

7. **Tanimoto-terskel for crystal anchoring** — Hvilken cutoff for
   "biologically plausible"? *Default: Tanimoto ≥ 0.3 (IFP) og
   pocket-RMSD ≤ 3.0 Å som soft flags, ikke harde gates.*

8. **Antall PLACER-modes** — 50 vs 100 vs 200?
   *Default: 100 for tuning, 200 for produksjon (jf. defaults.yaml).*

9. **Predictive model — feature selection** — Skal vi bruke alle
   IFP-bits + geometri, eller bare cluster-signaturer?
   *Default: cluster-signaturer + topp-20 IFP-bits (etter variansfilter).*

10. **Reproduserbarhet — random seeds** — Skal HDBSCAN bruke fast seed?
    (HDBSCAN er deterministisk for gitt input, men ProLIF/MDAnalysis
    kan ha floating-point-variasjon.)
    *Default: sett `numpy.random.seed(42)` kun for predictive_models.*
