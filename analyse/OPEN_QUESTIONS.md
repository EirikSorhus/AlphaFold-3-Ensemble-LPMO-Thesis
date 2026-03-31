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

4. **Pre-QC aktiv-sete terskel** — Hvilken hard cutoff for
   `min_cu_ligand_distance` for gate foran PoseBusters/Privateer?
   *Default: <= 10.0 A (hard gate), med Cu-C1/C4 fortsatt logget som metrikker.*

5. **Kristallstrukturer for anchoring** — Hvilke PDB-koder skal brukes
   som referanse? Trenger vi ligand-bundet + apo for alle LPMO-familier?
   *Default: bruk AA9-referanser fra litteraturen (4EIS, 5ACF, etc.).*

6. **CBM-varianter (DEL A / DEL B)** — Skal begge CBM-deletions kjøres
   for alle systemer, eller bare for CBM-bærende LPMOer?
   *Default: kun for systemer der full-length har annotert CBM.*

7. **DP-scope for aktiv implementasjon** — Skal andre DP enn 4/6/8 inn i
   samme hovedpipeline, eller beholdes de som egne sideanalyser?
   *Default: hovedpipeline = DP4/DP6/DP8 (9 uavhengige delanalyser).* 

8. **Tanimoto-terskel for crystal anchoring** — Hvilken cutoff for
   "biologically plausible"? *Default: Tanimoto >= 0.3 (IFP) og
   pocket-RMSD <= 3.0 A som soft flags, ikke harde gates.*

9. **R-modellvalg for cluster-rader** — Hvilken primarmodell skal brukes
   i R for regioselektivitet (glmnet vs glmer)?
   *Default: penalized logistisk regresjon (glmnet), mixed model som sensitivitet.*

10. **EC 1.14.99.- ikke-AA17 mapping** — Hvilken endelig tekstetikett og
    hvilket standardsubstrat for "xylan ol"-tilfeller?
    *Default: substrate_class=`xylan_or_other`, regio_class=`unknown`, aktivitet=`xylan_like_oxidative`.*

11. **Geometri-planaritet** — Operasjonelle planaritetskrav mangler forelopig.
    Dette ma spesifiseres (metode + terskler) for endelig analyse og rapportering.
    *Status: avventer definisjon fra prosjektleder.*

12. **Statusoppsummering av avklarte punkter**
    - AVKLART: hovedanalyse aggregerer ikke cluster -> enzym.
    - AVKLART: tuning skjer etter hovedanalyse (valgfritt, hvis tid).
    - AVKLART: pre-QC gate foran PoseBusters/Privateer er obligatorisk.
    - AVKLART: beregnede numeriske metrikker beholdes i output.
    - AVKLART: PoseBusters kjøres via SIF-container (`/cluster/projects/nn1003k/prog/posebusters/`).
    - AVKLART: CIF→PDB konvertering bruker PDBFixer (ikke Biopython), adaptert fra PoseBench.
    - AVKLART: Crystal anchoring bruker PyMOL `pair_fit` for optimal lokal superposisjon.
    - AVKLART: Pipeline har gått fra pseudokode til steg-for-steg implementasjon (2026-03-26).

13. **Substrat-recognition residues for alignment** — Hvordan identifisere
    surface residues involvert i substratgjenkjenning for PyMOL `pair_fit`?
    Foretrukket: litteratursøk for kjente LPMO-substrat-bindende residuer.
    Fallback: alle protein-residuer innen en cutoff (f.eks. 5 Å) fra ligand i predikert struktur.
    Merk: proximity-basert utvalg kan gi ulike residuer mellom prediksjonsmodeller.
    *Status: må avklares per LPMO-familie. Litteraturbasert er best men krever manuelt arbeid.*
