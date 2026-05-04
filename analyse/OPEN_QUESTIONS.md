# OPEN QUESTIONS

Spørsmål som må avklares under implementasjon. Hvert punkt har et
foreslått default-valg slik at arbeidet kan fortsette uten blokkering.

**⚠️ OPPDATERT 2026-04-21: Konflikter A–D avklart av bruker. Se bunnen av filen.**

---

1. ~~**PLACER GPU-krav**~~ — **AVKLART 2026-04-21: PLACER er fjernet fra analysen helt. Steg 1 er avviklet.**

2. ~~**Privateer CLI-versjon**~~ — **AVKLART 2026-04-30.**
   Privateer kjøres via SIF: `/cluster/projects/nn1003k/prog/privateer/privateer.sif`.
   Pipeline-policy: Privateer kjøres via `apptainer run --cleanenv` mot SIF (ikke via `analyse_env`/PATH).
   Avklart outputkontrakt: cluster-buildet gir ikke stabil JSON på stdout; den operative parse-kilden er `validation_data-privateer` skrevet av `-mode ccp4i2`.
   `qc/privateer_runner.py` bygger bind-aware SIF-kjøring, parser `validation_data-privateer`, henter versjon via `-list`, og beholder rå stdout/stderr kun ved feil eller eksplisitt debug-flag.
   Gjenstående arbeid er ikke outputformat-avklaring, men full hard-QC-verifikasjon på ekte poser.

3. **Reduce-versjon** — MolProbity Reduce vs AmberTools reduce?
   *Default: AmberTools `reduce` (mest tilgjengelig via conda).*

4. ~~**Pre-QC aktiv-sete terskel**~~ — **AVKLART 2026-04-21.**
   Hard cutoff for `min_cu_ligand_distance` foran PoseBusters/Privateer er
   `<= 10.0 A` (permissiv pre-QC gate). `<= 8.0 A` beholdes kun som mulig
   rapporterings-/soft-flag terskel ved behov.

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
   - AVKLART (2026-04-21): pre-QC aktiv-sete hard gate er `min_cu_ligand_distance <= 10.0 A`.
    - AVKLART (2026-04-21): PLACER er fjernet fra analysen helt.
    - AVKLART (2026-04-21): AF3-kjøringer bruker `num_diffusion_samples=5` → 75 poser per system.
    - AVKLART (2026-04-21): Kanonisk term er "protein" / `protein_id`.
    - AVKLART (2026-04-21): Fem separate pose-tabeller (ingen samlet `pose_table.tsv`).

13. **Substrat-recognition residues for alignment** — Hvordan identifisere
    surface residues involvert i substratgjenkjenning for PyMOL `pair_fit`?
    Foretrukket: litteratursøk for kjente LPMO-substrat-bindende residuer.
    Fallback: alle protein-residuer innen en cutoff (f.eks. 5 Å) fra ligand i predikert struktur.
    Merk: proximity-basert utvalg kan gi ulike residuer mellom prediksjonsmodeller.
    *Status: må avklares per LPMO-familie. Litteraturbasert er best men krever manuelt arbeid.*

---

## Nye konflikter identifisert 2026-04-21 (krever eksplisitt beslutning)

14. ~~**⚠️ KONFLIKT A — Antall prøver per system (KRITISK)**~~ — **AVKLART 2026-04-21.**

    AF3-kjøringene er fullført med `num_diffusion_samples=5` → `15 seeds × 5 samples = 75 poses` per protein–ligand betingelse. Ingen nye kjøringer med 10 samples planlegges.
    Alle plandokumenter er oppdatert til å reflektere 75 poser.

15. ~~**⚠️ KONFLIKT B — PLACER-rolle og terminologi**~~ — **AVKLART 2026-04-21.**
    PLACER fjernes fra analysen helt. Ingen kode, ingen tabeller, ingen skårer. Steg 17 avviklet.

16. ~~**⚠️ KONFLIKT C — Terminologi: "protein" vs "enzym"**~~ — **AVKLART 2026-04-21.**
    Kanonisk term er "protein" / `protein_id`. Ny kode bruker `protein_id`. Eksisterende kode
    endres ved neste refaktor (ikke blokkerende).

17. ~~**⚠️ KONFLIKT D — Output-tabellstruktur: samlet vs splittet pose-tabell**~~ — **AVKLART 2026-04-21.**
    Går over til fem separate tabeller. Ingen samlet `pose_table.tsv`. Kanoniske tabeller:
    - `pose_manifest.tsv`
    - `pose_confidence.tsv`
    - `pose_ifp_table.tsv`
    - `pose_residue_contact_table.tsv`
    - `pose_geometry.tsv`

18. **PoseBusters hard-fail policy** — Hvilke PoseBusters-feil skal telle som hard QC fail vs soft flag/pass?
   Bakgrunn: Etter korrigert PoseBusters-kontrakt for kombinerte AF3-eksporter
   (auto-splitt til ligand `mol_pred` + protein `mol_cond` i `dock`-modus)
   forsvant de tidligere falske real-case feilene `all_atoms_connected` og
   `internal_steric_clash`. I siste 3-pose PoseBusters-validering er eneste
   gjenstående PB-fail `minimum_distance_to_protein`.
   *Default: behold dagens konservative klassifisering midlertidig, men avklar
   eksplisitt om `minimum_distance_to_protein` og andre pocket-/distance-relaterte
   PoseBusters-feil skal gi `hard_fail`, `soft_flag` eller kun rapporteres.*

19. **ProLIF interaction-type pruning** — Skal endelig utvalg av ProLIF-interaksjonstyper
   bestemmes bare ved enkel sparsity-/nyttevurdering, eller er det verdt å lage en
   form for statistisk analyse/test for dette?
   *Default: ikke blokker på dette. Behold enkel manuell vurdering senere som første
   steg. Eventuell statistisk analyse er lav prioritet, ikke spesielt viktig akkurat nå,
   og kan ta tid å designe og teste på en meningsfull måte.*
