# OPEN QUESTIONS

Spørsmål som må avklares under implementasjon. Hvert punkt har et
foreslått default-valg slik at arbeidet kan fortsette uten blokkering.

**OPPDATERT 2026-05-22:** ProLIF-spørsmålet om sparse eksplisitte H-bonds er
avklart ved implementasjon av implicit-H ProLIF. Aktivt råsett er
`ImplicitHBAcceptor`, `ImplicitHBDonor`, `VdWContact`; main clustering bruker
kun implicit H-bonds. VdW-only/low-specific-contact-diagnostikk beholdes.

**⚠️ OPPDATERT 2026-04-21: Konflikter A–D avklart av bruker. Se bunnen av filen.**

---

1. ~~**PLACER GPU-krav**~~ — **AVKLART 2026-04-21: PLACER er fjernet fra analysen helt. Steg 1 er avviklet.**

2. ~~**Privateer CLI-versjon**~~ — **AVKLART 2026-04-30.**
   Privateer kjøres via SIF: `/cluster/projects/nn1003k/prog/privateer/privateer.sif`.
   Pipeline-policy: Privateer kjøres via `apptainer run --cleanenv` mot SIF (ikke via `analyse_env`/PATH).
   Avklart outputkontrakt: cluster-buildet gir ikke stabil JSON på stdout; den operative parse-kilden er `validation_data-privateer` skrevet av `-mode ccp4i2`.
   `qc/privateer_runner.py` bygger bind-aware SIF-kjøring, parser `validation_data-privateer`, henter versjon via `-list`, og beholder rå stdout/stderr kun ved feil eller eksplisitt debug-flag.
   Gjenstående arbeid er ikke outputformat-avklaring, men full hard-QC-verifikasjon på ekte poser.

3. ~~**Reduce-versjon**~~ — **AVKLART 2026-05-23.**
   Endelig valg: AmberTools `reduce`.

4. ~~**Pre-QC aktiv-sete terskel**~~ — **AVKLART 2026-04-21.**
   Hard cutoff for `min_cu_ligand_distance` foran PoseBusters/Privateer er
   `<= 10.0 A` (permissiv pre-QC gate). `<= 8.0 A` beholdes kun som mulig
   rapporterings-/soft-flag terskel ved behov.

5. ~~**Kristallstrukturer for anchoring**~~ — **AVKLART 2026-05-23.**
   Autoritativt referansesett for aktiv analyse fryses til
   `input_data/pdb_structure_data.csv` + filer under `crystal_structures/`, med
   kuratert metadata skrevet til `metadata/crystal_reference_list.tsv`.
   Nye/endrede referanser skjer kun som eksplisitt dataoppdatering i disse
   filene (ikke via kodeendring).

6. ~~**CBM-varianter (DEL A / DEL B)**~~ — **AVKLART 2026-05-23.**
   Endelig scope: kun for systemer der full-length har annotert CBM.

7. ~~**DP-scope for aktiv implementasjon**~~ — **AVKLART 2026-05-23.**
   Endelig scope: kun DP4/DP6/DP8 i hovedpipeline (9 uavhengige delanalyser).

8. ~~**Tanimoto-terskel for crystal anchoring**~~ — **AVKLART 2026-05-23.**
    Endelig biologisk tolkning i aktiv analyse:
    - Sammenligning er bare IFP-tolkbar når `ifp_comparison_eligible=True`.
    - `ifp_tanimoto >= 0.50` tolkes som moderat/stottende crystal-IFP-overlapp.
    - `0.30 <= ifp_tanimoto < 0.50` tolkes som svak/stoyutsatt stotte.
    - `ifp_tanimoto < 0.30` tolkes som lav overlapp.
    - `local_pocket_rmsd < 2.5 A` beholdes som binding-region-stotteflagg,
       ikke alene som biologisk bevis.
    - VdW-dominans vurderes via `crystal_ifp_diagnostic_summary.tsv`:
       hvis `contact_eligible_fraction < 0.50` eller `n_vdw_only > n_contact_eligible`,
       nedgraderes crystal-IFP-tolkning til "geometry-first, IFP-limited".
   Historisk kontekst: tidligere statusnotater (2026-05-08/16/20) er nå
   overstyrt av denne avklarte policyen.

9. ~~**Predictive modellvalg og prediksjonsvariabler**~~ — **AVKLART
   2026-05-22 for aktiv implementasjon.**
   Primær backend er compact Python/sklearn logistic regression med L2,
   `class_weight="balanced"`, fold-lokal numerisk preprocessing og grouped CV
   på `protein_id`. Prediktorsettene følger de reviderte kompakte
   aktivitetsplanene. `run_analysis_core` kan kjøre Stage 14 direkte fra den
   produserte `condition_table.tsv` når `production.predictive.enabled=true` og
   proteinmetadata er satt. Gjenstående arbeid er full-run tolkning etter at
   `predictive_summary.json["validation"]` viser nok rader, begge targetklasser
   og evaluerbare folds, ikke valg av backend.

10. ~~**EC 1.14.99.- ikke-AA17 mapping**~~ — **AVKLART 2026-05-23.**
   Endelig mapping for "xylan ol"-tilfeller:
   `substrate_class=xylan_or_other`, `regio_class=unknown`,
   `aktivitet=xylan_like_oxidative`.

11. **Geometri-planaritet** — Avklart 2026-05-22: separat planaritetsgate er
    ikke del av aktiv analysekontrakt. Sluttbrukerflaten bruker de implementerte
    downstream-feltene `sugar_face_orientation`, `ring_normal_vs_brace_normal`,
    `Cu_oxyl_H_C1/C4_angle`, `oxyl_H_*_distance` og `geometry_status_C1/C4`.
    Nye planaritetskrav skal bare legges til som en ny eksplisitt feature hvis
    biologisk terskel og metode bestemmes senere.

12. **Statusoppsummering av avklarte punkter**
    - AVKLART: hovedanalyse aggregerer ikke cluster -> enzym.
    - AVKLART: tuning skjer etter hovedanalyse (valgfritt, hvis tid).
    - AVKLART: pre-QC gate foran PoseBusters/Privateer er obligatorisk.
    - AVKLART: beregnede numeriske metrikker beholdes i output.
    - AVKLART: PoseBusters kjøres via SIF-container (`/cluster/projects/nn1003k/prog/posebusters/`).
    - AVKLART: CIF→PDB konvertering bruker PDBFixer (ikke Biopython), adaptert fra PoseBench.
   - AVKLART: dagens operative crystal-anchoring-RMSD bruker lokal gemmi/numpy Kabsch-superposisjon på delte pocket C-alpha-atomer; PyMOL `pair_fit` er ikke operativ backend per i dag.
    - AVKLART: Pipeline har gått fra pseudokode til steg-for-steg implementasjon (2026-03-26).
   - AVKLART (2026-04-21): pre-QC aktiv-sete hard gate er `min_cu_ligand_distance <= 10.0 A`.
    - AVKLART (2026-04-21): PLACER er fjernet fra analysen helt.
    - AVKLART (2026-04-21): AF3-kjøringer bruker `num_diffusion_samples=5` → 75 poser per system.
    - AVKLART (2026-04-21): Kanonisk term er "protein" / `protein_id`.
    - AVKLART (2026-04-21): Fem separate pose-tabeller (ingen samlet `pose_table.tsv`).
      - AVKLART (2026-05-23): Primær clusteringmetode er HDBSCAN Jaccard på
         contact-eligible IFP rows med `min_cluster_size=5`, `min_samples=null`
         og `cluster_selection_method=eom`. Begrunnelsen er bedre balanse mellom
         cluster recovery, noise og cluster-granularitet enn både HDBSCAN
         `min_cluster_size=3` og agglomerative `distance_threshold=0.55`,
         `min_cluster_size=5`. Sensitivitet: HDBSCAN `min_cluster_size=3`,
         agglomerative `distance_threshold=0.55`, `min_cluster_size=5`, og
         HDBSCAN `min_cluster_size=10` som konservativ negativ kontroll.
    - AVKLART (2026-05-22): Primær predictive backend er compact sklearn
      logistic regression, kjørt som valgfri Stage 14 fra main analysis output
      når metadata er konfigurert.
      - AVKLART (2026-05-23): Reduce-versjon er AmberTools `reduce`.
      - AVKLART (2026-05-23): CBM DEL A/DEL B kjøres kun for systemer der
         full-length har annotert CBM.
      - AVKLART (2026-05-23): DP-scope i hovedpipeline er kun DP4/DP6/DP8.
      - AVKLART (2026-05-23): EC 1.14.99.- ikke-AA17 "xylan ol" mapping er
         `substrate_class=xylan_or_other`, `regio_class=unknown`,
         `aktivitet=xylan_like_oxidative`.
      - AVKLART (2026-05-23): Endelig ProLIF interaction-type utvalg i aktiv
         analyse er implicit H-bond og VdW-interaksjoner (implementert i kode).

13. ~~**Substrat-recognition/pocket-residuer for crystal anchoring**~~ —
   **AVKLART 2026-05-23.**
   Aktiv analyse fryser dagens operative pocket-heuristikk (5 A fra ligand/Cu,
   med apo-sekvensprojeksjon og residunavn-normalisering) som beslutningsgrunnlag.
   Familievis residue-kilde dokumenteres eksplisitt i
   `metadata/alignment_residue_definitions.tsv`.
   I denne release er kilden satt til `proximity_fallback` for alle familier.

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

18. ~~**PoseBusters hard-fail policy**~~ — **AVKLART 2026-05-23.**
   Bakgrunn: Etter korrigert PoseBusters-kontrakt for kombinerte AF3-eksporter
   (auto-splitt til ligand `mol_pred` + protein `mol_cond` i `dock`-modus)
   forsvant de tidligere falske real-case feilene `all_atoms_connected` og
   `internal_steric_clash`. I siste 3-pose PoseBusters-validering er eneste
   gjenstående PB-fail `minimum_distance_to_protein`.
   Viktig presisering: i PoseBusters `dock`-mode er
   `protein-ligand_maximum_distance` far-away-testen fra
   `posebusters.modules.intermolecular_distance.check_intermolecular_distance()`
   med standard `max_distance=5.0 A` og `search_distance=6.0 A`.
   `minimum_distance_to_protein` er derimot den renamed `no_clashes`-utgangen,
   ikke selve avstandsterskelen.
   Endelig severity-policy i aktiv analyse:
   - behold PoseBusters `dock`-defaults uendret,
   - behold `minimum_distance_to_protein` som hard-fail (`posebusters_critical:*`),
   - behold kun eksplisitte soft-typer som `posebusters_soft:*`,
   - ukjente PoseBusters-feil forblir konservativt hard-fail.

19. ~~**ProLIF interaction-type pruning**~~ — **AVKLART 2026-05-23.**
   Endelig utvalg i aktiv analyse: implicit H-bond og VdW-interaksjoner.
   Dette er implementert i koden.
