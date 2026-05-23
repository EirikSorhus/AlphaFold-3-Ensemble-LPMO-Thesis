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

3. **Reduce-versjon** — MolProbity Reduce vs AmberTools reduce?
   *Default: AmberTools `reduce` (mest tilgjengelig via conda).*

4. ~~**Pre-QC aktiv-sete terskel**~~ — **AVKLART 2026-04-21.**
   Hard cutoff for `min_cu_ligand_distance` foran PoseBusters/Privateer er
   `<= 10.0 A` (permissiv pre-QC gate). `<= 8.0 A` beholdes kun som mulig
   rapporterings-/soft-flag terskel ved behov.

5. **Kristallstrukturer for anchoring** — Hvilke PDB-koder skal brukes
   som referanse på tvers av familier? Dagens operative referansesett kommer
   fra `input_data/pdb_structure_data.csv` + filer under `crystal_structures/`.
   Gjenstående spørsmål er om listen skal kurateres/utvides videre per familie,
   ikke hvordan dagens kode velger referanser.
   *Default: bruk dagens CSV + `crystal_structures/` som autoritativt
   referansesett. Eventuell utvidelse skjer som eksplisitt dataoppdatering.*

6. **CBM-varianter (DEL A / DEL B)** — Skal begge CBM-deletions kjøres
   for alle systemer, eller bare for CBM-bærende LPMOer?
   *Default: kun for systemer der full-length har annotert CBM.*

7. **DP-scope for aktiv implementasjon** — Skal andre DP enn 4/6/8 inn i
   samme hovedpipeline, eller beholdes de som egne sideanalyser?
   *Default: hovedpipeline = DP4/DP6/DP8 (9 uavhengige delanalyser).* 

8. **Tanimoto-terskel for crystal anchoring** — Hvilken cutoff for
   "biologically plausible"? *Default: Tanimoto >= 0.3 (IFP), og dagens kode
   bruker `pocket_rmsd < 2.5 A` som soft flag; endelig biologisk cutoff må
   fortsatt signeres eksplisitt.*
   Status 2026-05-08: dagens integrerte analysis-core smoke-tester fullforer,
   men den testede crystal-anchoring-kjoringen ender forelopig med
   `comparison_count = 0`. Det ser ogsa ut som de testede crystal-IFP-ene kan
   vaere dominert av VdW-interaksjoner. Etter senere crystal-prep-hardening og
   standalone real-data-validering regnes dette ikke lenger som et prima facie
   teknisk prep-problem; hvis crystal-IFP i praksis blir tom eller
   non-comparable er den mest sannsynlige forklaringen for tiden for svak eller
   for lite spesifikk biologisk-strukturell kontakt under den delte non-vdW
   contact-eligibility-regelen. Dette skal fortsatt kvantifiseres bredere, men
   skal ikke lenger default-tolkes som kjent Cu-/normaliseringskontaminasjon.
   Status 2026-05-16: standalone real-data-harnessen for `A0A0S2GKZ1` gir nå
   ikke-tomme sammenligninger mot `5ACI` og `7PXW` og lave pocket-RMSD-er.
   Det som fortsatt mangler er en integrert produksjonskjøring der en faktisk
   medoid-backed betingelse med crystal-referanser når crystal-anchoring-steget,
   samt en eksplisitt vurdering av om crystal-IFP-ene blir for VdW-dominerte.
   Status 2026-05-20: produksjonskoden sammenligner nå alle beholdte medoids,
   bruker top-level AF3 model CIF som hard-QC-gated fallback for no-cluster
   conditions, blokkerer IFP-Tanimoto når crystal-IFP ikke passerer non-vdW
   contact-eligibility, og skriver `crystal_ifp_diagnostic_summary.tsv` for å
   kvantifisere VdW-/low-specific-contact-problemet. Spørsmålet som gjenstår er
   biologisk cutoff/tolkning, ikke selve comparability-gaten. Arbeidshypotesen
   er nå at crystal-IFP som fortsatt blir tom eller non-comparable oftest
   reflekterer reell biologisk-strukturell kontaktsparsomhet i de preparerte
   deposited ligandene, ikke en kjent implementasjonsfeil i crystal-IFP-laget.

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

10. **EC 1.14.99.- ikke-AA17 mapping** — Hvilken endelig tekstetikett og
    hvilket standardsubstrat for "xylan ol"-tilfeller?
    *Default: substrate_class=`xylan_or_other`, regio_class=`unknown`, aktivitet=`xylan_like_oxidative`.*

11. **Geometri-planaritet** — Avklart 2026-05-22: separat planaritetsgate er
    ikke del av aktiv analysekontrakt. Sluttbrukerflaten bruker de implementerte
    downstream-feltene `sugar_face_orientation`, `ring_normal_vs_brace_normal`,
    `attack_angle_C1/C4`, `oxyl_H_*_distance` og `geometry_status_C1/C4`.
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

13. **Substrat-recognition/pocket-residuer for crystal anchoring** — Skal dagens
   operative pocket-heuristikk beholdes, eller erstattes/utvides med
   familie-spesifikke litteraturresiduer?
   Dagens implementasjon bruker proteinrester innen 5 A fra ligand eller Cu i
   holo-referanser. For apo-referanser projiseres pocket fra representant/medoid
   over på crystal-sekvensen med residunavn-normalisering (f.eks. `HIC -> HIS`).
   Foretrukket videre arbeid: litteratursøk for kjente LPMO-substrat-bindende
   residuer per familie, og beslutning om disse skal erstatte eller bare annotere
   dagens proximity-baserte pocket.
   *Status: dagens heuristic fungerer operativt, men familie-spesifikk
   litteraturforankring er fortsatt uavklart og krever manuelt arbeid.*

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
   Viktig presisering: i PoseBusters `dock`-mode er
   `protein-ligand_maximum_distance` far-away-testen fra
   `posebusters.modules.intermolecular_distance.check_intermolecular_distance()`
   med standard `max_distance=5.0 A` og `search_distance=6.0 A`.
   `minimum_distance_to_protein` er derimot den renamed `no_clashes`-utgangen,
   ikke selve avstandsterskelen.
   *Default: behold PoseBusters sine innebygde `dock`-defaults uendret.
   Gjenstående avklaring er bare severity-policyen for pocket-/distance-relaterte
   PoseBusters-feil: om de skal gi `hard_fail`, `soft_flag` eller kun rapporteres.*

19. **ProLIF interaction-type pruning** — Skal endelig utvalg av ProLIF-interaksjonstyper
   bestemmes bare ved enkel sparsity-/nyttevurdering, eller er det verdt å lage en
   form for statistisk analyse/test for dette?
   *Default: ikke blokker på dette. Behold enkel manuell vurdering senere som første
   steg. Eventuell statistisk analyse er lav prioritet, ikke spesielt viktig akkurat nå,
   og kan ta tid å designe og teste på en meningsfull måte.*
