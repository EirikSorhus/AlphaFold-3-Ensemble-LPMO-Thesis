> **⚠️ ARCHIVED — LEGACY PSEUDOCODE**
> This file is a legacy pseudocode document. It predates the PLACER removal (2026-04-21) and the AF3-only simplification. It still references cross-model clustering and a PLACER validation step (STEP G) which are now obsolete. For the current pipeline specification, see [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) (primary) and [MASTERPLAN.md](MASTERPLAN.md).

# PSEUDOKODE - CLUSTER-FIRST ANALYSEPIPELINE

Status:
- Kun pseudokode. Ingen installasjon, ingen kjoreklar kode.
- Strukturprediksjon antas ferdig kjort separat.
- **ARKIVERT**: Refererer til PLACER og cross-model clustering som er fjernet. Se gjeldende plandokumenter.

Prioritet:
1. Siste kommentarer i chat
2. plan_implementation_spec.txt
3. plan_analyse.txt

Relaterte pseudokodefiler:
- PSEUDOKODE_MODUL_FOR_MODUL.md (detalj per modul)
- PSEUDOKODE_STATISTIKK.md (statistikkdel i separat flyt)

## 1) Hovedorkestrering

PROCEDURE KjorAnalyse(config):
    ASSERT config.predictions_exist == TRUE

    subanalyser <- [
        ("chitin", 4), ("chitin", 6), ("chitin", 8),
        ("cellulose", 4), ("cellulose", 6), ("cellulose", 8),
        ("starch", 4), ("starch", 6), ("starch", 8)
    ]

    FOR hver (substrate_class, DP) i subanalyser:
        KjorSubanalyse(config, substrate_class, DP)

    ByggSamlerapporter(config)

    IF config.optional_post_analysis_tuning == TRUE:
        KjorValgfriTuningEtterAnalyse(config)


## 2) Subanalyse (cluster er primarenhet)

PROCEDURE KjorSubanalyse(config, substrate_class, DP):
    pose_table <- tom tabell

    systems <- FinnSystemer(config, substrate_class, DP)

    FOR hvert system i systems:
        poses <- LastPredikertePoser(system)

        # STEP A: hard ingest/normalisering
        normalized_poses <- []
        FOR hver pose i poses:
            normalized <- NormaliserMMCIF(pose)
            IF normalized.atom_mapping_coverage < 1.0:
                LoggDropp(pose, "atom_mapping_lt_100")
                LagreNumeriskeMetrikker(pose, normalized.metrics)
                CONTINUE
            normalized_poses <- normalized_poses + [normalized]

        # STEP B: pre-QC active-site proximity
        preqc_pass <- []
        FOR hver pose i normalized_poses:
            proximity <- BeregnActiveSiteProximity(pose)
            LagreNumeriskeMetrikker(pose, proximity.metrics)

            IF proximity.min_cu_ligand_distance > config.pre_qc_active_site_max_a:
                LoggDropp(pose, "ligand_too_far_from_active_site")
                CONTINUE

            preqc_pass <- preqc_pass + [pose]

        # STEP C: hard QC
        qc_pass <- []
        FOR hver p in preqc_pass:
            pb <- KjorPoseBusters(p)
            priv <- KjorPrivateer(p)
            geom_gate <- KjorCuHisGate(p)

            LagreNumeriskeMetrikker(p, pb.metrics)
            LagreNumeriskeMetrikker(p, priv.metrics)
            LagreNumeriskeMetrikker(p, geom_gate.metrics)

            IF pb.critical_fail == TRUE:
                LoggDropp(p, "posebusters_critical")
                CONTINUE
            IF priv.severe_fail == TRUE:
                LoggDropp(p, "privateer_severe")
                CONTINUE
            IF geom_gate.pass == FALSE:
                LoggDropp(p, "cu_his_gate_fail")
                CONTINUE

            qc_pass <- qc_pass + [p]

        # STEP D: Branch A (IFP)
        FOR hver p in qc_pass:
            # Viktig: bruker AF3-poser direkte, ingen PLACER-raffinering
            # Ingen Cu-reposition, ingen virtuell oxyl/H i IFP-branch
            ifp <- KjorProLIF(p)
            p.ifp_vector <- ifp.vector

        # STEP E: Branch B (geometry)
        FOR hver p in qc_pass:
            g <- BeregnGeometriMedVirtuelleObjekter(p)
            # behold alle tall, ogsaa for borderline geometri
            p.geometry <- g

        # STEP F: Clustering (kun IFP)
        within_model_clusters <- ClusterWithinModel(qc_pass, metric="jaccard")
        final_clusters <- ClusterCrossModel(within_model_clusters, metric="jaccard")

        # STEP G: PLACER uavhengig validering (post-clustering)
        # PLACER-koordinater brukes IKKE videre. Kun skårer (prmsd, rmsd, plddt).
        FOR hver cluster i final_clusters:
            medoid <- cluster.medoid
            glycam_file <- HentGLYCAMLigandFil(system.substrate_type, system.DP)
            placer_scores <- KjorPLACERValidering(medoid.pdb, glycam_file, n_samples=50)
            cluster.placer_prmsd <- Median(placer_scores.prmsd)
            cluster.placer_rmsd <- Median(placer_scores.rmsd)
            cluster.placer_plddt <- Median(placer_scores.plddt)

        # STEP H: Cluster-annotering med geometri/QC/support/PLACER-skårer
        annotated_clusters <- []
        FOR hver cluster i final_clusters:
            ann <- AnnoterCluster(cluster, qc_pass)
            # bruk median/IQR for geometri
            # medoid brukes kun representasjon/visualisering
            annotated_clusters <- annotated_clusters + [ann]

        # STEP I: tabeller
        pose_table <- pose_table + ByggPoseRader(system, poses, normalized_poses, qc_pass)
        cluster_table <- ByggClusterTable(system, annotated_clusters)  # inkluderer placer-skårer
        predictive_cluster_table <- ByggPredictiveClusterTable(system, annotated_clusters)

        # STEP J: valgfri crystal anchoring
        IF FinnCrystal(system) == TRUE:
            crystal_anchor_rows <- KjorCrystalAnchoring(system, annotated_clusters)
        ELSE:
            crystal_anchor_rows <- []

        # STEP K: CBM sammenligning (for relevante enzymer)
        IF HarCBM(system.enzyme_id) == TRUE:
            cbm_rows <- KjorCBMParanalyse(system, annotated_clusters)
        ELSE:
            cbm_rows <- []

        LagreOutputPerSystem(
            system,
            pose_table,
            cluster_table,
            predictive_cluster_table,
            crystal_anchor_rows,
            cbm_rows
        )

    # Deskriptiv/prediktiv statistikk foretrekkes i R
    KjorRDeskriptivStatistikk(substrate_class, DP)
    KjorRPrediktiveModeller(substrate_class, DP)


## 3) Cluster-basert prediktiv analyse (hovedregel)

PROCEDURE KjorRPrediktiveModeller(substrate_class, DP):
    data <- LesTSV("predictive_cluster_table.tsv")

    # HOVEDREGEL:
    # - En rad per cluster
    # - Ikke aggreger til en rad per enzym i hovedmodell

    folds <- LagGroupedCVFolds(data, group_col="enzyme_id", stratify_cols=["experimental_regio_label", "substrate_class"])

    model <- FitPenalizedLogistic(
        data,
        target="experimental_regio_label",
        features=[
            "occupancy",
            "cluster_weight",
            "median_Cu_C1",
            "median_Cu_C4",
            "median_oxyl_H_C1",
            "median_oxyl_H_C4",
            "ifp_features_selected",
            "cross_model_support",
            "convergence_support",
            "family",
            "cbm_status"
        ],
        grouped_folds=folds
    )

    Evaluer(model, metrics=["balanced_accuracy", "macro_f1"])
    LagreModelRapport(model)


## 4) Enzymniva er kun sekundar/sensitivitet

PROCEDURE KjorSekundarEnzymOppsummering():
    cluster_data <- LesTSV("predictive_cluster_table.tsv")

    enzyme_summary <- AggregerMedOccupancyVekting(cluster_data)
    LagreTSV(enzyme_summary, "enzyme_summary_table.tsv")

    # Brukes kun i sensitivitetstester, ikke hovedmodell


## 5) EC -> aktivitet mapping (pseudokode)

FUNCTION MapECToActivity(ec_number, family):
    IF ec_number == "1.14.99.54":
        RETURN ("cellulose", "C1", "cellulose_C1_hydroxylating")

    IF ec_number == "1.14.99.56":
        RETURN ("cellulose", "C4", "cellulose_C4_dehydrogenating")

    IF ec_number == "1.14.99.53":
        RETURN ("chitin", "mixed", "chitin_C1_C4_mixed")

    IF ec_number == "1.14.99.55":
        RETURN ("starch", "C1", "starch_C1_hydroxylating")

    IF ec_number == "1.14.99.-":
        IF family inneholder "AA17":
            RETURN ("homogalacturonan", "C4", "homogalacturonan_C4_oxidation")
        ELSE:
            RETURN ("xylan_or_other", "unknown", "xylan_like_oxidative")

    RETURN ("unknown", "unknown", "unknown")


## 6) Valgfri tuning etter analyse

PROCEDURE KjorValgfriTuningEtterAnalyse(config):
    # Kjor kun hvis tid/ressurser
    # Ikke blokkering for baseline analyse

    FOR model_platform i ["AF3", "RF3", "Boltz2"]:
        tuning_jobs <- LagTuningJobber(model_platform, config)
        tuning_results <- KjorTuningJobber(tuning_jobs)
        EvaluerTuning(tuning_results)

    LagreTuningSammendrag()


## 7) Ikke-slett regel for numeriske metrikker

PROCEDURE LagreNumeriskeMetrikker(pose, metrics):
    # Dette kjores alltid, ogsaa ved hard fail
    # Ingen beregnede tall slettes fra logg/tabeller
    AppendTilPoseMetricsLog(pose.id, metrics)
