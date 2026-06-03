> **⚠️ ARCHIVED — LEGACY PSEUDOCODE**
> This file is a legacy module-level pseudocode document. It contains a PLACER package section (section 4) which is now obsolete (PLACER removed 2026-04-21), and cross-model clustering logic that is not applicable for the AF3-only pipeline. For the current specification, see [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) (primary) and [MASTERPLAN.md](MASTERPLAN.md).

# PSEUDOKODE - MODUL FOR MODUL

Status:
- Kun pseudokode.
- Ingen installasjon.
- Ingen kjoreklar implementasjon.
- **ARKIVERT**: Seksjon 4 (PLACER) og cross-model clustering er foreldet. Se gjeldende plandokumenter.

Designregler som gjelder alle moduler:
- Strukturprediksjon er separat steg utenfor analysemodulene.
- Cluster er primarenhet i hovedanalyse.
- Pre-QC active-site proximity skal kjores for PoseBusters/Privateer.
- Beregnede numeriske metrikker skal ikke slettes, selv ved hard fail.
- Tuning er valgfri etteranalyse.

## 0) Top-level CLI

### MODULE: src/lpmo_pipeline/cli.py

PROCEDURE MainCLI(args):
    Parse subcommand
    IF subcommand == run:
        CmdRun(args)
    ELSE IF subcommand == tune:
        CmdTune(args)
    ELSE:
        PrintHelp()

PROCEDURE CmdRun(args):
    config <- LoadConfig(args.config)
    manifest <- InitManifest(mode="production")
    KjorAnalysePipeline(config, args.output)
    WriteManifest(manifest)

PROCEDURE CmdTune(args):
    config <- LoadConfig(args.config)
    manifest <- InitManifest(mode="tune_optional")
    KjorValgfriTuningEtterAnalyse(config, args.output)
    WriteManifest(manifest)

## 1) IO package

### MODULE: src/lpmo_pipeline/io/mmcif_ingest.py

PROCEDURE IngestMMCIF(path):
    structure <- ParseMMCIF(path)
    ValidateRequiredCategories(structure)
    RETURN {
        structure,
        ingest_ok,
        ingest_metrics
    }

### MODULE: src/lpmo_pipeline/io/normalize_mmcif.py

PROCEDURE NormalizeMMCIF(structure, config):
    AssignChainSchema(structure, protein="A", glycans=["B","C","D"], metal="E")
    EnsureChemCompBondComplete(structure)
    EnsureStructConnForGlycosidicAndCu(structure)
    confidence <- ExtractConfidence(structure)
    RETURN {
        normalized_structure,
        confidence,
        normalization_metrics
    }

### MODULE: src/lpmo_pipeline/io/ccd_lookup.py

FUNCTION ValidateCCDMonosaccharides(glycan_residues, whitelist):
    recognized <- []
    unrecognized <- []
    FOR residue IN glycan_residues:
        IF residue.comp_id IN whitelist:
            recognized <- recognized + [residue]
        ELSE:
            unrecognized <- unrecognized + [residue]
    RETURN {
        recognized,
        unrecognized,
        recognition_rate
    }

### MODULE: src/lpmo_pipeline/io/analysis_export.py

PROCEDURE ExportAnalysisArtifacts(normalized_structure):
    complex_pdb <- ExportNonProtonatedComplexPDB(normalized_structure)
    ligand_pdb <- ExportNonProtonatedLigandOnlyPDB(normalized_structure)
    RETURN {
        complex_for_prolif_pdb,
        ligand_only_for_prolif_pdb,
        analysis_export_report
    }

## 2) Mapping package

### MODULE: src/lpmo_pipeline/mapping/cross_model_atom_mapping.py

PROCEDURE BuildCrossModelAtomMap(reference_structure, target_structure):
    FOR atom_ref IN reference_structure:
        candidates <- FindCandidatesBy(
            element,
            residue_ccd,
            local_bond_graph,
            three_d_proximity
        )
        best <- ResolveBestCandidate(candidates)
        AppendMapping(atom_ref, best)
    coverage <- ComputeMappingCoverage()
    RETURN {
        atom_map,
        coverage,
        mapping_log
    }

### MODULE: src/lpmo_pipeline/mapping/rename_atoms.py

PROCEDURE ApplyAtomMap(structure, atom_map):
    RenameAtoms(structure, atom_map)
    ValidateNoUnmappedAtoms()
    RETURN {
        renamed_structure,
        rename_log
    }

## 3) QC package

### MODULE: src/lpmo_pipeline/qc/active_site_proximity.py

PROCEDURE CheckActiveSiteProximity(refined_pose, threshold):
    cu <- FindCu(refined_pose)
    min_cu_ligand <- ComputeMinDistance(cu, ligand_heavy_atoms)
    min_cu_c1 <- ComputeMinDistance(cu, candidate_C1_atoms)
    min_cu_c4 <- ComputeMinDistance(cu, candidate_C4_atoms)

    pass <- (min_cu_ligand <= threshold)

    RETURN {
        pass,
        min_cu_ligand,
        min_cu_c1,
        min_cu_c4,
        nearest_atoms,
        warnings
    }

### MODULE: src/lpmo_pipeline/qc/posebusters_runner.py

PROCEDURE RunPoseBusters(pose):
    raw <- ExecutePoseBusters(pose)
    classified <- ClassifyCriticalVsSoft(raw)
    RETURN {
        critical_fail,
        soft_flags,
        pb_metrics
    }

### MODULE: src/lpmo_pipeline/qc/privateer_runner.py

PROCEDURE RunPrivateer(privateer_input_pose):
    raw <- ExecutePrivateer(privateer_input_pose)
    severe <- DetectSeverePrivateerFailures(raw)
    RETURN {
        severe_fail,
        privateer_metrics,
        per_residue_privateer
    }

### MODULE: src/lpmo_pipeline/qc/custom_geometry_checks.py

PROCEDURE RunCuHisGate(pose):
    cu_his <- MeasureCuHisDistances(pose)
    pass <- AllBetween(cu_his, 1.9, 2.6)
    RETURN {
        pass,
        cu_his,
        geometry_gate_metrics
    }

### MODULE: src/lpmo_pipeline/qc/gates.py

PROCEDURE EvaluateHardGates(gate_inputs):
    hard_results <- []
    hard_results <- hard_results + CheckAtomMappingCoverage(gate_inputs.mapping_coverage)
    hard_results <- hard_results + CheckActiveSiteProximity(gate_inputs.min_cu_ligand)
    hard_results <- hard_results + CheckPrivateerRecognition(gate_inputs.privateer_rate)
    hard_results <- hard_results + CheckPoseBustersCritical(gate_inputs.pb_critical)
    hard_results <- hard_results + CheckCuHisRange(gate_inputs.cu_his)
    RETURN hard_results

### MODULE: src/lpmo_pipeline/qc/qc_report.py

PROCEDURE BuildQCVerdict(pose_id, pre_qc, pb, privateer, cu_his):
    metrics <- MergeMetrics(pre_qc, pb, privateer, cu_his)
    IF pre_qc.pass == FALSE:
        status <- "dropped"
    ELSE IF pb.critical_fail == TRUE:
        status <- "dropped"
    ELSE IF privateer.severe_fail == TRUE:
        status <- "dropped"
    ELSE IF cu_his.pass == FALSE:
        status <- "dropped"
    ELSE IF AnySoftFlags(pre_qc, pb, privateer, cu_his):
        status <- "flagged"
    ELSE:
        status <- "passed"

    RETURN {
        pose_id,
        status,
        drop_reasons,
        warnings,
        metrics
    }

## 4) PLACER package (Independent Validation, Post-Clustering)

### MODULE: src/lpmo_pipeline/placer/run_placer.py

PROCEDURE RunPLACERValidation(cluster_medoids, placer_config, glycam_ligand_files):
    # PLACER brukes som uavhengig validering, IKKE refinement.
    # Koordinater fra PLACER brukes IKKE videre i analysen.
    # Kun skårer (prmsd, rmsd, plddt) ekstraheres.
    validation_results <- []
    FOR medoid IN cluster_medoids:
        ligand_file <- glycam_ligand_files[medoid.substrate_type]
        placer_output <- ExecutePLACER(medoid.pdb, ligand_file, n_samples=50)
        scores <- ExtractScores(placer_output)  # prmsd, rmsd, plddt
        validation_results <- validation_results + [{medoid.cluster_id, scores}]
    RETURN validation_results

### MODULE: src/lpmo_pipeline/placer/rank_ensemble.py

# NOTE: rank_ensemble er ikke lenger brukt for å rangere poser i hovedanalysen.
# Beholdes kun for eventuell intern PLACER-analyse av validerings-ensemblet.

PROCEDURE SummarizePLACERValidation(validation_results):
    FOR result IN validation_results:
        median_prmsd <- Median(result.scores.prmsd)
        median_rmsd <- Median(result.scores.rmsd)
        median_plddt <- Median(result.scores.plddt)
    RETURN placer_summary_per_cluster

## 5) Analysis package

### MODULE: src/lpmo_pipeline/analysis/prolif_ifp.py

PROCEDURE ComputeIFP(post_qc_pose):
    # Viktig: bruker AF3-poser direkte, ikke PLACER-raffinerte
    # Ingen engineered geometry objekter
    ifp <- RunProLIF(post_qc_pose)
    RETURN {
        ifp_vector,
        ifp_feature_names,
        ifp_status
    }

### MODULE: src/lpmo_pipeline/analysis/mdanalysis_metrics.py

PROCEDURE ComputeGeometryBranch(post_qc_pose):
    his_brace <- IdentifyHisBrace(post_qc_pose)
    brace_plane <- BuildBracePlane(his_brace)
    cu_repositioned <- RepositionCuForGeometryOnly(post_qc_pose, brace_plane)
    virtual_oxyl <- PlaceVirtualOxyl(cu_repositioned, brace_plane)
    virtual_h_c1 <- PlaceVirtualHAtC1(post_qc_pose)
    virtual_h_c4 <- PlaceVirtualHAtC4(post_qc_pose)

    metrics <- {
        oxyl_H_C1_distance,
        oxyl_H_C4_distance,
        Cu_C1_distance,
        Cu_C4_distance,
        Cu_oxyl_H_C1_angle,
        Cu_oxyl_H_C4_angle,
        sugar_face_orientation,
        optional_ring_normal_vs_brace_normal
    }
    RETURN metrics

### MODULE: src/lpmo_pipeline/analysis/clustering_hdbscan.py

PROCEDURE ClusterWithinModel(pose_rows):
    groups <- GroupBy(pose_rows, keys=[enzyme_id, substrate_class, DP, model_platform])
    FOR g IN groups:
        g.clusters <- RunHDBSCAN(g.ifp_vectors, metric="jaccard")
    RETURN groups

PROCEDURE ClusterCrossModel(within_model_result):
    systems <- GroupBy(within_model_result, keys=[enzyme_id, substrate_class, DP])
    FOR s IN systems:
        reps <- SelectRepresentatives(s)
        s.final_clusters <- RunHDBSCAN(reps.ifp_vectors, metric="jaccard")
    RETURN systems

### MODULE: src/lpmo_pipeline/analysis/cluster_signatures.py

PROCEDURE BuildClusterSignatures(final_clusters, pose_rows):
    signatures <- []
    FOR cluster IN final_clusters:
        pose_subset <- FilterByCluster(pose_rows, cluster.id)
        summary <- {
            occupancy,
            cluster_size,
            medoid_pose_id,
            model_platform_support,
            convergence_support,
            qc_pass_fraction,
            invalid_geometry_fraction,
            plausible_geometry_fraction,
            scorable_geometry_fraction,
            median_oxyl_H_C1,
            IQR_oxyl_H_C1,
            median_oxyl_H_C4,
            IQR_oxyl_H_C4,
            median_Cu_C1,
            IQR_Cu_C1,
            median_Cu_C4,
            IQR_Cu_C4,
            median_Cu_oxyl_H_C1_angle,
            median_Cu_oxyl_H_C4_angle,
            median_face_orientation,
            ifp_enrichment_features
        }
        signatures <- signatures + [summary]
    RETURN signatures

### MODULE: src/lpmo_pipeline/analysis/activity_mapping.py

FUNCTION MapECToActivity(ec_number, family):
    IF ec_number == "1.14.99.54": RETURN ("cellulose","C1","cellulose_C1_hydroxylating")
    IF ec_number == "1.14.99.56": RETURN ("cellulose","C4","cellulose_C4_dehydrogenating")
    IF ec_number == "1.14.99.53": RETURN ("chitin","mixed","chitin_C1_C4_mixed")
    IF ec_number == "1.14.99.55": RETURN ("starch","C1","starch_C1_hydroxylating")
    IF ec_number == "1.14.99.-" AND Contains(family, "AA17"):
        RETURN ("homogalacturonan","C4","homogalacturonan_C4_oxidation")
    IF ec_number == "1.14.99.-":
        RETURN ("xylan_or_other","unknown","xylan_like_oxidative")
    RETURN ("unknown","unknown","unknown")

PROCEDURE BuildPredictiveClusterRows(cluster_signatures, metadata):
    rows <- []
    FOR c IN cluster_signatures:
        labels <- LookupLabels(metadata, c.enzyme_id)
        row <- {
            analysis_id,
            enzyme_id,
            family,
            cbm_status,
            substrate_class,
            DP,
            cluster_id,
            occupancy,
            qc_factor,
            support_factor,
            cluster_weight = occupancy * qc_factor * support_factor,
            cross_model_support,
            convergence_support,
            invalid_geometry_fraction,
            plausible_geometry_fraction,
            scorable_geometry_fraction,
            median_oxyl_H_C1,
            median_oxyl_H_C4,
            median_Cu_C1,
            median_Cu_C4,
            median_Cu_oxyl_H_C1_angle,
            median_Cu_oxyl_H_C4_angle,
            face_orientation_summary,
            ifp_features_selected,
            experimental_regio_label,
            experimental_ligand_specificity_label
        }
        rows <- rows + [row]
    RETURN rows

PROCEDURE BuildSecondaryEnzymeSummary(predictive_cluster_rows):
    # Kun sensitivitet, ikke hovedanalyse
    summary <- OccupancyWeightedAggregationByEnzyme(predictive_cluster_rows)
    RETURN summary

### MODULE: src/lpmo_pipeline/analysis/predictive_models.py

PROCEDURE RunPrimaryClusterModels(predictive_cluster_rows):
    folds <- GroupedCVByEnzyme(predictive_cluster_rows)
    model <- FitPrimaryModelClusterRows(predictive_cluster_rows, folds)
    metrics <- Evaluate(model, [balanced_accuracy, macro_f1])
    RETURN {model, metrics}

### MODULE: src/lpmo_pipeline/analysis/crystal_anchoring.py

PROCEDURE RunCrystalAnchoringIfAvailable(system, cluster_signatures):
    IF CrystalAvailable(system) == FALSE:
        RETURN EmptyTable()

    # 1) Prepare crystal reference in the same canonical chain schema as AF3
    selected_reference <- SelectProteinLigandCuSite(system.crystal)
    selected_reference_cif <- WriteSelectedReferenceCif(
        selected_reference,
        keep_metadata = [_entity, _chem_comp, _chem_comp_bond, _struct_conn],
        chain_schema = {protein: A, glycan: B/C/D, cu: E}
    )
    IF selected_reference.has_ligand:
        normalized_reference_cif <- NormalizeMMCIF(selected_reference_cif)
        crystal_artifacts <- ExportAnalysisArtifacts(normalized_reference_cif)
        # complex_for_prolif keeps Cu for geometry/RMSD; ligand-only PDB is glycan-only for ProLIF.
        crystal_ifp <- ComputeProLIF(
            crystal_artifacts.complex_for_prolif_pdb,
            crystal_artifacts.ligand_only_for_prolif_pdb
        )

    # 2) Identify alignment atoms
    his_brace_residues <- FindHistidineBrace(system.predicted)
    cu_near_residues <- FindResiduesNearCu(system.predicted, cutoff=3.5)

    # 3) Substrate-recognition residues
    IF LiteratureResiduesAvailable(system.enzyme_family):
        surface_residues <- LoadLiteratureResidues(system.enzyme_family)
    ELSE:
        surface_residues <- IdentifyPocketResiduesByProximity(
            system.predicted, ligand_chain, protein_chain, cutoff=5.0
        )
        LogWarning("Using proximity-based residues, not literature-based")

    alignment_residues <- Union(his_brace_residues, cu_near_residues, surface_residues)

    # 4) Current implementation: local Kabsch superposition on shared pocket C-alpha atoms
    rmsd_after_alignment <- LocalKabschPocketRMSD(
        system.predicted, selected_reference, alignment_residues
    )

    # 5) Metrics after optimized alignment
    metrics <- {
        pocket_rmsd: rmsd_after_alignment,
        alignment_type: "optimized_local_kabsch",
        alignment_residue_count: Count(alignment_residues),
        residue_source: "literature" OR "proximity",
        ligand_or_proximal_segment_rmsd,
        ifp_similarity_if_comparable,
        per_residue_deviations
    }
    RETURN metrics

### MODULE: src/lpmo_pipeline/analysis/cbm_variant.py

PROCEDURE BuildCBMVariants(system):
    domain_only <- PrepareDomainOnly(system)
    full_length <- PrepareFullLength(system)
    RETURN {domain_only, full_length}

### MODULE: src/lpmo_pipeline/analysis/cbm_comparison.py

PROCEDURE CompareCBMPaired(domain_rows, full_rows):
    paired <- PairByEnzymeAndSubanalysis(domain_rows, full_rows)
    comparison <- {
        delta_C1_occupancy,
        delta_C4_occupancy,
        delta_convergence,
        delta_invalid_geometry_fraction,
        delta_selected_ifp_features,
        cbm_ligand_distance,
        cbm_active_site_distance
    }
    RETURN comparison

## 6) Report package

### MODULE: src/lpmo_pipeline/report/build_metrics_csv.py

PROCEDURE BuildMetricsCSV(all_pose_rows, all_cluster_rows):
    flat <- FlattenToMetricsRows(all_pose_rows, all_cluster_rows)
    WriteCSV(flat, "metrics.csv")

### MODULE: src/lpmo_pipeline/report/build_summary_json.py

PROCEDURE BuildSummaryJSON(all_outputs):
    summary <- {
        dataset_stats,
        qc_stats,
        cluster_stats,
        geometry_stats,
        modeling_stats,
        failure_counts
    }
    WriteJSON(summary, "summary.json")

### MODULE: src/lpmo_pipeline/report/build_report_html.py

PROCEDURE BuildHTMLReport(all_outputs):
    figures <- BuildFigures(all_outputs)
    tables <- BuildTables(all_outputs)
    html <- RenderTemplate(figures, tables)
    WriteFile(html, "report.html")

## 7) Utils package

### MODULE: src/lpmo_pipeline/utils/data_models.py

PROCEDURE DefineCanonicalSchemas():
    DefinePoseRowSchema()
    DefineClusterRowSchema()
    DefinePredictiveClusterRowSchema()
    DefineEnzymeSummarySchema()

### MODULE: src/lpmo_pipeline/utils/exceptions.py

PROCEDURE DefinePipelineExceptions():
    DefineHardGateExceptionTypes()
    DefineRecoverableWarningTypes()

### MODULE: src/lpmo_pipeline/utils/hashing.py

FUNCTION ComputeInputHashes(paths):
    RETURN sha256_per_input

### MODULE: src/lpmo_pipeline/utils/logging.py

PROCEDURE InitStructuredLogging(run_id):
    SetupJsonAndConsoleLogs(run_id)

### MODULE: src/lpmo_pipeline/utils/manifest.py

PROCEDURE BuildManifest(config, tools, inputs, gates):
    manifest <- {
        pipeline_version,
        timestamp,
        config_hash,
        tool_versions,
        input_hashes,
        gate_outcomes,
        prediction_artifact_refs
    }
    WriteJSON(manifest, "run_manifest.json")

### MODULE: src/lpmo_pipeline/utils/parallel.py

PROCEDURE RunParallelJobs(job_list, n_jobs):
    SplitJobs(job_list, n_jobs)
    CollectResultsWithFailureMetadata()

### MODULE: src/lpmo_pipeline/utils/paths.py

FUNCTION ResolveRunPaths(base_output, analysis_id):
    RETURN {
        raw_predictions_dir,
        qc_dir,
        ifp_dir,
        geometry_dir,
        clustering_dir,
        report_dir
    }

## 8) Tuning package (valgfri etteranalyse)

### MODULE: src/lpmo_pipeline/tuning/parameter_grid.py

PROCEDURE BuildParameterGrid(model_platform, config):
    RETURN list_of_parameter_sets

### MODULE: src/lpmo_pipeline/tuning/sweep_runner.py

PROCEDURE RunSweepJobs(parameter_sets, tuning_subset):
    FOR params IN parameter_sets:
        result <- RunComparableMiniAnalysis(tuning_subset, params)
        StoreSweepResult(result)
    RETURN all_sweep_results

### MODULE: src/lpmo_pipeline/tuning/tune_orchestrator.py

PROCEDURE RunOptionalPostAnalysisTuning(config):
    FOR model_platform IN ["AF3","RF3","Boltz2"]:
        grid <- BuildParameterGrid(model_platform, config)
        sweep_results <- RunSweepJobs(grid, config.tuning_subset)
        SummarizePlatform(sweep_results)
    RETURN tuning_summary

### MODULE: src/lpmo_pipeline/tuning/tune_summarize.py

PROCEDURE SummarizeTuningResults(all_platform_results):
    rank <- RankBy(
        qc_pass_rate,
        outlier_rate_inverse,
        cluster_stability,
        crystal_sanity
    )
    WriteSummary(rank)

## 9) End-to-end pseudoflyt

PROCEDURE KjorAnalysePipeline(config, output_dir):
    InitRun(output_dir)

    FOR each subanalysis in 9 independent analyses:
        KjorSubanalyseClusterFirst(subanalysis)

    KjorStatistikkDelIR(output_dir)
    ByggSamletRapport(output_dir)

    IF config.optional_post_analysis_tuning == TRUE:
        KjorValgfriTuningEtterAnalyse(config)
