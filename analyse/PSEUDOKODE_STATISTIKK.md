> **⚠️ ARCHIVED — LEGACY PSEUDOCODE**
> This file is legacy statistics pseudocode. It references `pose_table.tsv` as a single merged table (now replaced by five separate pose-level tables) and `cross_model_support` (not applicable for AF3-only pipeline). Core statistical methodology is still conceptually valid. For the current specification, see [AF3_LPMO_pipeline_detailed_plan.md](AF3_LPMO_pipeline_detailed_plan.md) stages 13–14 and [MASTERPLAN.md](MASTERPLAN.md) stages 13–14.

# PSEUDOKODE - STATISTIKKDEL (R-ORIENTERT)

Status:
- Kun pseudokode for statistikk.
- Ingen kjoreklar kode.
- **ARKIVERT**: Refererer til `pose_table.tsv` (enkeltfil) og `cross_model_support` som er foreldet. Se gjeldende plandokumenter.

Hovedprinsipp:
- Cluster-rader er primardata i hovedmodeller.
- Enzym-oppsummering er kun sensitivitet.
- Grouped CV pa enzymniva.

## 1) Inputs

INPUT FILES:
- predictive_cluster_table.tsv
- cluster_table.tsv
- pose_table.tsv
- cbm_comparison_table.tsv (valgfri)
- crystal_anchor_table.tsv (valgfri)
- metadata med labels (regio, familie, CBM, EC-avledet aktivitet)

FOR each subanalysis in {
    chitin_DP4, chitin_DP6, chitin_DP8,
    cellulose_DP4, cellulose_DP6, cellulose_DP8,
    starch_DP4, starch_DP6, starch_DP8
}:
    RunStatistikkSubanalyse(subanalysis)

## 2) Datapreparering

PROCEDURE RunStatistikkSubanalyse(subanalysis):
    cluster_data <- LoadTSV("predictive_cluster_table.tsv", filter=subanalysis)
    pose_data <- LoadTSV("pose_table.tsv", filter=subanalysis)

    ASSERT OneRowPerCluster(cluster_data)
    ASSERT RequiredColumnsExist(cluster_data)

    cluster_data <- AddDerivedFields(cluster_data):
        - occupancy_capped
        - cluster_weight_norm
        - support_bin
        - family_factor
        - cbm_factor

    # Ikke slett numeriske felt
    cluster_data <- PreserveAllNumericColumns(cluster_data)

    SaveQCDataSnapshot(cluster_data)

## 3) Deskriptiv statistikk (cluster-niva)

PROCEDURE KjorDeskriptiv(cluster_data):
    desc_counts <- ComputeCountsAndProportions(
        by=[enzyme_id, family, substrate_class, DP]
    )

    desc_geometry <- ComputeMedianIQR(
        columns=[
            median_oxyl_H_C1,
            median_oxyl_H_C4,
            median_Cu_C1,
            median_Cu_C4,
            median_Cu_oxyl_H_C1_angle,
            median_Cu_oxyl_H_C4_angle
        ],
        by=[cluster_type_annotation, family, cbm_status]
    )

    desc_ifp <- ComputeIFPFeatureEnrichment(cluster_data)

    tests_cat <- RunCategoricalTests(
        method=[fisher_or_chi_square],
        adjust="BH_FDR"
    )

    ExportDescriptiveTables(desc_counts, desc_geometry, desc_ifp, tests_cat)
    ExportDescriptivePlots()

## 4) Primarmodell: regioselektivitet (cluster-rader)

PROCEDURE FitPrimaryRegioModel(cluster_data):
    # Target kan vaere binar C1/C4 eller multinomial C1/C4/mixed
    data <- FilterKnownRegioLabels(cluster_data)

    folds <- BuildGroupedStratifiedFolds(
        data,
        group_col="enzyme_id",
        stratify_cols=["experimental_regio_label", "family"]
    )

    features <- [
        occupancy,
        cluster_weight,
        median_oxyl_H_C1,
        median_oxyl_H_C4,
        median_Cu_C1,
        median_Cu_C4,
        median_Cu_oxyl_H_C1_angle,
        median_Cu_oxyl_H_C4_angle,
        convergence_support,
        cross_model_support,
        plausible_geometry_fraction,
        ifp_features_selected,
        family,
        cbm_status
    ]

    model <- FitPenalizedLogistic(
        data,
        target="experimental_regio_label",
        features=features,
        grouped_folds=folds,
        class_weights="balanced"
    )

    perf <- EvaluateCV(
        model,
        metrics=[balanced_accuracy, macro_f1]
    )

    importance <- ExtractFeatureImportance(model)

    SavePrimaryModelResults(perf, importance)

## 5) Alternativmodeller (sekundar)

PROCEDURE FitAlternativeModels(cluster_data):
    model_mixed <- FitMixedEffectsLogistic(
        formula="regio ~ features + family + cbm_status + (1|enzyme_id)",
        data=cluster_data
    )

    model_multinomial <- FitMultinomialIfMixedClassPresent(cluster_data)

    model_rf <- FitRandomForestSensitivity(cluster_data)

    CompareWithPrimary(model_mixed, model_multinomial, model_rf)
    SaveAlternativeModelResults()

## 6) Ligandspesifisitet pa tvers av subanalyser

PROCEDURE RunLigandSpecificityAcrossSubanalyses(all_subanalysis_results):
    # Ikke bland IFP-space ved modelltrening;
    # sammenligning skjer etter at hver subanalyse er ferdig analysert.

    per_enzyme_profiles <- BuildPerEnzymeProfiles(all_subanalysis_results):
        - plausible_cluster_existence
        - dominant_cluster_signatures
        - best_geometry_signatures

    paired_tests <- RunPairedComparisons(per_enzyme_profiles)
    mixed_models <- RunRepeatedMeasuresModels(per_enzyme_profiles)

    SaveLigandSpecificityResults(paired_tests, mixed_models)

## 7) CBM-analyse (parret)

PROCEDURE RunCBMStats(cbm_comparison_table):
    paired <- LoadCBMComparisonRows(cbm_comparison_table)

    stats <- {
        wilcoxon_delta_C1,
        wilcoxon_delta_C4,
        wilcoxon_delta_convergence,
        wilcoxon_delta_invalid_geometry
    }

    IF sample_size_sufficient:
        glmm <- FitGLMMForCBMEffects(paired)
        SaveGLMM(glmm)

    SaveCBMStats(stats)

## 8) Crystal anchoring som sanity check

PROCEDURE RunCrystalSanity(crystal_anchor_table):
    IF table_is_empty:
        RETURN "no_crystal_data"

    summary <- SummarizeCrystalMetrics(
        columns=[pocket_rmsd, ligand_or_segment_rmsd, ifp_similarity]
    )

    # Ikke hard valideringsregel i statistikkdelen
    FlagPotentialMismatches(summary)
    SaveCrystalSanity(summary)

## 9) Sensitivitetsanalyser

PROCEDURE RunSensitivityAnalyses(cluster_data):
    s1_top_cluster_only <- SelectTopClusterPerEnzymeLigand(cluster_data)
    s1_results <- RefitPrimaryModel(s1_top_cluster_only)

    s2_enzyme_summary <- BuildOccupancyWeightedEnzymeSummary(cluster_data)
    s2_results <- FitSecondaryEnzymeModel(s2_enzyme_summary)

    s3_high_confidence <- FilterBySupportAndPlausibility(cluster_data)
    s3_results <- RefitPrimaryModel(s3_high_confidence)

    s4_medoid_only <- KeepMedoidRepresentation(cluster_data)
    s4_results <- RefitPrimaryModel(s4_medoid_only)

    s5_strict_geometry <- ApplyStricterGeometryThresholds(cluster_data)
    s5_results <- RefitPrimaryModel(s5_strict_geometry)

    CompareSensitivityWithPrimary([
        s1_results,
        s2_results,
        s3_results,
        s4_results,
        s5_results
    ])

    SaveSensitivityReport()

## 10) EC-avledede labels i statistikken

PROCEDURE MergeECLabels(cluster_data, ec_mapping_table):
    merged <- JoinByEnzymeID(cluster_data, ec_mapping_table)

    ValidateECMappingCompleteness(merged)

    # mapping_rule beholdes som audit-spor
    ASSERT ColumnExists(merged, "mapping_rule")

    RETURN merged

## 11) Outputkontrakt for statistikk

PROCEDURE ExportStatisticsOutputs(subanalysis):
    WriteTable("08_descriptive/descriptive_summary.tsv")
    WriteTable("09_predictive/cv_metrics.tsv")
    WriteTable("09_predictive/feature_importance.tsv")
    WriteTable("09_predictive/model_coefficients.tsv")
    WriteTable("09_predictive/sensitivity_summary.tsv")
    WriteTable("11_cbm/cbm_stats.tsv")
    WriteTable("10_crystal_anchor/crystal_sanity.tsv")
    WriteFigureSet("figures/")

## 12) Masterrunner for statistikkdelen

PROCEDURE RunAllStatistics():
    all_subanalysis_results <- []

    FOR each subanalysis in nine_independent_subanalyses:
        data <- RunStatistikkSubanalyse(subanalysis)
        desc <- KjorDeskriptiv(data)
        primary <- FitPrimaryRegioModel(data)
        alternatives <- FitAlternativeModels(data)
        sensitivity <- RunSensitivityAnalyses(data)
        ExportStatisticsOutputs(subanalysis)

        all_subanalysis_results <- all_subanalysis_results + [
            {subanalysis, desc, primary, alternatives, sensitivity}
        ]

    RunLigandSpecificityAcrossSubanalyses(all_subanalysis_results)
    BuildGlobalStatisticsSummary(all_subanalysis_results)
