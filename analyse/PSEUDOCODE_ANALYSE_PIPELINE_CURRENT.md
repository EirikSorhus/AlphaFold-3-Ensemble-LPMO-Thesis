# Current Analyse Pipeline Pseudocode

Source of truth: current code in `src/lpmo_pipeline` and `scripts/run_analysis_full_pipeline.slurm.sh`.

## Current Scope / Excluded Legacy Logic

- Current production analysis is AF3-oriented and reads prediction artifacts
  from the structure-pipeline work roots.
- There is no active PLACER step in the current analyse pipeline.
- There is no active cross-model clustering step in the current analyse pipeline.
- Stale README/codewalkthrough files were not used as source of truth.

## Module Map

```text
CLI / config / runtime
    cli.py                       command routing and argument resolution
    config.py                    runtime paths, defaults, run-config merging
    utils/manifest.py            manifests, tool versions, gate outcomes

IO / preparation
    io/discovery.py              prediction artifact discovery
    io/normalize_mmcif.py        CIF normalization, atom maps, reports
    io/cif_to_pdb.py             PoseBusters PDB export
    io/analysis_export.py        complex/ligand PDB export for ProLIF
    io/gemmi_compat.py           gemmi compatibility wrapper

QC
    qc/hard_qc_orchestrator.py   hard QC orchestration
    qc/posebusters_runner.py     PoseBusters checks
    qc/privateer_runner.py       Privateer input/checks
    qc/*gates*, qc_report.py     gate evaluation and qc_report.json

Analysis core
    analysis/analysis_orchestrator.py  production control path
    analysis/mdanalysis_metrics.py     pose geometry
    analysis/prolif_ifp.py             ProLIF IFP and eligibility
    analysis/convergence_metrics.py    condition convergence

Clustering / signatures / residue importance
    clustering_hdbscan.py        primary condition-wise Jaccard/HDBSCAN
    clustering_agglomerative.py  orthogonal sensitivity clustering
    clustering_pilot.py          optional pilot matrices/sensitivity outputs
    cluster_signatures.py        cluster/IFP/residue signature outputs
    residue_importance.py        residue scores and regio deltas
    residue_region_annotation.py core/non-core annotations
    residue_contact_extraction.py pose-residue contact table

Crystal anchoring
    analysis/crystal_anchoring.py  representative-vs-crystal screening

Reports / postprocess
    condition_summary.py         condition/protein summary tables
    predictive_postprocess.py    predictive modeling
    cbm_comparison.py            full-length vs domain-only paired analysis
    family_enrichment_postprocess.py AA9/AA10 residue enrichment
    report/*                     metrics.csv, summary.json, report.html

SLURM production wrapper
    scripts/run_analysis_full_pipeline.slurm.sh
        shard configs, arrays, merge, optional global postprocess
```

## 1. PROGRAM lpmo analyse pipeline

```text
PROGRAM lpmo analyse pipeline

INPUTS:
    runtime_paths.yaml
    production config YAML or hybrid run-config YAML
    AF3 prediction work roots
    optional protein metadata TSV
    optional catalytic-core FASTA
    optional reused prepared cases / QC / IFP / analysis artifacts

PRIMARY SINGLE-RUN OUTPUTS:
    run_manifest.json
    analysis_core_summary.json
    pose accounting: pose_manifest.tsv, pose_confidence.tsv,
        structure_index.tsv, qc_attrition_table.tsv
    QC/geometry/IFP: qc_report.json, pose_geometry.tsv,
        pose_ifp_table.tsv, pose_residue_contact_table.tsv
    convergence: pose_convergence.tsv, condition_convergence_summary.tsv
    clustering: cluster_assignments.tsv, medoid_manifest.tsv,
        condition_cluster_summary.tsv, cluster_table.tsv,
        cluster_ifp_signature.tsv, cluster_residue_signature.tsv,
        cluster_signatures.json
    residue importance: protein_condition_residue_scores.tsv,
        protein_residue_regio_delta.tsv, condition_patch_summary.tsv,
        protein_patch_summary.tsv
    crystal anchoring: crystal_anchor_table.tsv, crystal_geometry_table.tsv,
        crystal_ifp_diagnostic_summary.tsv
    report surfaces: metrics.csv, summary.json, report.html,
        condition_table.tsv, protein_summary_table.tsv
    optional predictive outputs

FULL ARRAY OUTPUTS:
    output/_plan/{shard_plan.json,generated_configs,shard_indexes,slurm_scripts}
    output/{domain_only,full_length}/shards/*
    output/{domain_only,full_length}/merged/*.tsv
    output/combined/*.tsv
    output/merge_summary.json
    output/postprocess/*
```

## 2. CLI Command Routing

```text
PROCEDURE MainCLI(args):
    Create argparse parser
    Register subcommands:
        run
        discover
        tune
        predictive
        cbm-paired
        family-enrichment

    Parse args

    IF command == "run":
        RETURN CmdRun(args)
    ELSE IF command == "discover":
        RETURN CmdDiscover(args)
    ELSE IF command == "tune":
        RETURN CmdTune(args)
    ELSE IF command == "predictive":
        RETURN CmdPredictive(args)
    ELSE IF command == "cbm-paired":
        RETURN CmdCbmPaired(args)
    ELSE IF command == "family-enrichment":
        RETURN CmdFamilyEnrichment(args)
    ELSE:
        Print help
        RETURN 1
```

```text
PROCEDURE ResolveCommandArgs(args, command_name):
    IF args.run_config exists:
        run_config = LoadPipelineRunConfig(args.run_config)
        command_config = Merge(run_config.base, run_config.commands[command_name])
    ELSE:
        command_config = empty mapping

    FOR each required argument:
        Use explicit CLI arg if provided
        ELSE use matching command_config value
        ELSE raise missing-required error

    FOR each optional argument:
        Use explicit CLI arg if provided
        ELSE use matching command_config value
        ELSE use command default

    RETURN resolved values
```

## 3. Runtime / Config Resolution

```text
PROCEDURE LoadRuntimePathsConfig(optional_config_path):
    raw_path =
        optional_config_path
        OR environment LPMO_PIPELINE_RUNTIME_PATHS_CONFIG
        OR configs/runtime_paths.yaml

    Resolve project root relative to config location
    Resolve external tool settings:
        apptainer executable
        Privateer container candidates
        PoseBusters container candidates
    Resolve runtime settings:
        default Privateer mode
        Python executable
    Resolve package assets:
        thresholds config
        ProLIF feature config
        defaults config
        geometry rules
        residue rules
        CV hierarchy
        QC report schema

    Cache and return RuntimePathsConfig
```

```text
PROCEDURE LoadProductionOptions(config, output_dir, del_variant, n_jobs):
    production = config.production
    construct_type = production.construct_type OR "domain_only"
    work_root = production.work_roots[construct_type] OR production.work_root
    REQUIRE work_root

    Read discovery controls:
        af3_only
        latest_only
        max_cases
        include_targets
        include_proteins

    Read construct / metadata controls:
        run_id
        protein_metadata_path
        run_posebusters
        run_privateer

    Read clustering feature policy:
        main_include_interaction_types
        extra_include_interaction_types
        excluded_interaction_types
        rare feature thresholds

    Read primary clustering policy:
        minimum_clusterable_n
        insufficient_clusterable_signal_label
        HDBSCAN min_cluster_size
        HDBSCAN min_samples
        HDBSCAN cluster_selection_method
        agglomerative sensitivity parameters

    IF clustering_pilot.enabled:
        Build ClusteringPilotOptions

    IF predictive.enabled:
        REQUIRE predictive.protein_metadata_path
        Build PredictiveAnalysisOptions

    Read reuse controls:
        reuse_prepared_cases_from
        reuse_hard_qc_report_path
        reuse_ifp_table_path
        reuse_analysis_artifacts_from
        allow_reuse_fallback_generation

    RETURN ProductionRunOptions
```

## 4. Single `run_analysis_core` Production Flow

```text
PROCEDURE RunAnalysisCore(config, output_dir, del_variant, n_jobs):
    run_start = now
    options = LoadProductionOptions(config, output_dir, del_variant, n_jobs)
    Create output_dir

    discovery_errors, pose_inputs, discovery_summary =
        DiscoverPoseInputs(options)

    fallback_errors, top_model_fallback_by_condition =
        DiscoverTopModelFallbackInputs(options)

    analysis_summary = Initialize summary with:
        run id, output dir, del variant, work root
        discovery options
        QC options
        reuse options and counters
        predictive options
        primary clustering policy
        primary clustering feature policy
        discovery summary and errors
        top-model crystal fallback candidates

    cases, prepared_poses = PreparePoseCases(
        pose_inputs,
        output_dir,
        n_jobs,
        reuse_prepared_cases_from
    )

    Record preparation cases and reuse counters

    IF no prepared_poses:
        Write always-on pose/condition summary surfaces
        Write analysis_core_summary.json with failure reason
        RETURN AnalysisCoreResult(success=False)

    qc_report = RunOrLoadHardQC(prepared_poses, options)
    Write and schema-validate qc_report.json

    Load ProLIF feature config
    Load contact eligibility rule
    Load protein region definitions

    FOR each prepared pose:
        Attach QC verdict and condition_id
        IF QC verdict is not dropped:
            Compute pose geometry metrics
            Export or reuse complex/ligand PDBs for ProLIF
            Register pose for condition-level IFP analysis
        ELSE:
            Mark analysis skipped
        Build metrics record and summary geometry row

    FOR each condition with QC-passing poses:
        RunPerConditionAnalysis(condition)

    IF clustering pilot enabled:
        Write pilot matrices, prevalence tables, and sensitivity outputs

    Write IFP, residue-contact, convergence, clustering, cluster-signature,
    residue-importance, crystal-anchoring, geometry, metrics, summary, and
    report outputs when inputs exist

    Always write:
        pose_manifest.tsv
        pose_confidence.tsv
        structure_index.tsv
        qc_attrition_table.tsv
        condition_table.tsv
        protein_summary_table.tsv

    IF production.predictive.enabled:
        Run predictive postprocess on just-produced condition_table.tsv

    Write analysis_core_summary.json

    RETURN AnalysisCoreResult with paths and stage-completion flags
```

## 5. Pose Discovery and Preparation

```text
PROCEDURE DiscoverPoseInputs(options):
    manifest = DiscoverWorkRoot(
        work_root = options.work_root,
        af3_only = options.af3_only,
        latest_only = options.latest_only
    )

    pose_inputs = []
    errors = []

    FOR each target, model, run_id, uniprot_target in manifest:
        IF include_targets is set and target not selected:
            CONTINUE
        FOR each sample in uniprot_target.samples:
            FOR each CIF path in sample.cif_paths:
                Extract protein_id and ligand_id
                IF include_proteins is set and protein_id not selected:
                    CONTINUE
                Infer source run id from CIF path
                Select matching confidence JSON if available
                Build PoseInputRecord:
                    cif path
                    pose id
                    protein id
                    ligand id
                    model
                    source/discovered run id
                    confidence JSON path
                    run status
                    seed/sample
                Append to pose_inputs
                Stop at max_cases if configured

    RETURN errors, pose_inputs, discovery_summary
```

```text
PROCEDURE PreparePoseCases(pose_inputs, output_dir, n_jobs, reuse_root):
    IF reuse_root is set:
        Load reusable case metadata by pose id

    IF n_jobs <= 1:
        FOR pose in pose_inputs:
            Prepare one pose case
    ELSE:
        Prepare pose cases in ProcessPoolExecutor

    FOR each pose:
        IF reusable prepared case exists:
            Validate required normalized CIF, PoseBusters PDB, Privateer input
            Rebuild PreparedPose from reused artifacts
        ELSE:
            case_dir = output_dir / "prepared_cases" / safe pose id
            Normalize original CIF
            Convert normalized CIF to PoseBusters PDB
            Prepare Privateer input CIF
            Read normalized structure with gemmi
            Build PreparedPose

        Record status, paths, warnings, and errors in case summary

    RETURN case summaries, prepared poses
```

## 6. Hard QC and Attrition

```text
PROCEDURE RunOrLoadHardQC(prepared_poses, options):
    IF options.reuse_hard_qc_report_path exists:
        Load previous QC report
        Require verdict for every prepared pose
        Return report scoped to this run

    hard_qc_inputs = []
    FOR prepared_pose in prepared_poses:
        hard_qc_inputs += HardQCInput(
            pose_id,
            normalized_cif,
            posebusters_pdb,
            privateer_input_cif,
            case_dir
        )

    report = RunHardQC(
        hard_qc_inputs,
        run_id = options.run_id,
        run_posebusters = options.run_posebusters,
        run_privateer = options.run_privateer,
        max_workers = options.n_jobs
    )

    RETURN report
```

```text
PROCEDURE BuildQCAttritionTable(pose_inputs, case_by_pose_id):
    FOR each original pose input:
        Read preparation status
        Read QC verdict if present
        Count and group drop reasons and warnings
        Include condition id, protein, ligand, model, seed/sample, paths
    Write qc_attrition_table.tsv
```

## 7. Per-Pose Geometry / Export / IFP Setup

```text
PROCEDURE AnalyzePreparedPose(prepared_pose, verdict):
    pose = prepared_pose.pose
    condition_id = ConditionIdForPose(pose, construct_type)

    IF verdict.status == "dropped":
        case.analysis_status = "skipped_dropped"
        case.geometry_metrics_status = "skipped"
        Record skipped timing events if enabled
        RETURN

    TRY:
        metrics = ComputePoseMetricsFromStructure(
            prepared_pose.structure,
            pose id, model, protein id, ligand id
        )
        case.analysis_status = "analyzed"
        case.geometry_metrics_status = "ok"
    EXCEPT:
        metrics = GeometryNotComputableMetrics(pose)
        case.analysis_status = "analyzed"
        case.geometry_metrics_status = "flagged"
        case.analysis_flags += "geometry_metrics_error"

    TRY:
        IF reusable analysis export exists:
            complex_pdb, ligand_pdb = reused paths
        ELSE:
            ExportAnalysisArtifacts(normalized_cif, analysis_export_dir)
    EXCEPT or export blockers:
        IF reusable IFP exists:
            Allow condition-level use of reused IFP
        ELSE:
            Create failed IFP result for this pose

    IF export ok or reusable IFP exists:
        Add pose entry to ifp_pose_inputs_by_condition[condition_id]

    Build pose metrics record:
        ids, model, seed, QC verdict, geometry, del variant,
        substrate class, DP, CBM present flag
```

## 8. Per-Condition Convergence, IFP, Eligibility, Clustering

```text
PROCEDURE RunPerConditionAnalysis(condition_id, condition_pose_inputs):
    clustering_result = empty result
    production_clustering_status = "no_ifp_inputs"
    formal_clustering_allowed = false

    convergence_inputs = condition entries with complex PDB files
    IF convergence_inputs exist:
        ComputeConditionConvergence(convergence_inputs)
        Update per-pose convergence metrics in case summaries

    reused_entries = entries with reusable IFP result
    compute_entries = entries needing ProLIF computation

    IF compute_entries exist:
        computed_batch = ComputeIFPBatch(
            complex PDB and ligand PDB for each entry,
            protein id,
            ligand id,
            model,
            max_workers = options.n_jobs
        )

    batch = BuildIFPBatchFromReusedAndComputedResults(
        condition_pose_inputs,
        reused_entries,
        computed_results
    )

    Write full condition IFP matrix:
        ifp_matrices/<condition>/ifp_matrix.csv

    region_annotations = BuildRegionAnnotationsForIFPResults(
        batch results,
        protein id,
        construct type,
        region definitions
    )

    core_batch = Keep only core-region IFP features when region definitions apply

    FOR each result in batch:
        Update case with IFP status, error, and contact counts

    clusterable_indices = []
    FOR each result in core_batch:
        IF result status is ok or zero_contacts:
            eligibility = EvaluateContactEligibility(result, rule)
            Store eligibility and case metrics
        IF result status is ok AND eligibility.eligible:
            clusterable_indices += result index

    IF clustering_pilot enabled:
        Save PilotConditionIFP for later pilot output generation

    IF clusterable_indices exist:
        main_matrix = BuildMainClusteringMatrix(
            batch,
            selected clusterable pose ids,
            clustering feature options,
            construct type,
            region definitions
        )
        Write:
            ifp_matrices/<condition>/main_clustering_ifp_matrix.csv

        production_clustering_status, formal_clustering_allowed =
            ClassifyPrimaryClusteringInput(
                n_clusterable poses,
                n selected main features,
                minimum_clusterable_n
            )

        IF formal_clustering_allowed:
            clustering_result = HDBSCANClusterer.cluster(
                Jaccard matrix,
                cluster pose ids
            )
            Compute Jaccard distance matrix
            Append cluster assignment rows
            Append medoid rows
            Update cluster id in case summaries and metrics records
    ELSE:
        production_clustering_status =
            ClassifyPrimaryClusteringInput(0, 0)

    Choose crystal representatives from clustering result and medoids
    IF no representatives:
        Try hard-QC-passing top-model fallback for condition
    Run crystal anchoring for selected representatives

    Build condition cluster summary:
        QC pass count
        IFP success count
        contact eligibility fractions
        clustering status
        formal clustering flag
        noise fraction
        cluster occupancy / entropy / gini
```

## 9. Cluster Signatures and Residue Importance

```text
PROCEDURE BuildClusterAnnotationOutputs():
    IF no QC-passing conditions:
        Skip cluster annotation outputs

    Build pose metadata by pose id:
        protein id
        ligand id
        construct type
        substrate class
        DP
        condition id

    Build condition metadata by condition id
    Read confidence summaries by pose id
    Collect convergence rows by pose id

    cluster_signature_tables = BuildClusterSignatureTables(
        cluster assignments,
        medoid rows,
        geometry rows by pose id,
        all IFP results,
        residue contact rows,
        pose metadata,
        condition metadata,
        confidence rows,
        convergence rows
    )

    Write:
        cluster_ifp_signature.tsv
        cluster_residue_signature.tsv
        cluster_signatures.json
        cluster_table.tsv

    residue_importance_outputs = ComputeResidueImportanceOutputs(
        residue signature rows,
        cluster summaries,
        cluster IFP signature rows,
        observed residue contact rows,
        condition metadata
    )

    Write:
        protein_condition_residue_scores.tsv
        protein_residue_regio_delta.tsv
        condition_patch_summary.tsv
        protein_patch_summary.tsv
```

## 10. Crystal Anchoring and Fallback Representatives

```text
PROCEDURE ChooseCrystalRepresentatives(clustering_result, cluster_pose_ids):
    representatives = []

    FOR each cluster id with medoid:
        medoid_pose_id = clustering_result.medoids[cluster id]
        IF medoid pose has prepared data:
            Add CrystalRepresentativeCandidate:
                role = medoid / cluster representative
                pose id
                normalized CIF
                complex PDB
                ligand PDB
                IFP result

    RETURN representatives
```

```text
PROCEDURE PrepareTopModelFallback(condition_id):
    fallback_pose = top_model_fallback_by_condition[condition_id]
    IF no fallback pose:
        RETURN none

    Prepare or reuse fallback normalized CIF, PDB exports, and QC result
    IF fallback passes hard QC:
        Build CrystalRepresentativeCandidate with role = top_model_fallback
    ELSE:
        Record fallback not usable

    RETURN candidate or none
```

```text
PROCEDURE RunCrystalAnchoring(representatives):
    FOR each representative:
        crystal_output_dir =
            output/crystal_anchoring/<condition>/<representative_pose>

        TRY:
            crystal_report = RunCrystalReferenceScreen(
                representative CIF,
                protein id,
                ligand id,
                representative pose id,
                output dir,
                optional normalized CIF,
                optional complex PDB,
                optional ligand PDB,
                optional representative IFP result
            )
            Write crystal_reference_screen.json

            Extract best IFP Tanimoto and best pocket RMSD
            Update representative case and metrics record

            FOR each crystal comparison:
                Append crystal_anchor_table row
                Append crystal_geometry_table row
        EXCEPT:
            Append crystal anchoring error row
            Mark representative case with error

    Write:
        crystal_anchor_table.tsv
        crystal_geometry_table.tsv
        crystal_ifp_diagnostic_summary.tsv
```

## 11. Output Tables, Summary JSON, HTML Report, Manifest Gates

```text
PROCEDURE WriteCoreOutputs():
    IF all_ifp_results exist:
        Write pose_ifp_table.tsv
        Build and write pose_residue_contact_table.tsv

    IF pose_convergence_metrics exist:
        Write pose_convergence.tsv
        Write condition_convergence_summary.tsv

    IF condition cluster summaries exist:
        Write cluster_assignments.tsv
        Write medoid_manifest.tsv
        Write condition_cluster_summary.tsv
        Build cluster signatures and residue importance outputs

    Write crystal anchoring outputs
    Write pose_geometry.tsv
    Write metrics.csv

    summary_json = BuildSummary(
        run id,
        mode = production,
        del variant,
        QC report payload,
        cluster results,
        geometry stats,
        crystal reports,
        pipeline version
    )
    Write summary.json
    Build report.html from summary.json and metrics.csv

    Always write run accounting outputs:
        pose_manifest.tsv
        pose_confidence.tsv
        structure_index.tsv
        qc_attrition_table.tsv

    condition_table = BuildConditionTableRows(
        qc attrition,
        condition cluster summary,
        condition convergence summary,
        condition patch summary,
        pose confidence,
        cluster table
    )
    Write condition_table.tsv

    protein_summary_table = BuildProteinSummaryRows(condition_table)
    Write protein_summary_table.tsv

    IF predictive enabled in production config:
        RunPredictivePostprocess(condition_table, metadata, output_dir)

    Write analysis_core_summary.json
```

```text
PROCEDURE CmdRun(args):
    Resolve mode, config path, output dir, del branch, n_jobs
    REQUIRE mode == production
    REQUIRE del branch is del_a or del_b

    Create output dir
    manifest = ManifestBuilder(mode="production", version="2.1")
    Record git info and tool versions
    Load config YAML and record config hash
    IF tuning_reference exists in config:
        Add tuning reference to manifest

    Record production_pipeline_started = true

    TRY:
        result = RunAnalysisCore(config, output, del branch, n_jobs)
        Record manifest gates:
            analysis_core_completed = result.success
            hard_qc_completed = qc_report_path exists
            geometry_stage_completed = pose_geometry path exists
            ifp_stage_completed = pose_ifp_table path exists
            clustering_stage_completed = condition_cluster_summary path exists
            cluster_annotation_stage_completed = result flag
            crystal_anchoring_stage_completed = result flag
            predictive_stage_completed if predictive enabled
        Write run_manifest.json

        IF result.success is false:
            Print error and RETURN 1
        Print key output paths
        RETURN 0
    EXCEPT:
        Record failed analysis and crystal gates
        Write run_manifest.json
        RETURN 1
```

## 12. Full SLURM Array Workflow

```text
PROCEDURE RunAnalysisFullPipelineSlurm(args):
    Parse wrapper args:
        runtime paths
        base/core/full configs
        output dir
        metadata
        core FASTA
        proteins per shard
        max parallel shards
        n_jobs per shard
        predictive task
        random state
        shard and summary resources
        postprocess toggle
        family alignment dir
        MAFFT executable
        account / partition
        submission metadata path
        dry-run flag

    Validate required args and files
    Export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG
    Resolve Python executable from runtime_paths.yaml

    Create:
        output/_plan
        output/_plan/generated_configs
        output/_plan/shard_indexes
        output/_plan/slurm_scripts
        output/logs

    BuildShardPlan(core_config, full_config)
    Submit domain_only and full_length shard arrays
    Submit dependent collect/postprocess job after successful arrays
    Optionally write submission metadata JSON
    Print submitted job ids and shard plan path
```

```text
PROCEDURE BuildShardPlan(core_config, full_config):
    entries = [
        ("domain_only", "del_a", core_config),
        ("full_length", "del_b", full_config)
    ]

    FOR each construct_type, del_branch, source_config:
        Load YAML config
        Set production.construct_type

        Discover proteins under selected construct work root:
            work_roots[construct_type] OR work_root
            include targets if configured
            scan target/af3/input/*_<target>_data.json

        IF include_proteins configured:
            Filter configured proteins by discovered proteins
            Warn for configured proteins not found
            Refuse unfiltered submit if no proteins discovered
        ELSE:
            Use discovered proteins

        Split proteins into shards of proteins_per_shard

        FOR each shard:
            Copy base config
            Set construct_type
            Set run_id suffix with construct and shard id
            Set include_proteins to shard protein list
            Write generated shard config YAML
            Register task:
                task index
                construct type
                del branch
                shard id
                config path
                shard output dir
                selected proteins

        Write construct shard index TSV

    Write output/_plan/shard_plan.json
```

```text
PROCEDURE ShardArrayWorker(index_tsv, repo_root, runtime_paths, n_jobs_per_shard):
    task_id = SLURM_ARRAY_TASK_ID
    Lookup shard id, del branch, config path, output dir in index TSV
    Export PYTHONPATH and runtime paths env var
    Resolve Python executable
    Create shard output dir and shard log

    Run:
        python -m lpmo_pipeline.cli run
            --mode production
            --config shard config
            --output shard output
            --del del_a or del_b
            --n-jobs n_jobs_per_shard

    Require shard output cluster_table.tsv
    Log completion
```

```text
PROCEDURE CollectAndPostprocess(
    repo_root,
    runtime_paths,
    output_dir,
    shard_plan_json,
    run_global_postprocess,
    protein_metadata,
    core_fasta,
    predictive_task,
    random_state,
    family_alignment_dir,
    mafft_executable
):
    Load shard_plan.json

    Per-construct surfaces to merge:
        condition_table.tsv
        cluster_table.tsv
        protein_condition_residue_scores.tsv
        protein_residue_regio_delta.tsv

    FOR each construct:
        FOR each surface:
            Read surface from each shard if present
            Add _construct_type and _shard_id columns
            Write output/<construct>/merged/<surface>

    FOR each combined surface:
        Read domain_only/merged and full_length/merged versions if present
        Write output/combined/<surface>

    Write merge_summary.json
    REQUIRE output/combined/condition_table.tsv

    IF run_global_postprocess == true:
        Run predictive over combined condition table
        Run cbm-paired over combined condition table and optional cluster table
        IF combined residue score and delta tables exist:
            Run family-enrichment over combined Stage 16b tables
        ELSE:
            Warn and skip family-enrichment
```

## 13. Optional Commands

```text
PROCEDURE CmdDiscover(args):
    DiscoverWorkRoot(work_root, af3_only, latest_only)
    Print summary, write manifest JSON, or print full JSON
```

```text
PROCEDURE CmdPredictive(args):
    Resolve condition table, protein metadata, output dir, task, folds, seed
    Create output dir
    RunPredictivePostprocess(...)
    Print summary, modeling table, metrics, prediction paths
```

```text
PROCEDURE CmdCbmPaired(args):
    Resolve condition table, optional cluster/signature/metadata tables,
    output dir, random_state
    RunCbmPairedAnalysis(...)
    Print summary and table paths
```

```text
PROCEDURE CmdFamilyEnrichment(args):
    Resolve residue score table, regio-delta table, metadata, core FASTA,
    output dir, optional alignment dir, families, substrates, MAFFT executable
    RunFamilyEnrichmentPostprocess(...)
    Print summary and output table paths
```

```text
PROCEDURE CmdTune(args):
    Resolve model, tuning config, output dir, n_jobs
    Create output dir
    manifest = ManifestBuilder(mode="tune", version="2.1")
    Record git info and tool versions
    Load tuning config and record hash

    TRY:
        RunTuning(model, tuning config, test_cases=[], output dir, config, n_jobs)
        Record tuning_completed = true
        Write run_manifest.json
        Print best params and summary paths
        RETURN 0
    EXCEPT:
        Record tuning_completed = false
        Write run_manifest.json
        RETURN 1
```
