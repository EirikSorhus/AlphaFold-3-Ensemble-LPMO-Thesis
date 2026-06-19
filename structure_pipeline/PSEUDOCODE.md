# Structure Pipeline Pseudocode

Source of truth: `CODE_WALKTHROUGH.md`, `README.md`, and current code under
`src/structure_pipeline`.

Note: `cli.py`, `cases.py`, tests, and docs reference
`src/structure_pipeline/manifest/*`, but that package is not present in this
checkout. Manifest pseudocode follows the documented API and current call sites.

## Module Map

```text
cli.py                  Typer commands, orchestration, submit flow
config.py               Pydantic config, path expansion, output dirs
manifest/proteins.py    FASTA -> proteins.csv
manifest/ligands.py     CIF scan or CCD list -> ligands.csv
manifest/msa.py         MSA matching -> msa.csv
cases.py                proteins x ligands x models -> cases.csv
resume.py               DONE.ok/status.jsonl resume and status
executors/*             SLURM/local script submission
runners/af3.py          AF3 CPU MSA stage + GPU inference stage
runners/boltz.py        Boltz YAML batch per ligand
runners/rf3.py          RF3 JSON batch per ligand
runners/ligand_utils.py CU ligand helpers
runners/oligo.py        AF3 oligosaccharide helpers
```

## 1. CLI Program

```text
PROGRAM structure-pipeline

COMMANDS:
    version       print package version
    init-config   write example YAML config
    manifest      normalize raw inputs into CSV manifests
    validate      check input/model paths
    run           submit pending prediction jobs
    status        summarize work directory status

PROCEDURE Main():
    Create Typer app
    Dispatch command
```

## 2. Config Loading

```text
PROCEDURE LoadConfig(optional_path):
    config_path =
        optional_path
        OR env STRUCTURE_PIPELINE_CONFIG
        OR config/pipeline.yaml
        OR ~/structure_pipeline.yaml
        OR fail

    Read YAML into PipelineConfig:
        paths: AF3/Boltz/RF3 containers, weights, databases, checkpoints
        slurm: account, partitions, time, CPU, GPU, memory, tmp base
        af3/boltz/rf3: model-specific runtime parameters
        inputs: FASTA, ligand dir, optional MSA dir/squashfs, optional oligos
        outputs: work_dir, manifest_dir, results_dir
        models_enabled: af3, boltz, rf3 subset

    Expand environment variables
    Resolve relative input/output paths
    Create output directories
    Validate enabled model names
    RETURN cfg
```

## 3. Manifest Command

```text
PROCEDURE Manifest(config, force, ligand_ccd_list):
    cfg = LoadConfig(config)
    manifest_dir = cfg.outputs.manifest_dir

    IF proteins.csv exists AND NOT force:
        proteins = LoadProteinsManifest(proteins.csv)
    ELSE:
        proteins = ParseFasta(cfg.inputs.proteins_fasta)
        WriteProteinsManifest(proteins, proteins.csv)

    IF ligands.csv exists AND NOT force:
        ligands = LoadLigandsManifest(ligands.csv)
    ELSE IF ligand_ccd_list:
        ligands = LoadCcdList(ligand_ccd_list)
        WriteLigandsManifest(ligands, ligands.csv)
    ELSE:
        ligands = ScanLigands(cfg.inputs.ligands_dir)
        WriteLigandsManifest(ligands, ligands.csv)

    IF msa.csv exists AND NOT force:
        msa_records = LoadMsaManifest(msa.csv)
        missing_msa = proteins without an MSARecord
    ELSE IF cfg.inputs.msa_dir OR cfg.inputs.msa_sqsh_files:
        msa_records, missing_msa = MatchMsa(proteins, msa_dir, msa_sqsh_files)
        WriteMsaManifest(msa_records, msa.csv)
    ELSE:
        msa_records = []
        missing_msa = all protein IDs

    IF cases.csv exists AND NOT force:
        cases = LoadCasesManifest(cases.csv)
    ELSE:
        cases = GenerateCases(proteins, ligands, msa_records,
                              cfg.models_enabled, missing_msa)
        WriteCasesManifest(cases, cases.csv)

    Print counts by ligand category
    Print model-wise pending and skipped case counts
```

## 4. Input Normalization

```text
PROCEDURE ParseFasta(fasta_path):
    FOR each FASTA entry:
        Parse header as pipe-separated, legacy underscore, UniProt, or fallback
        protein_id = first UniProt ID before semicolon
        sequence = uppercase joined sequence lines
        Yield ProteinRecord with IDs, organism, annotation, sequence,
            sequence hash, length, raw header

PROCEDURE ScanLigands(ligands_dir):
    Find sorted CIF files in root, category dirs, and category/cif dirs
    Skip hidden dirs and boltz_ccd_lib
    FOR each CIF:
        ccd_code = _chem_comp.id OR data_ block OR filename stem
        Yield LigandRecord(L001..., category, ccd_code, filename,
            resolved cif_path, cif SHA256)

PROCEDURE LoadCcdList(path):
    Read non-empty, non-comment CCD lines
    Deduplicate in input order
    Yield LigandRecord records with category ccd_list and empty CIF fields

PROCEDURE MatchMsa(proteins, msa_dir, msa_sqsh_files):
    index = UniProt-like ID -> MSA location
    FOR each .a3m/.a3m.gz in msa_dir:
        extract IDs from query header and filename
        index each ID to full file path
    FOR each squashfs file:
        list contents with unsquashfs -l
        extract IDs from contained MSA filenames
        index each ID to filename plus sqsh_file path

    FOR each protein:
        candidate IDs = protein_id plus all semicolon-separated UniProt IDs
        matches = indexed entries for candidate IDs
        IF no matches:
            add protein_id to missing_msa
        ELSE:
            choose deterministic first match if ambiguous
            create MSARecord with path, source, matched_id, method, sqsh_file

    RETURN msa_records, missing_msa
```

## 5. Case Planning

```text
PROCEDURE GenerateCases(proteins, ligands, msa_records, models, missing_msa):
    msa_by_protein = map protein_id -> MSARecord

    FOR each protein, ligand, model:
        IF model == af3:
            msa_source = native
            msa_path = ""
            status = pending
        ELSE IF protein_id is in missing_msa:
            msa_source = mmseqs
            msa_path = ""
            status = skipped
            error = "Missing MSA"
        ELSE:
            msa_source = mmseqs
            msa_path = matched MSA path if present
            status = pending

        case_id = SHA256(protein_id, ligand_id, model, msa_source)[0:12]
        Yield Case(case_id, protein_id, ligand_id, ligand_ccd_code,
                   model, msa_source, msa_path, status, error)

PROCEDURE GroupCasesByLigandModel(cases):
    FOR each non-skipped case:
        grouped[(case.ligand_ccd_code, case.model)].append(case)
    RETURN grouped
```

## 6. Validate Command

```text
PROCEDURE Validate(config):
    cfg = LoadConfig(config)
    Check cfg.inputs.proteins_fasta and cfg.inputs.ligands_dir
    Check optional cfg.inputs.msa_dir and each cfg.inputs.msa_sqsh_file

    FOR each enabled model:
        runner = matching AF3Runner, BoltzRunner, or RF3Runner
        errors = runner.ValidatePaths()
        Print missing container, weights, database, or checkpoint paths

    Exit with error if any check failed, otherwise print success
```

## 7. Run Command

```text
PROCEDURE Run(config, dry_run, resume, model_filter, protein_filter, local,
              overrides, ligand_ccd_list, oligo_definitions):
    cfg = LoadConfig(config)
    Apply CLI overrides to cfg.af3, cfg.boltz, cfg.rf3

    Require cases.csv
    cases = LoadCasesManifest(cases.csv)
    proteins = LoadProteinsManifest(proteins.csv)
    protein_sequences = map protein_id -> sequence
    msa_by_protein = LoadMsaManifest(msa.csv) if present else {}

    ligand_cif_paths = {}
    IF ligands.csv exists AND ligand_ccd_list is not provided:
        ligand_cif_paths = ligand_id -> cif_path for ligands with CIF paths

    Filter cases by optional model and protein filters
    pending_cases = StatusTracker(work_dir).GetPendingCases(cases, resume)
    IF no pending cases:
        print all completed
        return

    grouped = GroupCasesByLigandModel(pending_cases)
    Print planned ligand/model groups

    IF dry_run:
        collect model paths, input paths, output paths, CIF paths
        Print ComputeBindMounts(paths)
        return

    oligo_registry = empty OR OligoRegistry.FromYaml(configured oligo path)
    executor = LocalExecutor if local else SlurmExecutor
    runners = AF3Runner(cfg, oligo_registry), BoltzRunner(cfg), RF3Runner(cfg)

    SubmitAF3MsaJobs(grouped)
    FOR each ligand/model group:
        IF model == af3:
            SubmitAF3Inference(group)
        ELSE IF model == boltz:
            SubmitBoltz(group)
        ELSE IF model == rf3:
            SubmitRF3(group)
```

## 8. AF3 Flow

```text
PROCEDURE SubmitAF3MsaJobs(grouped):
    af3_msa_jobs = protein_id -> job/work_dir
    first_ligand_for_protein = first AF3 ligand seen for each protein

    FOR each AF3 case:
        IF protein already has MSA job:
            continue
        msa_work_dir = work_dir/af3_msa/protein_id
        IF msa_work_dir/DONE.ok exists:
            remember completed MSA with no dependency
            continue
        script = AF3Runner.BuildMsaSlurmScript(protein_id, sequence,
            first_ligand_for_protein, msa_work_dir, optional ligand CIF)
        job = executor.Submit(script, msa_work_dir, af3_msa_<protein_id>)
        remember job and msa_work_dir

PROCEDURE SubmitAF3Inference(group):
    work_dir = cfg.outputs.work_dir/ligand_ccd/af3
    msa_work_dirs = protein_id -> af3_msa work_dir
    dependency = afterok:<MSA job IDs> for SLURM-submitted MSA jobs
    script = AF3Runner.BuildSlurmScript(group, sequences, work_dir,
        af3_<ligand>_inf, msa_work_dirs, ligand_cif_paths)
    executor.Submit(script, work_dir, job_name, dependency)

PROCEDURE AF3Runner scripts:
    MSA script:
        build AF3 JSON with protein A, CU, main ligand, optional oligo bonds/CIF
        run AF3 CPU with --run_inference=false
        copy *_data.json to input/produced_msa.json
        touch DONE.ok
    Inference script:
        patch produced_msa.json for each protein/ligand
        run AF3 GPU with configured diffusion samples and optional recycles
        symlink latest and touch DONE.ok
```

## 9. Boltz and RF3 Flow

```text
PROCEDURE SubmitBoltz(group):
    Collect group protein sequences and MSA paths
    script = BoltzRunner.BuildSlurmScript(group, sequences, work_dir,
        job_name, msa_paths, optional msa_sqsh_file)
    executor.Submit(script, work_dir, job_name)

PROCEDURE BoltzRunner script:
    Write one YAML per protein under work/input
    Include sequence, optional MSA, CU plus main ligand, optional affinity
    Run apptainer boltz predict with configured parameters
    touch DONE.ok

PROCEDURE SubmitRF3(group):
    Collect group protein sequences, MSA paths, optional ligand CIF paths
    script = RF3Runner.BuildSlurmScript(group, sequences, work_dir,
        job_name, msa_paths, optional msa_sqsh_file, ligand_cif_paths)
    executor.Submit(script, work_dir, job_name)

PROCEDURE RF3Runner script:
    Write one JSON array with one example per protein
    Include sequence, optional MSA, CU, and main ligand CCD or CIF
    Set Foundry/RF3 environment and run RF3 fold with configured parameters
    symlink latest and touch DONE.ok
```

## 10. Executors and Status

```text
PROCEDURE SlurmExecutor.Submit(script, work_dir, job_name, dependency):
    Write executable work_dir/job_name.sh
    IF dry_run:
        return SubmittedJob(dry_run_<job_name>)
    Run sbatch with optional --dependency
    Parse "Submitted batch job <id>"
    Return SubmittedJob

PROCEDURE StatusTracker.GetPendingCases(cases, resume):
    IF resume is false:
        return all non-skipped cases
    FOR each non-skipped case:
        skip if work_dir/ligand/model/DONE.ok exists
        skip if status.jsonl has completed event for case_id
        otherwise keep pending
    RETURN pending cases

PROCEDURE Status(config):
    cfg = LoadConfig(config)
    Walk cfg.outputs.work_dir by ligand/model directories
    Count completed when DONE.ok exists
    Count pending when status.jsonl exists but DONE.ok is absent
    Print model, completed, pending, failed table

END PROGRAM
```
