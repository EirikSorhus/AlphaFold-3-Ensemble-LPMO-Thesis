# Pipeline flowcharts for Masteroppgave_clean

This document provides four Mermaid flowcharts describing the main processing
stages in `Masteroppgave_clean`:

- `pipe_test/scripts/module_2_new`: protein metadata and sequence enrichment
- `structure_pipeline`: manifest-driven structure prediction on HPC
- `analyse`: downstream analysis of predicted AF3 structures
- an integrated overview showing how the modules are connected

The diagrams are based on the project documentation and current code
walkthroughs in `module_2_new`, `structure_pipeline`, and `analyse`.

## 1. Metadata and Sequence Enrichment

```mermaid
flowchart TD
  A[Start main_driver.py] --> B[Parse command-line arguments]
  B --> C{Select input mode}

  C -->|fasta| D1[Read FASTA input]
  C -->|cazy| D2[Read local CAZy file or download CAZy family]
  C -->|list| D3[Read plain UniProt accession list]
  C -->|characterized| D4[Read curated semicolon-delimited CSV]

  D1 --> E1[Parse FASTA headers and identify UniProt accessions]
  D1 --> E2{Unknown headers?}
  E2 -->|yes, optional| E3[Identify sequences by NCBI BLAST]
  E2 -->|no| F
  E3 --> F

  D2 --> E4[Parse CAZy NCBI accessions and JGI groups]
  E4 --> E5[Batch fetch NCBI or GenBank records through UniProt]
  E4 --> E6[Resolve JGI groups using UniProt queries]
  E5 --> F
  E6 --> F

  D3 --> F[Batch fetch records from UniProt]

  D4 --> E7[Separate UniProt and GenBank or RefSeq identifiers]
  E7 --> E8[Fetch UniProt accessions directly]
  E7 --> E9[Map GenBank or RefSeq identifiers to UniProt cross-references]
  E9 --> E10{Mapping successful?}
  E10 -->|no, optional fallback| E11[Fetch sequence directly from NCBI efetch]
  E10 -->|yes| F
  E8 --> F
  E11 --> F

  F --> G[Retrieve InterPro domain annotations]
  G --> H[Extract UniProt sequence features]
  H --> I[Annotate signal peptide, transmembrane regions, LPMO core, and H1]
  I --> J[Select domains and CBMs using source priority and overlap filtering]

  J --> K{Curated CSV mode?}
  K -->|yes| L[Deduplicate identical sequences and assign co-accessions]
  K -->|no| M[Retain one metadata row per resolved accession]
  L --> N[Detect CSV rows that resolve to multiple distinct sequences]
  M --> O
  N --> O

  O[Collect successful records and failed identifiers] --> P[Write metadata_expanded_TIMESTAMP.tsv]
  O --> Q[Write all_sequences_TIMESTAMP.fasta]
  O --> R[Write failed_ids_TIMESTAMP.txt]
  O --> S[Write run_metadata_TIMESTAMP.json]

  P --> T[Curated input for protein selection and structure prediction]
  Q --> T
  R --> U[Failure list for manual review]
  S --> V[Run statistics and provenance]
```

## 2. Structure Prediction Pipeline

```mermaid
flowchart TD
  A[Start structure-pipeline CLI] --> B{Command}

  B -->|init-config| C[Write example YAML configuration]
  B -->|manifest| D[Load PipelineConfig]
  B -->|validate| E[Validate inputs, containers, databases, and checkpoints]
  B -->|run| F[Load manifests and execution plan]
  B -->|status| G[Summarize work directory from DONE.ok and status.jsonl]

  D --> H[Parse protein FASTA]
  H --> I[Write or load proteins.csv]
  D --> J{Ligand source}
  J -->|ligand directory| K[Scan CIF files and extract CCD codes]
  J -->|CCD list| L[Read deduplicated CCD list]
  K --> M[Write or load ligands.csv]
  L --> M

  D --> N{MSA source available?}
  N -->|MSA directory or squashfs| O[Match MSA files to protein identifiers]
  N -->|none| P[Mark all proteins as missing MSA for non-AF3 models]
  O --> Q[Write or load msa.csv]
  P --> R
  Q --> R[Generate prediction cases]

  R --> S[Enumerate protein by ligand by model combinations]
  S --> T{Model}
  T -->|AF3| U[Create pending case with native MSA generation]
  T -->|Boltz or RF3 with MSA| V[Create pending case with matched MSA]
  T -->|Boltz or RF3 without MSA| W[Skip case with missing MSA]

  U --> X[Write cases.csv]
  V --> X
  W --> X
  X --> Y[Manifest files define an explicit execution plan]

  E --> Z[Run model-specific path validation]
  Z --> AA{All paths valid?}
  AA -->|no| AB[Stop before GPU or SLURM submission]
  AA -->|yes| AC[Proceed to execution]

  F --> AD[Filter cases by model, protein, and resume state]
  AD --> AE[Group cases by ligand and model]
  AE --> AF{Model group}

  AF -->|AF3| AG[Submit CPU MSA job per protein]
  AG --> AH[Patch AF3 data JSON with ligand-specific inputs]
  AH --> AI[Submit GPU inference per ligand with afterok dependency]

  AF -->|Boltz-2| AJ[Build one YAML input per protein]
  AJ --> AK[Submit one batch job per ligand]

  AF -->|RF3| AL[Build JSON array of protein and ligand components]
  AL --> AM[Submit one Foundry/RF3 batch job per ligand]

  AI --> AN[Write model outputs under work_dir/runs/JOBID]
  AK --> AN
  AM --> AN
  AN --> AO[Update latest symlink, DONE.ok, and confidence outputs]
  AO --> AP[Predicted CIF files and model confidence metrics]
  AP --> AQ[Input for downstream analysis]
```

## 3. Downstream Analysis Pipeline

```mermaid
flowchart TD
  A[Start lpmo-pipeline run] --> B[Parse CLI arguments: mode, config, output, del, n_jobs]
  B --> C[Load production YAML and runtime_paths.yaml]
  C --> D[Build ProductionRunOptions]
  D --> E[Create output root and run_manifest.json]

  E --> F[Discover prediction outputs in structure_pipeline work root]
  F --> G[Select target, model, run, protein, seed, and sample CIF files]
  G --> H[Apply af3_only, latest_only, target, protein, and case filters]
  H --> I[Build PoseInputRecord list]
  F --> J[Build separate AF3 top-model fallback map for crystal anchoring]

  I --> K[Prepare one analysis case per pose]
  K --> L[Normalize mmCIF to canonical chain layout]
  L --> M[Assign protein A, glycans B-D, and metal E]
  M --> N[Validate atom mapping and glycan CCD codes]
  N --> O[Export PoseBusters-compatible PDB]
  N --> P[Prepare Privateer input CIF]
  O --> Q[PreparedPose]
  P --> Q

  Q --> R[Hard quality control]
  R --> S[Check active-site proximity and Cu geometry]
  S --> T[Run PoseBusters in dock mode]
  T --> U[Run Privateer glycan validation]
  U --> V{QC verdict}
  V -->|dropped| W[Exclude pose from downstream analysis]
  V -->|passed or flagged| X[Retain pose for downstream analysis]

  X --> Y[Compute downstream geometry metrics]
  X --> Z[Export non-protonated ProLIF analysis structures]
  Z --> AA[Compute convergence metrics per condition]
  Z --> AB[Compute ProLIF interaction fingerprints]
  AB --> AC[Extract residue-level contacts]

  AA --> AD[Assess clustering eligibility]
  AC --> AD
  AD --> AE[Cluster poses using HDBSCAN with Jaccard distance]
  AE --> AF[Assign clusters and identify medoids]
  AF --> AG[Annotate cluster signatures]
  AG --> AH[Calculate residue importance and patch scores]
  AH --> AI[Build condition and protein summaries]

  AI --> AJ[Perform crystal anchoring]
  J -. fallback for no-cluster conditions .-> AJ
  AJ --> AK[Compare retained medoids with crystal references]

  AK --> AL[Write reports and analysis tables]
  W --> AL
  AL --> AM[pose_manifest.tsv, pose_confidence.tsv, structure_index.tsv]
  AL --> AN[qc_attrition_table.tsv and QC reports]
  AL --> AO[pose_geometry.tsv, pose_ifp_table.tsv, and residue contacts]
  AL --> AP[cluster_table.tsv, condition_table.tsv, protein_summary_table.tsv]
  AL --> AQ[metrics.csv, summary.json, report.html, analysis_core_summary.json]

  B -. side branch .-> AR[lpmo-pipeline discover]
  B -. side branch .-> AS[lpmo-pipeline tune]
  AE -. side branch .-> AT[Clustering pilot and parameter sensitivity analysis]
```

## 4. Integrated Workflow

```mermaid
flowchart TD
  A[Public and curated protein sources] --> B

  subgraph M2[Part 1: Protein set construction]
    B[Input records from FASTA, CAZy, accession lists, or curated tables]
    C[Resolve identifiers and retrieve sequence metadata]
    D[Annotate domains, signal peptides, LPMO core regions, H1, and CBMs]
    E[Export curated protein metadata and sequence sets]
    B --> C --> D --> E
  end

  E --> F[Selected protein and ligand set]

  subgraph SP[Part 2: Structure prediction]
    G[Define proteins, ligands, model settings, and MSA resources]
    H[Build explicit manifests and prediction cases]
    I[Validate HPC environment and model dependencies]
    J[Run structure prediction jobs]
    K[Collect predicted structures and confidence metrics]
    G --> H --> I --> J --> K
  end

  F --> G
  K --> L[Predicted protein-ligand structures]

  subgraph AN[Part 3: Downstream structural analysis]
    M[Discover prediction outputs and prepare analysis cases]
    N[Normalize structures and perform hard quality control]
    O[Compute geometric, interaction, and convergence descriptors]
    P[Cluster retained poses and identify representative medoids]
    Q[Compare predictions with crystal references]
    R[Summarize results for figures, tables, and interpretation]
    M --> N --> O --> P --> Q --> R
  end

  L --> M
  R --> S[Thesis figures, statistical summaries, and biological interpretation]

  classDef module2 fill:#e8f4ff,stroke:#2f6f9f,stroke-width:1.5px,color:#111827
  classDef structure fill:#ecfdf3,stroke:#2f7d4f,stroke-width:1.5px,color:#111827
  classDef analysis fill:#fff4e6,stroke:#b36b00,stroke-width:1.5px,color:#111827
  classDef bridge fill:#f7f7f7,stroke:#6b7280,stroke-width:1px,color:#111827

  class B,C,D,E module2
  class G,H,I,J,K structure
  class M,N,O,P,Q,R analysis
  class A,F,L,S bridge

  style M2 fill:#e8f4ff,stroke:#2f6f9f,stroke-width:2px
  style SP fill:#ecfdf3,stroke:#2f7d4f,stroke-width:2px
  style AN fill:#fff4e6,stroke:#b36b00,stroke-width:2px
```
