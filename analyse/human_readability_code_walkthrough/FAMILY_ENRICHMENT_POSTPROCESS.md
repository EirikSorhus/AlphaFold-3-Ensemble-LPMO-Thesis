# Family Enrichment Postprocess

## Purpose

This document explains the optional AA9/AA10 within-family residue comparison
that now lives outside the main production run.

The implementation is intentionally postprocess-only:

- it starts from existing Stage 16b outputs
- it does not alter `run_analysis_core()`
- nothing in the main analysis/reporting path depends on it

The entrypoint is `lpmo-pipeline family-enrichment`, which calls
`run_family_enrichment_postprocess()` in
`src/lpmo_pipeline/analysis/family_enrichment_postprocess.py`.

## Inputs

The postprocess consumes four required inputs:

- `protein_condition_residue_scores.tsv`
- `protein_residue_regio_delta.tsv`
- `input_data/metadata_final_ec_fixed.tsv`
- `input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta`

Optional input:

- `--alignment-dir` containing `AA9.aligned.fasta` and/or `AA10.aligned.fasta`
- `--substrates`, the predicted ligand classes to summarize in the
  substrate-enrichment table. The default is `cellulose chitin`.

These files are not inferred from the production config. The command receives
them explicitly.

## Active scope

The current implementation is deliberately narrow.

- Only `AA9` and `AA10` are processed.
- Only `domain_only` Stage 16b rows are retained.
- Full-length rows are ignored.

The construct filtering happens by parsing `condition_id` and selecting the
second token. This is intentional because the first version assumes that domain
only residue numbers already correspond to catalytic-core-local numbering.

## Sequence source and mapping

The catalytic-core FASTA is the canonical sequence source.

The parser in `analysis/family_alignment.py`:

- reads `UniProtIDs|...|...` FASTA headers
- expands semicolon-separated accession lists
- resolves metadata rows through `UniProt_ID`, `CoAccessions`, and
  `Cat_CoAccessions`

The family aggregation unit is `Cat_Seq_Group` when present. This matters
because several proteins or UniProt accessions may share the exact same
catalytic-core sequence.

The pipeline therefore keeps two levels at once:

- per-protein provenance rows in `family_aligned_residue_table.tsv`
- deduplicated family-level summary statistics over unique catalytic sequence
  groups in `family_residue_enrichment.tsv`
- deduplicated target-vs-other substrate residue summaries in
  `family_substrate_residue_enrichment.tsv`
- deduplicated right-vs-wrong predicted ligand summaries for unambiguous
  cellulose-only/chitin-only proteins in
  `family_wrong_ligand_residue_enrichment.tsv`

## Alignment strategy

The postprocess resolves one alignment per family.

Preferred path:

- use a precomputed `<family>.aligned.fasta` from `--alignment-dir`

Fallback path:

- write one sequence per catalytic sequence group to a temporary FASTA
- run `mafft --localpair --maxiterate 1000`
- read the aligned FASTA back and map ungapped residue positions to 1-based
  alignment columns

Alignment record IDs are the catalytic sequence-group identifiers, not raw
protein IDs. This is what lets the summary layer deduplicate identical
catalytic cores without losing the per-protein provenance rows.

## Residue mapping

The postprocess does not align every residue in the metadata sequence space.
Instead it only maps residues that already occur in Stage 16b.

For each family:

1. Collect the union of residues present in the condition-level Stage 16b table
   and the protein-level regio-delta table.
2. Resolve the protein to its catalytic sequence group.
3. Convert the ungapped catalytic-core residue number to the aligned column.
4. Build explicit `residue_alignment_rows` for
   `compute_family_residue_enrichment_outputs()`.

These alignment rows now carry explicit provenance fields:

- `family_aggregation_id`
- `construct_type`
- `catalytic_core_residue_index`
- `alignment_residue`

## Outputs

The command writes under `08_family_residue_enrichment/`:

- `family_aligned_residue_table.tsv`
- `family_residue_enrichment.tsv`
- `family_substrate_residue_enrichment.tsv`
- `family_wrong_ligand_residue_enrichment.tsv`
- `family_alignment_manifest.tsv`
- `family_enrichment_summary.json`

`family_substrate_residue_enrichment.tsv` uses `substrate_class` from the
condition-level Stage 16b rows as the predicted/input ligand type. For each
requested target substrate and aligned residue position, it reports target and
non-target condition counts, protein counts, family aggregation-unit counts,
target-minus-other contact-score deltas, and `substrate_group_status`. Statuses
flag missing target rows, target-only/no-contrast cases, small groups, and
strongly unbalanced groups. The table is descriptive; it is not treated as a
confirmatory statistical test.

`family_wrong_ligand_residue_enrichment.tsv` is stricter. It derives active
cellulose/chitin labels from explicit metadata labels when present, otherwise
from EC numbers (`1.14.99.54`/`.56` for cellulose and `1.14.99.53` for
chitin). It keeps only proteins with exactly one active substrate in the
cellulose/chitin pair and compares that right predicted ligand against the
opposite predicted ligand within the same protein/alignment position.
Cellulose+chitin dual-active proteins and unknown-activity proteins are
excluded from the contrast and counted in `wrong_ligand_group_status`.

`family_alignment_manifest.tsv` is the bridge between the alignment input and
the computed residue outputs. It records which family aggregation IDs were used,
which proteins belong to them, the alignment source, and the aligned FASTA path.

`family_enrichment_summary.json` records:

- processed families
- skipped families and skip reasons
- unresolved metadata-to-sequence resolutions
- output file paths

## Skip semantics

This layer is intentionally soft-failing because it is not prioritized and may
be cut.

The command skips a family instead of failing the whole job when:

- no AA9/AA10 members remain after filtering
- fewer than 2 catalytic sequence groups remain in the family
- no precomputed alignment exists and `mafft` is unavailable
- MAFFT fails
- the alignment is structurally inconsistent with the canonical catalytic-core
  sequences
- no Stage 16b residues can be mapped into alignment columns

The overall command still completes and writes a summary JSON in these cases.

## Validation surfaces

Focused automated validation now exists in:

- `tests/test_family_alignment.py`
- `tests/test_family_residue_enrichment.py`
- `tests/test_family_enrichment_postprocess.py`
- `tests/test_cli_run.py`

Cluster-backed real-data validation is available through:

- `tests/run_tests_scripts/test_family_enrichment_validation.sh`
- `tests/run_tests_scripts/run_family_enrichment_validation.py`

The `sbatch` wrapper runs the focused family pytest suite first, then executes
the Python harness against the staged domain-only shard root and requires at
least one processed family. The underlying Python harness merges Stage 16b
outputs from that staged shard root, runs the optional postprocess, and writes
a standalone validation summary.
