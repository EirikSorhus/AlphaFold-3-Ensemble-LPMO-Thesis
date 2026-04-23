---
name: r-bioinformatics-analyse
description: "Default R coding workflow for analyse/ bioinformatics tasks: implement or refactor R analysis scripts, statistical workflows, plotting/report artifacts, and R-Python handoffs for lpmo_pipeline outputs. Use this as the default for R work in analyse/, with tidyverse + ggplot2 + readr/dplyr as first-choice packages, and CSV-first artifact handoff to Python unless another format is explicitly requested or better justified. Do not use for package installation tasks or for re-implementing structure prediction runner internals."
argument-hint: "Describe task + target path + expected input/output artifacts + plotting/statistical constraints"
---

# R Bioinformatics Coding for analyse/

Use this skill as the default for R development in the unfinished analyse pipeline when implementing, refactoring, or reviewing statistical analysis, plotting logic, and R-generated artifacts that must work with Python-produced artifacts.

## Goals

- Deliver robust, reproducible R code for bioinformatics analysis, statistics, and visualization.
- Keep changes small, explicit, and compatible with existing `analyse/` artifact contracts.
- Preserve separation of concerns:
  - `analyse/` owns analysis/validation/report artifacts.
  - `structure_pipeline/` remains the source of truth for prediction runner logic.
- Ensure R outputs can be consumed by Python steps (and vice versa) without schema drift.

## Non-Goals

- Do not use this workflow for package installation or environment provisioning.
- Do not re-implement AF3/RF3/Boltz prediction-running internals in `analyse/`.
- Do not hardcode site-specific HPC paths in source code.

## When to Use

- Use this as the default skill for R coding tied to `analyse/` pipeline artifacts.
- Prefer `tidyverse` + `ggplot2` + `readr`/`dplyr` as the baseline package stack.
- Add or update R scripts/notebooks that:
  - read normalized CSV/TSV/JSON/parquet outputs from `analyse/`
  - compute statistical summaries, comparisons, and QC-oriented derived metrics
  - generate publication-ready figures and tables for reports/thesis material
  - produce machine-readable outputs that are fed back into Python-driven pipeline steps
- Improve reproducibility, statistical clarity, and figure consistency across runs.
- Diagnose/fix contract mismatches between R outputs and Python consumers.

## Inputs to Collect First

Before coding, identify:

1. Target path(s) and ownership of the step in the analyse workflow.
2. Required artifact contract(s), including filenames, schema, and delimiters.
3. Statistical method expectations (test/model assumptions, correction methods, effect sizes).
4. Plot constraints (required aesthetics, facets, scales, labels, units, themes, output formats).
5. Minimum regression checks needed for changed behavior and R-Python compatibility.

## Procedure

1. Locate governing docs and current implementation context.
2. Confirm input artifact contracts from Python-producing steps (columns, types, units, IDs).
3. Define clear boundaries:
   - ingest/validate inputs
   - transform/analyze
   - emit deterministic figures/tables/artifacts
4. Implement the smallest viable change first.
5. Keep statistical logic explicit:
   - state assumptions
   - encode multiple-testing correction when needed
   - include effect sizes/confidence intervals where relevant
6. Build plotting code that is deterministic and reusable:
   - stable factor ordering
   - explicit palette and theme
   - explicit axis/unit labeling
   - fixed output dimensions/DPI for reproducibility
7. Emit outputs with stable schemas and predictable paths for downstream Python consumers.
  - Default to CSV for handoff artifacts.
  - Use alternative formats (for example, parquet/JSON) only when explicitly requested or clearly more beneficial.
8. Add focused checks/tests for:
   - happy path
   - malformed/missing input fields
   - boundary cases in statistical thresholds
   - schema/contract conformance for cross-language handoff
9. Re-check that prediction-runner internals are not reimplemented in `analyse/`.
10. Avoid heavy full-pipeline runs unless explicitly requested; prefer targeted checks and artifact-contract validation.
11. If a step requires manual visual/statistical review, stop at a clear handoff and wait for user confirmation.

## Design Rules for Reusability

- Keep data ingestion, statistical logic, and plotting as separate functions/sections.
- Prefer pure functions that accept data frames and return data frames/plot objects.
- Centralize plotting defaults (theme/palette/labels) to avoid drift across figures.
- Use explicit joins/keys for merges; fail early on duplicate/missing IDs.
- Keep output writers deterministic (column order, sort order, locale-safe formatting).
- Surface actionable error messages with file/column/context details.

## R Coding Standards

- R only for this skill's code-generation path, while preserving interoperability with Python outputs.
- Use clear object/function names and keep functions small.
- Prefer tidy, explicit pipelines over side-effect-heavy scripts.
- Avoid hidden global state; pass dependencies/configs explicitly.
- Keep public function/script interfaces stable unless breaking change is required by task.

## Bioinformatics-Specific Checklist

- Biological entity identifiers are preserved and traceable across joins.
- Units are explicit and consistent before statistical comparisons.
- Grouping variables/factors reflect experimental design and are not inferred implicitly.
- QC-related thresholds are represented explicitly and documented in outputs.
- Figure captions/labels align with computed metrics and artifact field names.
- Cross-artifact traceability is maintained for report/manifests when applicable.

## R-Python Interoperability Checklist

- IDs and key columns exactly match Python-side expectations.
- Data types are normalized for handoff (no ambiguous logical/factor encodings).
- Missing values use a documented representation understood by both languages.
- Output schemas are versioned/validated when contract-sensitive.
- CSV remains the primary interchange format unless task constraints call for another format.
- Round-trip spot checks confirm Python can parse R-produced artifacts.

## Test and Validation Strategy

- Prefer lightweight fixture-based checks for edge cases and schema conformance.
- Validate statistical outputs against known small examples where practical.
- Validate plot generation for expected facets/scales/labels and non-empty outputs.
- Add contract checks for artifacts consumed by Python modules.
- For tasks requiring manual interpretation, provide a clear review checklist and await confirmation.

## Definition of Done

- R code changes are implemented in the correct analyse-related paths.
- Statistical and plotting behavior is explicit, reproducible, and documented in code/output naming.
- Produced artifacts conform to expected schema/contract requirements.
- Cross-language (R-Python) compatibility checks were performed at an appropriate level.
- Validation/check results were reviewed, and failures were fixed or reported with impact.
- Any manual review task is paused at a clear handoff and awaits explicit user confirmation.

## Anti-Patterns to Avoid

- Blending ingestion, modeling, and plotting into one opaque script block.
- Implicit factor ordering or locale-dependent output formatting.
- Undocumented statistical choices (tests, corrections, thresholds).
- Emitting figures/tables without stable naming or schema guarantees.
- Running heavy full workflows when targeted validation is sufficient.

## Adaptation Knobs

Use these knobs to adapt this skill across tasks:

- Strictness mode:
  - conservative: minimal edits + compatibility first
  - progress: implement pending analysis/report logic end-to-end
- Validation depth:
  - schema-only
  - schema + statistical sanity checks
  - schema + cross-artifact and cross-language consistency
- Output focus:
  - stats-first
  - plot-first
  - report-first

## Prompt Pattern for Invocation

Provide this information when invoking the skill:

- Task objective
- Target path(s)
- Expected input and output artifacts
- Required statistical/plot constraints
- Compatibility constraints for Python consumers
- Validation expectations

Example:
"Implement an R analysis script that consumes analyse QC tables, computes groupwise effect sizes with adjusted p-values, emits a normalized summary CSV for Python ingestion, and generates a publication-ready ggplot2 figure with deterministic ordering and labels."