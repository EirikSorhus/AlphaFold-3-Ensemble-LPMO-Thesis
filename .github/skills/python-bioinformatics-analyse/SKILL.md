---
name: python-bioinformatics-analyse
description: "Default Python coding workflow for almost all analyse/ tasks: implement or refactor parsers, validators, adapters, QC gates, artifact builders, manifest/report logic, schema contracts, and tests in analyse/src/lpmo_pipeline and analyse/tests. Use for mmCIF, Gemmi, PoseBusters, Privateer, MDAnalysis, ProLIF, HDBSCAN, and manifest consistency work. Do not use for installing programs/dependencies or R coding. Avoid re-implementing structure prediction runner internals."
argument-hint: "Describe task + module path + expected input/output artifacts + constraints"
---

# Python Bioinformatics Coding for analyse/

Use this skill for Python development in the unfinished analyse pipeline when implementing, refactoring, or reviewing modules in `analyse/src/lpmo_pipeline/` and related tests/configs.

## Goals

- Deliver robust, testable Python code for bioinformatics pipeline steps.
- Keep changes small, explicit, and compatible with existing I/O contracts.
- Preserve separation of concerns:
  - `analyse/` orchestrates analysis and validation.
  - `structure_pipeline/` remains the source of truth for prediction runner logic.
- Prefer reusable patterns that work across multiple submodules (IO, QC, analysis, reporting, tuning).

## Non-Goals

- Do not produce R code in this workflow.
- Do not re-implement AF3/RF3/Boltz prediction-running internals in `analyse/`.
- Do not hardcode site-specific HPC paths in source code.

## When to Use

- Use this as the default skill for Python work in `analyse/` unless the task is package/tool installation or R coding.
- Add or update Python modules in:
  - `analyse/src/lpmo_pipeline/io/`
  - `analyse/src/lpmo_pipeline/mapping/`
  - `analyse/src/lpmo_pipeline/qc/`
  - `analyse/src/lpmo_pipeline/analysis/`
  - `analyse/src/lpmo_pipeline/report/`
  - `analyse/src/lpmo_pipeline/tuning/`
  - `analyse/src/lpmo_pipeline/utils/`
- Implement validators, runners, parsers, and artifact builders.
- Add or adjust tests and fixture-driven checks in `analyse/tests/`.
- Improve typing, error handling, logging, and schema validation.
- Diagnose/fix artifact contract mismatches across JSON/CSV/manifest outputs.

## Inputs to Collect First

Before coding, identify:

1. Target module and step ownership (IO/QC/analysis/report/tuning/utils).
2. Required artifact contract(s), including filenames and schema requirements.
3. Hard gates and thresholds relevant to the step.
4. Whether behavior must be config-driven via YAML instead of constants.
5. Minimum regression tests needed for the changed behavior.

## Procedure

1. Locate governing docs and current implementation context.
2. Define clear function boundaries:
   - Parse/validate input
   - Transform/analyze
   - Emit artifacts with deterministic structure
3. Implement smallest viable change first.
4. Use typed interfaces and explicit dataclasses/TypedDicts for structured records.
5. Route thresholds and tunables through config objects, not literals.
6. Add focused tests for:
   - happy path
   - malformed input
   - gate boundary cases
  - schema or contract conformance
  - manifest reconciliation for artifact IDs/paths/checksums when relevant
7. Keep logs meaningful for pipeline debugging:
   - step start/end
   - key parameters
   - artifact output paths
   - gate results
8. Re-check that changes do not move prediction-runner logic into `analyse/`.
9. Avoid running heavy pipeline jobs interactively; validate with targeted/unit/contract tests unless the user explicitly asks for a full run.
10. If a step requires manual control/inspection, stop after preparing outputs and wait for user confirmation before treating the task as finished.

## Design Rules for Reusability

- Write pure helper functions where possible.
- Keep external tool wrappers thin and isolated.
- Separate data extraction from decision logic (gates/thresholds).
- Prefer composable records (dict/dataclass) over many positional tuples.
- Make file-format adapters explicit (mmCIF in, table/json out).
- Fail with actionable errors that include context (file, step, expected field).

## Python Coding Standards

- Python only for this skill.
- Follow existing project style (type hints, small functions, clear names).
- Use `snake_case` for functions/modules and `PascalCase` for classes.
- Avoid hidden side effects in utility functions.
- Keep public function signatures stable unless the task requires a breaking change.

## Bioinformatics-Specific Checklist

- mmCIF remains the canonical internal structure format.
- CCD/monosaccharide assumptions are validated before downstream QC.
- Atom mapping logic is topology/geometry aware, not name-only.
- QC outputs distinguish hard-fail vs soft-flag outcomes.
- Geometric metrics include units and stable field names.
- Clustering inputs are reproducible and traceable in manifest/report artifacts.

## Test Strategy

- Prefer unit tests with small synthetic fixtures for edge-case coverage.
- Add contract tests when artifact structure is part of pipeline guarantees.
- For external-tool-dependent paths, isolate parsing/business logic so most tests run without those tools.
- Use deterministic expected outputs for JSON/CSV where practical.
- When manifests are produced or consumed, test that artifact IDs, output paths, and checksums reconcile with manifest entries.
- For changes needing manual validation, provide a clear checklist and wait for user confirmation before marking completion.

## Definition of Done

- Code changes are implemented in the correct `analyse/` module(s) with separation of concerns preserved.
- Automated validation has been run at an appropriate level (at minimum targeted tests for changed behavior).
- Test results were reviewed, and failures were either fixed or explicitly reported with impact.
- Produced artifacts (when applicable) conform to schema/contract requirements.
- If manifests are in scope, artifact IDs/paths/checksums reconcile with manifest entries and this is covered by tests where relevant.
- Logs and outputs were inspected for obvious regressions in gate outcomes and key metrics.
- Any task requiring manual control/review is paused at a clear handoff point and awaits explicit user confirmation before being considered finished.

## Anti-Patterns to Avoid

- Mixing parse/extract logic with gate decision logic in one opaque block.
- Name-only atom mapping when topology/geometry context is required.
- Hardcoding thresholds in business logic instead of config/schema-driven values.
- Running full heavy pipeline jobs interactively when targeted tests are sufficient.
- Writing partial artifacts without validating required fields and manifest linkage.

## Adaptation Knobs

Use these knobs to adapt this skill across tasks:

- Strictness mode:
  - conservative: minimal edits + compatibility first
  - progress: implement pending module behavior end-to-end
- Validation depth:
  - schema-only
  - schema + numerical/range checks
  - schema + cross-artifact consistency
- Output focus:
  - QC-first
  - analysis-first
  - reporting-first

## Prompt Pattern for Invocation

Provide this information when invoking the skill:

- Task objective
- Target module(s)
- Expected input and output artifacts
- Required constraints (gates, schema, backward compatibility)
- Test expectations

Example:
"Implement Python parser updates in qc/privateer_runner.py to emit normalized QC fields, enforce hard-gate mapping, and add tests for malformed tool output and pass/fail boundaries."
