---
description: "Plan and start the next implementation step in the analyse pipeline. Reads status tables, asks clarifying questions, proposes priorities, and begins implementation."
name: "Next Implementation Step"
argument-hint: "Optional: specific module or step to focus on (e.g. 'stage 13b', 'residue_contact_extraction')"
agent: "agent"
---

You are helping implement the next step in the LPMO analyse pipeline (`Masteroppgave/analyse/`).

## Step 1 — Read governing documents

Read the following files to understand current status and priorities **before doing anything else**:

1. [IMPLEMENTATION_PLAYBOOK.md](../../analyse/IMPLEMENTATION_PLAYBOOK.md) — Status table (✅/⚠️/❌) and ordered implementation steps
2. [AF3_LPMO_pipeline_detailed_plan.md](../../analyse/AF3_LPMO_pipeline_detailed_plan.md) — Primary governing source for pipeline design and artefact contracts
3. [OPEN_QUESTIONS.md](../../analyse/OPEN_QUESTIONS.md) — Unresolved decisions that must not be guessed
4. [MASTERPLAN.md](../../analyse/MASTERPLAN.md) — Scope invariants, output table specs, conflict priority

Also scan the relevant source modules and test files under `analyse/src/lpmo_pipeline/` and `analyse/tests/` for the steps you are about to work on.

## Step 2 — Assess status

From the status table in IMPLEMENTATION_PLAYBOOK.md, identify:

- All steps marked ❌ (not started)
- All steps marked ⚠️ (started but not verified on real data)
- Which of these have blockers recorded in OPEN_QUESTIONS.md

Produce a **short numbered list** of the top candidates for the next step, with one sentence of rationale each (dependency order, blocking others, complexity).

## Step 3 — Ask clarifying questions before proceeding

Before writing any code, ask the user:

1. **Which step to implement first** — Show the ranked candidate list from Step 2 and ask which to start with. If `$input` was provided, use that as the suggested default.
2. **Real data availability** — Does real AF3 data exist and is it accessible for manual verification after the step is done? (Relevant for steps marked ⚠️.)
3. **Specific constraints or comments** — Are there any design decisions, edge cases, performance concerns, or open questions from OPEN_QUESTIONS.md that need resolving before starting?
4. **Test strategy** — Should the implementation include new narrow contract tests, extend existing ones, or rely on a manual smoke-test on real data? What test data or fixtures are available?

Do **not** start implementation until the user has answered.

## Step 4 — Plan the implementation

Once the user has answered, produce a concise implementation plan:

- List the files to create or modify
- Describe the key logic for each file (≤3 bullet points each)
- Identify any schema or config changes needed in `schemas/` or `configs/`
- Identify which existing tests to update and what new tests to add
- Flag any hard rules from `copilot.md` or `analyse.instructions.md` that apply

Present the plan and ask: **"Does this plan look right, or do you want to adjust anything before I start?"**

Wait for confirmation.

## Step 5 — Implement

After confirmation:

- Implement one module or closely related artefact step at a time
- Follow all hard rules in `copilot.md` and `analyse.instructions.md`
- Do **not** change locked parameters, thresholds, or stage ordering without explicit user approval
- Do **not** implement PLACER (removed from pipeline, 2026-04-21)
- Keep prediction logic in `structure_pipeline/`; analysis code in `analyse/`
- Add or update tests in `analyse/tests/` near changed behaviour
- Update the status row in IMPLEMENTATION_PLAYBOOK.md when a step is complete

After each module is done, summarise:
- What was implemented
- How to verify it (command or test to run)
- What the next natural step would be

---

**Governing document priority (conflict resolution):**
1. `AF3_LPMO_pipeline_detailed_plan.md` — primary
2. User comments in this chat — override all documents
3. `MASTERPLAN.md` and `IMPLEMENTATION_PLAYBOOK.md`
4. Legacy specs (`plan_implementation_spec.txt`, `plan_analyse.txt`) — context only

**Hard invariants (never violate):**
- Main analysis: AF3 only; RF3 and Boltz-2 are excluded
- 75 poses per protein–ligand condition (15 seeds × 5 diffusion samples); runs completed
- Primary analysis unit: cluster (do not collapse to protein/enzyme as primary unit)
- IFP clustering uses IFP features only; geometry is linked after clustering
- mmCIF is the master format
- Locked thresholds live in `configs/thresholds.yaml` (`geometry_plausibility.locked = true`)
- Large changes to stage order or which analyses run require explicit user approval
