#!/usr/bin/env bash
#SBATCH --job-name=test_full_pipeline_3proteins
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --output=logs/test_full_pipeline_3proteins_%j.out
#SBATCH --error=logs/test_full_pipeline_3proteins_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G

# ---------------------------------------------------------------------------
# End-to-end pipeline smoke test for 3 proteins with full-length predictions.
#
# Proteins: B2ADG1, Q7SHI8, D6Y7U3
# Constructs: domain-only (del_a) AND full-length (del_b)
#
# Stages:
#   1. pytest unit tests
#   2. run (domain-only / del_a)
#   3. run (full-length / del_b)
#   4. merge condition_table + cluster_table from both constructs
#   5. predictive postprocess (on domain-only condition_table)
#   6. cbm-paired postprocess (on merged condition_table)
#   7. family-enrichment postprocess (on domain-only outputs)
#
# Usage (from the analyse/ project root):
#   sbatch scripts/test_full_pipeline_3proteins.slurm.sh [--output <dir>]
# ---------------------------------------------------------------------------

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "[ERROR] Submit with sbatch, not bash directly." >&2
  echo "[INFO]  sbatch scripts/test_full_pipeline_3proteins.slurm.sh" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT="/cluster/work/projects/nn1003k/eirik/Masteroppgave_clean/analyse"
PYTHON_BIN="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
RUNTIME_PATHS="$PROJECT_ROOT/configs/runtime_paths.yaml"
PROTEIN_METADATA="$PROJECT_ROOT/input_data/metadata_final_ec_fixed.tsv"
CORE_FASTA="$PROJECT_ROOT/input_data/lpmo_core_domain_2026-03-14_06-52-15_deduplicated.fasta"

OUTPUT_DIR="$PROJECT_ROOT/tests/tests_results/test_full_pipeline_3proteins_${SLURM_JOB_ID}"

# Parse optional --output override
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) OUTPUT_DIR="$2"; shift 2 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
  esac
done

N_JOBS="${SLURM_CPUS_PER_TASK:-9}"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for f in "$RUNTIME_PATHS" "$PROTEIN_METADATA" "$CORE_FASTA"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Required file not found: $f" >&2
    exit 1
  fi
done

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Python not executable: $PYTHON_BIN" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT/src:${PYTHONPATH:-}"
export LPMO_PIPELINE_RUNTIME_PATHS_CONFIG="$RUNTIME_PATHS"
export XDG_CACHE_HOME="$OUTPUT_DIR/cache"

mkdir -p logs
mkdir -p "$OUTPUT_DIR"
mkdir -p "$XDG_CACHE_HOME"

CORE_OUTPUT="$OUTPUT_DIR/core"
FL_OUTPUT="$OUTPUT_DIR/full_length"
COMBINED_OUTPUT="$OUTPUT_DIR/combined"
PREDICTIVE_OUTPUT="$OUTPUT_DIR/predictive"
CBM_OUTPUT="$OUTPUT_DIR/cbm_paired"
FAMILY_OUTPUT="$OUTPUT_DIR/family_enrichment"

echo "[INFO] Project root:    $PROJECT_ROOT"
echo "[INFO] Python:          $PYTHON_BIN"
echo "[INFO] Output root:     $OUTPUT_DIR"
echo "[INFO] n_jobs:          $N_JOBS"

# ---------------------------------------------------------------------------
# Stage 0: Unit tests
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 0] Running unit tests..."
"$PYTHON_BIN" -m pytest tests/ -q --tb=short 2>&1 | tee "$OUTPUT_DIR/pytest_output.txt"
echo "[STAGE 0] Unit tests passed"

# ---------------------------------------------------------------------------
# Stage 1: domain-only production run (del_a)
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 1] Running domain-only (del_a) analysis..."
"$PYTHON_BIN" -m lpmo_pipeline.cli run \
  --mode production \
  --config configs/test_3proteins_core.yaml \
  --output "$CORE_OUTPUT" \
  --del del_a \
  --n-jobs "$N_JOBS"

CORE_CONDITION_TABLE="$CORE_OUTPUT/condition_table.tsv"
CORE_CLUSTER_TABLE="$CORE_OUTPUT/cluster_table.tsv"

if [[ ! -f "$CORE_CONDITION_TABLE" ]]; then
  echo "[ERROR] domain-only condition_table.tsv missing after run." >&2
  exit 1
fi
echo "[STAGE 1] domain-only run complete: $CORE_CONDITION_TABLE"

# ---------------------------------------------------------------------------
# Stage 2: full-length production run (del_b)
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 2] Running full-length (del_b) analysis..."
"$PYTHON_BIN" -m lpmo_pipeline.cli run \
  --mode production \
  --config configs/test_3proteins_full_length.yaml \
  --output "$FL_OUTPUT" \
  --del del_b \
  --n-jobs "$N_JOBS"

FL_CONDITION_TABLE="$FL_OUTPUT/condition_table.tsv"
FL_CLUSTER_TABLE="$FL_OUTPUT/cluster_table.tsv"

if [[ ! -f "$FL_CONDITION_TABLE" ]]; then
  echo "[ERROR] full-length condition_table.tsv missing after run." >&2
  exit 1
fi
echo "[STAGE 2] full-length run complete: $FL_CONDITION_TABLE"

# ---------------------------------------------------------------------------
# Stage 3: Merge condition_table and cluster_table from both constructs
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 3] Merging condition_table and cluster_table..."
mkdir -p "$COMBINED_OUTPUT"

COMBINED_CONDITION_TABLE="$COMBINED_OUTPUT/condition_table.tsv"
COMBINED_CLUSTER_TABLE="$COMBINED_OUTPUT/cluster_table.tsv"

"$PYTHON_BIN" - <<PYEOF
import csv, sys
from pathlib import Path

def read_tsv(p):
    p = Path(p)
    if not p.exists():
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(rows, p):
    if not rows:
        print(f"[WARN] No rows to write: {p}", file=sys.stderr)
        return
    fieldnames = list(dict.fromkeys(k for row in rows for k in row))
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})
    print(f"[INFO] Wrote {len(rows)} rows -> {p}")

core_cond  = read_tsv("$CORE_CONDITION_TABLE")
fl_cond    = read_tsv("$FL_CONDITION_TABLE")
write_tsv(core_cond + fl_cond, "$COMBINED_CONDITION_TABLE")

core_clust = read_tsv("$CORE_CLUSTER_TABLE") if Path("$CORE_CLUSTER_TABLE").exists() else []
fl_clust   = read_tsv("$FL_CLUSTER_TABLE") if Path("$FL_CLUSTER_TABLE").exists() else []
write_tsv(core_clust + fl_clust, "$COMBINED_CLUSTER_TABLE")
PYEOF

echo "[STAGE 3] Merge complete"

# ---------------------------------------------------------------------------
# Stage 4: Predictive postprocess (domain-only condition_table only)
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 4] Running predictive postprocess..."
mkdir -p "$PREDICTIVE_OUTPUT"

"$PYTHON_BIN" -m lpmo_pipeline.cli predictive \
  --condition-table "$CORE_CONDITION_TABLE" \
  --protein-metadata "$PROTEIN_METADATA" \
  --output "$PREDICTIVE_OUTPUT" \
  --task all \
  --n-folds 5 \
  --random-state 42

PREDICTIVE_SUMMARY="$PREDICTIVE_OUTPUT/predictive_summary.json"
if [[ ! -f "$PREDICTIVE_SUMMARY" ]]; then
  echo "[WARN] predictive_summary.json not produced (may be expected with only 3 proteins)."
else
  echo "[STAGE 4] Predictive complete: $PREDICTIVE_SUMMARY"
fi

# ---------------------------------------------------------------------------
# Stage 5: CBM paired postprocess (merged condition_table)
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 5] Running cbm-paired postprocess..."
mkdir -p "$CBM_OUTPUT"

CBM_ARGS=(
  --condition-table "$COMBINED_CONDITION_TABLE"
  --output "$CBM_OUTPUT"
  --protein-metadata "$PROTEIN_METADATA"
  --random-state 42
)
if [[ -f "$COMBINED_CLUSTER_TABLE" ]]; then
  CBM_ARGS+=(--cluster-table "$COMBINED_CLUSTER_TABLE")
else
  echo "[WARN] No combined cluster_table.tsv – running cbm-paired without it."
fi

"$PYTHON_BIN" -m lpmo_pipeline.cli cbm-paired "${CBM_ARGS[@]}"

CBM_SUMMARY="$CBM_OUTPUT/cbm_paired_analysis_summary.json"
if [[ ! -f "$CBM_SUMMARY" ]]; then
  echo "[WARN] cbm_paired_analysis_summary.json not produced."
else
  echo "[STAGE 5] CBM paired complete: $CBM_SUMMARY"
fi

# ---------------------------------------------------------------------------
# Stage 6: Family enrichment postprocess (domain-only outputs)
# ---------------------------------------------------------------------------
echo ""
echo "[STAGE 6] Running family-enrichment postprocess..."
mkdir -p "$FAMILY_OUTPUT"

RESIDUE_SCORES="$CORE_OUTPUT/protein_condition_residue_scores.tsv"
REGIO_DELTA="$CORE_OUTPUT/protein_residue_regio_delta.tsv"

if [[ ! -f "$RESIDUE_SCORES" || ! -f "$REGIO_DELTA" ]]; then
  echo "[WARN] Residue score files missing from core output – skipping family-enrichment."
  echo "[WARN]   Expected: $RESIDUE_SCORES"
  echo "[WARN]   Expected: $REGIO_DELTA"
else
  "$PYTHON_BIN" -m lpmo_pipeline.cli family-enrichment \
    --protein-condition-residue-scores "$RESIDUE_SCORES" \
    --protein-residue-regio-delta "$REGIO_DELTA" \
    --protein-metadata "$PROTEIN_METADATA" \
    --core-fasta "$CORE_FASTA" \
    --output "$FAMILY_OUTPUT" \
    --families AA9 AA10

  echo "[STAGE 6] Family enrichment complete: $FAMILY_OUTPUT"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "[INFO] Test pipeline completed successfully"
echo "[INFO] Output root:              $OUTPUT_DIR"
echo "[INFO] core results:             $CORE_OUTPUT"
echo "[INFO] full-length results:      $FL_OUTPUT"
echo "[INFO] combined tables:          $COMBINED_OUTPUT"
echo "[INFO] predictive results:       $PREDICTIVE_OUTPUT"
echo "[INFO] cbm-paired results:       $CBM_OUTPUT"
echo "[INFO] family-enrichment:        $FAMILY_OUTPUT"
echo "================================================================"
