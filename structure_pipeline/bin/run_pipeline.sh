#!/bin/bash
#SBATCH --job-name=structure-pipeline
#SBATCH --account=nn1003k
#SBATCH --time=00:30:00
#SBATCH --partition=small
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# ============================================================================
# Structure Prediction Pipeline - Main Runner Script
# ============================================================================
# This script runs the manifest generation and submits prediction jobs.
#
# Usage:
#   sbatch bin/run_pipeline.sh [options]
#
# Options are passed directly to structure-pipeline run:
#   --model af3|boltz|rf3       Run only specific model
#   --protein PROTEIN_ID        Run only specific protein
#   --dry-run                   Show what would be submitted (no actual jobs)
#   --no-resume                 Rerun completed jobs
#   --ligand-ccd-list FILE.txt  Use plain CCD codes from a text file
#   --oligo-definitions FILE    YAML with oligo prefix→monomer+bond mappings
#
# Examples:
#   sbatch bin/run_pipeline.sh                             # Run all models, all proteins
#   sbatch bin/run_pipeline.sh --model rf3                 # Run only RF3
#   sbatch bin/run_pipeline.sh --dry-run                   # Preview without submitting
#   sbatch bin/run_pipeline.sh --protein A0A1Y2N3J1        # Single protein
#   sbatch bin/run_pipeline.sh --ligand-ccd-list ccd.txt   # CCD-list mode
#   sbatch bin/run_pipeline.sh --oligo-definitions config/oligo_definitions.yaml
#
# ============================================================================

set -euo pipefail

# Project directory (hardcoded for SLURM compatibility)
PROJECT_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline"
cd "$PROJECT_DIR"

# Activate conda environment (fallback to PATH if activate is missing)
ENV_ROOT="/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env"
if [[ -f "$ENV_ROOT/bin/activate" ]]; then
	source "$ENV_ROOT/bin/activate"
else
	export PATH="$ENV_ROOT/bin:$PATH"
	echo "WARNING: $ENV_ROOT/bin/activate not found; using PATH only" >&2
fi

echo "============================================"
echo "Structure Prediction Pipeline"
echo "============================================"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Project directory: $PROJECT_DIR"
echo ""

# ── Parse flags that must go to BOTH manifest and run ──
MANIFEST_EXTRA=()
RUN_EXTRA=()
i=1
while [[ $i -le $# ]]; do
	case "${!i}" in
		--ligand-ccd-list)
			MANIFEST_EXTRA+=("${!i}")
			((i++))
			if [[ $i -le $# ]]; then
				MANIFEST_EXTRA+=("${!i}")
			fi
			;;
		--oligo-definitions)
			# Only passed to 'run', not 'manifest'
			RUN_EXTRA+=("${!i}")
			((i++))
			if [[ $i -le $# ]]; then
				RUN_EXTRA+=("${!i}")
			fi
			;;
	esac
	((i++))
done

# Step 1: Generate/update manifests
echo "Step 1: Generating manifests..."
echo "---"
python -m structure_pipeline.cli manifest --config config/pipeline.yaml "${MANIFEST_EXTRA[@]}"
echo ""

# Step 2: Validate configuration
echo "Step 2: Validating configuration..."
echo "---"
python -m structure_pipeline.cli validate --config config/pipeline.yaml
echo ""

# Step 3: Submit jobs
echo "Step 3: Submitting prediction jobs..."
echo "---"
# Pass through any command-line arguments to the run command
python -m structure_pipeline.cli run --config config/pipeline.yaml "${RUN_EXTRA[@]}" "$@"

echo ""
echo "============================================"
echo "Pipeline orchestration complete"
echo "Check SLURM queue with: squeue -u $USER"
echo "============================================"
