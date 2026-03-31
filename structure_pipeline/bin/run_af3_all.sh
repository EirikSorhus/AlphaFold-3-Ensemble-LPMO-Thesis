#!/bin/bash
#SBATCH --job-name=af3-pipeline
#SBATCH --account=nn1003k
#SBATCH --time=00:30:00
#SBATCH --partition=small
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/af3_pipeline_%j.out
#SBATCH --error=logs/af3_pipeline_%j.err

# ============================================================================
# AF3 Pipeline Runner - All Proteins and Ligands
# ============================================================================
# Runs AlphaFold 3 for all proteins and ligands with fixed parameters:
#   seeds               = 1-15
#   num_diffusion_samples = 5
#   num_recycles        = 10
#
# Usage:
#   sbatch bin/run_af3_all.sh [options]
#
# Options are passed through to structure-pipeline run:
#   --dry-run                   Show what would be submitted (no actual jobs)
#   --no-resume                 Rerun already completed jobs
#   --protein PROTEIN_ID        Limit to a single protein
#   --ligand-ccd-list FILE.txt  Use plain CCD codes from a text file
#   --oligo-definitions FILE    YAML with oligo prefix→monomer+bond mappings
#
# Examples:
#   sbatch bin/run_af3_all.sh
#   sbatch bin/run_af3_all.sh --dry-run
#   sbatch bin/run_af3_all.sh --protein A0A1Y2N3J1
#   sbatch bin/run_af3_all.sh --ligand-ccd-list ccd.txt
#
# ============================================================================

set -euo pipefail

# ── Fixed AF3 parameters ───────────────────────────────────────────────────
AF3_SEEDS=15
AF3_NUM_DIFFUSION_SAMPLES=5
AF3_NUM_RECYCLES=10

# ── Project setup ──────────────────────────────────────────────────────────
PROJECT_DIR="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline"
cd "$PROJECT_DIR"

ENV_ROOT="/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env"
if [[ -f "$ENV_ROOT/bin/activate" ]]; then
	source "$ENV_ROOT/bin/activate"
else
	export PATH="$ENV_ROOT/bin:$PATH"
	echo "WARNING: $ENV_ROOT/bin/activate not found; using PATH only" >&2
fi

echo "============================================"
echo "AF3 Pipeline - All Proteins & Ligands"
echo "============================================"
echo "Date:                  $(date)"
echo "Node:                  $(hostname)"
echo "Seeds:                 1-${AF3_SEEDS}"
echo "num_diffusion_samples: ${AF3_NUM_DIFFUSION_SAMPLES}"
echo "num_recycles:          ${AF3_NUM_RECYCLES}"
echo ""

mkdir -p logs

# ── Parse flags that must go to BOTH manifest and run ──────────────────────
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
python -m structure_pipeline.cli manifest \
	--config config/pipeline.yaml \
	--force \
	"${MANIFEST_EXTRA[@]}"
echo ""

# Step 2: Validate configuration
echo "Step 2: Validating configuration..."
echo "---"
python -m structure_pipeline.cli validate --config config/pipeline.yaml || {
	echo "WARNING: Validation reported issues (see above). Continuing anyway..." >&2
}
echo ""

# Step 3: Submit AF3 jobs
echo "Step 3: Submitting AF3 jobs..."
echo "---"
python -m structure_pipeline.cli run \
	--config config/pipeline.yaml \
	--model af3 \
	--af3-seeds "${AF3_SEEDS}" \
	--af3-diffusion-samples "${AF3_NUM_DIFFUSION_SAMPLES}" \
	--af3-num-recycles "${AF3_NUM_RECYCLES}" \
	"${RUN_EXTRA[@]}" \
	"$@"

echo ""
echo "============================================"
echo "AF3 pipeline submission complete"
echo "Check SLURM queue with: squeue -u $USER"
echo "============================================"
