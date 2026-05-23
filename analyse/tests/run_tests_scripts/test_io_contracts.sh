#!/bin/bash
#SBATCH --job-name=run_tests_io_contracts
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=run_tests_io_contracts_%j.log

set -euo pipefail
# Activate conda environment
export PATH="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin:$PATH"

# Paths
project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
cd "$project_root"
export PYTHONPATH="$project_root/src:${PYTHONPATH:-}"

test_cif="/cluster/work/projects/nn1003k/eirik/Masteroppgave/structure_pipeline/work/NAG8/af3/runs/408004/A0A0A1ED04_NAG8/seed-1_sample-0/A0A0A1ED04_NAG8_seed-1_sample-0_model.cif"
out_dir="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/tests_results/"

if [[ ! -f "$test_cif" ]]; then
	echo "[ERROR] test_cif does not exist: $test_cif"
	exit 1
fi

mkdir -p "$out_dir"

job_suffix="${SLURM_JOB_ID:-manual}"
run_dir="$out_dir/io_contracts_${job_suffix}"
mkdir -p "$run_dir"

echo "[INFO] Running CCD lookup contract tests"
pytest tests/test_io_contracts.py::TestCCDLookup -v

echo "[INFO] Running real-CIF normalize contract check on: $test_cif"
export TEST_CIF="$test_cif"
export RUN_DIR="$run_dir"
python - <<'PY'
import json
import os
import re
import sys
from pathlib import Path

from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner

test_cif = Path(os.environ["TEST_CIF"])
run_dir = Path(os.environ["RUN_DIR"])
normalize_out = run_dir / "normalize_real_cif"
normalize_out.mkdir(parents=True, exist_ok=True)

ok, normalized_path = NormalizeMMCIFRunner(test_cif, normalize_out).run()
report_path = normalize_out / "normalize_report.json"

if not report_path.exists():
	print(f"[ERROR] Missing normalize report: {report_path}")
	sys.exit(1)

report = json.loads(report_path.read_text())
ccd = report.get("ccd_validation", {})

print("[INFO] Normalize result summary")
print(json.dumps(
	{
		"ok": ok,
		"normalized_path": str(normalized_path) if normalized_path else None,
		"all_valid": ccd.get("all_valid"),
		"invalid_comp_ids": ccd.get("invalid_comp_ids"),
		"observed_comp_ids": ccd.get("observed_comp_ids"),
	},
	indent=2,
))

if ok is not True:
	print("[ERROR] Normalization failed for real test CIF")
	sys.exit(1)

if ccd.get("all_valid") is not True:
	print("[ERROR] CCD validation did not pass for real test CIF")
	sys.exit(1)

if not normalized_path or not Path(normalized_path).exists():
	print("[ERROR] normalized.cif was not produced")
	sys.exit(1)

print("[INFO] Running explicit negative CCD test with invalid comp_id CEL6")
observed = ccd.get("observed_comp_ids") or []
if len(observed) == 0:
	print("[ERROR] Real test CIF has no glycan comp_id entries to mutate for negative test")
	sys.exit(1)

source_comp = str(observed[0])
invalid_cif = run_dir / "invalid_glycan_input.cif"
invalid_out = run_dir / "normalize_invalid_cif"
invalid_out.mkdir(parents=True, exist_ok=True)

raw_text = test_cif.read_text()
mutated_text, n_repl = re.subn(rf"\b{re.escape(source_comp)}\b", "CEL6", raw_text)
if n_repl == 0:
	print(f"[ERROR] Could not mutate comp_id {source_comp} to CEL6 in test CIF")
	sys.exit(1)

invalid_cif.write_text(mutated_text)
ok_invalid, normalized_invalid_path = NormalizeMMCIFRunner(invalid_cif, invalid_out).run()
invalid_report_path = invalid_out / "normalize_report.json"

if not invalid_report_path.exists():
	print(f"[ERROR] Missing normalize report for invalid CIF: {invalid_report_path}")
	sys.exit(1)

invalid_report = json.loads(invalid_report_path.read_text())
invalid_ccd = invalid_report.get("ccd_validation", {})
invalid_comp_ids = invalid_ccd.get("invalid_comp_ids") or []

print("[INFO] Invalid-case normalize summary")
print(json.dumps(
	{
		"ok": ok_invalid,
		"normalized_path": str(normalized_invalid_path) if normalized_invalid_path else None,
		"all_valid": invalid_ccd.get("all_valid"),
		"invalid_comp_ids": invalid_comp_ids,
		"mutated_from": source_comp,
		"replacement_count": n_repl,
	},
	indent=2,
))

if ok_invalid is not False:
	print("[ERROR] Negative CCD test expected normalization failure but got success")
	sys.exit(1)

if invalid_ccd.get("all_valid") is not False:
	print("[ERROR] Negative CCD test expected all_valid=false")
	sys.exit(1)

if "CEL6" not in invalid_comp_ids:
	print("[ERROR] Negative CCD test expected CEL6 in invalid_comp_ids")
	sys.exit(1)

if normalized_invalid_path is not None:
	print("[ERROR] Negative CCD test should not produce normalized.cif path")
	sys.exit(1)
PY

echo "[INFO] IO contract checks completed successfully"
echo "[INFO] Results directory: $run_dir"
