#!/bin/bash
#SBATCH --job-name=run_tests_protonation_contracts
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_protonation_contracts_%j.log

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
run_dir="$out_dir/protonation_contracts_${job_suffix}"
mkdir -p "$run_dir"

echo "[INFO] Running protonation/export contract generation on real AF3 CIF"
export TEST_CIF="$test_cif"
export RUN_DIR="$run_dir"
python - <<'PY'
import json
import os
import sys
from pathlib import Path

from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner
from lpmo_pipeline.io.protonate_export import protonate_and_export, _count_atom_lines


def fail(msg: str) -> None:
    print(f"[ERROR] {msg}")
    sys.exit(1)


test_cif = Path(os.environ["TEST_CIF"])
run_dir = Path(os.environ["RUN_DIR"])

# -------------------------------------------------------------------------
# Step 1 — normalize
# -------------------------------------------------------------------------
normalize_out = run_dir / "normalize_for_protonation"
normalize_out.mkdir(parents=True, exist_ok=True)

ok_norm, normalized_path = NormalizeMMCIFRunner(test_cif, normalize_out).run()
if ok_norm is not True:
    fail("Normalization failed before protonation/export")
if not normalized_path:
    fail("Normalization did not return normalized.cif path")

normalized_path = Path(normalized_path)
if not normalized_path.exists():
    fail(f"normalized.cif missing: {normalized_path}")

# -------------------------------------------------------------------------
# Step 2 — protonate/export
# -------------------------------------------------------------------------
protonation_out = run_dir / "protonated"
protonation_out.mkdir(parents=True, exist_ok=True)

ok_prot, report = protonate_and_export(normalized_path, protonation_out)
report_path = protonation_out / "protonation_report.json"
if not report_path.exists():
    fail(f"Missing protonation_report.json: {report_path}")

data = json.loads(report_path.read_text())

for_posebusters = protonation_out / "for_posebusters.pdb"
complex_h = protonation_out / "complex_H.pdb"
ligand_mol2 = protonation_out / "ligand_for_prolif.mol2"

if not for_posebusters.exists():
    fail(f"Missing for_posebusters.pdb: {for_posebusters}")
if not complex_h.exists():
    fail(f"Missing complex_H.pdb: {complex_h}")
if not ligand_mol2.exists():
    fail(f"Missing ligand_for_prolif.mol2: {ligand_mol2}")

# -------------------------------------------------------------------------
# Check: cif_to_pdb backend
# -------------------------------------------------------------------------
cif_to_pdb_report_path = protonation_out / "cif_to_pdb_report.json"
if not cif_to_pdb_report_path.exists():
    fail(f"Missing cif_to_pdb_report.json: {cif_to_pdb_report_path}")

cif_report = json.loads(cif_to_pdb_report_path.read_text())
cif_backend = cif_report.get("backend", "unknown")
print(f"[INFO] cif_to_pdb backend: {cif_backend}")
if cif_report.get("backend_fallback_reason"):
    print(f"[WARN] cif_to_pdb fallback reason: {cif_report['backend_fallback_reason']}")

# -------------------------------------------------------------------------
# Check: complex_H backend — must not be "none" (no silent copy allowed)
# -------------------------------------------------------------------------
complex_h_backend = data.get("complex_h_backend", "none")
print(f"[INFO] complex_H backend: {complex_h_backend}")
if complex_h_backend == "none":
    fail(
        f"complex_H.pdb protonation failed (backend=none). "
        f"warnings={data.get('warnings')}"
    )
if data.get("blockers"):
    fail(f"Protonation blockers: {data['blockers']}")

# -------------------------------------------------------------------------
# Check: complex_H.pdb has more atoms than for_posebusters.pdb (hydrogens added)
# -------------------------------------------------------------------------
n_base = _count_atom_lines(for_posebusters)
n_h = _count_atom_lines(complex_h)
print(f"[INFO] for_posebusters.pdb atoms: {n_base}")
print(f"[INFO] complex_H.pdb atoms:       {n_h}")
if n_h <= n_base:
    fail(
        f"complex_H.pdb must have more atoms than for_posebusters.pdb "
        f"(base={n_base}, complex_h={n_h}, backend={complex_h_backend}). "
        "Hydrogens were not added."
    )

# -------------------------------------------------------------------------
# Check: MOL2 has required TRIPOS sections
# -------------------------------------------------------------------------
mol2_text = ligand_mol2.read_text(errors="ignore")
if "@<TRIPOS>ATOM" not in mol2_text:
    fail("ligand_for_prolif.mol2 missing @<TRIPOS>ATOM section")
if "@<TRIPOS>BOND" not in mol2_text:
    fail("ligand_for_prolif.mol2 missing @<TRIPOS>BOND section")

# -------------------------------------------------------------------------
# Check: MOL2 is ligand-only (no protein residue labels)
# -------------------------------------------------------------------------
PROTEIN_RESIDUES = {
    "HIS", "GLY", "TRP", "SER", "ALA", "LEU", "ILE", "VAL", "PHE",
    "PRO", "MET", "CYS", "ASN", "ASP", "GLN", "GLU", "LYS", "ARG",
    "THR", "TYR",
}
protein_hits = [r for r in PROTEIN_RESIDUES if r in mol2_text]
if protein_hits:
    fail(f"ligand_for_prolif.mol2 contains protein residue labels: {protein_hits}")

# -------------------------------------------------------------------------
# Check: MOL2 must be materially smaller than the full complex PDB
# -------------------------------------------------------------------------
n_mol2_atom_lines = sum(
    1 for line in mol2_text.splitlines()
    if line.strip() and not line.startswith(("@", "#"))
)
if n_mol2_atom_lines >= n_base:
    print(
        f"[WARN] MOL2 non-header lines ({n_mol2_atom_lines}) >= PDB atoms ({n_base}). "
        "Ligand may not be correctly isolated."
    )

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
summary = {
    "ok_normalize": ok_norm,
    "ok_protonate_export": ok_prot,
    "test_cif": str(test_cif),
    "normalized_cif": str(normalized_path),
    "cif_to_pdb_backend": cif_backend,
    "cif_to_pdb_fallback_reason": cif_report.get("backend_fallback_reason", ""),
    "for_posebusters_pdb": str(for_posebusters),
    "for_posebusters_atom_count": n_base,
    "complex_h_pdb": str(complex_h),
    "complex_h_atom_count": n_h,
    "complex_h_backend": complex_h_backend,
    "hydrogen_atoms_added": n_h - n_base,
    "ligand_for_prolif_mol2": str(ligand_mol2),
    "ligand_mol2_size_bytes": ligand_mol2.stat().st_size,
    "protonation_report": str(report_path),
    "blockers": data.get("blockers"),
    "warnings": data.get("warnings"),
    "checks_passed": {
        "cif_to_pdb_backend_known": cif_backend != "unknown",
        "complex_h_backend_not_none": complex_h_backend != "none",
        "hydrogens_added": n_h > n_base,
        "mol2_tripos_atom": "@<TRIPOS>ATOM" in mol2_text,
        "mol2_tripos_bond": "@<TRIPOS>BOND" in mol2_text,
        "mol2_no_protein_residues": len(protein_hits) == 0,
    },
}

summary_path = run_dir / "protonation_contract_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))

print("[INFO] Protonation contract summary")
print(json.dumps(summary, indent=2))
print("[INFO] All checks passed.")
PY

echo "[INFO] Protonation contract generation completed successfully"
echo "[INFO] Results directory: $run_dir"