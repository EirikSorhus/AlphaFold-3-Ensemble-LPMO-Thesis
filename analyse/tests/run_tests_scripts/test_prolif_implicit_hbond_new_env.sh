#!/bin/bash
#SBATCH --job-name=run_tests_prolif_implicit_hbond_new_env
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse/tests/run_tests_scripts/logs/run_tests_prolif_implicit_hbond_new_env_%j.log

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "[ERROR] This script must be submitted with sbatch, not run directly on the login node."
    echo "[INFO] Example: sbatch tests/run_tests_scripts/test_prolif_implicit_hbond_new_env.sh"
    exit 1
fi

project_root="/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse"
env_python="/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env/bin/python"
cd "$project_root"

run_dir="$project_root/tests/tests_results/prolif_implicit_hbond_new_env_${SLURM_JOB_ID}"
mkdir -p "$run_dir/cache"

export XDG_CACHE_HOME="$run_dir/cache"

echo "[INFO] Using Python: $env_python"
echo "[INFO] Running focused implicit H-bond probe matching tests/test_prolif_ifp.py"
RUN_DIR="$run_dir" "$env_python" - <<'PY'
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import MDAnalysis as mda
import prolif as plf
from prolif import Molecule
from prolif.io.protein_helper import ProteinHelper
from rdkit import Chem
from rdkit.Chem.rdDetermineBonds import DetermineConnectivity

project_root = Path('/cluster/work/projects/nn1003k/eirik/Masteroppgave/analyse')
run_dir = Path(Path.cwd().joinpath(Path(__import__('os').environ['RUN_DIR'])).resolve())

probe_dirs = {
    '7PXW': project_root / 'tests' / 'tests_results' / 'crystal_fil_test_1156504' / 'cases' / '0001_A0A0S2GKZ1_CEL4' / 'crystal_anchoring' / 'references' / '7PXW_A0A0S2GKZ1' / 'protonated',
    '5ACI': project_root / 'tests' / 'tests_results' / 'crystal_fil_test_1156504' / 'cases' / '0002_A0A0S2GKZ1_CEL6' / 'crystal_anchoring' / 'references' / '5ACI_A0A0S2GKZ1' / 'protonated',
}


def run_probe(complex_pdb: Path, ligand_pdb: Path):
    with tempfile.TemporaryDirectory() as tmp_dir:
        protein_pdb = Path(tmp_dir) / 'protein_only.pdb'
        universe = mda.Universe(str(complex_pdb))
        protein_atoms = universe.select_atoms('protein and chainID A')
        if protein_atoms.n_atoms == 0:
            raise RuntimeError(f'Protein selection produced zero atoms for {complex_pdb}')
        protein_atoms.write(str(protein_pdb))

        protein_rdkit = Chem.MolFromPDBFile(str(protein_pdb), removeHs=False)
        if protein_rdkit is None:
            raise RuntimeError(f'RDKit failed to load protein PDB: {protein_pdb}')
        protein_mol = Molecule.from_rdkit(protein_rdkit)
        protein_mol = ProteinHelper().standardize_protein(protein_mol)

        ligand_rdkit = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
        if ligand_rdkit is None:
            raise RuntimeError(f'RDKit failed to load ligand PDB: {ligand_pdb}')
        try:
            DetermineConnectivity(ligand_rdkit, useHueckel=True)
        except TypeError:
            DetermineConnectivity(ligand_rdkit)
        for atom in ligand_rdkit.GetAtoms():
            atom.SetNoImplicit(False)
        ligand = Molecule.from_rdkit(ligand_rdkit)

        fingerprint = plf.Fingerprint(
            interactions=['ImplicitHBAcceptor', 'ImplicitHBDonor', 'VdWContact'],
            count=True,
        )
        fingerprint.run_from_iterable([ligand], protein_mol)
        dataframe = fingerprint.to_dataframe().T

    active_contacts = int(dataframe.to_numpy().sum()) if not dataframe.empty else 0
    interaction_types = sorted({str(value) for value in dataframe.index.get_level_values(-1)})
    return active_contacts, interaction_types, dataframe


summary: dict[str, object] = {
    'job_id': __import__('os').environ.get('SLURM_JOB_ID', ''),
    'prolif_version': getattr(plf, '__version__', 'unknown'),
    'cases': {},
}

total_active_contacts = 0
observed_interactions: set[str] = set()

for label, probe_dir in probe_dirs.items():
    complex_pdb = probe_dir / 'for_posebusters.pdb'
    ligand_pdb = probe_dir / 'ligand_only_for_prolif.pdb'
    if not complex_pdb.exists() or not ligand_pdb.exists():
        raise FileNotFoundError(f'Missing probe artifacts for {label}: {complex_pdb} / {ligand_pdb}')

    active_contacts, interaction_types, dataframe = run_probe(complex_pdb, ligand_pdb)
    if active_contacts < 1:
        raise AssertionError(f'No implicit H-bond contacts detected for {label}')
    if not set(interaction_types) <= {'ImplicitHBAcceptor', 'ImplicitHBDonor', 'VdWContact'}:
        raise AssertionError(f'Unexpected interaction types for {label}: {interaction_types}')

    detail_tsv = run_dir / f'{label}_implicit_hbond_output.tsv'
    if dataframe.empty:
        detail_tsv.write_text('ligand\tprotein\tinteraction\tcount\n')
    else:
        detail_frame = dataframe.rename_axis(index=['ligand', 'protein', 'interaction'])
        detail_frame.columns = ['count']
        detail_frame.reset_index().to_csv(detail_tsv, sep='\t', index=False)

    total_active_contacts += active_contacts
    observed_interactions.update(interaction_types)
    summary['cases'][label] = {
        'complex_pdb': str(complex_pdb),
        'ligand_pdb': str(ligand_pdb),
        'active_contacts': active_contacts,
        'interaction_types': interaction_types,
        'detail_tsv': str(detail_tsv),
    }

if total_active_contacts < len(probe_dirs):
    raise AssertionError('Total implicit H-bond contacts lower than expected across probes')
if not observed_interactions:
    raise AssertionError('No implicit H-bond interaction types were observed')

summary['total_active_contacts'] = total_active_contacts
summary['observed_interactions'] = sorted(observed_interactions)

summary_path = run_dir / 'implicit_hbond_probe_summary.json'
summary_path.write_text(json.dumps(summary, indent=2) + '\n')

print(json.dumps(summary, indent=2))
print(f'[INFO] Summary written to {summary_path}')
PY

echo "[INFO] Focused implicit H-bond probe completed"
echo "[INFO] Result dir: $run_dir"
