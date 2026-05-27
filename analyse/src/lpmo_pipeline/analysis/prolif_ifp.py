# src/lpmo_pipeline/analysis/prolif_ifp.py
"""
Responsibility: Generate interaction fingerprints (IFP) using ProLIF.
Input:  non-protonated complex PDB + ligand-only PDB
Output: IFP matrix (binary: residues × interaction types)

RQ relevance:
  RQ1: IFP signatures distinguish C1 vs C4 binding modes
  RQ2: IFP per substrate type → specificity patterns
  RQ3: IFP_LPMO vs IFP_CBM in DEL B
"""
from __future__ import annotations

import json
import logging
import re
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from lpmo_pipeline.config import load_defaults_config, load_runtime_paths_config

logger = logging.getLogger(__name__)

_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
_PROLIF_CONFIG_PATH = _RUNTIME_PATHS.pipeline_assets.prolif_features_config
_DEFAULT_INTERACTION_TYPES = [
    "ImplicitHBAcceptor",
    "ImplicitHBDonor",
    "VdWContact",
]
_INTERACTION_COUNT_COLUMNS = [
    ("ImplicitHBAcceptor", "n_implicit_hbond_acceptor"),
    ("ImplicitHBDonor", "n_implicit_hbond_donor"),
    ("VdWContact", "n_vdw_contact"),
]
DEFAULT_PROTEIN_CHAIN = str(_CHAIN_SCHEMA.get("protein") or "A")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))
_DEFAULT_LIGAND_CHAIN = DEFAULT_GLYCAN_CHAINS[0] if DEFAULT_GLYCAN_CHAINS else "B"
_FEATURE_SEPARATOR = "|"


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------
@dataclass
class IFPResult:
    """Interaction fingerprint result for a single pose."""

    pose_id: str = ""
    status: str = "ok"
    error: str | None = None
    n_residues: int = 0
    n_interaction_types: int = 0
    residue_names: list[str] = field(default_factory=list)
    interaction_types: list[str] = field(default_factory=list)
    feature_names: list[str] = field(default_factory=list)
    fingerprint: list[list[int]] = field(default_factory=list)  # n_residues × n_types
    flat_bitvector: list[int] = field(default_factory=list)     # Flattened for clustering
    n_total_contacts: int = 0
    interaction_counts: dict[str, int] = field(default_factory=dict)
    interaction_occurrence_counts: dict[str, int] = field(default_factory=dict)

    def to_row(self) -> dict[str, Any]:
        """Serialize the per-pose IFP result to the planned TSV row contract."""
        row = {
            "pose_id": self.pose_id,
            "ifp_generation_status": self.status,
            "ifp_vector": json.dumps(self.flat_bitvector),
            "ifp_feature_names": json.dumps(self.feature_names),
            "n_total_contacts": self.n_total_contacts,
            "ifp_interaction_counts": json.dumps(self.interaction_counts, sort_keys=True),
            "ifp_interaction_occurrence_counts": json.dumps(
                self.interaction_occurrence_counts or self.interaction_counts,
                sort_keys=True,
            ),
            "ifp_error": self.error or "",
        }
        for interaction_name, column_name in _INTERACTION_COUNT_COLUMNS:
            row[column_name] = int(
                (self.interaction_occurrence_counts or self.interaction_counts).get(interaction_name, 0)
            )
        return row


@dataclass
class IFPBatch:
    """IFP results for multiple poses (same protein×ligand)."""

    protein_id: str = ""
    ligand_id: str = ""
    model: str = ""
    results: list[IFPResult] = field(default_factory=list)
    # Matrix: n_poses × n_features (for clustering input)
    matrix: list[list[int]] = field(default_factory=list)
    feature_names: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class ContactEligibilityRule:
    """Minimum non-vdW signal required for formal binding-mode clustering."""

    min_non_vdw_interactions: int
    min_non_vdw_contact_residues: int


@dataclass(frozen=True)
class ContactEligibility:
    """Derived clustering eligibility summary for one per-pose IFP result."""

    pose_id: str
    eligible: bool
    exclusion_class: str | None
    n_total_interactions: int
    n_vdw_interactions: int
    n_non_vdw_interactions: int
    n_non_vdw_contact_residues: int


def load_prolif_features_config(config_path: Path | None = None) -> dict[str, Any]:
    """Load the configured ProLIF interaction set and selections."""
    import yaml

    path = config_path or _PROLIF_CONFIG_PATH
    return yaml.safe_load(path.read_text()) or {}


def load_contact_eligibility_rule(
    config_path: Path | None = None,
    *,
    rule_name: str = "main_rule",
) -> ContactEligibilityRule:
    """Load one named contact-eligibility rule from prolif_features.yaml."""
    config = load_prolif_features_config(config_path)
    eligibility_config = config.get("contact_eligibility") or {}

    if rule_name == "main_rule":
        rule_config = eligibility_config.get("main_rule") or {}
    else:
        rule_config = ((eligibility_config.get("sensitivity_rules") or {}).get(rule_name) or {})

    return ContactEligibilityRule(
        min_non_vdw_interactions=int(rule_config.get("min_non_vdw_interactions", 2)),
        min_non_vdw_contact_residues=int(rule_config.get("min_non_vdw_contact_residues", 1)),
    )


def parse_ifp_feature_name(feature_name: str) -> tuple[str, str, str]:
    """Split a flattened IFP feature label into ligand, protein, and interaction."""
    parts = feature_name.split(_FEATURE_SEPARATOR)
    if len(parts) != 3:
        raise ValueError(f"Malformed IFP feature name: {feature_name}")
    ligand_residue, protein_residue, interaction_name = parts
    return ligand_residue, protein_residue, interaction_name


def evaluate_contact_eligibility(
    result: IFPResult,
    rule: ContactEligibilityRule,
    *,
    vdw_interaction_name: str = "VdWContact",
) -> ContactEligibility:
    """Classify whether a pose has enough specific non-vdW signal for clustering."""
    n_vdw_interactions = int(result.interaction_counts.get(vdw_interaction_name, 0))
    n_non_vdw_interactions = int(result.n_total_contacts) - n_vdw_interactions

    non_vdw_contact_residues: set[str] = set()
    for feature_name, value in zip(result.feature_names, result.flat_bitvector, strict=True):
        if not int(value):
            continue
        _, protein_residue, interaction_name = parse_ifp_feature_name(feature_name)
        if interaction_name == vdw_interaction_name:
            continue
        non_vdw_contact_residues.add(protein_residue)

    n_non_vdw_contact_residues = len(non_vdw_contact_residues)

    if result.n_total_contacts == 0:
        exclusion_class = "null_ifp"
        eligible = False
    elif n_non_vdw_interactions == 0 and n_vdw_interactions > 0:
        exclusion_class = "vdw_only"
        eligible = False
    elif (
        n_non_vdw_interactions < rule.min_non_vdw_interactions
        or n_non_vdw_contact_residues < rule.min_non_vdw_contact_residues
    ):
        exclusion_class = "low_specific_contact"
        eligible = False
    else:
        exclusion_class = None
        eligible = True

    return ContactEligibility(
        pose_id=result.pose_id,
        eligible=eligible,
        exclusion_class=exclusion_class,
        n_total_interactions=int(result.n_total_contacts),
        n_vdw_interactions=n_vdw_interactions,
        n_non_vdw_interactions=n_non_vdw_interactions,
        n_non_vdw_contact_residues=n_non_vdw_contact_residues,
    )


def _resolve_interaction_types(
    interaction_types: list[str] | None,
    config: dict[str, Any],
) -> list[str]:
    if interaction_types:
        return list(interaction_types)

    configured = [
        str(name)
        for name in config.get("active_interaction_types", ())
        if str(name).strip()
    ]
    return configured or list(_DEFAULT_INTERACTION_TYPES)


def _resolve_protein_selection(config: dict[str, Any], protein_chain: str) -> str:
    selection = str(
        (config.get("selection") or {}).get("protein_selection")
        or f"chainID {protein_chain}"
    ).strip()
    if "protein" not in selection.lower():
        return f"protein and ({selection})"
    return selection


def _resolve_ligand_chain(config: dict[str, Any]) -> str:
    selection = config.get("selection") or {}
    ligand_chains = selection.get("ligand_chains") or []
    if ligand_chains:
        ligand_chain = str(ligand_chains[0]).strip()
        if ligand_chain:
            return ligand_chain

    ligand_selection = str(selection.get("ligand_selection") or "").strip()
    match = re.search(r"chainID\s+([A-Za-z0-9])", ligand_selection)
    if match:
        return match.group(1)

    return _DEFAULT_LIGAND_CHAIN


def _standardize_protein_molecule(protein_rdkit: Any) -> Any:
    import prolif as plf
    from prolif.io.protein_helper import ProteinHelper

    helper = ProteinHelper()
    try:
        standardized = helper.standardize_protein(protein_rdkit)
        return plf.Molecule.from_rdkit(standardized)
    except Exception:
        protein = plf.Molecule.from_rdkit(protein_rdkit)
        return helper.standardize_protein(protein)


def _load_protein_molecule_from_complex_pdb(
    complex_pdb: Path,
    protein_selection: str,
) -> Any:
    import MDAnalysis as mda
    from rdkit import Chem
    import tempfile

    universe = mda.Universe(str(complex_pdb))
    protein_atoms = universe.select_atoms(protein_selection)
    if protein_atoms.n_atoms == 0:
        raise ValueError("Protein selection produced zero atoms")

    with tempfile.TemporaryDirectory() as tmp_dir:
        protein_pdb = Path(tmp_dir) / "protein_only.pdb"
        protein_atoms.write(str(protein_pdb))
        protein_rdkit = Chem.MolFromPDBFile(str(protein_pdb), removeHs=False)
        if protein_rdkit is None:
            raise ValueError(f"RDKit could not parse protein PDB extracted from {complex_pdb}")
        return _standardize_protein_molecule(protein_rdkit)


def _load_ligand_pdb_molecule(ligand_pdb: Path) -> Any:
    import prolif as plf
    from rdkit import Chem
    from rdkit.Chem.rdDetermineBonds import DetermineConnectivity

    ligand_rdkit = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
    if ligand_rdkit is None:
        raise ValueError(f"RDKit could not parse ligand PDB: {ligand_pdb}")
    try:
        DetermineConnectivity(ligand_rdkit, useHueckel=True)
    except TypeError:
        DetermineConnectivity(ligand_rdkit)
    for atom in ligand_rdkit.GetAtoms():
        atom.SetNoImplicit(False)
    return plf.Molecule.from_rdkit(ligand_rdkit)


def _empty_result(
    pose_id: str,
    interaction_types: list[str],
    *,
    status: str,
    error: str | None = None,
) -> IFPResult:
    return IFPResult(
        pose_id=pose_id,
        status=status,
        error=error,
        n_residues=0,
        n_interaction_types=len(interaction_types),
        residue_names=[],
        interaction_types=list(interaction_types),
        feature_names=[],
        fingerprint=[],
        flat_bitvector=[],
        n_total_contacts=0,
        interaction_counts={name: 0 for name in interaction_types},
        interaction_occurrence_counts={name: 0 for name in interaction_types},
    )


def _result_from_dataframe(
    pose_id: str,
    interaction_types: list[str],
    dataframe: Any,
) -> IFPResult:
    if getattr(dataframe, "shape", (0, 0))[1] == 0:
        return _empty_result(pose_id, interaction_types, status="zero_contacts")

    interaction_order = list(dict.fromkeys(interaction_types))
    residue_rows: dict[tuple[str, str], dict[str, int]] = {}
    residue_names: list[str] = []

    for column, value in zip(dataframe.columns.tolist(), dataframe.iloc[0].tolist()):
        ligand_label, protein_residue, interaction_name = tuple(map(str, column))
        residue_key = (ligand_label, protein_residue)
        if residue_key not in residue_rows:
            residue_rows[residue_key] = {}
            residue_names.append(
                f"{ligand_label}{_FEATURE_SEPARATOR}{protein_residue}"
            )
        if interaction_name not in interaction_order:
            interaction_order.append(interaction_name)
        residue_rows[residue_key][interaction_name] = int(value)

    feature_names: list[str] = []
    fingerprint: list[list[int]] = []
    flat_bitvector: list[int] = []
    interaction_counts = {name: 0 for name in interaction_order}
    interaction_occurrence_counts = {name: 0 for name in interaction_order}

    for residue_name in residue_names:
        ligand_label, protein_residue = residue_name.split(_FEATURE_SEPARATOR, maxsplit=1)
        row: list[int] = []
        residue_values = residue_rows[(ligand_label, protein_residue)]
        for interaction_name in interaction_order:
            occurrence_count = int(residue_values.get(interaction_name, 0))
            bit = 1 if occurrence_count > 0 else 0
            row.append(bit)
            feature_names.append(
                f"{ligand_label}{_FEATURE_SEPARATOR}{protein_residue}{_FEATURE_SEPARATOR}{interaction_name}"
            )
            flat_bitvector.append(bit)
            interaction_counts[interaction_name] += bit
            interaction_occurrence_counts[interaction_name] += occurrence_count
        fingerprint.append(row)

    return IFPResult(
        pose_id=pose_id,
        status="ok",
        n_residues=len(residue_names),
        n_interaction_types=len(interaction_order),
        residue_names=residue_names,
        interaction_types=interaction_order,
        feature_names=feature_names,
        fingerprint=fingerprint,
        flat_bitvector=flat_bitvector,
        n_total_contacts=sum(flat_bitvector),
        interaction_counts=interaction_counts,
        interaction_occurrence_counts=interaction_occurrence_counts,
    )


def _feature_value_map(result: IFPResult) -> dict[str, int]:
    return {
        feature_name: int(value)
        for feature_name, value in zip(result.feature_names, result.flat_bitvector)
    }


def filter_ifp_result_features(
    result: IFPResult,
    allowed_feature_names: set[str],
) -> IFPResult:
    """Return a copy of one IFP result restricted to selected flattened features."""

    feature_values = [
        (feature_name, int(value))
        for feature_name, value in zip(result.feature_names, result.flat_bitvector, strict=True)
        if feature_name in allowed_feature_names
    ]
    feature_names = [feature_name for feature_name, _ in feature_values]
    flat_bitvector = [value for _, value in feature_values]
    interaction_types = list(dict.fromkeys(parse_ifp_feature_name(feature_name)[2] for feature_name in feature_names))
    interaction_counts = {name: 0 for name in interaction_types}
    for feature_name, value in feature_values:
        interaction_counts[parse_ifp_feature_name(feature_name)[2]] += int(value)
    return IFPResult(
        pose_id=result.pose_id,
        status=result.status,
        error=result.error,
        n_residues=len({parse_ifp_feature_name(feature_name)[1] for feature_name in feature_names}),
        n_interaction_types=len(interaction_types),
        residue_names=list(dict.fromkeys(
            f"{parse_ifp_feature_name(feature_name)[0]}{_FEATURE_SEPARATOR}{parse_ifp_feature_name(feature_name)[1]}"
            for feature_name in feature_names
        )),
        interaction_types=interaction_types,
        feature_names=feature_names,
        fingerprint=[],
        flat_bitvector=flat_bitvector,
        n_total_contacts=sum(flat_bitvector),
        interaction_counts=interaction_counts,
        interaction_occurrence_counts=interaction_counts,
    )


def filter_ifp_batch_features(
    batch: IFPBatch,
    allowed_feature_names: set[str],
) -> IFPBatch:
    """Return an IFP batch restricted to selected flattened features."""

    results = [filter_ifp_result_features(result, allowed_feature_names) for result in batch.results]
    feature_names = [feature_name for feature_name in batch.feature_names if feature_name in allowed_feature_names]
    matrix = []
    for result in results:
        feature_map = _feature_value_map(result)
        matrix.append([int(feature_map.get(feature_name, 0)) for feature_name in feature_names])
    return IFPBatch(
        protein_id=batch.protein_id,
        ligand_id=batch.ligand_id,
        model=batch.model,
        results=results,
        matrix=matrix,
        feature_names=feature_names,
    )


# ---------------------------------------------------------------------------
# IFP generation
# ---------------------------------------------------------------------------
def compute_ifp_single(
    complex_pdb: Path,
    ligand_pdb: Path | None = None,
    pose_id: str = "",
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    interaction_types: list[str] | None = None,
    config_path: Path | None = None,
) -> IFPResult:
    """Compute ProLIF IFP for a single pose.

    Args:
        complex_pdb: Path to non-protonated complex PDB.
        ligand_pdb: Path to ligand-only PDB.
        pose_id: Identifier for this pose.
        protein_chain: Protein chain for IFP computation.
        interaction_types: Which interactions to compute.
        config_path: Optional override for the ProLIF feature config.

    Returns:
        IFPResult with binary fingerprint matrix and aligned feature names.
    """
    config = load_prolif_features_config(config_path)
    resolved_interaction_types = _resolve_interaction_types(interaction_types, config)
    ligand_input = Path(ligand_pdb) if ligand_pdb is not None else None

    logger.info("Computing IFP for pose %s: %s + %s", pose_id, complex_pdb, ligand_input)

    if ligand_input is None or not complex_pdb.exists() or not ligand_input.exists():
        missing = []
        if not complex_pdb.exists():
            missing.append(str(complex_pdb))
        if ligand_input is None:
            missing.append("<ligand_pdb>")
        elif not ligand_input.exists():
            missing.append(str(ligand_input))
        return _empty_result(
            pose_id,
            resolved_interaction_types,
            status="input_missing",
            error=f"Missing input artifact(s): {', '.join(missing)}",
        )

    try:
        import prolif as plf

        protein = _load_protein_molecule_from_complex_pdb(
            complex_pdb,
            _resolve_protein_selection(config, protein_chain),
        )
        ligand = _load_ligand_pdb_molecule(ligand_input)
        fingerprint = plf.Fingerprint(resolved_interaction_types, count=True)
        try:
            fingerprint.run_from_iterable([ligand], protein, progress=False)
        except TypeError:
            fingerprint.run_from_iterable([ligand], protein)
        dataframe = fingerprint.to_dataframe()
        result = _result_from_dataframe(pose_id, resolved_interaction_types, dataframe)
    except Exception as exc:
        logger.exception("ProLIF failed for pose %s", pose_id)
        return _empty_result(
            pose_id,
            resolved_interaction_types,
            status="prolif_error",
            error=f"{exc.__class__.__name__}: {exc}",
        )

    logger.info(
        "IFP for %s: status=%s, %d residues × %d types = %d features, %d active",
        pose_id,
        result.status,
        result.n_residues,
        result.n_interaction_types,
        len(result.flat_bitvector),
        result.n_total_contacts,
    )
    return result


def _compute_ifp_batch_item(args: tuple[dict[str, Path], str]) -> IFPResult:
    pd, protein_chain = args
    return compute_ifp_single(
        complex_pdb=pd["complex_pdb"],
        ligand_pdb=pd["ligand_pdb"],
        pose_id=pd["pose_id"],
        protein_chain=protein_chain,
    )


def compute_ifp_batch(
    pose_data: list[dict[str, Path]],
    protein_id: str = "",
    ligand_id: str = "",
    model: str = "",
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    max_workers: int | None = None,
) -> IFPBatch:
    """Compute IFP for all poses in a protein×ligand×model combination.

    Args:
        pose_data: List of {"pose_id": str, "complex_pdb": Path, "ligand_pdb": Path}
        protein_id: For metadata.
        ligand_id: For metadata.
        model: For metadata.
        protein_chain: Protein chain for IFP.

    Returns:
        IFPBatch with aligned matrix ready for clustering.
    """
    worker_count = max(1, int(max_workers or 1))
    if worker_count <= 1 or len(pose_data) <= 1:
        results = [_compute_ifp_batch_item((pd, protein_chain)) for pd in pose_data]
    else:
        with ProcessPoolExecutor(max_workers=min(worker_count, len(pose_data))) as executor:
            results = list(
                executor.map(
                    _compute_ifp_batch_item,
                    [(pd, protein_chain) for pd in pose_data],
                )
            )

    # Align all fingerprints to the same feature set
    all_feature_names: list[str] = []
    if results:
        feature_set: set[str] = set()
        for r in results:
            feature_set.update(r.feature_names)
        all_feature_names = sorted(feature_set)

    matrix: list[list[int]] = []
    for r in results:
        feature_map = _feature_value_map(r)
        row = [int(feature_map.get(feature_name, 0)) for feature_name in all_feature_names]
        matrix.append(row)

    batch = IFPBatch(
        protein_id=protein_id,
        ligand_id=ligand_id,
        model=model,
        results=results,
        matrix=matrix,
        feature_names=all_feature_names,
    )

    logger.info(
        "IFP batch: %s × %s × %s → %d poses × %d features",
        protein_id, ligand_id, model, len(results), len(all_feature_names),
    )
    return batch


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_ifp_matrix(batch: IFPBatch, output_path: Path) -> None:
    """Write IFP matrix to CSV (columns = features, rows = poses)."""
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pose_id"] + batch.feature_names)
        for r, row in zip(batch.results, batch.matrix):
            writer.writerow([r.pose_id] + row)
    logger.info("Wrote IFP matrix (%d×%d) to %s",
                len(batch.matrix), len(batch.feature_names), output_path)


def write_pose_ifp_table(results: list[IFPResult], output_path: Path) -> None:
    """Write the canonical per-pose IFP table as TSV."""
    import csv

    fieldnames = [
        "pose_id",
        "ifp_generation_status",
        "ifp_vector",
        "ifp_feature_names",
        "n_total_contacts",
        *[column_name for _, column_name in _INTERACTION_COUNT_COLUMNS],
        "ifp_interaction_counts",
        "ifp_interaction_occurrence_counts",
        "ifp_error",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow(result.to_row())

    logger.info("Wrote pose IFP table (%d rows) to %s", len(results), output_path)


def compute_tanimoto_similarity(bitvec_a: list[int], bitvec_b: list[int]) -> float:
    """Tanimoto similarity between two binary fingerprints.

    T(A,B) = |A∩B| / |A∪B|
    """
    a = np.array(bitvec_a, dtype=bool)
    b = np.array(bitvec_b, dtype=bool)
    intersection = np.sum(a & b)
    union = np.sum(a | b)
    if union == 0:
        return 0.0
    return float(intersection / union)
