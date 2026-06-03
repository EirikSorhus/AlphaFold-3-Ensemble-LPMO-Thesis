# src/lpmo_pipeline/analysis/crystal_anchoring.py
"""Crystal reference preparation and comparison for real deposited structures.

This module now covers the first working crystal-anchoring slice:
    - resolve all crystal references for one protein from input_data/pdb_structure_data.csv
    - ignore the non-authoritative CSV columns Oligo_Activity and Comment
    - prepare multichain deposited crystal mmCIFs by selecting the relevant
        protein copy, the associated carbohydrate ligand, and Cu
    - fall back from preferred chain A to another protein chain if A is not the
        ligand-bound copy in a multimeric crystal
    - compare representative-pose IFPs against crystal IFPs using feature-aligned
        Tanimoto instead of assuming identical raw bitvector layouts
    - compute a first working local pocket RMSD via ligand-proximal C-alpha
        superposition on the representative-pose-vs-crystal slice

PyMOL pair_fit hardening/parity still remains for later follow-up.
"""

import json
import logging
import csv
import math
import re
from difflib import SequenceMatcher
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from lpmo_pipeline.analysis.prolif_ifp import (
    ContactEligibility,
    IFPBatch,
    IFPResult,
    compute_ifp_single,
    compute_tanimoto_similarity,
    evaluate_contact_eligibility,
    load_contact_eligibility_rule,
    write_ifp_matrix,
    write_pose_ifp_table,
)
from lpmo_pipeline.analysis.mdanalysis_metrics import (
    POSE_GEOMETRY_COLUMNS,
    PoseGeometryMetrics,
    compute_pose_metrics,
)
from lpmo_pipeline.config import load_defaults_config
from lpmo_pipeline.io.analysis_export import export_analysis_artifacts
from lpmo_pipeline.io.gemmi_compat import gemmi
from lpmo_pipeline.io.normalize_mmcif import NormalizeMMCIFRunner

logger = logging.getLogger(__name__)

_DEFAULTS_CONFIG = load_defaults_config()
_CHAIN_SCHEMA = _DEFAULTS_CONFIG.get("chain_schema") or {}
DEFAULT_PROTEIN_CHAIN = str(_CHAIN_SCHEMA.get("protein") or "A")
DEFAULT_GLYCAN_CHAINS = tuple(str(chain) for chain in (_CHAIN_SCHEMA.get("glycans") or ["B", "C", "D"]))
DEFAULT_LIGAND_CHAIN = DEFAULT_GLYCAN_CHAINS[0] if DEFAULT_GLYCAN_CHAINS else "B"
DEFAULT_CU_CHAIN = str(_CHAIN_SCHEMA.get("metal") or "E")

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_INDEX_CSV = _PROJECT_ROOT / "input_data" / "pdb_structure_data.csv"
DEFAULT_CRYSTAL_ROOT = _PROJECT_ROOT / "crystal_structures"


CRYSTAL_SIM_MODERATE_THRESHOLD = 0.5  # Empirical; flag below this
POCKET_RMSD_THRESHOLD = 2.5  # Å
DEFAULT_LIGAND_DISTANCE_CUTOFF_A = 6.0
DEFAULT_COPPER_DISTANCE_CUTOFF_A = 4.0

_PROTEIN_RESIDUE_NAMES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "MSE", "SEC", "HIC", "PYL",
    }
)
_COPPER_MONOMERS = frozenset({"CU", "CU1", "CU2"})
_NONPOLY_EXCLUDED_MONOMERS = frozenset(
    {
        "HOH", "WAT", "DOD", "CL", "NA", "K", "CA", "MG", "ZN", "MN",
        "SO4", "PO4", "ACT", "FMT", "EOH", "PEG", "GOL", "MPD", "EDO",
        "MES", "TRS",
    }
)
_PROTEIN_RESIDUE_ALIASES = {
    "HIC": "HIS",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "MSE": "MET",
}
_NORMALIZED_PROTEIN_RESIDUE_NAMES = frozenset(
    _PROTEIN_RESIDUE_ALIASES.get(name, name) for name in _PROTEIN_RESIDUE_NAMES
)


@dataclass(frozen=True)
class CrystalReferenceRecord:
    """One crystal reference row resolved to an on-disk mmCIF file."""

    family: str
    protein_id: str
    pdb_code: str
    source_cif: Path
    carbohydrate_ligands: str = ""
    dp: int | None = None
    resolution: str = ""
    preferred_author_chain: str = DEFAULT_PROTEIN_CHAIN
    expected_ligand: bool = False


@dataclass(frozen=True)
class CrystalReferenceSelection:
    """Chosen crystal binding-site view for one deposited structure."""

    source_cif: Path
    preferred_protein_chain: str
    selected_protein_chain: str
    ligand_chain_ids: tuple[str, ...] = ()
    copper_chain_ids: tuple[str, ...] = ()
    expected_ligand: bool = False
    preferred_chain_has_ligand: bool = False
    selected_chain_has_ligand: bool = False
    used_fallback_protein_chain: bool = False
    selected_ligand_distance_a: float | None = None


@dataclass(frozen=True)
class PreparedCrystalReference:
    """Prepared subset/reference artifacts for one crystal structure."""

    record: CrystalReferenceRecord
    selection: CrystalReferenceSelection
    subset_cif: Path
    normalized_cif: Path | None = None
    complex_pdb: Path | None = None
    ligand_pdb: Path | None = None
    ifp_artifact_dir: Path | None = None
    ifp_result_json: Path | None = None
    pose_ifp_table_tsv: Path | None = None
    ifp_matrix_csv: Path | None = None
    ifp_result: IFPResult | None = None
    ifp_contact_eligibility: ContactEligibility | None = None
    geometry_metrics: PoseGeometryMetrics | None = None
    status: str = "prepared"
    error: str | None = None


@dataclass
class CrystalReferencePoseComparison:
    """Representative-pose comparison against one crystal reference."""

    pdb_code: str
    source_cif: str
    prepared_subset_cif: str = ""
    prepared_normalized_cif: str = ""
    selected_protein_chain: str = ""
    ligand_chain_ids: list[str] = field(default_factory=list)
    copper_chain_ids: list[str] = field(default_factory=list)
    preferred_protein_chain: str = DEFAULT_PROTEIN_CHAIN
    preferred_chain_has_ligand: bool = False
    used_fallback_protein_chain: bool = False
    representative_pose_id: str = ""
    representative_pose_ifp_status: str = ""
    representative_pose_ifp_result_json: str = ""
    representative_pose_pose_ifp_table_tsv: str = ""
    representative_pose_ifp_matrix_csv: str = ""
    crystal_ifp_status: str = ""
    crystal_ifp_result_json: str = ""
    crystal_pose_ifp_table_tsv: str = ""
    crystal_ifp_matrix_csv: str = ""
    crystal_ifp_contact_eligible: bool = False
    crystal_ifp_exclusion_class: str = ""
    crystal_n_vdw_interactions: int = 0
    crystal_n_non_vdw_interactions: int = 0
    crystal_n_non_vdw_contact_residues: int = 0
    ifp_comparison_eligible: bool = False
    crystal_geometry: dict[str, Any] = field(default_factory=dict)
    pocket_residues: list[int] = field(default_factory=list)
    pocket_rmsd: float | None = None
    pocket_rmsd_below_threshold: bool = False
    ifp_tanimoto: float | None = None
    status: str = "not_compared"
    error: str | None = None


CRYSTAL_GEOMETRY_COLUMNS = [
    "protein_id",
    "pdb_code",
    "source_cif",
    "prepared_subset_cif",
    "selected_protein_chain",
    "ligand_chain_ids",
    "copper_chain_ids",
    "geometry_status",
    *[
        column
        for column in POSE_GEOMETRY_COLUMNS
        if column not in {"pose_id", "model", "protein_id", "ligand_id"}
    ],
]


@dataclass
class CrystalReferenceScreenReport:
    """All crystal comparisons for one representative pose."""

    protein_id: str
    representative_pose_id: str
    representative_pose_cif: str
    ligand_id: str = ""
    representative_pose_ifp_result_json: str = ""
    representative_pose_pose_ifp_table_tsv: str = ""
    representative_pose_ifp_matrix_csv: str = ""
    comparisons: list[CrystalReferencePoseComparison] = field(default_factory=list)


@dataclass
class CrystalComparisonResult:
    """Comparison of a single cluster against crystal."""

    cluster_id: int
    medoid_pose_id: str = ""
    ifp_tanimoto: float = 0.0
    pocket_rmsd: float | None = None
    ifp_above_threshold: bool = False
    pocket_rmsd_below_threshold: bool = False


@dataclass
class CrystalAnchoringReport:
    """Full crystal anchoring report for a protein×ligand pair."""

    protein_id: str = ""
    ligand_id: str = ""
    crystal_pdb: str = ""
    has_crystal_ligand: bool = False
    comparisons: list[CrystalComparisonResult] = field(default_factory=list)
    best_cluster_id: int = -1
    best_tanimoto: float = 0.0


def _normalize_reference_protein_id(value: str) -> str:
    tokens = [token.strip() for token in str(value).split() if token.strip()]
    return "__".join(tokens)


def _parse_pdb_code(value: str) -> tuple[str, str | None]:
    text = str(value).strip()
    match = re.match(r"^([A-Za-z0-9]{4})(?:\[([^\]]+)\])?$", text)
    if not match:
        raise ValueError(f"Malformed PDB field in crystal reference index: {value!r}")
    return match.group(1).upper(), match.group(2)


def _parse_optional_int(value: str) -> int | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def load_crystal_reference_records(
    protein_id: str,
    *,
    reference_index_csv: Path = DEFAULT_REFERENCE_INDEX_CSV,
    crystal_root: Path = DEFAULT_CRYSTAL_ROOT,
) -> list[CrystalReferenceRecord]:
    """Load every crystal reference row for one protein.

    The CSV is authoritative for protein→PDB membership, while the columns
    Oligo_Activity and Comment are intentionally ignored.
    """
    protein_key = _normalize_reference_protein_id(protein_id)
    records: list[CrystalReferenceRecord] = []

    with reference_index_csv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        for row in reader:
            row_protein_id = _normalize_reference_protein_id(row.get("Uniprot_ID", ""))
            if row_protein_id != protein_key:
                continue

            pdb_code, chain_hint = _parse_pdb_code(str(row.get("PDB", "") or ""))
            crystal_dir = crystal_root / protein_key
            ligand_cif = crystal_dir / f"{pdb_code}_{protein_key}_ligand.cif"
            apo_cif = crystal_dir / f"{pdb_code}_{protein_key}_no_ligand.cif"

            if ligand_cif.exists():
                source_cif = ligand_cif
            elif apo_cif.exists():
                source_cif = apo_cif
            else:
                logger.warning(
                    "Skipping crystal reference %s for %s because no matching mmCIF exists in %s",
                    pdb_code,
                    protein_key,
                    crystal_dir,
                )
                continue

            carbohydrate_ligands = str(row.get("Carbohydrate_Ligands", "") or "").strip()
            records.append(
                CrystalReferenceRecord(
                    family=str(row.get("Family", "") or "").strip(),
                    protein_id=protein_key,
                    pdb_code=pdb_code,
                    source_cif=source_cif.resolve(),
                    carbohydrate_ligands=carbohydrate_ligands,
                    dp=_parse_optional_int(str(row.get("DP", "") or "")),
                    resolution=str(row.get("Resolution", "") or "").strip(),
                    preferred_author_chain=str(chain_hint or DEFAULT_PROTEIN_CHAIN).strip() or DEFAULT_PROTEIN_CHAIN,
                    expected_ligand=source_cif.name.endswith("_ligand.cif") or bool(carbohydrate_ligands),
                )
            )

    records.sort(key=lambda record: (record.pdb_code, str(record.source_cif)))
    return records


def _parse_ligand_id_components(ligand_id: str) -> tuple[str, int | None]:
    text = str(ligand_id).strip().upper()
    if not text:
        return "", None
    match = re.match(r"^(?P<ligand>.*?)(?P<dp>\d+)$", text)
    if match is None:
        return text, None
    return match.group("ligand"), int(match.group("dp"))


def filter_crystal_reference_records(
    records: Iterable[CrystalReferenceRecord],
    *,
    ligand_id: str = "",
    require_expected_ligand: bool = True,
) -> list[CrystalReferenceRecord]:
    """Return the subset of crystal references matching one ligand condition.

    The ligand condition uses the AF3 naming convention like ``CEL4`` or ``NAG6``.
    Crystal rows are matched on both carbohydrate ligand label and DP so different
    holo references for the same protein stay separable.
    """
    ligand_code, dp = _parse_ligand_id_components(ligand_id)
    filtered: list[CrystalReferenceRecord] = []

    for record in records:
        if require_expected_ligand and not record.expected_ligand:
            continue
        if ligand_code and str(record.carbohydrate_ligands).strip().upper() != ligand_code:
            continue
        if dp is not None and record.dp != dp:
            continue
        filtered.append(record)

    filtered.sort(key=lambda record: (record.pdb_code, str(record.source_cif)))
    return filtered


def _unique_nonempty(values: Iterable[str]) -> tuple[str, ...]:
    ordered: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value).strip()
        if not text or text in {"?", "."} or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return tuple(ordered)


def _read_reference_block(source_cif: Path) -> Any:
    return gemmi.cif.read(str(source_cif)).sole_block()


def _find_branch_ligand_chain_ids(block: Any) -> tuple[str, ...]:
    return _unique_nonempty(block.find_values("_pdbx_branch_scheme.asym_id"))


def _resolve_structure_chain_id(
    candidates: Iterable[Any],
    available_chain_ids: set[str],
) -> str | None:
    normalized_candidates = _unique_nonempty(str(candidate) for candidate in candidates)
    for chain_id in normalized_candidates:
        if chain_id in available_chain_ids:
            return chain_id
    return normalized_candidates[0] if normalized_candidates else None


def _find_branch_ligand_chain_mappings(
    block: Any,
    available_chain_ids: set[str],
) -> dict[str, str]:
    table = block.find(
        "_pdbx_branch_scheme.",
        ["asym_id", "pdb_asym_id", "auth_asym_id"],
    )
    if not table:
        return {}

    chain_mappings: dict[str, str] = {}
    for source_chain_id, pdb_chain_id, auth_chain_id in table:
        source_chain_id = str(source_chain_id).strip()
        if not source_chain_id or source_chain_id in chain_mappings:
            continue
        structure_chain_id = _resolve_structure_chain_id(
            (source_chain_id, auth_chain_id, pdb_chain_id),
            available_chain_ids,
        )
        if structure_chain_id:
            chain_mappings[source_chain_id] = structure_chain_id
    return chain_mappings


def _find_nonpoly_chain_ids_for_monomers(block: Any, allowed_monomers: set[str]) -> tuple[str, ...]:
    table = block.find("_pdbx_nonpoly_scheme.", ["asym_id", "mon_id"])
    if not table:
        return ()
    return _unique_nonempty(
        asym_id
        for asym_id, mon_id in table
        if str(mon_id).strip().upper() in allowed_monomers
    )


def _find_nonpoly_carbohydrate_chain_ids(block: Any) -> tuple[str, ...]:
    table = block.find("_pdbx_nonpoly_scheme.", ["asym_id", "mon_id"])
    if not table:
        return ()
    return _unique_nonempty(
        asym_id
        for asym_id, mon_id in table
        if str(mon_id).strip().upper() not in _NONPOLY_EXCLUDED_MONOMERS
        and str(mon_id).strip().upper() not in _COPPER_MONOMERS
    )


def _normalize_comp_id(value: Any) -> str:
    return str(value).strip().upper()


def _is_protein_comp_id(value: Any) -> bool:
    return _normalize_comp_id(value) in _PROTEIN_RESIDUE_NAMES


def _is_copper_comp_id(value: Any) -> bool:
    return _normalize_comp_id(value) in _COPPER_MONOMERS


def _is_ligand_comp_id(value: Any) -> bool:
    comp_id = _normalize_comp_id(value)
    return (
        bool(comp_id)
        and comp_id not in _PROTEIN_RESIDUE_NAMES
        and comp_id not in _NONPOLY_EXCLUDED_MONOMERS
        and comp_id not in _COPPER_MONOMERS
    )


def _atom_element_name(atom: Any) -> str:
    element = getattr(atom, "element", None)
    name = getattr(element, "name", "") if element is not None else ""
    text = str(name).strip()
    if text:
        return text.upper()
    return str(getattr(atom, "name", "") or "").strip()[:1].upper()


def _atom_xyz(atom: Any) -> tuple[float, float, float]:
    pos = getattr(atom, "pos", None)
    if pos is not None:
        return float(pos.x), float(pos.y), float(pos.z)
    return (
        float(getattr(atom, "x", 0.0)),
        float(getattr(atom, "y", 0.0)),
        float(getattr(atom, "z", 0.0)),
    )


def _atom_b_iso(atom: Any) -> float:
    return float(getattr(atom, "b_iso", 0.0))


def _chain_looks_like_protein(chain: Any) -> bool:
    protein_like = 0
    total = 0
    for residue in chain:
        total += 1
        if str(getattr(residue, "name", "") or "").strip().upper() in _PROTEIN_RESIDUE_NAMES:
            protein_like += 1

    if protein_like >= 20:
        return True
    return protein_like >= 3 and protein_like == total


def _find_chain(model: Any, chain_id: str) -> Any | None:
    for chain in model:
        if str(getattr(chain, "name", "") or "") == chain_id:
            return chain
    return None


def _iter_heavy_atoms(chain: Any) -> Iterable[Any]:
    for residue in chain:
        for atom in residue:
            if _atom_element_name(atom) == "H":
                continue
            yield atom


def _min_chain_distance(chain_a: Any, chain_b: Any) -> float | None:
    atoms_a = list(_iter_heavy_atoms(chain_a))
    atoms_b = list(_iter_heavy_atoms(chain_b))
    if not atoms_a or not atoms_b:
        return None

    best = float("inf")
    for atom_a in atoms_a:
        ax, ay, az = _atom_xyz(atom_a)
        for atom_b in atoms_b:
            bx, by, bz = _atom_xyz(atom_b)
            distance = math.dist((ax, ay, az), (bx, by, bz))
            if distance < best:
                best = distance
    return None if not math.isfinite(best) else best


def _select_copper_chain_ids(block: Any) -> tuple[str, ...]:
    return _find_nonpoly_chain_ids_for_monomers(block, set(_COPPER_MONOMERS))


def _struct_conn_partner_chain_for_copper_site(
    row: Iterable[Any],
    tags: list[str],
    selected_protein_chain: str,
) -> str | None:
    values = list(row)
    index_by_tag = {tag: index for index, tag in enumerate(tags)}

    def _value(tag: str) -> str:
        index = index_by_tag.get(tag)
        if index is None:
            return ""
        return str(values[index]).strip()

    partners: list[dict[str, str]] = []
    for partner_index in (1, 2, 3):
        label_asym = _value(f"_struct_conn.ptnr{partner_index}_label_asym_id")
        label_comp = _value(f"_struct_conn.ptnr{partner_index}_label_comp_id")
        auth_asym = _value(f"_struct_conn.ptnr{partner_index}_auth_asym_id")
        auth_comp = _value(f"_struct_conn.ptnr{partner_index}_auth_comp_id")
        if not any((label_asym, label_comp, auth_asym, auth_comp)):
            continue
        partners.append(
            {
                "label_asym_id": label_asym,
                "label_comp_id": label_comp,
                "auth_asym_id": auth_asym,
                "auth_comp_id": auth_comp,
            }
        )

    copper_partners = [
        partner
        for partner in partners
        if _is_copper_comp_id(partner["label_comp_id"] or partner["auth_comp_id"])
    ]
    if not copper_partners:
        return None

    protein_matches = [
        partner
        for partner in partners
        if not _is_copper_comp_id(partner["label_comp_id"] or partner["auth_comp_id"])
        and (
            partner["label_asym_id"] == selected_protein_chain
            or partner["auth_asym_id"] == selected_protein_chain
        )
    ]
    if not protein_matches:
        return None

    return copper_partners[0]["label_asym_id"] or copper_partners[0]["auth_asym_id"] or None


def _select_copper_chain_ids_for_site(
    block: Any,
    model: Any,
    selected_protein_chain: str,
) -> tuple[str, ...]:
    """Select Cu chains for the chosen protein copy, not every deposited Cu."""
    struct_conn_tags, struct_conn_rows = _loop_rows_for_tag(block, "_struct_conn.id")
    if struct_conn_rows:
        chain_ids = _unique_nonempty(
            chain_id
            for row in struct_conn_rows
            for chain_id in [
                _struct_conn_partner_chain_for_copper_site(
                    row,
                    struct_conn_tags,
                    selected_protein_chain,
                )
            ]
            if chain_id
        )
        if chain_ids:
            return chain_ids

    nonpoly_tags, nonpoly_rows = _loop_rows_for_tag(
        block,
        "_pdbx_nonpoly_scheme.asym_id",
        fallback_prefix="_pdbx_nonpoly_scheme.",
        fallback_columns=[
            "asym_id",
            "entity_id",
            "mon_id",
            "ndb_seq_num",
            "pdb_seq_num",
            "auth_seq_num",
            "pdb_mon_id",
            "auth_mon_id",
            "pdb_strand_id",
            "pdb_ins_code",
        ],
    )
    if nonpoly_rows:
        chain_ids = _unique_nonempty(
            _row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.asym_id")
            for row in nonpoly_rows
            if _is_copper_comp_id(_row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.mon_id"))
            and selected_protein_chain
            in {
                _row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.auth_asym_id"),
                _row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.pdb_strand_id"),
            }
        )
        if chain_ids:
            return chain_ids

    protein_chain = _find_chain(model, selected_protein_chain)
    copper_chain_ids = _select_copper_chain_ids(block)
    if protein_chain is None or not copper_chain_ids:
        return copper_chain_ids[:1]

    distances: list[tuple[str, float]] = []
    for copper_chain_id in copper_chain_ids:
        copper_chain = _find_chain(model, copper_chain_id)
        if copper_chain is None:
            continue
        distance = _min_chain_distance(protein_chain, copper_chain)
        if distance is not None:
            distances.append((copper_chain_id, distance))
    if distances:
        return (min(distances, key=lambda item: item[1])[0],)
    return copper_chain_ids[:1]


def select_crystal_reference_site(
    source_cif: Path,
    *,
    preferred_protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    expected_ligand: bool | None = None,
    ligand_distance_cutoff_a: float = DEFAULT_LIGAND_DISTANCE_CUTOFF_A,
) -> CrystalReferenceSelection:
    """Pick the protein copy and ligand chains to use from one crystal mmCIF."""
    structure = gemmi.read_structure(str(source_cif))
    if len(structure) == 0:
        raise ValueError(f"Crystal reference has no models: {source_cif}")
    model = structure[0]
    block = _read_reference_block(source_cif)
    available_model_chain_ids = {
        str(getattr(chain, "name", "") or "")
        for chain in model
        if str(getattr(chain, "name", "") or "")
    }

    protein_chain_ids = tuple(
        str(getattr(chain, "name", "") or "")
        for chain in model
        if _chain_looks_like_protein(chain)
    )
    if not protein_chain_ids:
        raise ValueError(f"No protein-like chains found in crystal reference: {source_cif}")

    ligand_chain_mappings = _find_branch_ligand_chain_mappings(block, available_model_chain_ids)
    ligand_chain_ids = tuple(ligand_chain_mappings)
    if not ligand_chain_ids:
        ligand_chain_ids = _find_nonpoly_carbohydrate_chain_ids(block)
        ligand_chain_mappings = {chain_id: chain_id for chain_id in ligand_chain_ids}
    expected_ligand = bool(ligand_chain_ids) if expected_ligand is None else bool(expected_ligand)

    ligand_distances_by_protein: dict[str, list[tuple[str, float]]] = {}
    for protein_chain_id in protein_chain_ids:
        protein_chain = _find_chain(model, protein_chain_id)
        if protein_chain is None:
            continue
        distances: list[tuple[str, float]] = []
        for ligand_chain_id in ligand_chain_ids:
            ligand_chain = _find_chain(
                model,
                ligand_chain_mappings.get(ligand_chain_id, ligand_chain_id),
            )
            if ligand_chain is None:
                continue
            distance = _min_chain_distance(protein_chain, ligand_chain)
            if distance is not None:
                distances.append((ligand_chain_id, distance))
        ligand_distances_by_protein[protein_chain_id] = distances

    selected_protein_chain = (
        preferred_protein_chain if preferred_protein_chain in protein_chain_ids else protein_chain_ids[0]
    )
    preferred_distances = ligand_distances_by_protein.get(preferred_protein_chain, [])
    preferred_chain_has_ligand = any(
        distance <= ligand_distance_cutoff_a for _, distance in preferred_distances
    )

    if expected_ligand and not preferred_chain_has_ligand:
        fallback_candidates = [
            protein_chain_id
            for protein_chain_id, distances in ligand_distances_by_protein.items()
            if any(distance <= ligand_distance_cutoff_a for _, distance in distances)
        ]
        if fallback_candidates:
            selected_protein_chain = min(
                fallback_candidates,
                key=lambda protein_chain_id: min(
                    distance
                    for _, distance in ligand_distances_by_protein[protein_chain_id]
                    if distance <= ligand_distance_cutoff_a
                ),
            )

    selected_distances = ligand_distances_by_protein.get(selected_protein_chain, [])
    selected_ligand_chain_ids = tuple(
        ligand_chain_id
        for ligand_chain_id, distance in selected_distances
        if distance <= ligand_distance_cutoff_a
    )
    selected_ligand_distance_a = min(
        (
            distance
            for _, distance in selected_distances
            if distance <= ligand_distance_cutoff_a
        ),
        default=None,
    )

    return CrystalReferenceSelection(
        source_cif=source_cif.resolve(),
        preferred_protein_chain=preferred_protein_chain,
        selected_protein_chain=selected_protein_chain,
        ligand_chain_ids=selected_ligand_chain_ids,
        copper_chain_ids=_select_copper_chain_ids_for_site(
            block,
            model,
            selected_protein_chain,
        ),
        expected_ligand=expected_ligand,
        preferred_chain_has_ligand=preferred_chain_has_ligand,
        selected_chain_has_ligand=bool(selected_ligand_chain_ids),
        used_fallback_protein_chain=selected_protein_chain != preferred_protein_chain,
        selected_ligand_distance_a=selected_ligand_distance_a,
    )


def _format_mmcif_value(value: Any) -> str:
    text = str(value)
    if "\n" in text or "\r" in text:
        text = " ".join(
            line.strip().strip(";")
            for line in text.replace("\r", "\n").split("\n")
            if line.strip().strip(";")
        )
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        text = text[1:-1]
    if text == "":
        return "?"
    if any(character.isspace() for character in text) or any(character in text for character in ("'", '"', "#")):
        return "'{}'".format(text.replace("'", "''"))
    return text


def _write_optional_filtered_loop(
    lines: list[str],
    block: Any,
    prefix: str,
    columns: list[str],
    chain_id_column: str,
    allowed_chain_ids: set[str],
    *,
    chain_id_map: dict[str, str] | None = None,
    remapped_chain_columns: tuple[str, ...] = (),
) -> None:
    table = block.find(prefix, columns)
    if not table:
        return

    chain_id_index = columns.index(chain_id_column)
    filtered_rows = [
        row
        for row in table
        if str(row[chain_id_index]).strip() in allowed_chain_ids
    ]
    if not filtered_rows:
        return

    lines.append("loop_")
    lines.extend([f"{prefix}{column}" for column in columns])
    for row in filtered_rows:
        row_values = list(row)
        if chain_id_map:
            for column_name in remapped_chain_columns:
                column_index = columns.index(column_name)
                source_chain_id = str(row_values[column_index]).strip()
                row_values[column_index] = chain_id_map.get(source_chain_id, row_values[column_index])
        lines.append(" ".join(_format_mmcif_value(value) for value in row_values))
    lines.append("#")


def _iter_remapped_chain_ids() -> Iterable[str]:
    used_defaults: set[str] = set()
    for chain_id in DEFAULT_GLYCAN_CHAINS:
        if chain_id and chain_id not in used_defaults:
            used_defaults.add(chain_id)
            yield chain_id

    for codepoint in range(ord("A"), ord("Z") + 1):
        chain_id = chr(codepoint)
        if chain_id == DEFAULT_PROTEIN_CHAIN or chain_id in used_defaults:
            continue
        yield chain_id


def _build_selected_chain_id_map(selection: CrystalReferenceSelection) -> dict[str, str]:
    chain_id_map = {selection.selected_protein_chain: DEFAULT_PROTEIN_CHAIN}
    used_chain_ids = {DEFAULT_PROTEIN_CHAIN}

    for source_chain_id, mapped_chain_id in zip(
        selection.ligand_chain_ids,
        DEFAULT_GLYCAN_CHAINS,
        strict=False,
    ):
        if source_chain_id in chain_id_map:
            continue
        chain_id_map[source_chain_id] = mapped_chain_id
        used_chain_ids.add(mapped_chain_id)

    extra_chain_ids = _iter_remapped_chain_ids()
    for source_chain_id in selection.ligand_chain_ids[len(DEFAULT_GLYCAN_CHAINS):]:
        if source_chain_id in chain_id_map:
            continue
        mapped_chain_id = next(
            candidate for candidate in extra_chain_ids if candidate not in used_chain_ids
        )
        chain_id_map[source_chain_id] = mapped_chain_id
        used_chain_ids.add(mapped_chain_id)

    for index, source_chain_id in enumerate(selection.copper_chain_ids):
        if source_chain_id in chain_id_map:
            continue
        if index == 0 and DEFAULT_CU_CHAIN not in used_chain_ids:
            mapped_chain_id = DEFAULT_CU_CHAIN
        else:
            mapped_chain_id = next(
                candidate for candidate in extra_chain_ids if candidate not in used_chain_ids
            )
        chain_id_map[source_chain_id] = mapped_chain_id
        used_chain_ids.add(mapped_chain_id)

    return chain_id_map


def _loop_rows_for_tag(
    block: Any,
    required_tag: str,
    *,
    fallback_prefix: str | None = None,
    fallback_columns: list[str] | None = None,
) -> tuple[list[str], list[list[str]]]:
    find_loop = getattr(block, "find_loop", None)
    if callable(find_loop):
        try:
            column = find_loop(required_tag)
        except Exception:
            column = None
        if column:
            loop = column.get_loop()
            tags = [str(tag) for tag in loop.tags]
            width = int(loop.width())
            length = int(loop.length())
            rows = [
                [str(loop.values[row_index * width + column_index]) for column_index in range(width)]
                for row_index in range(length)
            ]
            return tags, rows

    if fallback_prefix and fallback_columns:
        rows = block.find(fallback_prefix, fallback_columns)
        if rows:
            return [f"{fallback_prefix}{column}" for column in fallback_columns], [list(row) for row in rows]
    return [], []


def _write_loop(lines: list[str], tags: list[str], rows: Iterable[Iterable[Any]]) -> None:
    row_values = [list(row) for row in rows]
    if not tags or not row_values:
        return
    lines.append("loop_")
    lines.extend(tags)
    for row in row_values:
        lines.append(" ".join(_format_mmcif_value(value) for value in row))
    lines.append("#")


def _row_value(tags: list[str], row: list[Any], tag: str) -> str:
    try:
        return str(row[tags.index(tag)]).strip()
    except ValueError:
        return ""


def _set_row_value(tags: list[str], row: list[Any], tag: str, value: str) -> None:
    try:
        row[tags.index(tag)] = value
    except ValueError:
        return


def _remapped_chain_for_atom_row(
    row: list[Any],
    tags: list[str],
    selection: CrystalReferenceSelection,
    chain_id_map: dict[str, str],
) -> str | None:
    chain_id = _row_value(tags, row, "_atom_site.label_asym_id")
    comp_id = _row_value(tags, row, "_atom_site.label_comp_id")
    if chain_id == selection.selected_protein_chain and _is_protein_comp_id(comp_id):
        return chain_id_map[selection.selected_protein_chain]
    if chain_id in selection.ligand_chain_ids and _is_ligand_comp_id(comp_id):
        return chain_id_map[chain_id]
    if chain_id in selection.copper_chain_ids and _is_copper_comp_id(comp_id):
        return chain_id_map[chain_id]
    return None


def _struct_conn_partner_allowed(
    chain_id: str,
    comp_id: str,
    selection: CrystalReferenceSelection,
) -> bool:
    if not chain_id or chain_id in {"?", "."}:
        return True
    if chain_id == selection.selected_protein_chain:
        return _is_protein_comp_id(comp_id)
    if chain_id in selection.ligand_chain_ids:
        return _is_ligand_comp_id(comp_id)
    if chain_id in selection.copper_chain_ids:
        return _is_copper_comp_id(comp_id)
    return False


def _filter_and_remap_struct_conn_rows(
    tags: list[str],
    rows: list[list[str]],
    selection: CrystalReferenceSelection,
    chain_id_map: dict[str, str],
) -> list[list[str]]:
    filtered_rows: list[list[str]] = []
    for source_row in rows:
        row = list(source_row)
        include = False
        allowed = True
        partner_destinations: dict[int, str] = {}

        for partner_index in (1, 2, 3):
            label_chain_tag = f"_struct_conn.ptnr{partner_index}_label_asym_id"
            label_comp_tag = f"_struct_conn.ptnr{partner_index}_label_comp_id"
            auth_chain_tag = f"_struct_conn.ptnr{partner_index}_auth_asym_id"
            auth_comp_tag = f"_struct_conn.ptnr{partner_index}_auth_comp_id"
            chain_id = _row_value(tags, row, label_chain_tag)
            comp_id = _row_value(tags, row, label_comp_tag) or _row_value(tags, row, auth_comp_tag)
            if not chain_id:
                chain_id = _row_value(tags, row, auth_chain_tag)
            if not chain_id and not comp_id:
                continue
            if chain_id in chain_id_map:
                include = True
                partner_destinations[partner_index] = chain_id_map[chain_id]
            if not _struct_conn_partner_allowed(chain_id, comp_id, selection):
                allowed = False
                break

        if not include or not allowed:
            continue

        for partner_index, destination_chain in partner_destinations.items():
            _set_row_value(tags, row, f"_struct_conn.ptnr{partner_index}_label_asym_id", destination_chain)
            _set_row_value(tags, row, f"_struct_conn.ptnr{partner_index}_auth_asym_id", destination_chain)
        filtered_rows.append(row)
    return filtered_rows


def _write_selected_reference_cif(
    source_cif: Path,
    selection: CrystalReferenceSelection,
    output_path: Path,
) -> Path:
    """Write a metadata-preserving mmCIF for the selected crystal view."""
    block = _read_reference_block(source_cif)
    chain_id_map = _build_selected_chain_id_map(selection)
    included_chain_ids = {
        selection.selected_protein_chain,
        *selection.ligand_chain_ids,
        *selection.copper_chain_ids,
    }
    atom_columns = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_model_num",
    ]
    atom_tags = [f"_atom_site.{column}" for column in atom_columns]
    atom_rows = [list(row) for row in block.find("_atom_site.", atom_columns)]
    if not atom_rows:
        raise ValueError(f"No _atom_site loop found in crystal reference: {source_cif}")

    filtered_atom_rows: list[list[str]] = []
    retained_entity_ids: set[str] = set()
    retained_comp_ids: set[str] = set()
    for atom_row in atom_rows:
        row = list(atom_row)
        destination_chain = _remapped_chain_for_atom_row(row, atom_tags, selection, chain_id_map)
        if destination_chain is None:
            continue
        _set_row_value(atom_tags, row, "_atom_site.label_asym_id", destination_chain)
        _set_row_value(atom_tags, row, "_atom_site.auth_asym_id", destination_chain)
        retained_entity_ids.add(_row_value(atom_tags, row, "_atom_site.label_entity_id"))
        retained_comp_ids.add(_normalize_comp_id(_row_value(atom_tags, row, "_atom_site.label_comp_id")))
        filtered_atom_rows.append(row)
    if not filtered_atom_rows:
        raise ValueError(
            f"Selected crystal subset for {source_cif} contained no atom_site rows for chains {sorted(included_chain_ids)}"
        )
    retained_entity_ids.discard("")
    retained_entity_ids.discard("?")
    retained_entity_ids.discard(".")
    retained_comp_ids.discard("")
    retained_comp_ids.discard("?")
    retained_comp_ids.discard(".")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"data_{output_path.stem}", f"_entry.id {output_path.stem}", "#"]

    entity_tags, entity_rows = _loop_rows_for_tag(
        block,
        "_entity.id",
        fallback_prefix="_entity.",
        fallback_columns=[
            "id",
            "type",
            "src_method",
            "pdbx_description",
            "formula_weight",
            "pdbx_number_of_molecules",
            "pdbx_ec",
            "pdbx_mutation",
            "pdbx_fragment",
            "details",
        ],
    )
    _write_loop(
        lines,
        entity_tags,
        [
            row
            for row in entity_rows
            if _row_value(entity_tags, row, "_entity.id") in retained_entity_ids
        ],
    )

    struct_asym_tags, struct_asym_rows = _loop_rows_for_tag(
        block,
        "_struct_asym.id",
        fallback_prefix="_struct_asym.",
        fallback_columns=["id", "pdbx_blank_PDB_chainid_flag", "pdbx_modified", "entity_id", "details"],
    )
    remapped_struct_asym_rows: list[list[str]] = []
    for source_row in struct_asym_rows:
        row = list(source_row)
        source_chain_id = _row_value(struct_asym_tags, row, "_struct_asym.id")
        if source_chain_id not in included_chain_ids:
            continue
        _set_row_value(struct_asym_tags, row, "_struct_asym.id", chain_id_map[source_chain_id])
        remapped_struct_asym_rows.append(row)
    _write_loop(lines, struct_asym_tags, remapped_struct_asym_rows)

    chem_comp_tags, chem_comp_rows = _loop_rows_for_tag(
        block,
        "_chem_comp.id",
        fallback_prefix="_chem_comp.",
        fallback_columns=[
            "id",
            "type",
            "mon_nstd_flag",
            "name",
            "pdbx_synonyms",
            "formula",
            "formula_weight",
        ],
    )
    _write_loop(
        lines,
        chem_comp_tags,
        [
            row
            for row in chem_comp_rows
            if _normalize_comp_id(_row_value(chem_comp_tags, row, "_chem_comp.id")) in retained_comp_ids
        ],
    )

    chem_comp_bond_tags, chem_comp_bond_rows = _loop_rows_for_tag(
        block,
        "_chem_comp_bond.comp_id",
        fallback_prefix="_chem_comp_bond.",
        fallback_columns=[
            "comp_id",
            "atom_id_1",
            "atom_id_2",
            "value_order",
            "pdbx_aromatic_flag",
            "pdbx_stereo_config",
            "pdbx_ordinal",
        ],
    )
    _write_loop(
        lines,
        chem_comp_bond_tags,
        [
            row
            for row in chem_comp_bond_rows
            if _normalize_comp_id(_row_value(chem_comp_bond_tags, row, "_chem_comp_bond.comp_id")) in retained_comp_ids
        ],
    )

    _write_loop(lines, atom_tags, filtered_atom_rows)

    branch_columns = [
            "asym_id",
            "entity_id",
            "mon_id",
            "num",
            "pdb_asym_id",
            "pdb_mon_id",
            "pdb_seq_num",
            "auth_asym_id",
            "auth_mon_id",
            "auth_seq_num",
            "hetero",
    ]
    branch_tags, branch_rows = _loop_rows_for_tag(
        block,
        "_pdbx_branch_scheme.asym_id",
        fallback_prefix="_pdbx_branch_scheme.",
        fallback_columns=branch_columns,
    )
    remapped_branch_rows: list[list[str]] = []
    for source_row in branch_rows:
        row = list(source_row)
        source_chain_id = _row_value(branch_tags, row, "_pdbx_branch_scheme.asym_id")
        if source_chain_id not in selection.ligand_chain_ids:
            continue
        destination_chain = chain_id_map[source_chain_id]
        _set_row_value(branch_tags, row, "_pdbx_branch_scheme.asym_id", destination_chain)
        _set_row_value(branch_tags, row, "_pdbx_branch_scheme.pdb_asym_id", destination_chain)
        _set_row_value(branch_tags, row, "_pdbx_branch_scheme.auth_asym_id", destination_chain)
        remapped_branch_rows.append(row)
    _write_loop(lines, branch_tags, remapped_branch_rows)

    nonpoly_columns = [
            "asym_id",
            "entity_id",
            "mon_id",
            "ndb_seq_num",
            "pdb_seq_num",
            "auth_seq_num",
            "pdb_mon_id",
            "auth_mon_id",
            "pdb_strand_id",
            "pdb_ins_code",
    ]
    nonpoly_tags, nonpoly_rows = _loop_rows_for_tag(
        block,
        "_pdbx_nonpoly_scheme.asym_id",
        fallback_prefix="_pdbx_nonpoly_scheme.",
        fallback_columns=nonpoly_columns,
    )
    remapped_nonpoly_rows: list[list[str]] = []
    for source_row in nonpoly_rows:
        row = list(source_row)
        source_chain_id = _row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.asym_id")
        mon_id = _row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.mon_id")
        if source_chain_id not in selection.copper_chain_ids or not _is_copper_comp_id(mon_id):
            continue
        destination_chain = chain_id_map[source_chain_id]
        _set_row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.asym_id", destination_chain)
        _set_row_value(nonpoly_tags, row, "_pdbx_nonpoly_scheme.pdb_strand_id", destination_chain)
        remapped_nonpoly_rows.append(row)
    _write_loop(lines, nonpoly_tags, remapped_nonpoly_rows)

    struct_conn_type_tags, struct_conn_type_rows = _loop_rows_for_tag(
        block,
        "_struct_conn_type.id",
        fallback_prefix="_struct_conn_type.",
        fallback_columns=["id", "criteria", "reference"],
    )
    _write_loop(lines, struct_conn_type_tags, struct_conn_type_rows)

    struct_conn_tags, struct_conn_rows = _loop_rows_for_tag(
        block,
        "_struct_conn.id",
    )
    _write_loop(
        lines,
        struct_conn_tags,
        _filter_and_remap_struct_conn_rows(
            struct_conn_tags,
            struct_conn_rows,
            selection,
            chain_id_map,
        ),
    )

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def compute_feature_aligned_tanimoto(
    feature_names_a: list[str],
    bitvec_a: list[int],
    feature_names_b: list[str],
    bitvec_b: list[int],
) -> float:
    """Tanimoto similarity after aligning onto the union feature set."""
    if len(feature_names_a) != len(bitvec_a):
        raise ValueError("feature_names_a and bitvec_a must have the same length")
    if len(feature_names_b) != len(bitvec_b):
        raise ValueError("feature_names_b and bitvec_b must have the same length")

    union_feature_names = sorted(set(feature_names_a) | set(feature_names_b))
    feature_map_a = {
        feature_name: int(value)
        for feature_name, value in zip(feature_names_a, bitvec_a)
    }
    feature_map_b = {
        feature_name: int(value)
        for feature_name, value in zip(feature_names_b, bitvec_b)
    }
    aligned_a = [feature_map_a.get(feature_name, 0) for feature_name in union_feature_names]
    aligned_b = [feature_map_b.get(feature_name, 0) for feature_name in union_feature_names]
    return compute_tanimoto_similarity(aligned_a, aligned_b)


def _write_ifp_artifacts(result: IFPResult, output_dir: Path) -> tuple[Path, Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    pose_ifp_table_tsv = output_dir / "pose_ifp_table.tsv"
    ifp_matrix_csv = output_dir / "ifp_matrix.csv"
    ifp_result_json = output_dir / "ifp_result.json"

    write_pose_ifp_table([result], pose_ifp_table_tsv)
    write_ifp_matrix(
        IFPBatch(
            results=[result],
            matrix=[list(result.flat_bitvector)],
            feature_names=list(result.feature_names),
        ),
        ifp_matrix_csv,
    )

    active_feature_names = [
        feature_name
        for feature_name, value in zip(result.feature_names, result.flat_bitvector, strict=True)
        if int(value)
    ]
    data = {
        "pose_id": result.pose_id,
        "status": result.status,
        "error": result.error,
        "n_residues": result.n_residues,
        "n_interaction_types": result.n_interaction_types,
        "residue_names": result.residue_names,
        "n_total_contacts": result.n_total_contacts,
        "interaction_counts": result.interaction_counts,
        "active_feature_names": active_feature_names,
    }
    ifp_result_json.write_text(json.dumps(data, indent=2))
    return output_dir, ifp_result_json, pose_ifp_table_tsv, ifp_matrix_csv


def _contact_eligibility_for_ifp(result: IFPResult) -> ContactEligibility:
    return evaluate_contact_eligibility(result, load_contact_eligibility_rule())


def _compute_crystal_geometry_metrics(
    prepared: PreparedCrystalReference,
) -> PoseGeometryMetrics | None:
    if prepared.complex_pdb is None:
        return None

    chain_id_map = _build_selected_chain_id_map(prepared.selection)
    remapped_copper_chains = [
        chain_id_map[source_chain_id]
        for source_chain_id in prepared.selection.copper_chain_ids
        if source_chain_id in chain_id_map
    ]
    remapped_ligand_chains = [
        chain_id_map[source_chain_id]
        for source_chain_id in prepared.selection.ligand_chain_ids
        if source_chain_id in chain_id_map
    ]
    cu_chain = remapped_copper_chains[0] if remapped_copper_chains else DEFAULT_CU_CHAIN
    glycan_chains = remapped_ligand_chains or [DEFAULT_LIGAND_CHAIN]

    return compute_pose_metrics(
        prepared.complex_pdb,
        pose_id=prepared.record.pdb_code,
        protein_chain=DEFAULT_PROTEIN_CHAIN,
        glycan_chains=glycan_chains,
        cu_chain=cu_chain,
        model="crystal",
        protein_id=prepared.record.protein_id,
        ligand_id=prepared.record.carbohydrate_ligands or prepared.record.pdb_code,
    )


def _crystal_geometry_payload(metrics: PoseGeometryMetrics | None) -> dict[str, Any]:
    if metrics is None:
        return {}
    return metrics.to_row()


def crystal_geometry_row(prepared: PreparedCrystalReference) -> dict[str, Any]:
    metrics = prepared.geometry_metrics
    geometry_status = (
        "computed"
        if metrics is not None
        else ("no_ligand_reference" if prepared.status == "prepared_no_ligand" else prepared.status)
    )
    row = {
        "protein_id": prepared.record.protein_id,
        "pdb_code": prepared.record.pdb_code,
        "source_cif": str(prepared.record.source_cif),
        "prepared_subset_cif": str(prepared.subset_cif),
        "selected_protein_chain": prepared.selection.selected_protein_chain,
        "ligand_chain_ids": json.dumps(list(prepared.selection.ligand_chain_ids)),
        "copper_chain_ids": json.dumps(list(prepared.selection.copper_chain_ids)),
        "geometry_status": geometry_status,
    }
    if metrics is None:
        row.update({column: "" for column in POSE_GEOMETRY_COLUMNS})
    else:
        row.update(metrics.to_row())
    return row


def write_crystal_geometry_table(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CRYSTAL_GEOMETRY_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in CRYSTAL_GEOMETRY_COLUMNS})


def prepare_crystal_reference(
    record: CrystalReferenceRecord,
    output_dir: Path,
    *,
    ligand_distance_cutoff_a: float = DEFAULT_LIGAND_DISTANCE_CUTOFF_A,
) -> PreparedCrystalReference:
    """Prepare one crystal reference through subset, normalization, export, and IFP."""
    selection = select_crystal_reference_site(
        record.source_cif,
        preferred_protein_chain=record.preferred_author_chain,
        expected_ligand=record.expected_ligand,
        ligand_distance_cutoff_a=ligand_distance_cutoff_a,
    )
    reference_dir = output_dir / f"{record.pdb_code}_{record.protein_id}"
    subset_cif = _write_selected_reference_cif(
        record.source_cif,
        selection,
        reference_dir / "selected_reference.cif",
    )
    prepared_input_cif = subset_cif.resolve()
    if not selection.selected_chain_has_ligand:
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            status="prepared_no_ligand",
        )

    try:
        normalize_ok, normalized_cif = NormalizeMMCIFRunner(
            prepared_input_cif,
            reference_dir / "normalize",
        ).run()
    except Exception as exc:
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            status="normalization_failed",
            error=f"{exc.__class__.__name__}: {exc}",
        )

    if not normalize_ok or normalized_cif is None:
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            status="normalization_failed",
            error="Crystal reference normalization failed",
        )

    prepared_input_cif = Path(normalized_cif).resolve()

    try:
        export_ok, export_report = export_analysis_artifacts(prepared_input_cif, reference_dir / "analysis_export")
    except Exception as exc:
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            status="analysis_export_failed",
            error=f"{exc.__class__.__name__}: {exc}",
        )

    if not export_ok or not export_report:
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            status="analysis_export_failed",
            error="Analysis export failed",
        )

    complex_pdb = Path(str(export_report.get("complex_for_prolif_pdb", ""))).resolve()
    ligand_pdb = Path(str(export_report.get("ligand_pdb", ""))).resolve()
    if not complex_pdb.exists() or not ligand_pdb.exists():
        return PreparedCrystalReference(
            record=record,
            selection=selection,
            subset_cif=subset_cif,
            normalized_cif=prepared_input_cif,
            complex_pdb=complex_pdb,
            ligand_pdb=ligand_pdb,
            status="analysis_export_artifacts_missing",
            error="Crystal analysis export artifacts were not written",
        )

    crystal_ifp = compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_pdb=ligand_pdb,
        pose_id=record.pdb_code,
        protein_chain=DEFAULT_PROTEIN_CHAIN,
    )
    crystal_eligibility = (
        _contact_eligibility_for_ifp(crystal_ifp)
        if crystal_ifp.status in {"ok", "zero_contacts"}
        else None
    )
    try:
        geometry_metrics = _compute_crystal_geometry_metrics(
            PreparedCrystalReference(
                record=record,
                selection=selection,
                subset_cif=subset_cif,
                normalized_cif=prepared_input_cif,
                complex_pdb=complex_pdb,
                ligand_pdb=ligand_pdb,
            )
        )
    except Exception as exc:
        logger.warning(
            "Crystal geometry failed for %s: %s: %s",
            record.pdb_code,
            exc.__class__.__name__,
            exc,
        )
        geometry_metrics = None
    ifp_artifact_dir, ifp_result_json, pose_ifp_table_tsv, ifp_matrix_csv = _write_ifp_artifacts(
        crystal_ifp,
        reference_dir / "ifp",
    )
    return PreparedCrystalReference(
        record=record,
        selection=selection,
        subset_cif=subset_cif,
        normalized_cif=prepared_input_cif,
        complex_pdb=complex_pdb,
        ligand_pdb=ligand_pdb,
        ifp_artifact_dir=ifp_artifact_dir,
        ifp_result_json=ifp_result_json,
        pose_ifp_table_tsv=pose_ifp_table_tsv,
        ifp_matrix_csv=ifp_matrix_csv,
        ifp_result=crystal_ifp,
        ifp_contact_eligibility=crystal_eligibility,
        geometry_metrics=geometry_metrics,
        status="prepared" if crystal_ifp.status == "ok" else "ifp_failed",
        error=crystal_ifp.error,
    )


def _prepare_representative_pose_ifp(
    representative_pose_cif: Path,
    output_dir: Path,
    *,
    representative_pose_id: str,
    representative_normalized_cif: Path | None = None,
    representative_complex_pdb: Path | None = None,
    representative_ligand_pdb: Path | None = None,
    representative_ifp_result: IFPResult | None = None,
) -> tuple[IFPResult | None, str | None, tuple[Path, Path, Path, Path] | None, Path | None]:
    complex_pdb = Path(representative_complex_pdb).resolve() if representative_complex_pdb else None
    ligand_pdb = Path(representative_ligand_pdb).resolve() if representative_ligand_pdb else None
    if representative_ifp_result is not None and complex_pdb is not None and complex_pdb.exists():
        ifp_artifacts = _write_ifp_artifacts(representative_ifp_result, output_dir / "ifp")
        return representative_ifp_result, representative_ifp_result.error, ifp_artifacts, complex_pdb

    if complex_pdb is None or ligand_pdb is None or not complex_pdb.exists() or not ligand_pdb.exists():
        normalized_cif: Path | None = None
        if representative_normalized_cif is not None and Path(representative_normalized_cif).is_file():
            normalized_cif = Path(representative_normalized_cif).resolve()
        else:
            try:
                normalize_ok, normalized_cif_raw = NormalizeMMCIFRunner(
                    representative_pose_cif,
                    output_dir / "normalize",
                ).run()
            except Exception as exc:
                return None, f"{exc.__class__.__name__}: {exc}", None, None

            if not normalize_ok or normalized_cif_raw is None:
                return None, "Normalization failed", None, None
            normalized_cif = Path(normalized_cif_raw)

        try:
            export_ok, export_report = export_analysis_artifacts(Path(normalized_cif), output_dir / "analysis_export")
        except Exception as exc:
            return None, f"{exc.__class__.__name__}: {exc}", None, None

        if not export_ok or not export_report:
            return None, "Analysis export failed", None, None

        complex_pdb = Path(str(export_report.get("complex_for_prolif_pdb", ""))).resolve()
        ligand_pdb = Path(str(export_report.get("ligand_pdb", ""))).resolve()

    if not complex_pdb.exists() or not ligand_pdb.exists():
        return None, "Representative pose analysis export artifacts missing", None, None

    pose_ifp = representative_ifp_result or compute_ifp_single(
        complex_pdb=complex_pdb,
        ligand_pdb=ligand_pdb,
        pose_id=representative_pose_id,
        protein_chain=DEFAULT_PROTEIN_CHAIN,
    )
    ifp_artifacts = _write_ifp_artifacts(pose_ifp, output_dir / "ifp")
    return pose_ifp, pose_ifp.error, ifp_artifacts, complex_pdb


def run_crystal_reference_screen(
    representative_pose_cif: Path,
    *,
    protein_id: str,
    output_dir: Path,
    representative_pose_id: str = "",
    ligand_id: str = "",
    reference_index_csv: Path = DEFAULT_REFERENCE_INDEX_CSV,
    crystal_root: Path = DEFAULT_CRYSTAL_ROOT,
    records: Iterable[CrystalReferenceRecord] | None = None,
    representative_normalized_cif: Path | None = None,
    representative_complex_pdb: Path | None = None,
    representative_ligand_pdb: Path | None = None,
    representative_ifp_result: IFPResult | None = None,
) -> CrystalReferenceScreenReport:
    """Compare one representative pose against all crystal references for a protein."""
    representative_pose_path = Path(representative_pose_cif).resolve()
    representative_pose_id = representative_pose_id or representative_pose_path.stem
    report = CrystalReferenceScreenReport(
        protein_id=_normalize_reference_protein_id(protein_id),
        representative_pose_id=representative_pose_id,
        representative_pose_cif=str(representative_pose_path),
        ligand_id=ligand_id,
    )

    pose_ifp, pose_error, representative_pose_ifp_artifacts, representative_pose_complex_pdb = _prepare_representative_pose_ifp(
        representative_pose_path,
        output_dir / "representative_pose",
        representative_pose_id=representative_pose_id,
        representative_normalized_cif=representative_normalized_cif,
        representative_complex_pdb=representative_complex_pdb,
        representative_ligand_pdb=representative_ligand_pdb,
        representative_ifp_result=representative_ifp_result,
    )
    if representative_pose_ifp_artifacts is not None:
        _, ifp_result_json, pose_ifp_table_tsv, ifp_matrix_csv = representative_pose_ifp_artifacts
        report.representative_pose_ifp_result_json = str(ifp_result_json)
        report.representative_pose_pose_ifp_table_tsv = str(pose_ifp_table_tsv)
        report.representative_pose_ifp_matrix_csv = str(ifp_matrix_csv)
    representative_pose_pocket_residues = []
    if representative_pose_complex_pdb is not None:
        representative_pose_pocket_residues = identify_pocket_residues_by_proximity(
            representative_pose_complex_pdb,
            ligand_chain=DEFAULT_LIGAND_CHAIN,
            protein_chain=DEFAULT_PROTEIN_CHAIN,
        )
    records_to_compare = list(records) if records is not None else load_crystal_reference_records(
        report.protein_id,
        reference_index_csv=reference_index_csv,
        crystal_root=crystal_root,
    )

    for record in records_to_compare:
        prepared = prepare_crystal_reference(record, output_dir / "references")
        comparison = CrystalReferencePoseComparison(
            pdb_code=record.pdb_code,
            source_cif=str(record.source_cif),
            prepared_subset_cif=str(prepared.subset_cif),
            prepared_normalized_cif=str(prepared.normalized_cif or ""),
            selected_protein_chain=prepared.selection.selected_protein_chain,
            ligand_chain_ids=list(prepared.selection.ligand_chain_ids),
            copper_chain_ids=list(prepared.selection.copper_chain_ids),
            preferred_protein_chain=prepared.selection.preferred_protein_chain,
            preferred_chain_has_ligand=prepared.selection.preferred_chain_has_ligand,
            used_fallback_protein_chain=prepared.selection.used_fallback_protein_chain,
            representative_pose_id=representative_pose_id,
            representative_pose_ifp_status=pose_ifp.status if pose_ifp is not None else "pose_ifp_failed",
            representative_pose_ifp_result_json=report.representative_pose_ifp_result_json,
            representative_pose_pose_ifp_table_tsv=report.representative_pose_pose_ifp_table_tsv,
            representative_pose_ifp_matrix_csv=report.representative_pose_ifp_matrix_csv,
            crystal_ifp_status=prepared.ifp_result.status if prepared.ifp_result is not None else prepared.status,
            crystal_ifp_result_json=str(prepared.ifp_result_json or ""),
            crystal_pose_ifp_table_tsv=str(prepared.pose_ifp_table_tsv or ""),
            crystal_ifp_matrix_csv=str(prepared.ifp_matrix_csv or ""),
            crystal_ifp_contact_eligible=(
                bool(prepared.ifp_contact_eligibility.eligible)
                if prepared.ifp_contact_eligibility is not None
                else False
            ),
            crystal_ifp_exclusion_class=(
                str(prepared.ifp_contact_eligibility.exclusion_class or "")
                if prepared.ifp_contact_eligibility is not None
                else ("no_ligand_reference" if prepared.status == "prepared_no_ligand" else "")
            ),
            crystal_n_vdw_interactions=(
                int(prepared.ifp_contact_eligibility.n_vdw_interactions)
                if prepared.ifp_contact_eligibility is not None
                else 0
            ),
            crystal_n_non_vdw_interactions=(
                int(prepared.ifp_contact_eligibility.n_non_vdw_interactions)
                if prepared.ifp_contact_eligibility is not None
                else 0
            ),
            crystal_n_non_vdw_contact_residues=(
                int(prepared.ifp_contact_eligibility.n_non_vdw_contact_residues)
                if prepared.ifp_contact_eligibility is not None
                else 0
            ),
            crystal_geometry=_crystal_geometry_payload(prepared.geometry_metrics),
            status=prepared.status,
            error=prepared.error,
        )

        if representative_pose_complex_pdb is not None:
            pocket_residue_pairs: list[tuple[int, int]] = []
            reference_structure_for_pocket = prepared.complex_pdb or prepared.normalized_cif or prepared.subset_cif

            if prepared.complex_pdb is not None:
                crystal_pocket_residues = identify_pocket_residues_by_proximity(
                    prepared.complex_pdb,
                    ligand_chain=DEFAULT_LIGAND_CHAIN,
                    protein_chain=DEFAULT_PROTEIN_CHAIN,
                )
                mapped_pairs = _map_residue_number_pairs_by_sequence(
                    prepared.complex_pdb,
                    representative_pose_complex_pdb,
                    crystal_pocket_residues,
                    chain_name=DEFAULT_PROTEIN_CHAIN,
                )
                pocket_residue_pairs = [
                    (representative_residue_number, crystal_residue_number)
                    for crystal_residue_number, representative_residue_number in mapped_pairs
                ]
            elif representative_pose_pocket_residues:
                pocket_residue_pairs = _map_residue_number_pairs_by_sequence(
                    representative_pose_complex_pdb,
                    reference_structure_for_pocket,
                    representative_pose_pocket_residues,
                    chain_name=DEFAULT_PROTEIN_CHAIN,
                )

            comparison.pocket_residues = [
                crystal_residue_number for _, crystal_residue_number in pocket_residue_pairs
            ]
            if pocket_residue_pairs:
                comparison.pocket_rmsd = _compute_pocket_rmsd_from_residue_pairs(
                    representative_pose_complex_pdb,
                    reference_structure_for_pocket,
                    pocket_residue_pairs,
                    protein_chain=DEFAULT_PROTEIN_CHAIN,
                )
                comparison.pocket_rmsd_below_threshold = (
                    comparison.pocket_rmsd is not None
                    and comparison.pocket_rmsd < POCKET_RMSD_THRESHOLD
                )

        if pose_ifp is None:
            comparison.status = "pose_ifp_failed"
            comparison.error = pose_error
        elif pose_ifp.status != "ok":
            comparison.status = "pose_ifp_failed"
            comparison.error = pose_ifp.error
        elif prepared.ifp_result is None or prepared.ifp_result.status != "ok":
            comparison.status = prepared.status
            comparison.error = prepared.error or (
                prepared.ifp_result.error if prepared.ifp_result is not None else None
            )
        elif not comparison.crystal_ifp_contact_eligible:
            comparison.status = "crystal_ifp_not_contact_eligible"
            comparison.error = comparison.crystal_ifp_exclusion_class or "crystal_ifp_not_contact_eligible"
        else:
            comparison.ifp_tanimoto = compute_feature_aligned_tanimoto(
                pose_ifp.feature_names,
                pose_ifp.flat_bitvector,
                prepared.ifp_result.feature_names,
                prepared.ifp_result.flat_bitvector,
            )
            comparison.ifp_comparison_eligible = True
            comparison.status = "ok"

        report.comparisons.append(comparison)

    return report


def write_crystal_reference_screen_report(
    report: CrystalReferenceScreenReport,
    output_path: Path,
) -> None:
    """Write representative-pose-vs-crystal comparisons to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "protein_id": report.protein_id,
        "ligand_id": report.ligand_id,
        "representative_pose_id": report.representative_pose_id,
        "representative_pose_cif": report.representative_pose_cif,
        "representative_pose_ifp_result_json": report.representative_pose_ifp_result_json,
        "representative_pose_pose_ifp_table_tsv": report.representative_pose_pose_ifp_table_tsv,
        "representative_pose_ifp_matrix_csv": report.representative_pose_ifp_matrix_csv,
        "comparisons": [
            {
                "pdb_code": comparison.pdb_code,
                "source_cif": comparison.source_cif,
                "prepared_subset_cif": comparison.prepared_subset_cif,
                "prepared_normalized_cif": comparison.prepared_normalized_cif,
                "selected_protein_chain": comparison.selected_protein_chain,
                "ligand_chain_ids": comparison.ligand_chain_ids,
                "copper_chain_ids": comparison.copper_chain_ids,
                "preferred_protein_chain": comparison.preferred_protein_chain,
                "preferred_chain_has_ligand": comparison.preferred_chain_has_ligand,
                "used_fallback_protein_chain": comparison.used_fallback_protein_chain,
                "representative_pose_id": comparison.representative_pose_id,
                "representative_pose_ifp_status": comparison.representative_pose_ifp_status,
                "representative_pose_ifp_result_json": comparison.representative_pose_ifp_result_json,
                "representative_pose_pose_ifp_table_tsv": comparison.representative_pose_pose_ifp_table_tsv,
                "representative_pose_ifp_matrix_csv": comparison.representative_pose_ifp_matrix_csv,
                "crystal_ifp_status": comparison.crystal_ifp_status,
                "crystal_ifp_result_json": comparison.crystal_ifp_result_json,
                "crystal_pose_ifp_table_tsv": comparison.crystal_pose_ifp_table_tsv,
                "crystal_ifp_matrix_csv": comparison.crystal_ifp_matrix_csv,
                "crystal_ifp_contact_eligible": comparison.crystal_ifp_contact_eligible,
                "crystal_ifp_exclusion_class": comparison.crystal_ifp_exclusion_class,
                "crystal_n_vdw_interactions": comparison.crystal_n_vdw_interactions,
                "crystal_n_non_vdw_interactions": comparison.crystal_n_non_vdw_interactions,
                "crystal_n_non_vdw_contact_residues": comparison.crystal_n_non_vdw_contact_residues,
                "ifp_comparison_eligible": comparison.ifp_comparison_eligible,
                "crystal_geometry": comparison.crystal_geometry,
                "pocket_residues": comparison.pocket_residues,
                "pocket_rmsd": comparison.pocket_rmsd,
                "pocket_rmsd_below_threshold": comparison.pocket_rmsd_below_threshold,
                "ifp_tanimoto": comparison.ifp_tanimoto,
                "status": comparison.status,
                "error": comparison.error,
            }
            for comparison in report.comparisons
        ],
    }
    output_path.write_text(json.dumps(data, indent=2))


def _comparison_tanimoto(
    pose_feature_names: list[str] | None,
    pose_bitvec: list[int],
    crystal_ifp_result: IFPResult,
) -> float:
    if pose_feature_names is None:
        if len(pose_bitvec) != len(crystal_ifp_result.flat_bitvector):
            raise ValueError(
                "Crystal-vs-pose IFP comparison requires feature names when the two bitvectors use different layouts"
            )
        return compute_tanimoto_similarity(pose_bitvec, crystal_ifp_result.flat_bitvector)

    return compute_feature_aligned_tanimoto(
        pose_feature_names,
        pose_bitvec,
        crystal_ifp_result.feature_names,
        crystal_ifp_result.flat_bitvector,
    )


def run_crystal_anchoring(
    crystal_pdb: Path,
    crystal_ligand_pdb: Path | None,
    cluster_medoid_data: list[dict[str, Any]],
    protein_id: str = "",
    ligand_id: str = "",
    pocket_residues: list[int] | None = None,
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
) -> CrystalAnchoringReport:
    """Compare predicted clusters to crystal structure.

    Args:
        crystal_pdb: Path to crystal complex PDB.
        crystal_ligand_pdb: Path to crystal ligand PDB (None if apo).
        cluster_medoid_data: List of {
            "cluster_id": int,
            "medoid_pose_id": str,
            "medoid_pdb": Path,
            "medoid_ligand_pdb": Path,
            "medoid_ifp": list[int],  # flat bitvector
        }
        protein_id: For metadata.
        ligand_id: For metadata.
        pocket_residues: Residue numbers for pocket RMSD.
        protein_chain: Protein chain.

    Returns:
        CrystalAnchoringReport.
    """
    report = CrystalAnchoringReport(
        protein_id=protein_id,
        ligand_id=ligand_id,
        crystal_pdb=str(crystal_pdb),
    )

    if crystal_ligand_pdb is not None and crystal_ligand_pdb.exists():
        report.has_crystal_ligand = True
        crystal_ifp_result = compute_ifp_single(
            complex_pdb=crystal_pdb,
            ligand_pdb=crystal_ligand_pdb,
            pose_id="crystal",
            protein_chain=protein_chain,
        )

        for cdata in cluster_medoid_data:
            tanimoto = _comparison_tanimoto(
                cdata.get("medoid_ifp_feature_names"),
                cdata["medoid_ifp"],
                crystal_ifp_result,
            )
            comparison = CrystalComparisonResult(
                cluster_id=cdata["cluster_id"],
                medoid_pose_id=cdata["medoid_pose_id"],
                ifp_tanimoto=tanimoto,
                ifp_above_threshold=tanimoto >= CRYSTAL_SIM_MODERATE_THRESHOLD,
            )

            if pocket_residues:
                rmsd = _compute_pocket_rmsd(
                    pred_pdb=cdata["medoid_pdb"],
                    crystal_pdb=crystal_pdb,
                    pocket_residues=pocket_residues,
                    protein_chain=protein_chain,
                )
                comparison.pocket_rmsd = rmsd
                comparison.pocket_rmsd_below_threshold = (
                    rmsd is not None and rmsd < POCKET_RMSD_THRESHOLD
                )

            report.comparisons.append(comparison)
    else:
        report.has_crystal_ligand = False
        if pocket_residues:
            for cdata in cluster_medoid_data:
                rmsd = _compute_pocket_rmsd(
                    pred_pdb=cdata["medoid_pdb"],
                    crystal_pdb=crystal_pdb,
                    pocket_residues=pocket_residues,
                    protein_chain=protein_chain,
                )
                comparison = CrystalComparisonResult(
                    cluster_id=cdata["cluster_id"],
                    medoid_pose_id=cdata["medoid_pose_id"],
                    pocket_rmsd=rmsd,
                    pocket_rmsd_below_threshold=(
                        rmsd is not None and rmsd < POCKET_RMSD_THRESHOLD
                    ),
                )
                report.comparisons.append(comparison)

    if report.comparisons:
        best = max(report.comparisons, key=lambda c: c.ifp_tanimoto)
        report.best_cluster_id = best.cluster_id
        report.best_tanimoto = best.ifp_tanimoto

    logger.info(
        "Crystal anchoring %s×%s: %d clusters compared, best Tanimoto=%.3f (cluster %d)",
        protein_id, ligand_id, len(report.comparisons),
        report.best_tanimoto, report.best_cluster_id,
    )
    return report


def _compute_pocket_rmsd(
    pred_pdb: Path,
    crystal_pdb: Path,
    pocket_residues: list[int],
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
) -> float | None:
    """Compute local pocket RMSD after optimal superposition on shared C-alpha atoms.

    This is the first working RMSD slice for crystal anchoring. It aligns the
    predicted pocket onto the crystal pocket using shared protein C-alpha atoms
    from `pocket_residues`, then reports the post-alignment RMSD.
    """
    residue_numbers = sorted({int(residue_number) for residue_number in pocket_residues})
    if not residue_numbers:
        return None

    return _compute_pocket_rmsd_from_residue_pairs(
        pred_pdb,
        crystal_pdb,
        [(residue_number, residue_number) for residue_number in residue_numbers],
        protein_chain=protein_chain,
    )


def _compute_pocket_rmsd_from_residue_pairs(
    pred_pdb: Path,
    crystal_pdb: Path,
    residue_number_pairs: list[tuple[int, int]],
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
) -> float | None:
    normalized_pairs = [
        (int(pred_residue_number), int(crystal_residue_number))
        for pred_residue_number, crystal_residue_number in residue_number_pairs
    ]
    if not normalized_pairs:
        return None

    pred_positions = _collect_ca_positions(
        pred_pdb,
        protein_chain,
        [pred_residue_number for pred_residue_number, _ in normalized_pairs],
    )
    crystal_positions = _collect_ca_positions(
        crystal_pdb,
        protein_chain,
        [crystal_residue_number for _, crystal_residue_number in normalized_pairs],
    )
    shared_pairs = [
        (pred_residue_number, crystal_residue_number)
        for pred_residue_number, crystal_residue_number in normalized_pairs
        if pred_residue_number in pred_positions and crystal_residue_number in crystal_positions
    ]
    if len(shared_pairs) < 2:
        return None

    pred_array = np.asarray(
        [pred_positions[pred_residue_number] for pred_residue_number, _ in shared_pairs],
        dtype=float,
    )
    crystal_array = np.asarray(
        [crystal_positions[crystal_residue_number] for _, crystal_residue_number in shared_pairs],
        dtype=float,
    )
    return _kabsch_rmsd(crystal_array, pred_array)


def write_crystal_anchoring_report(
    report: CrystalAnchoringReport,
    output_path: Path,
) -> None:
    """Write crystal anchoring report to JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "protein_id": report.protein_id,
        "ligand_id": report.ligand_id,
        "crystal_pdb": report.crystal_pdb,
        "has_crystal_ligand": report.has_crystal_ligand,
        "best_cluster_id": report.best_cluster_id,
        "best_tanimoto": report.best_tanimoto,
        "comparisons": [
            {
                "cluster_id": c.cluster_id,
                "medoid_pose_id": c.medoid_pose_id,
                "ifp_tanimoto": c.ifp_tanimoto,
                "pocket_rmsd": c.pocket_rmsd,
                "ifp_above_threshold": c.ifp_above_threshold,
                "pocket_rmsd_below_threshold": c.pocket_rmsd_below_threshold,
            }
            for c in report.comparisons
        ],
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info("Wrote crystal anchoring report to %s", output_path)


def identify_pocket_residues_by_proximity(
    complex_pdb: Path,
    ligand_chain: str = DEFAULT_LIGAND_CHAIN,
    protein_chain: str = DEFAULT_PROTEIN_CHAIN,
    cutoff_a: float = 5.0,
) -> list[int]:
    """Identify protein residues near the ligand or Cu by heavy-atom proximity."""
    structure = gemmi.read_structure(str(complex_pdb))
    if len(structure) == 0:
        return []
    model = structure[0]
    ligand_atoms = _collect_non_hydrogen_atoms(model, ligand_chain)
    copper_atoms = _collect_atoms_by_element(model, "CU")
    focus_atoms = ligand_atoms + copper_atoms
    if not focus_atoms:
        return []

    protein_chain_obj = _find_chain(model, protein_chain)
    if protein_chain_obj is None:
        return []

    cutoff_sq = float(cutoff_a) ** 2
    residue_numbers: set[int] = set()
    for residue in protein_chain_obj:
        protein_atoms = [atom for atom in residue if not _is_hydrogen(atom)]
        if not protein_atoms:
            continue
        if any(
            _squared_distance(focus_atom, protein_atom) <= cutoff_sq
            for focus_atom in focus_atoms
            for protein_atom in protein_atoms
        ):
            residue_numbers.add(int(residue.seqid.num))
    return sorted(residue_numbers)


def _find_chain(model: Any, chain_name: str) -> Any | None:
    for chain in model:
        if chain.name == chain_name:
            return chain
    return None


def _collect_non_hydrogen_atoms(model: Any, chain_name: str) -> list[Any]:
    chain = _find_chain(model, chain_name)
    if chain is None:
        return []
    return [atom for residue in chain for atom in residue if not _is_hydrogen(atom)]


def _collect_atoms_by_element(model: Any, element_name: str) -> list[Any]:
    normalized_element_name = str(element_name).strip().upper()
    return [
        atom
        for chain in model
        for residue in chain
        for atom in residue
        if _atom_element_name(atom) == normalized_element_name
    ]


def _atom_position(atom: Any) -> np.ndarray:
    return np.asarray([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def _squared_distance(atom_a: Any, atom_b: Any) -> float:
    delta = _atom_position(atom_a) - _atom_position(atom_b)
    return float(np.dot(delta, delta))


def _is_hydrogen(atom: Any) -> bool:
    return str(atom.element.name).strip().upper().startswith("H")


def _normalize_protein_residue_name(residue_name: str) -> str:
    normalized_name = str(residue_name or "").strip().upper()
    return _PROTEIN_RESIDUE_ALIASES.get(normalized_name, normalized_name)


def _collect_protein_residue_identities(
    structure_path: Path,
    chain_name: str,
) -> list[tuple[int, str]]:
    structure = gemmi.read_structure(str(structure_path))
    if len(structure) == 0:
        return []
    model = structure[0]
    chain = _find_chain(model, chain_name)
    if chain is None:
        return []

    residue_identities: list[tuple[int, str]] = []
    for residue in chain:
        normalized_name = _normalize_protein_residue_name(str(getattr(residue, "name", "") or ""))
        if normalized_name not in _NORMALIZED_PROTEIN_RESIDUE_NAMES:
            continue
        residue_identities.append((int(residue.seqid.num), normalized_name))
    return residue_identities


def _map_residue_number_pairs_by_sequence(
    source_structure_path: Path,
    target_structure_path: Path,
    residue_numbers: list[int],
    *,
    chain_name: str = DEFAULT_PROTEIN_CHAIN,
) -> list[tuple[int, int]]:
    source_residues = _collect_protein_residue_identities(source_structure_path, chain_name)
    target_residues = _collect_protein_residue_identities(target_structure_path, chain_name)
    if not source_residues or not target_residues:
        return []

    source_index_by_residue_number = {
        residue_number: index for index, (residue_number, _residue_name) in enumerate(source_residues)
    }
    requested_source_indices = [
        source_index_by_residue_number[residue_number]
        for residue_number in sorted({int(residue_number) for residue_number in residue_numbers})
        if int(residue_number) in source_index_by_residue_number
    ]
    if not requested_source_indices:
        return []

    source_names = [residue_name for _residue_number, residue_name in source_residues]
    target_names = [residue_name for _residue_number, residue_name in target_residues]
    source_to_target_index: dict[int, int] = {}

    if len(source_names) == len(target_names) and source_names == target_names:
        source_to_target_index = {index: index for index in range(len(source_names))}
    else:
        matcher = SequenceMatcher(a=source_names, b=target_names, autojunk=False)
        for block in matcher.get_matching_blocks():
            for offset in range(block.size):
                source_to_target_index[block.a + offset] = block.b + offset

    residue_pairs: list[tuple[int, int]] = []
    for source_index in requested_source_indices:
        target_index = source_to_target_index.get(source_index)
        if target_index is None:
            continue

        source_residue_number, source_residue_name = source_residues[source_index]
        target_residue_number, target_residue_name = target_residues[target_index]
        if source_residue_name != target_residue_name:
            continue
        residue_pairs.append((source_residue_number, target_residue_number))

    return residue_pairs


def _collect_ca_positions(
    structure_path: Path,
    chain_name: str,
    residue_numbers: list[int],
) -> dict[int, np.ndarray]:
    structure = gemmi.read_structure(str(structure_path))
    if len(structure) == 0:
        return {}
    model = structure[0]
    chain = _find_chain(model, chain_name)
    if chain is None:
        return {}

    residue_number_set = {int(residue_number) for residue_number in residue_numbers}
    positions: dict[int, np.ndarray] = {}
    for residue in chain:
        residue_number = int(residue.seqid.num)
        if residue_number not in residue_number_set:
            continue
        for atom in residue:
            if atom.name.strip() == "CA":
                positions[residue_number] = _atom_position(atom)
                break
    return positions


def _kabsch_rmsd(reference_positions: np.ndarray, mobile_positions: np.ndarray) -> float:
    if reference_positions.shape != mobile_positions.shape:
        raise ValueError(
            "Pocket RMSD atom count mismatch: "
            f"{reference_positions.shape} != {mobile_positions.shape}"
        )

    reference_centered = reference_positions - reference_positions.mean(axis=0)
    mobile_centered = mobile_positions - mobile_positions.mean(axis=0)
    covariance = mobile_centered.T @ reference_centered
    left, _singular_values, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    aligned_mobile = mobile_centered @ rotation
    squared = np.sum((aligned_mobile - reference_centered) ** 2, axis=1)
    return float(np.sqrt(np.mean(squared)))
