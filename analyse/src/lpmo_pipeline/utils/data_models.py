"""
LPMO Pipeline Data Models
Responsibility: Typed containers for all major data structures

Key entities:
  - RunConfig: full run configuration + tuning flags
  - PoseMetrics: per-pose geometry + QC + cluster assignment
  - ClusterSignature: cluster definition + members + geometry aggregate
  - ActivityFeature: protein-level aggregates (C1/C4 occupancy, IFP freqs)
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from enum import Enum
import json


class QCStatus(Enum):
    """Pass/fail/soft-flag outcomes per QC gate."""
    PASS = "pass"
    SOFT_FLAG = "soft_flag"
    HARD_FAIL = "hard_fail"
    SKIP = "skip"


class FailureReason(Enum):
    """Catalog of failure modes for logging."""
    ATOM_MAPPING_INCOMPLETE = "atom_mapping_incomplete"
    GLYKAN_NOT_CCD = "glykan_not_ccd"
    PLACER_NO_OUTPUT = "placer_no_output"
    POSEBUSTERS_CRITICAL = "posebusters_critical"
    PRIVATEER_ANOMER_FAIL = "privateer_anomer_fail"
    CU_HIS_OUT_OF_RANGE = "cu_his_out_of_range"
    INGEST_PARSE_FAIL = "ingest_parse_fail"
    UNKNOWN = "unknown"


@dataclass
class RunConfig:
    """Unified run configuration (tuning or production)."""
    mode: str  # "tune" or "production"
    model: str  # AF3, RF3, Boltz2
    protein_id: str
    ligand_id: str
    run_id: str
    
    # Model-specific parameters (will be dict or nested dataclass)
    model_params: Dict = field(default_factory=dict)
    
    # QC + pipeline thresholds
    cu_his_dist_min: float = 1.9
    cu_his_dist_max: float = 2.6
    cu_c_proximity_threshold_a: float = 7.0
    crystal_ifp_similarity_threshold: float = 0.3
    
    # HDBSCAN (locked after tuning)
    hdbscan_min_cluster_size: int = 10
    hdbscan_metric: str = "jaccard"
    
    # Paths
    input_protein_seq: Optional[str] = None
    input_template: Optional[str] = None
    input_ligand: Optional[str] = None
    output_dir: Optional[str] = None
    
    # Reproducibility
    seed: Optional[int] = None
    timestamp: Optional[str] = None
    
    def to_json(self) -> str:
        """Serialize to JSON string (for config hashing)."""
        return json.dumps(self.__dict__, default=str, sort_keys=True)


@dataclass
class AtomMap:
    """Cross-model atom mapping result."""
    old_atom_name: str
    new_atom_name: str
    element: str
    residue_ccd: str
    confidence: float  # 0.0–1.0 matching score
    reason: str  # e.g., "topology_match", "3d_proximity", "ccd_lookup"


@dataclass
class QCFlag:
    """Single QC check result."""
    check_name: str
    status: QCStatus
    message: str
    data: Optional[Dict] = None


@dataclass
class PoseMetrics:
    """All metrics for a single pose (after full analysis)."""
    run_id: str
    pose_id: int
    cluster_id: Optional[int] = None
    
    # QC
    qc_flags: List[QCFlag] = field(default_factory=list)
    hard_fail: bool = False
    failure_reason: Optional[FailureReason] = None
    
    # Geometry (3D structures)
    cu_c1_distance_a: Optional[float] = None  # Ångström
    cu_c4_distance_a: Optional[float] = None
    his_brace_angle_deg: Optional[float] = None
    substrate_planarity_score: Optional[float] = None  # 0–1 (1 = planar)
    
    # Confidence
    atom_b_iso_mean: Optional[float] = None
    residue_qa_metric_median: Optional[float] = None
    
    # IFP (interaction fingerprint)
    ifp_vector: Optional[List[int]] = None  # binary {0,1}
    ifp_residues: Optional[List[str]] = None  # residue IDs (for index mapping)
    
    # Clustering
    cluster_distance_to_medoid: Optional[float] = None  # Jaccard/Tanimoto
    
    # Crystal anchoring (optional)
    crystal_ifp_similarity: Optional[float] = None  # Tanimoto
    crystal_pocket_rmsd_a: Optional[float] = None
    
    def to_csv_row(self) -> Dict:
        """Flatten to CSV-friendly dict."""
        return {
            "run_id": self.run_id,
            "pose_id": self.pose_id,
            "cluster_id": self.cluster_id,
            "qc_hard_fail": self.hard_fail,
            "failure_reason": self.failure_reason.value if self.failure_reason else None,
            "cu_c1_dist_a": self.cu_c1_distance_a,
            "cu_c4_dist_a": self.cu_c4_distance_a,
            "his_angle_deg": self.his_brace_angle_deg,
            "planarity": self.substrate_planarity_score,
            "b_iso_mean": self.atom_b_iso_mean,
            "qa_metric_median": self.residue_qa_metric_median,
            "ifp_jaccard_to_medoid": self.cluster_distance_to_medoid,
            "crystal_ifp_sim": self.crystal_ifp_similarity,
            "crystal_pocket_rmsd_a": self.crystal_pocket_rmsd_a,
        }


@dataclass
class ClusterSignature:
    """Definition of a single cluster."""
    cluster_id: int
    protein_id: str
    ligand_id: str
    
    # Composition
    pose_ids: List[str] = field(default_factory=list)
    occupancy: float = 0.0  # fraction of total poses
    outlier_rate: Optional[float] = None  # fraction of outliers (label=-1)
    
    # Representative
    medoid_pose_id: Optional[str] = None
    medoid_geometry: Optional[Dict] = None  # {cu_c1, cu_c4, angle, ...}
    
    # Aggregate geometry
    cu_c1_mean_a: Optional[float] = None
    cu_c1_std_a: Optional[float] = None
    cu_c4_mean_a: Optional[float] = None
    cu_c4_std_a: Optional[float] = None
    his_angle_mean_deg: Optional[float] = None
    planarity_mean: Optional[float] = None
    
    # Dominant interactions (top 3–5)
    dominant_interactions: List[Tuple[str, float]] = field(default_factory=list)  # (interaction, frequency)
    
    # Mechanistic tag (inferred)
    mechanism_tag: Optional[str] = None  # e.g., "C1_like", "C4_like", "mixed"


@dataclass
class ActivityFeature:
    """Protein-level aggregated features (for RQ1–RQ4)."""
    protein_id: str
    activity_class: Optional[str] = None  # C1, C4, mixed (from annotation)
    substrate_specificity: Optional[str] = None  # chitin, cellulose, amylose, etc.
    
    # Cluster occupancies
    c1_like_occupancy: float = 0.0
    c4_like_occupancy: float = 0.0
    mixed_occupancy: float = 0.0
    
    # Signature interactions per substrate
    interaction_freq_by_substrate: Dict[str, Dict[str, float]] = field(default_factory=dict)
    # Example: {"chitin": {"Asp_His_contact": 0.8, "Tyr_stacking": 0.6}, ...}
    
    # CBM contribution (for DEL B)
    cbm_presence: bool = False
    cbm_ligand_distance_a: Optional[float] = None
    cbm_cu_distance_a: Optional[float] = None


@dataclass
class TuningResult:
    """Result of one parameter-sweep point."""
    parameter_hash: str
    model: str
    param_point: Dict  # {num_recycles, num_seeds, ...}
    
    # QC metrics
    posebusters_pass_rate: float
    privateer_pass_rate: float
    
    # Cluster stability
    outlier_rate: float
    n_clusters: int
    cluster_occupancy_std: float  # std dev of cluster sizes (lower = more stable)
    
    # Crystal sanity (if available)
    crystal_ifp_mean_similarity: Optional[float] = None
    crystal_pocket_rmsd_mean: Optional[float] = None
    
    # Decision
    rank: Optional[int] = None  # lower is better
    decision_notes: str = ""


@dataclass
class RunManifest:
    """Full reproducibility record written at end of every run."""
    pipeline_version: str
    timestamp: str
    mode: str  # tune | production
    
    # Reproducibility
    git_commit: Optional[str] = None
    config_hash: Optional[str] = None
    
    # Tool versions (dict of str → str)
    tool_versions: Dict[str, str] = field(default_factory=dict)
    
    # Seeds (all sources of randomness)
    seeds: Dict[str, any] = field(default_factory=dict)
    
    # Input checksums (sha256)
    input_checksums: Dict[str, str] = field(default_factory=dict)
    
    # Gate outcomes (pass/fail)
    gates_passed: Dict[str, bool] = field(default_factory=dict)
    
    # Tuning reference (if production mode)
    tuning_reference: Optional[str] = None


# Registry of CCD monosaccharides (whitelist for Privateer prep)
CCD_MONOSACCHARIDES = {
    "NAG",  # N-acetyl-D-glucosamine (chitin)
    "BGC",  # beta-D-glucose (cellulose)
    "GLC",  # D-glucose
    "MAN",  # D-mannose
    "GLA",  # D-glucuronic acid
    "ARA",  # L-arabinose
    "XYS",  # D-xylose
    "FUC",  # L-fucose
    "GAL",  # D-galactose
    "GAC",  # N-acetyl-D-galactosamine
    "SIA",  # neuraminic acid
    # Add more as needed
}
