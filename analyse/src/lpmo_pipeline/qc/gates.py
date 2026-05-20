"""
LPMO Pipeline: QC Gates & Threshold Enforcement
Responsibility: Centralized pass/fail decision logic for all gates

Hard gates (must pass for pose to continue):
  1. atom_mapping_coverage = 100%
    2. active_site_proximity <= threshold
    3. privateer_recognized_sugars = 100%
    4. no_critical_posebusters_errors = true
    5. cu_his_distance ∈ [1.9, 2.6] Å
  
Soft flags (kept but marked):
  - PoseBusters minor warnings
  - Outlier poses in clustering
  - Crystal similarity < 0.3 (flag for alternative mode)
"""

from dataclasses import dataclass
from typing import Dict, Optional, List
from pathlib import Path

from lpmo_pipeline.config import load_runtime_paths_config


_RUNTIME_PATHS = load_runtime_paths_config()
_DEFAULT_GATE_CONFIG_PATH = _RUNTIME_PATHS.pipeline_assets.thresholds_config


@dataclass
class GateConfig:
    """Configuration for all QC gates."""
    # Atom mapping
    atom_mapping_coverage_min: float = 1.0  # 100%

    # Pre-QC active-site proximity
    active_site_proximity_max_a: float = 10.0
    
    # Privateer
    privateer_recognized_min: float = 1.0  # 100%
    
    # PoseBusters (critical errors that cause hard fail)
    posebusters_critical_errors: List[str] = None  # e.g., ["steric_clash", "valence_error"]
    
    # Cu geometry hard gate
    cu_his_dist_min: float = 1.9  # Ångström
    cu_his_dist_max: float = 2.6

    # Cu geometry soft QC band
    cu_his_soft_min_a: float = 1.8
    cu_his_soft_max_a: float = 2.6

    cu_c_proximity_threshold_a: float = 7.0
    his_brace_max_search_a: float = 3.0
    
    # Crystal anchoring (soft threshold)
    crystal_ifp_similarity_soft_threshold: float = 0.3
    
    # Clustering
    hdbscan_min_cluster_size: int = 10
    hdbscan_min_pose_occupancy: float = 0.01  # 1% of poses in cluster
    
    def __post_init__(self):
        if self.posebusters_critical_errors is None:
            self.posebusters_critical_errors = [
                "sanitization",
                "all_atoms_connected",
                "no_radicals",
                "internal_steric_clash",
                "bond_lengths",
                "bond_angles",
                "tetrahedral_chirality",
                "protein-ligand_maximum_distance",
                "minimum_distance_to_protein",
                "volume_overlap_with_protein",
            ]


@dataclass
class GateResult:
    """Result of one gate check."""
    gate_name: str
    passed: bool
    message: str
    threshold: Optional[float] = None
    measured_value: Optional[float] = None
    data: Optional[Dict] = None


class QCGateChecker:
    """Evaluate pose against all QC gates."""
    
    def __init__(self, config: Optional[GateConfig] = None):
        self.config = config or load_gate_config_from_yaml(_DEFAULT_GATE_CONFIG_PATH)
        self.results: List[GateResult] = []
    
    def check_atom_mapping(self, coverage: float) -> GateResult:
        """Gate 1: atom_mapping_coverage."""
        passed = coverage >= self.config.atom_mapping_coverage_min
        result = GateResult(
            gate_name="atom_mapping_coverage",
            passed=passed,
            message=f"Coverage: {coverage:.1%}",
            threshold=self.config.atom_mapping_coverage_min,
            measured_value=coverage,
        )
        self.results.append(result)
        return result

    def check_active_site_proximity(self, min_distance_a: float) -> GateResult:
        """Gate 2: active-site proximity before PoseBusters/Privateer."""
        passed = min_distance_a <= self.config.active_site_proximity_max_a
        result = GateResult(
            gate_name="active_site_proximity_range",
            passed=passed,
            message=(
                f"Min Cu-ligand distance: {min_distance_a:.2f} Å "
                f"(max: {self.config.active_site_proximity_max_a:.2f} Å)"
            ),
            threshold=self.config.active_site_proximity_max_a,
            measured_value=min_distance_a,
        )
        self.results.append(result)
        return result
    
    def check_privateer_sugars(self, n_recognized: int, n_total: int) -> GateResult:
        """Gate 2: privateer_recognized_sugars."""
        if n_total == 0:
            # No sugars (protein only) → pass
            result = GateResult(
                gate_name="privateer_recognized_sugars",
                passed=True,
                message="No glycans in structure (pass)",
            )
        else:
            recognition_rate = n_recognized / n_total
            passed = recognition_rate >= self.config.privateer_recognized_min
            result = GateResult(
                gate_name="privateer_recognized_sugars",
                passed=passed,
                message=f"Recognition rate: {recognition_rate:.1%}",
                threshold=self.config.privateer_recognized_min,
                measured_value=recognition_rate,
                data={"recognized": n_recognized, "total": n_total}
            )
        
        self.results.append(result)
        return result
    
    def check_posebusters(self, errors: Dict[str, int]) -> GateResult:
        """Gate 3: no_critical_posebusters_errors."""
        critical_errors = {k: v for k, v in errors.items() 
                          if k in self.config.posebusters_critical_errors}
        
        passed = len(critical_errors) == 0
        result = GateResult(
            gate_name="no_critical_posebusters_errors",
            passed=passed,
            message=f"Critical errors: {list(critical_errors.keys())}",
            data={"all_errors": errors, "critical": critical_errors}
        )
        self.results.append(result)
        return result
    
    def check_cu_his_distance(self, distance_a: float) -> GateResult:
        """Gate 4: cu_his_distance in range."""
        passed = self.config.cu_his_dist_min <= distance_a <= self.config.cu_his_dist_max
        result = GateResult(
            gate_name="cu_his_distance_range",
            passed=passed,
            message=f"Distance: {distance_a:.2f} Å (range: {self.config.cu_his_dist_min}–{self.config.cu_his_dist_max})",
            threshold=f"[{self.config.cu_his_dist_min}, {self.config.cu_his_dist_max}]",
            measured_value=distance_a,
        )
        self.results.append(result)
        return result
    
    def soft_flag_crystal_similarity(self, similarity: Optional[float]) -> GateResult:
        """Soft flag: crystal_ifp_similarity (not hard fail)."""
        if similarity is None:
            result = GateResult(
                gate_name="crystal_ifp_similarity",
                passed=True,
                message="No crystal available (skip check)",
            )
        else:
            # Not a hard gate, just flag if low
            flag = similarity < self.config.crystal_ifp_similarity_soft_threshold
            result = GateResult(
                gate_name="crystal_ifp_similarity",
                passed=not flag,  # inverted: "pass" means similarity acceptable
                message=f"Similarity: {similarity:.2f} (threshold: {self.config.crystal_ifp_similarity_soft_threshold})",
                threshold=self.config.crystal_ifp_similarity_soft_threshold,
                measured_value=similarity,
            )
        
        self.results.append(result)
        return result
    
    def all_hard_gates_passed(self) -> bool:
        """Check if all hard gates passed."""
        hard_gate_names = [
            "atom_mapping_coverage",
            "active_site_proximity_range",
            "privateer_recognized_sugars",
            "no_critical_posebusters_errors",
            "cu_his_distance_range",
        ]
        
        hard_results = [r for r in self.results if r.gate_name in hard_gate_names]
        return all(r.passed for r in hard_results)
    
    def summary(self) -> Dict:
        """Generate summary of all gate checks."""
        hard_gates = [
            "atom_mapping_coverage",
            "active_site_proximity_range",
            "privateer_recognized_sugars",
            "no_critical_posebusters_errors",
            "cu_his_distance_range",
        ]
        
        summary = {
            "total_checks": len(self.results),
            "hard_gates_passed": self.all_hard_gates_passed(),
            "checks": [
                {
                    "gate": r.gate_name,
                    "passed": r.passed,
                    "message": r.message,
                    "is_hard": r.gate_name in hard_gates,
                    "measured_value": r.measured_value,
                    "threshold": r.threshold,
                }
                for r in self.results
            ]
        }
        
        return summary


# Load default gate config from YAML if available
def load_gate_config_from_yaml(config_path: Path) -> GateConfig:
    """Load gate thresholds from YAML config file."""
    import yaml
    try:
        with open(config_path) as f:
            data = yaml.safe_load(f)

        hard = data.get("gates") or data.get("hard_gates") or {}
        soft = data.get("soft_thresholds") or data.get("soft") or {}
        qc = data.get("qc") or {}
        geometry_rules = data.get("geometry_rules") or {}
        his_brace = geometry_rules.get("his_brace") or {}
        clustering = data.get("clustering") or data.get("hdbscan") or {}

        # Map YAML keys to GateConfig fields (supports legacy and current schema)
        return GateConfig(
            atom_mapping_coverage_min=hard.get("atom_mapping_coverage_min", 1.0),
            active_site_proximity_max_a=hard.get(
                "active_site_proximity_max_a",
                qc.get("pre_qc_active_site_max_a", 10.0),
            ),
            privateer_recognized_min=hard.get(
                "privateer_recognized_min",
                hard.get("privateer_recognized_sugars_min", 1.0),
            ),
            posebusters_critical_errors=data.get("posebusters_critical_errors"),
            cu_his_dist_min=hard.get("cu_his_dist_min", hard.get("cu_his_distance_min_a", 1.9)),
            cu_his_dist_max=hard.get("cu_his_dist_max", hard.get("cu_his_distance_max_a", 2.6)),
            cu_his_soft_min_a=qc.get("cu_his_min_a", hard.get("cu_his_distance_min_a", 1.9)),
            cu_his_soft_max_a=qc.get("cu_his_max_a", hard.get("cu_his_distance_max_a", 2.6)),
            cu_c_proximity_threshold_a=soft.get("cu_c_proximity_max_a", 7.0),
            his_brace_max_search_a=his_brace.get("max_search_dist_a", 3.0),
            crystal_ifp_similarity_soft_threshold=soft.get(
                "crystal_ifp_similarity_soft_threshold",
                soft.get("crystal_ifp_similarity_min", 0.3),
            ),
            hdbscan_min_cluster_size=clustering.get("min_cluster_size", 10),
        )
    except FileNotFoundError:
        return GateConfig()
