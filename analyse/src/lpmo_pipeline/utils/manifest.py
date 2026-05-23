"""
LPMO Pipeline Run Manifest & Reproducibility
Responsibility: Create + serialize run_manifest.json at end of every run

Manifest contents:
  - Pipeline version, timestamp, mode (tune/production)
  - Git commit hash, config hash
  - Tool versions (Gemmi, PLACER, PoseBusters, Privateer, MDAnalysis, ProLIF, HDBSCAN)
  - All seeds (AF3, RF3, Boltz2, HDBSCAN if applicable)
  - Input checksums (sha256 of config, protein seq, ligand)
  - Gates passed (atom_mapping=100%, privateer_recognized=100%, etc.)
  - Tuning reference (if production mode)
"""

import json
from pathlib import Path
from typing import Dict, Optional
from dataclasses import asdict
import hashlib
from datetime import datetime
import subprocess

from lpmo_pipeline.io.gemmi_compat import gemmi_version
from lpmo_pipeline.qc.privateer_runner import get_privateer_version
from lpmo_pipeline.utils.data_models import RunManifest


def _probe_gemmi() -> str:
    """Run a lightweight functional probe against the active Gemmi binding."""
    try:
        from lpmo_pipeline.io.gemmi_compat import gemmi

        if not hasattr(gemmi, "read_structure") or not hasattr(gemmi, "cif"):
            return "missing_required_api"

        doc = gemmi.cif.read_string("data_probe\n_entry.id probe\n")
        block = doc.sole_block()
        if block.find_value("_entry.id") != "probe":
            return "unexpected_probe_result"
    except Exception as exc:
        return f"{exc.__class__.__name__}: {exc}"

    return "ok"


def _probe_mdanalysis() -> str:
    """Run a minimal in-memory Universe probe for MDAnalysis."""
    try:
        import numpy as np
        import MDAnalysis as mda

        universe = mda.Universe.empty(
            1,
            n_residues=1,
            atom_resindex=np.array([0], dtype=int),
            trajectory=True,
        )
        if universe.atoms.n_atoms != 1 or universe.residues.n_residues != 1:
            return "unexpected_probe_result"
        if universe.trajectory.ts.positions.shape != (1, 3):
            return "unexpected_probe_result"
    except Exception as exc:
        return f"{exc.__class__.__name__}: {exc}"

    return "ok"


def _probe_prolif() -> str:
    """Confirm ProLIF can instantiate the configured fingerprint surface."""
    try:
        import prolif as plf

        expected = {"ImplicitHBAcceptor", "ImplicitHBDonor", "VdWContact"}
        fingerprint = plf.Fingerprint(sorted(expected), count=True)
        interaction_names = set(fingerprint.interactions.keys())
        if expected - interaction_names:
            return "unexpected_probe_result"
    except Exception as exc:
        return f"{exc.__class__.__name__}: {exc}"

    return "ok"


def _probe_hdbscan() -> str:
    """Run a tiny deterministic clustering job to catch broken native wheels."""
    try:
        import numpy as np
        import hdbscan

        probe_points = np.array(
            [
                [0.0, 0.0],
                [0.0, 0.01],
                [0.01, 0.0],
                [5.0, 5.0],
                [5.0, 5.01],
                [5.01, 5.0],
            ],
            dtype=float,
        )
        labels = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1).fit_predict(probe_points)
        if len(labels) != len(probe_points):
            return "unexpected_probe_result"
        if all(label == -1 for label in labels):
            return "unexpected_probe_result"
    except Exception as exc:
        return f"{exc.__class__.__name__}: {exc}"

    return "ok"


class ManifestBuilder:
    """Build run_manifest.json incrementally throughout execution."""
    
    def __init__(self, mode: str, version: str):
        """
        Args:
            mode: "tune" or "production"
            version: pipeline version string
        """
        self.mode = mode
        self.version = version
        self.manifest = RunManifest(
            pipeline_version=version,
            timestamp=datetime.utcnow().isoformat(),
            mode=mode,
        )
    
    def set_git_info(self):
        """Retrieve git commit hash (if in git repo)."""
        try:
            commit = subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                text=True,
                stderr=subprocess.DEVNULL
            ).strip()
            self.manifest.git_commit = commit[:8]  # short hash
        except Exception:
            self.manifest.git_commit = "unknown"
    
    def set_config_hash(self, config_dict: Dict):
        """Hash the full config for reproducibility."""
        config_json = json.dumps(config_dict, sort_keys=True, default=str)
        h = hashlib.sha256(config_json.encode()).hexdigest()
        self.manifest.config_hash = h
    
    def set_tool_versions(self, versions: Dict[str, str]):
        """Set external tool versions."""
        self.manifest.tool_versions = versions
    
    def add_seed(self, source: str, seed_value):
        """Record a seed (e.g., "af3_seeds", "rf3_seed", "hdbscan_seed")."""
        self.manifest.seeds[source] = seed_value
    
    def add_input_checksum(self, name: str, file_path: Path):
        """Compute SHA256 checksum of an input file."""
        sha256_hash = hashlib.sha256()
        try:
            with open(file_path, "rb") as f:
                for byte_block in iter(lambda: f.read(4096), b""):
                    sha256_hash.update(byte_block)
            self.manifest.input_checksums[name] = sha256_hash.hexdigest()
        except FileNotFoundError:
            self.manifest.input_checksums[name] = "file_not_found"
    
    def record_gate(self, gate_name: str, passed: bool):
        """Record if a pass/fail gate was satisfied."""
        self.manifest.gates_passed[gate_name] = passed
    
    def set_tuning_reference(self, ref_path: str):
        """If production mode, link to tuning result."""
        if self.mode == "production":
            self.manifest.tuning_reference = ref_path
    
    def to_dict(self) -> Dict:
        """Serialize to dict (for JSON writing)."""
        return asdict(self.manifest)
    
    def to_json_str(self) -> str:
        """Serialize to JSON string."""
        return json.dumps(self.to_dict(), indent=2, default=str)
    
    def write(self, output_path: Path):
        """Write manifest to file."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(self.to_json_str())
        print(f"Manifest written to {output_path}")


class ToolVersionFetcher:
    """Utility to query installed tool versions."""
    
    @staticmethod
    def get_versions() -> Dict[str, str]:
        """
        Query tool versions from their respective CLIs.
        This is a skeleton; implementations depend on tool-specific flags.
        """
        versions = {}
        
        # Gemmi (Python library version)
        try:
            versions["gemmi"] = gemmi_version()
        except Exception:
            versions["gemmi"] = "not_installed"
        
        # PLACER (external, query via CLI)
        try:
            result = subprocess.run(
                ["placer", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            versions["placer"] = result.stdout.strip() or "unknown_version"
        except Exception:
            versions["placer"] = "not_found"
        
        # PoseBusters
        try:
            import posebusters
            versions["posebusters"] = getattr(posebusters, "__version__", "unknown_version")
        except ImportError:
            versions["posebusters"] = "not_installed"
        
        # Privateer (external)
        try:
            versions["privateer"] = get_privateer_version()
        except Exception:
            versions["privateer"] = "not_found"
        
        # MDAnalysis
        try:
            import MDAnalysis
            versions["mdanalysis"] = getattr(MDAnalysis, "__version__", "unknown_version")
        except ImportError:
            versions["mdanalysis"] = "not_installed"
        
        # ProLIF
        try:
            import prolif
            versions["prolif"] = getattr(prolif, "__version__", "unknown_version")
        except ImportError:
            versions["prolif"] = "not_installed"
        
        # HDBSCAN
        try:
            import hdbscan
            versions["hdbscan"] = getattr(hdbscan, "__version__", "unknown_version")
        except ImportError:
            versions["hdbscan"] = "not_installed"
        
        return versions

    @staticmethod
    def get_analysis_dependency_statuses() -> Dict[str, str]:
        """Run lightweight functional probes for clustering-critical dependencies."""
        return {
            "gemmi": _probe_gemmi(),
            "mdanalysis": _probe_mdanalysis(),
            "prolif": _probe_prolif(),
            "hdbscan": _probe_hdbscan(),
        }
