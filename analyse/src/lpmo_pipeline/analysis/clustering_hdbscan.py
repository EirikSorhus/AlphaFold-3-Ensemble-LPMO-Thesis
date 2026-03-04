"""
LPMO Pipeline: Clustering with HDBSCAN
Responsibility: Cluster poses by IFP similarity, identify binding modes

Step 8 contract:
  Input: ifp_matrix.csv (all poses, single protein×ligand×model)
  Output: clusters.json (cluster assignments + medoids)
  
  Procedure:
    1. Phase 1: HDBSCAN within single run (compress noise)
    2. Phase 2: HDBSCAN cross-run per protein×ligand (aggregate medoid IFPs)
  
  HDBSCAN params LOCKED post-tuning (anti p-hack):
    - min_cluster_size: determined from tuning
    - metric: 'jaccard' (for binary IFP)
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
from dataclasses import dataclass, asdict
import numpy as np

try:
    import hdbscan
except ImportError:
    hdbscan = None

from lpmo_pipeline.utils.logging import StructuredLogger
from lpmo_pipeline.utils.data_models import ClusterSignature


@dataclass
class ClusteringResult:
    """Outcome of clustering procedure."""
    n_clusters: int
    n_outliers: int
    outlier_rate: float
    cluster_sizes: Dict[int, int]  # cluster_id → n_poses
    cluster_occupancy: Dict[int, float]  # cluster_id → fraction of total
    outlier_indices: List[int]
    cluster_labels: np.ndarray  # labels for all poses (-1 for outliers)


class HDBANSCANClusterer:
    """Cluster poses by IFP similarity (binary Jaccard distance)."""
    
    def __init__(self, min_cluster_size: int, metric: str = "jaccard", output_dir: Optional[Path] = None):
        """
        Args:
            min_cluster_size: minimum cluster size (locked from tuning)
            metric: distance metric ('jaccard' for binary IFP)
            output_dir: for logging
        """
        self.min_cluster_size = min_cluster_size
        self.metric = metric
        self.output_dir = Path(output_dir) if output_dir else Path(".")
        self.logger = StructuredLogger("clustering", self.output_dir)
        
        if hdbscan is None:
            raise ImportError("HDBSCAN not installed")
    
    def cluster(self, ifp_matrix: np.ndarray, pose_ids: List[str]) -> ClusteringResult:
        """
        Cluster poses using HDBSCAN.
        
        Args:
            ifp_matrix: shape (n_poses, n_residues × n_interaction_types)
                        binary (0/1) values
            pose_ids: list of pose identifiers
        
        Returns:
            ClusteringResult with labels, cluster info, etc.
        """
        self.logger.log_step_start("hdbscan_clustering", {
            "n_poses": ifp_matrix.shape[0],
            "ifp_dim": ifp_matrix.shape[1],
            "min_cluster_size": self.min_cluster_size,
            "metric": self.metric,
        })
        
        # Compute Jaccard distance matrix
        # For binary vectors: Jaccard(u, v) = 1 - (u·v) / (||u||² + ||v||² - u·v)
        distance_matrix = self._compute_jaccard_distances(ifp_matrix)
        
        # Run HDBSCAN (use precomputed distance + linkage)
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=self.min_cluster_size,
            min_samples=None,
            metric='precomputed',
            cluster_selection_epsilon=0.0,
            cluster_selection_method='eom',  # excess of mass
        )
        
        labels = clusterer.fit_predict(distance_matrix)  # -1 for outliers
        
        # Aggregate cluster info
        result = self._aggregate_clusters(labels, pose_ids)
        
        # Log results
        self.logger.log_step_end("hdbscan_clustering", "success", {
            "n_clusters": result.n_clusters,
            "n_outliers": result.n_outliers,
            "outlier_rate": result.outlier_rate,
        }, 0.0)
        
        return result
    
    def _compute_jaccard_distances(self, ifp_matrix: np.ndarray) -> np.ndarray:
        """
        Compute pairwise Jaccard distance for binary IFP vectors.
        
        Jaccard distance: 1 - (intersection / union)
        For binary: 1 - (u·v) / (||u||² + ||v||² - u·v)
        """
        n = ifp_matrix.shape[0]
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                u = ifp_matrix[i]
                v = ifp_matrix[j]
                
                intersection = np.dot(u, v)
                union = np.sum((u + v) > 0)  # at least one is 1
                
                if union == 0:
                    jaccard_sim = 1.0  # both empty vectors
                else:
                    jaccard_sim = intersection / union
                
                jaccard_dist = 1.0 - jaccard_sim
                distances[i, j] = jaccard_dist
                distances[j, i] = jaccard_dist
        
        return distances
    
    def _aggregate_clusters(self, labels: np.ndarray, pose_ids: List[str]) -> ClusteringResult:
        """
        Aggregate clustering results.
        """
        unique_labels = set(labels)
        n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
        n_outliers = np.sum(labels == -1)
        outlier_rate = n_outliers / len(labels) if len(labels) > 0 else 0.0
        
        cluster_sizes = {}
        cluster_occupancy = {}
        
        for label in unique_labels:
            if label == -1:  # outliers
                continue
            
            size = np.sum(labels == label)
            cluster_sizes[int(label)] = int(size)
            cluster_occupancy[int(label)] = size / len(labels)
        
        outlier_indices = np.where(labels == -1)[0].tolist()
        
        result = ClusteringResult(
            n_clusters=n_clusters,
            n_outliers=n_outliers,
            outlier_rate=outlier_rate,
            cluster_sizes=cluster_sizes,
            cluster_occupancy=cluster_occupancy,
            outlier_indices=outlier_indices,
            cluster_labels=labels,
        )
        
        return result
    
    def select_medoids(self, ifp_matrix: np.ndarray, labels: np.ndarray) -> Dict[int, int]:
        """
        For each cluster, select the medoid (most central pose).
        
        Args:
            ifp_matrix: shape (n_poses, n_features)
            labels: cluster labels (-1 for outliers)
        
        Returns:
            {cluster_id: medoid_index}
        """
        medoids = {}
        
        for label in set(labels):
            if label == -1:
                continue
            
            cluster_indices = np.where(labels == label)[0]
            cluster_data = ifp_matrix[cluster_indices]
            
            # Medoid = pose closest to cluster center (in Jaccard space)
            # Simplified: use centroid + find nearest
            # Real impl: compute all pairwise distances within cluster, find minimum-sum
            centroid = np.mean(cluster_data, axis=0)
            
            # Find pose with minimum euclidean distance to centroid
            distances_to_centroid = np.linalg.norm(cluster_data - centroid, axis=1)
            medoid_in_cluster = np.argmin(distances_to_centroid)
            medoid_global_index = cluster_indices[medoid_in_cluster]
            
            medoids[int(label)] = int(medoid_global_index)
        
        return medoids


class CrossRunClusterer:
    """Phase 2: Cluster aggregated IFPs across runs (protein×ligand)."""
    
    def __init__(self, min_cluster_size: int, output_dir: Optional[Path] = None):
        """
        Args:
            min_cluster_size: minimum cluster size (locked from tuning)
            output_dir: for logging
        """
        self.min_cluster_size = min_cluster_size
        self.output_dir = Path(output_dir) if output_dir else Path(".")
        self.logger = StructuredLogger("crossrun_clustering", self.output_dir)
    
    def aggregate_medoid_ifps(self, 
                              phase1_results: List[Tuple[str, ClusteringResult, np.ndarray]],
                              ) -> np.ndarray:
        """
        Aggregate medoid IFPs from phase-1 clustering.
        
        Args:
            phase1_results: list of (run_id, clustering_result, ifp_matrix)
        
        Returns:
            aggregated_ifp_matrix: shape (n_medoids, n_features)
        """
        medoid_ifps = []
        
        for run_id, clustering_result, ifp_matrix in phase1_results:
            # Identify medoids for this run
            for cluster_id, medoid_idx in clustering_result.medoids.items():
                medoid_ifp = ifp_matrix[medoid_idx]
                medoid_ifps.append(medoid_ifp)
        
        return np.array(medoid_ifps)
    
    def cluster_cross_run(self, aggregated_ifps: np.ndarray) -> ClusteringResult:
        """
        Cluster aggregated medoid IFPs across all runs.
        
        Args:
            aggregated_ifps: medoid IFP vectors
        
        Returns:
            ClusteringResult for cross-run binding modes
        """
        clusterer = HDBANSCANClusterer(
            min_cluster_size=self.min_cluster_size,
            output_dir=self.output_dir
        )
        
        pose_ids = [f"medoid_{i}" for i in range(len(aggregated_ifps))]
        return clusterer.cluster(aggregated_ifps, pose_ids)
