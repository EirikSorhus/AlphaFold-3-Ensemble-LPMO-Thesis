# scripts/module_2/config_improved.py
"""
Configuration file for Module 2: LPMO sequence fetching pipeline.

This module defines all constants, file paths, and configuration parameters
used by the LPMO fetching pipeline.
"""

import os
import sys
import json
from pathlib import Path
from typing import Optional, List
from datetime import datetime


# ============================================================================
# DEPENDENCY CHECK
# ============================================================================

def check_dependencies() -> bool:
    """
    Check if all required Python packages are installed.
    
    Returns:
        bool: True if all dependencies are available, False otherwise
    """
    required_packages = {
        'bioservices': 'bioservices',
        'Bio': 'biopython',
        'bs4': 'beautifulsoup4',
        'requests': 'requests'
    }
    
    missing = []
    for module_name, package_name in required_packages.items():
        try:
            __import__(module_name)
        except ImportError:
            missing.append(package_name)
    
    if missing:
        print("[CONFIG] ERROR: Missing required packages:")
        for pkg in missing:
            print(f"  - {pkg}")
        print("\nInstall with: pip install " + " ".join(missing))
        return False
    
    print("[CONFIG] All required packages are installed")
    return True


# ============================================================================
# PATH CONFIGURATION
# ============================================================================

# Determine project root (2 levels up from this file: module_2 -> scripts -> project_root)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RAW_DIR = PROJECT_ROOT / "data_raw"  # Legacy, kept for compatibility
DATA_SEQUENCES_DIR = PROJECT_ROOT / "data" / "sequences"
DATA_METADATA_DIR = PROJECT_ROOT / "data" / "metadata"

# Output file paths (legacy, not used in improved version)
UNIPROT_FASTA = DATA_RAW_DIR / "uniprot_LPMO_raw.fasta"
CAZY_FASTA = DATA_RAW_DIR / "cazy_LPMO_raw.fasta"
MERGED_FASTA = DATA_RAW_DIR / "LPMO_all_raw.fasta"


# ============================================================================
# DATABASE QUERY PARAMETERS
# ============================================================================

# Default CAZy families for LPMOs (Lytic Polysaccharide Monooxygenases)
DEFAULT_CAZY_FAMILIES = ["AA9", "AA10", "AA11", "AA13"]


def build_uniprot_query(families: List[str]) -> str:
    """
    Build UniProt query string for given protein families.
    
    Args:
        families: List of CAZy family names (e.g., ['AA9', 'AA10'])
        
    Returns:
        str: UniProt query string
    """
    if not families:
        raise ValueError("At least one family must be specified")
    
    if len(families) == 1:
        return f'(family:"{families[0]}")'
    
    family_queries = [f'family:"{fam}"' for fam in families]
    return f'({" OR ".join(family_queries)})'


def get_run_id() -> str:
    """
    Generate a unique run ID based on current timestamp.
    Format: YYYYMMDD_HHMMSS
    
    Returns:
        str: Run ID string
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def get_family_output_filename(family: str, run_id: str = None, source: str = "") -> str:
    """
    Generate output filename for a specific family.
    Format: lpmo_AA9_raw.fasta (stable naming without timestamp)
    
    Args:
        family: CAZy family name (e.g., 'AA9')
        run_id: Run ID from get_run_id() (ignored for compatibility)
        source: Optional source identifier (ignored for compatibility)
        
    Returns:
        str: Filename in format lpmo_{FAMILY}_raw.fasta
    """
    return f"lpmo_{family}_raw.fasta"


def get_merged_output_filename(run_id: str = None) -> str:
    """
    Generate merged output filename.
    Format: lpmo_all_raw.fasta (stable name)
    
    Args:
        run_id: Run ID from get_run_id() (ignored, kept for compatibility)
        
    Returns:
        str: Filename
    """
    return "lpmo_all_raw.fasta"


def get_metadata_filename(run_id: str = None) -> str:
    """
    Generate metadata filename for a run.
    Format: m2_run_metadata.json (stable name)
    
    Args:
        run_id: Run ID from get_run_id() (ignored, kept for compatibility)
        
    Returns:
        str: Filename
    """
    return "m2_run_metadata.json"


def save_metadata(
    run_id: str,
    families: List[str],
    uniprot_query: str,
    output_dir: Path,
    additional_info: Optional[dict] = None
) -> Path:
    """
    Save metadata file for a pipeline run.
    
    Args:
        run_id: Run ID for this run
        families: List of families processed
        uniprot_query: UniProt query string used
        output_dir: Directory to save metadata file
        additional_info: Optional dictionary with additional metadata
        
    Returns:
        Path: Path to created metadata file
    """
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "date": datetime.now().strftime("%Y-%m-%d"),
        "time": datetime.now().strftime("%H:%M:%S"),
        "families": families,
        "uniprot_query": uniprot_query,
        "num_families": len(families)
    }
    
    if additional_info:
        metadata.update(additional_info)
    
    metadata_file = output_dir / get_metadata_filename(run_id)
    
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    
    print(f"[CONFIG] Metadata saved to: {metadata_file}")
    return metadata_file


# ============================================================================
# API CONFIGURATION
# ============================================================================

# NCBI requires a valid email address for API access
# Priority: 1) Environment variable, 2) Parameter passed to functions
# See: https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
def get_ncbi_email(email: Optional[str] = None) -> str:
    """
    Get NCBI email from parameter, environment variable, or raise error.
    
    Args:
        email: Optional email address to use
        
    Returns:
        str: Valid email address
        
    Raises:
        ValueError: If no valid email is provided
    """
    if email and email != "din.epost@institusjon.no":
        return email
    
    env_email = os.getenv("NCBI_EMAIL")
    if env_email and env_email != "din.epost@institusjon.no":
        return env_email
    
    raise ValueError(
        "NCBI email not configured. Please either:\n"
        "  1) Set environment variable: export NCBI_EMAIL='your.email@institution.no'\n"
        "  2) Pass email parameter to fetch functions\n"
        "NCBI requires a valid email address for API access."
    )


# ============================================================================
# NETWORK CONFIGURATION
# ============================================================================

# HTTP request timeout (seconds)
REQUEST_TIMEOUT = 30

# NCBI rate limiting (requests per second)
# Without API key: max 3 requests/second
# With API key: max 10 requests/second
NCBI_REQUESTS_PER_SECOND = 3
NCBI_BATCH_SIZE = 100  # Number of sequences to fetch per request


# ============================================================================
# VALIDATION
# ============================================================================

def validate_config(families: Optional[List[str]] = None) -> bool:
    """
    Validate configuration settings.
    
    Args:
        families: Optional list of families to validate
    
    Returns:
        bool: True if configuration is valid
    """
    issues = []
    
    # Check that PROJECT_ROOT exists
    if not PROJECT_ROOT.exists():
        issues.append(f"Project root does not exist: {PROJECT_ROOT}")
    
    # Check families if provided
    if families is not None:
        if not families:
            issues.append("Families list is empty")
        else:
            # Validate each family name
            for fam in families:
                if not fam or not isinstance(fam, str):
                    issues.append(f"Invalid family name: {fam}")
    
    if issues:
        print("[CONFIG] Configuration validation failed:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    
    print("[CONFIG] Configuration is valid")
    return True


# ============================================================================
# INITIALIZATION
# ============================================================================

def ensure_data_directory():
    """Ensure all data directories exist."""
    DATA_RAW_DIR.mkdir(parents=True, exist_ok=True)
    DATA_SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)
    DATA_METADATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[CONFIG] Data directories ready:")
    print(f"[CONFIG]   - Sequences: {DATA_SEQUENCES_DIR}")
    print(f"[CONFIG]   - Metadata: {DATA_METADATA_DIR}")


if __name__ == "__main__":
    print("=== Module 2 Configuration ===")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data directory: {DATA_RAW_DIR}")
    print(f"Default CAZy families: {DEFAULT_CAZY_FAMILIES}")
    print(f"Example UniProt query: {build_uniprot_query(DEFAULT_CAZY_FAMILIES)}")
    print()
    
    # Run checks
    deps_ok = check_dependencies()
    config_ok = validate_config(DEFAULT_CAZY_FAMILIES)
    
    # Test helper functions
    print("\n=== Testing helper functions ===")
    test_run_id = get_run_id()
    print(f"Run ID: {test_run_id}")
    print(f"Family file example: {get_family_output_filename('AA9', test_run_id, 'uniprot')}")
    print(f"Merged file: {get_merged_output_filename(test_run_id)}")
    print(f"Metadata file: {get_metadata_filename(test_run_id)}")
    
    if deps_ok and config_ok:
        print("\n✓ Configuration is ready")
        sys.exit(0)
    else:
        print("\n✗ Configuration has issues")
        sys.exit(1)
