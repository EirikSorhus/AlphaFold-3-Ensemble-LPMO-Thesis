"""
Taxonomy information extraction from CAZy data files.

This module extracts organism and kingdom information from CAZy .txt data files
and associates it with GenBank/protein IDs.

NAMING CONVENTION (intentionally inverted for consistency with pipeline):
  - Column 2 (Kingdom in CAZy) → stored as "organism" in metadata
    Examples: Eukaryota, Bacteria, Archaea
  
  - Column 3 (Organism/Species in CAZy) → stored as "species" in metadata
    Examples: Alternaria alternata, Saccharomyces cerevisiae
    
This inversion is deliberate to align with downstream module expectations.
For details, see module_2 README or TAXONOMY_MAPPING.md
"""

from typing import Dict, Tuple, Optional
from pathlib import Path


def parse_cazy_taxonomy_from_file(cazy_txt_path: Path, family: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse CAZy .txt data file to extract kingdom and organism for each protein ID.
    
    Args:
        cazy_txt_path: Path to CAZy data file (e.g., data/domains/dbcan/AA13.txt)
        family: Family name (e.g., 'AA13') for validation
        
    Returns:
        dict: Mapping of protein_id → (kingdom, species)
        Example: {"RYN72100.1": ("Eukaryota", "Alternaria alternata")}
    """
    taxonomy_map: Dict[str, Tuple[str, str]] = {}
    
    if not cazy_txt_path.exists():
        print(f"[TAXONOMY] WARNING: CAZy file not found: {cazy_txt_path}")
        return taxonomy_map
    
    try:
        with open(cazy_txt_path, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                # Skip header
                if line_num == 1 or line.startswith("#"):
                    continue
                
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                
                # Column structure: Family | Kingdom | Organism | Protein_ID | Source | ...
                cazy_family = parts[0].strip()
                kingdom = parts[1].strip()
                organism = parts[2].strip()
                protein_id = parts[3].strip()
                
                # Skip if not matching family
                if cazy_family != family:
                    continue
                
                # Store mapping (avoid duplicates, keep first occurrence)
                if protein_id and protein_id not in taxonomy_map:
                    taxonomy_map[protein_id] = (kingdom, organism)
        
        print(f"[TAXONOMY] Extracted taxonomy for {len(taxonomy_map)} proteins from {cazy_txt_path.name}")
        return taxonomy_map
        
    except Exception as e:
        print(f"[TAXONOMY] ERROR: Failed to parse CAZy file {cazy_txt_path}: {e}")
        return taxonomy_map


def get_taxonomy_for_protein(
    protein_id: str,
    taxonomy_map: Dict[str, Tuple[str, str]]
) -> Tuple[Optional[str], Optional[str]]:
    """
    Look up kingdom and species for a protein ID.
    
    Args:
        protein_id: GenBank protein ID (e.g., "RYN72100.1")
        taxonomy_map: Taxonomy mapping from parse_cazy_taxonomy_from_file()
        
    Returns:
        tuple: (kingdom, species) or (None, None) if not found
    """
    if protein_id in taxonomy_map:
        return taxonomy_map[protein_id]
    return None, None


if __name__ == "__main__":
    print("=== Taxonomy Info Module ===")
    # Test: parse a CAZy file
    test_file = Path("data/domains/dbcan/AA13.txt")
    if test_file.exists():
        tax_map = parse_cazy_taxonomy_from_file(test_file, "AA13")
        print(f"Found {len(tax_map)} entries")
        if tax_map:
            first_id, (kingdom, species) = list(tax_map.items())[0]
            print(f"Example: {first_id} → kingdom={kingdom}, species={species}")
    else:
        print(f"Test file not found: {test_file}")
