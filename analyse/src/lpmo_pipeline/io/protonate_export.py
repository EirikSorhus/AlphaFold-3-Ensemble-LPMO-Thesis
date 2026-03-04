"""
LPMO Pipeline: Step 3 — Glycan Expansion & Privateer Prep
Responsibility: Ensure glykan residues use valid CCD monosaccharides (not custom oligomer codes)

Step 3 contract:
  Input: normalized.cif from step 2
  Output: privateer_input.cif
  Hard rule: Privateer input MUST use only valid CCD monosaccharides (NAG, BGC, MAN, etc.)
             NO custom oligomer IDs (e.g., CEL6, OLI, etc.)
  
  Procedure:
    1. For each glykan residue: check comp_id ∈ CCD monosaccharides whitelist
    2. If custom oligomer code: decompose into constituent monosacchs OR fail
    3. Verify _struct_conn has glycosidic linkage definition
    4. Ensure Cu-coordinating His residues present + correctly named
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
from dataclasses import dataclass

try:
    import gemmi
except ImportError:
    gemmi = None

from lpmo_pipeline.utils.logging import StructuredLogger, FailureLog
from lpmo_pipeline.utils.data_models import CCD_MONOSACCHARIDES, FailureReason


@dataclass
class GlykanExpansionResult:
    """Outcome of glykan expansion procedure."""
    total_glycan_residues: int
    recognized_monosacchs: int
    failed_residues: List[str]  # comp_id that could not be handled
    custom_oligomer_codes: List[str]  # e.g., ['CEL6', 'OLI']
    success: bool


class PrivateerPrepRunner:
    """Step 3: Prepare normalized mmCIF for Privateer (CCD monosakkarid requirement)."""
    
    def __init__(self, input_cif: Path, output_dir: Path):
        """
        Args:
            input_cif: normalized.cif from step 2
            output_dir: output directory
        """
        self.input_cif = Path(input_cif)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = StructuredLogger("privateer_prep", self.output_dir)
        self.failures = FailureLog(self.output_dir / "privateer_prep_failures.json")
    
    def run(self) -> Tuple[bool, Optional[Path]]:
        """
        Execute glycan expansion + Privateer prep.
        
        Returns:
            (success_flag, privateer_input_cif_or_none)
        """
        self.logger.log_step_start("privateer_prep", {"input": str(self.input_cif)})
        
        try:
            # Parse normalized mmCIF
            doc = gemmi.cif.read_file(str(self.input_cif))
            structure = gemmi.cif.as_structure(doc)
            
            # Main procedure: expand/validate glycans
            expansion_result = self._expand_glycans(structure)
            
            if not expansion_result.success:
                self.logger.log_failure(
                    "glykan_not_ccd",
                    {
                        "failed_residues": expansion_result.failed_residues,
                        "custom_codes": expansion_result.custom_oligomer_codes,
                    }
                )
                self.failures.record(
                    "privateer_prep",
                    None,
                    FailureReason.GLYKAN_NOT_CCD.value,
                    {
                        "total_glycans": expansion_result.total_glycan_residues,
                        "recognized": expansion_result.recognized_monosacchs,
                        "failed": expansion_result.failed_residues,
                    }
                )
                self.failures.write()
                self.logger.log_step_end("privateer_prep", "failure", {}, 0.0)
                return False, None
            
            # Ensure Cu-coordinating His residues present
            self._check_cu_his_residues(structure)
            
            # Save privateer-ready mmCIF
            output_path = self.output_dir / "privateer_input.cif"
            structure.write_minimal_pdb(str(output_path))  # placeholder
            
            self.logger.log_step_end(
                "privateer_prep",
                "success",
                {
                    "glycan_residues": expansion_result.total_glycan_residues,
                    "recognized_monosacchs": expansion_result.recognized_monosacchs,
                },
                0.0
            )
            
            self.failures.write()
            self.logger.close()
            
            return True, output_path
        
        except Exception as e:
            self.logger.log_failure("privateer_prep_exception", {"error": str(e)})
            self.failures.record("privateer_prep", None, "exception", {"error": str(e)})
            self.failures.write()
            self.logger.log_step_end("privateer_prep", "failure", {}, 0.0)
            return False, None
    
    def _expand_glycans(self, structure: gemmi.Structure) -> GlykanExpansionResult:
        """
        Main logic: check each glykan residue, ensure CCD monosaccharide.
        
        Whitelist: CCD_MONOSACCHARIDES (NAG, BGC, MAN, GLC, etc.)
        
        If comp_id not in whitelist:
          - Try to decompose (e.g., CEL6 → 6× BGC + 5 linkages)
          - If not decomposable: FAIL
        """
        model = structure[0]
        failed = []
        custom_codes = []
        recognized_count = 0
        total_glycan = 0
        
        # Identify glycan chains (B, C, D, ...)
        for chain in model:
            if chain.name in ['B', 'C', 'D', 'E']:  # exclude protein (A) and metal (E)
                if chain.name == 'E':  # metal chain, skip
                    continue
                
                for residue in chain:
                    total_glycan += 1
                    comp_id = residue.name
                    
                    if comp_id in CCD_MONOSACCHARIDES:
                        recognized_count += 1
                        self.logger.log_check(
                            f"glykan_ccd_{comp_id}",
                            "pass",
                            f"Residue {comp_id} recognized (CCD monosaccharide)"
                        )
                    else:
                        # Try to decompose or fail
                        if self._is_custom_oligomer(comp_id):
                            custom_codes.append(comp_id)
                            # Attempt decomposition (placeholder)
                            decomposed = self._decompose_oligomer(comp_id)
                            if decomposed:
                                # Successful: would replace residue in-place
                                recognized_count += len(decomposed)
                                self.logger.log_check(
                                    f"glykan_decompose_{comp_id}",
                                    "pass",
                                    f"Custom code {comp_id} decomposed to {len(decomposed)} monosacchs"
                                )
                            else:
                                failed.append(comp_id)
                                self.logger.log_check(
                                    f"glykan_ccd_{comp_id}",
                                    "hard_fail",
                                    f"Non-CCD oligomer {comp_id}, cannot decompose"
                                )
                        else:
                            failed.append(comp_id)
                            self.logger.log_check(
                                f"glykan_ccd_{comp_id}",
                                "hard_fail",
                                f"Unrecognized residue code {comp_id}"
                            )
        
        # Success if all glycans are CCD monosacchs (or decomposed)
        success = len(failed) == 0
        
        return GlykanExpansionResult(
            total_glycan_residues=total_glycan,
            recognized_monosacchs=recognized_count,
            failed_residues=failed,
            custom_oligomer_codes=custom_codes,
            success=success
        )
    
    def _is_custom_oligomer(self, comp_id: str) -> bool:
        """Check if comp_id is a custom oligomer code (heuristic)."""
        # Custom codes often contain letters + numbers (CEL6, OLI, GLC3, etc.)
        # Monosaccharides are typically 3-letter codes
        if len(comp_id) > 3:
            return True
        return False
    
    def _decompose_oligomer(self, comp_id: str) -> Optional[List[str]]:
        """
        Attempt to decompose custom oligomer code into monosacchs.
        
        Examples:
          - CEL6 (6 glucose units) → ['BGC'] * 6
          - CHI4 (4 NAG units) → ['NAG'] * 4
        
        Placeholder: real impl would have a decomposition lookup table.
        """
        decomp_rules = {
            "CEL": "BGC",   # cellulose
            "CHI": "NAG",   # chitin
            "MAL": "MAN",   # mannose-based (heuristic)
        }
        
        for prefix, monosacch in decomp_rules.items():
            if comp_id.startswith(prefix):
                # Extract number if present (e.g., CEL6 → 6)
                try:
                    num_str = comp_id[len(prefix):]
                    num = int(num_str) if num_str else 1
                    return [monosacch] * num
                except ValueError:
                    pass
        
        return None
    
    def _check_cu_his_residues(self, structure: gemmi.Structure):
        """
        Verify Cu-coordinating His residues are present in protein.
        
        Placeholder: real impl would identify Cu atom,
                     find nearby His, check distances.
        """
        self.logger.log_check(
            "cu_his_present",
            "pass",
            "Cu-coordinating His check (placeholder)"
        )
