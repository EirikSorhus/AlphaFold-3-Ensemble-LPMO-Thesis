# src/lpmo_pipeline/utils/exceptions.py
"""
Responsibility: Custom exception hierarchy for the LPMO pipeline.
All hard-fail gate violations raise a GateFailure subclass.
Soft flags use flag objects, not exceptions.
"""
from __future__ import annotations


class PipelineError(Exception):
    """Base for all pipeline errors."""

    def __init__(self, message: str, run_id: str = "", step: str = ""):
        self.run_id = run_id
        self.step = step
        super().__init__(message)


# --------------------------------------------------------------------------- #
# Gate failures (hard-fail → pose is dropped)
# --------------------------------------------------------------------------- #
class GateFailure(PipelineError):
    """A mandatory pass/fail gate was not met."""


class AtomMappingCoverageError(GateFailure):
    """atom_mapping_coverage < 100 %."""


class PrivateerCCDError(GateFailure):
    """Glycan residue not a valid CCD monosaccharide."""


class PlacerZeroPosesError(GateFailure):
    """PLACER returned zero poses."""


class PoseBustersCriticalError(GateFailure):
    """Critical PoseBusters failure (steric clash, chem valence)."""


class PrivateerAnomerError(GateFailure):
    """Privateer anomer / ring pucker failure."""


class CuHisDistanceError(GateFailure):
    """Cu–His distance outside 1.9–2.6 Å."""


# --------------------------------------------------------------------------- #
# Non-fatal issues (logged, but processing continues)
# --------------------------------------------------------------------------- #
class SoftFlag(PipelineError):
    """Pose kept but flagged for downstream attention."""


class PoseBustersSoftWarning(SoftFlag):
    """Minor PoseBusters warning (cosmetic, non-structural)."""


class CrystalSimilarityLow(SoftFlag):
    """Crystal IFP similarity below empirical threshold."""


# --------------------------------------------------------------------------- #
# Infrastructure errors
# --------------------------------------------------------------------------- #
class IngestError(PipelineError):
    """mmCIF could not be parsed or required categories missing."""


class ProtonationError(PipelineError):
    """Reduce / OpenBabel protonation step failed."""
