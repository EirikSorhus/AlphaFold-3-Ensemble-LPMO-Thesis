"""Oligosaccharide definitions and AF3 builder helpers.

An *OligoRegistry* maps oligosaccharide prefixes (e.g. ``NAG``, ``CEL``,
``STA``) to their monomer CCD code and the bond-atom pair that links
successive residues.  The registry is loaded from a user-supplied YAML
file and is consumed exclusively by the AF3 runner to build multi-residue
``ccdCodes`` lists and ``bondedAtomPairs``.

Example YAML::

    oligo_definitions:
      NAG:
        monomer: NAG
        bond_atom_pair: ["C1", "O4"]
      CEL:
        monomer: GLC
        bond_atom_pair: ["C1", "O4"]
      STA:
        monomer: BGC
        bond_atom_pair: ["C1", "O4"]
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

# Regex for parsing oligo CCD codes like "NAG6", "CEL12", "STA3".
# Group 1 = alphabetical prefix, Group 2 = numeric chain length.
_OLIGO_RE = re.compile(r"^([A-Za-z]+)(\d+)$")


@dataclass(frozen=True)
class OligoSpec:
    """Definition of one oligosaccharide type.

    Attributes:
        prefix: Oligosaccharide prefix (e.g. ``"NAG"``).
        monomer: CCD code of the repeating monomer (e.g. ``"NAG"``).
        bond_atom_pair: Tuple *(atom_in_i, atom_in_i_plus_1)* that
            describes the covalent bond between consecutive residues.
    """

    prefix: str
    monomer: str
    bond_atom_pair: tuple[str, str]

    def __post_init__(self) -> None:
        if not self.prefix:
            raise ValueError("OligoSpec.prefix must be non-empty")
        if not self.monomer:
            raise ValueError("OligoSpec.monomer must be non-empty")
        if len(self.bond_atom_pair) != 2:
            raise ValueError(
                f"bond_atom_pair must have exactly 2 elements, got {len(self.bond_atom_pair)}"
            )
        if not self.bond_atom_pair[0] or not self.bond_atom_pair[1]:
            raise ValueError("bond_atom_pair atom names must be non-empty strings")


@dataclass(frozen=True)
class OligoRegistry:
    """Registry of known oligosaccharide definitions.

    An empty registry (no specs) is a valid no-op: ``parse()`` will
    always return ``None`` and existing CIF-based behaviour is
    preserved.
    """

    specs: dict[str, OligoSpec] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def parse(self, ccd_code: str) -> tuple[str, int, OligoSpec] | None:
        """Try to interpret *ccd_code* as an oligosaccharide.

        Returns ``(prefix, n, spec)`` when the code matches a known
        prefix followed by a numeric length ≥ 2, otherwise ``None``.
        """
        m = _OLIGO_RE.match(ccd_code)
        if m is None:
            return None

        prefix = m.group(1).upper()
        n = int(m.group(2))

        spec = self.specs.get(prefix)
        if spec is None:
            return None

        if n < 2:
            logger.warning(
                "Oligo '%s' has chain length %d (<2) — treating as regular ligand",
                ccd_code, n,
            )
            return None

        return prefix, n, spec

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: Path) -> "OligoRegistry":
        """Load an ``OligoRegistry`` from a YAML file.

        Expected format::

            oligo_definitions:
              NAG:
                monomer: NAG
                bond_atom_pair: ["C1", "O4"]
        """
        with open(path) as fh:
            raw = yaml.safe_load(fh)

        if raw is None:
            raise ValueError(f"Empty YAML file: {path}")

        defs: dict[str, Any] = raw.get("oligo_definitions", {})
        if not defs:
            raise ValueError(
                f"No 'oligo_definitions' key (or empty) in {path}"
            )

        specs: dict[str, OligoSpec] = {}
        for prefix_raw, body in defs.items():
            prefix = str(prefix_raw).upper()
            monomer = body.get("monomer")
            bap = body.get("bond_atom_pair")

            if monomer is None:
                raise ValueError(
                    f"Missing 'monomer' for oligo prefix '{prefix}' in {path}"
                )
            if bap is None or len(bap) != 2:
                raise ValueError(
                    f"'bond_atom_pair' for '{prefix}' must be a 2-element list in {path}"
                )

            specs[prefix] = OligoSpec(
                prefix=prefix,
                monomer=str(monomer).upper(),
                bond_atom_pair=(str(bap[0]), str(bap[1])),
            )

        logger.info(
            "Loaded oligo registry with %d definitions: %s",
            len(specs),
            ", ".join(sorted(specs)),
        )
        return cls(specs=specs)

    @classmethod
    def empty(cls) -> "OligoRegistry":
        """Return an empty (no-op) registry."""
        return cls()


# ------------------------------------------------------------------
# AF3 JSON builders for oligosaccharides
# ------------------------------------------------------------------


def build_af3_oligo_ligand_entry(
    chain_id: str,
    spec: OligoSpec,
    n: int,
) -> dict[str, Any]:
    """Build a single AF3 ``ligand`` sequence entry for an oligosaccharide.

    Returns::

        {"ligand": {"id": "<chain_id>", "ccdCodes": ["GLC", "GLC", ...]}}
    """
    return {
        "ligand": {
            "id": chain_id,
            "ccdCodes": [spec.monomer] * n,
        }
    }


def build_af3_bonded_atom_pairs(
    chain_id: str,
    spec: OligoSpec,
    n: int,
) -> list[list[list[Any]]]:
    """Generate ``bondedAtomPairs`` for *n* residues linked end-to-end.

    Each bond connects residue *i* → *i+1* via the spec's
    ``bond_atom_pair``.  Residue indices are **1-based** as required by
    AF3.

    Returns a list of bond entries, each of the form::

        [[chain_id, res_i, atom_a], [chain_id, res_i+1, atom_b]]
    """
    atom_a, atom_b = spec.bond_atom_pair
    pairs: list[list[list[Any]]] = []
    for i in range(1, n):  # i = 1 … n-1
        pairs.append([
            [chain_id, i, atom_a],
            [chain_id, i + 1, atom_b],
        ])
    return pairs
