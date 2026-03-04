# src/lpmo_pipeline/tuning/parameter_grid.py
"""
Responsibility: Load and expand parameter grids from tuning YAML configs.
Input:  configs/tuning_{af3,rf3,boltz2}.yaml
Output: List of parameter point dicts, each representing one run config
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)


@dataclass
class ParamPoint:
    """A single parameter configuration to evaluate."""

    model: str
    phase: str  # "refinement" | "diversity"
    params: dict[str, Any] = field(default_factory=dict)
    held_fixed: dict[str, Any] = field(default_factory=dict)
    param_hash: str = ""  # Filled by hashing


def load_param_grid(config_path: Path) -> list[ParamPoint]:
    """Load and expand a tuning YAML config into ParamPoint list.

    Expected YAML structure:
    ```yaml
    model: AF3
    phases:
      refinement:
        vary:
          num_recycles: [10, 15, 20]
        hold:
          num_diffusion_samples: 5
          num_seeds: 10
      diversity:
        vary:
          num_seeds: [20, 50, 100]
        hold:
          num_recycles: "{best_from_refinement}"
          num_diffusion_samples: 5
    ```

    Args:
        config_path: Path to tuning YAML.

    Returns:
        List of ParamPoint (one per grid point, across all phases).
    """
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    model = cfg["model"]
    phases = cfg.get("phases", {})
    points: list[ParamPoint] = []

    for phase_name, phase_cfg in phases.items():
        vary = phase_cfg.get("vary", {})
        hold = phase_cfg.get("hold", {})

        # Expand grid: cartesian product of all "vary" keys
        vary_keys = sorted(vary.keys())
        vary_values = [vary[k] for k in vary_keys]

        for combo in product(*vary_values):
            params = dict(zip(vary_keys, combo))
            point = ParamPoint(
                model=model,
                phase=phase_name,
                params=params,
                held_fixed=hold.copy(),
            )
            points.append(point)

    logger.info(
        "Loaded param grid from %s: %s model, %d phases, %d total points",
        config_path, model, len(phases), len(points),
    )
    return points


def resolve_dependencies(
    points: list[ParamPoint],
    best_params: dict[str, Any],
) -> list[ParamPoint]:
    """Replace '{best_from_*}' placeholders with actual values.

    Args:
        points: List of ParamPoints (some may have placeholder values).
        best_params: {param_name: best_value} from earlier phases.

    Returns:
        Updated list with placeholders resolved.
    """
    for point in points:
        for key, val in point.held_fixed.items():
            if isinstance(val, str) and val.startswith("{best_from_"):
                # Extract the referenced parameter name
                ref_key = key  # The held param name IS the param to look up
                if ref_key in best_params:
                    point.held_fixed[key] = best_params[ref_key]
                    logger.debug("Resolved %s=%s for %s", key, best_params[ref_key], point.phase)
                else:
                    logger.warning("Cannot resolve %s: '%s' not in best_params", key, val)

    return points
