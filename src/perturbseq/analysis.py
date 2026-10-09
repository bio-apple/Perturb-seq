"""Compat facade — re-exports from ``perturbation`` + ``statistics``.

Owns: nothing (no logic). Prefer importing domain modules directly.
Does NOT own: Mixscape/E-distance (``perturbation``) or DE/design (``statistics``).
"""

from perturbseq.perturbation import (
    build_perturbation_space,
    cluster_perturbations,
    run_edistance,
    run_mixscape,
)
from perturbseq.statistics import exploratory_wilcoxon, run_deseq2_or_wilcoxon

__all__ = [
    "build_perturbation_space",
    "cluster_perturbations",
    "exploratory_wilcoxon",
    "run_deseq2_or_wilcoxon",
    "run_edistance",
    "run_mixscape",
]
