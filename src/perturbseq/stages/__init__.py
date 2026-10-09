"""Pipeline stage bodies (1→5). Orchestration and resume live in ``pipeline`` / ``repro``.

Each submodule owns one stage's wiring; domain algorithms stay in sibling packages
(``io``, ``preprocessing``, ``perturbation``, ``statistics``, ``robustness``, ``report``).
"""

from perturbseq.stages.input_validation import stage_input_validation
from perturbseq.stages.perturbation_modeling import stage_perturbation_modeling
from perturbseq.stages.preprocessing_qc import stage_preprocessing_qc
from perturbseq.stages.report import stage_report
from perturbseq.stages.robustness import stage_robustness
from perturbseq.stages.statistical_inference import stage_statistical_inference

__all__ = [
    "stage_input_validation",
    "stage_perturbation_modeling",
    "stage_preprocessing_qc",
    "stage_report",
    "stage_robustness",
    "stage_statistical_inference",
]
