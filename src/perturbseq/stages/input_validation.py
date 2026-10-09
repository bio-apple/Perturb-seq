"""Stage 1 — input & validation.

Owns: load/validate DRAGEN MEX + guide reference + assignments into AnnData.
Does NOT own: guide annotation, QC, or downstream modeling (see later stages / domain modules).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from anndata import AnnData

from perturbseq.io import load_dragen_sample, validate_sample_inputs
from perturbseq.stages._helpers import _mark_stage

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def stage_input_validation(config: PipelineConfig, report: dict) -> tuple[AnnData, dict[str, Path]]:
    """Stage 1: DRAGEN MEX + guide reference + assignments."""
    print(f"[1/5] Input & validation: {config.sample_id} from {config.input_dir}", flush=True)
    files = validate_sample_inputs(config.input_dir, config.sample_id)
    rna, crispr, assignments, feature_ref = load_dragen_sample(config.input_dir, config.sample_id)
    print(f"Loaded {rna.n_obs} cells × {rna.n_vars} genes", flush=True)
    report["n_cells_loaded"] = int(rna.n_obs)
    report["n_genes_loaded"] = int(rna.n_vars)
    report["n_crispr_features"] = int(crispr.n_vars)
    report["input_files"] = {k: str(v) for k, v in files.items()}
    _mark_stage(report, "1_input_validation")
    # Carry assignments/ref/CRISPR via temporary attrs for stage 2 (avoid widening return API).
    rna.uns["_perturbseq_assignments"] = assignments
    rna.uns["_perturbseq_feature_ref"] = feature_ref
    rna.uns["_perturbseq_crispr"] = crispr
    return rna, files
