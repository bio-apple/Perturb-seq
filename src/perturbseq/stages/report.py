"""Stage 5 — reproducible report.

Owns: finalizing H5AD / JSON / HTML / provenance via the report module.
Does NOT own: report writers themselves (``report``) or upstream analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from anndata import AnnData

from perturbseq.report import finalize_outputs
from perturbseq.stages._helpers import _mark_stage

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def stage_report(
    rna: AnnData,
    config: PipelineConfig,
    report: dict,
    *,
    input_files: dict[str, Path],
) -> dict:
    """Stage 5: H5AD, JSON/HTML reports, provenance manifest."""
    print("[5/5] Reproducible report", flush=True)
    report = finalize_outputs(
        rna,
        report,
        config.output_dir,
        config.sample_id,
        inputs={k: str(v) for k, v in input_files.items()},
    )
    _mark_stage(report, "5_report")
    return report
