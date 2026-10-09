"""Stage 4 — robustness & evidence integration.

Owns: calling evidence integration and marking the stage complete.
Does NOT own: confidence-flag logic (``robustness`` domain module).
"""

from __future__ import annotations

import pandas as pd

from perturbseq.robustness import integrate_evidence
from perturbseq.stages._helpers import _mark_stage


def stage_robustness(
    report: dict,
    *,
    guide_summary: dict,
    consistency_df: pd.DataFrame,
    edistances: pd.DataFrame,
    mixscape_ok: bool,
    n_de_contrasts: int,
    design: dict,
) -> None:
    """Stage 4: sensitivity / effect consistency / confidence flags."""
    print("[4/5] Robustness & evidence integration", flush=True)
    integrate_evidence(
        report,
        guide_summary=guide_summary,
        consistency_df=consistency_df,
        edistances=edistances,
        mixscape_ok=mixscape_ok,
        n_de_contrasts=n_de_contrasts,
        design=design,
    )
    _mark_stage(report, "4_robustness")
