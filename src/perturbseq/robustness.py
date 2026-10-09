"""Domain: sensitivity / consistency signals and integrated confidence flags.

Owns: deriving confidence flags from guide QC, Mixscape success, DE, design.
Does NOT own: running Mixscape/DE themselves (``perturbation``, ``statistics``) or
writing the HTML/JSON report (``report``). Stage 4 only calls into this module.
"""

from __future__ import annotations

from typing import Any

import pandas as pd


def effect_consistency_flags(
    guide_summary: dict | None = None,
    consistency_df: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Derive confidence flags from guide-level consistency / QC summaries."""
    flags: dict[str, Any] = {
        "guides_mostly_consistent": None,
        "n_inconsistent_genes": 0,
        "n_warnings": 0,
        "confidence": "unknown",
    }
    if guide_summary:
        n_consistent = int(guide_summary.get("n_genes_consistent_guides") or 0)
        n_inconsistent = int(guide_summary.get("n_genes_inconsistent_guides") or 0)
        n_unknown = int(guide_summary.get("n_genes_consistency_unknown") or 0)
        flags["n_inconsistent_genes"] = n_inconsistent
        flags["n_warnings"] = int(guide_summary.get("n_warnings") or 0)
        scored = n_consistent + n_inconsistent
        if scored == 0 and n_unknown == 0:
            flags["confidence"] = "insufficient_guide_data"
        elif n_inconsistent == 0 and n_consistent > 0:
            flags["guides_mostly_consistent"] = True
            flags["confidence"] = "high" if flags["n_warnings"] == 0 else "moderate"
        elif n_inconsistent > 0 and n_inconsistent >= n_consistent:
            flags["guides_mostly_consistent"] = False
            flags["confidence"] = "low"
        else:
            flags["guides_mostly_consistent"] = True
            flags["confidence"] = "moderate"
    elif consistency_df is not None and not consistency_df.empty:
        n_inconsistent = int(consistency_df["guides_consistent"].eq(False).sum())
        flags["n_inconsistent_genes"] = n_inconsistent
        flags["guides_mostly_consistent"] = n_inconsistent == 0
        flags["confidence"] = "high" if n_inconsistent == 0 else "low"
    return flags


def sensitivity_summary(
    *,
    mixscape_ok: bool,
    edistances: pd.DataFrame | None,
    n_de_contrasts: int,
    design: dict | None = None,
    guide_flags: dict | None = None,
) -> dict[str, Any]:
    """Lightweight stability snapshot across modeling and inference tracks."""
    ed = edistances if edistances is not None else pd.DataFrame()
    top_ed = None
    if not ed.empty and "edistance" in ed.columns:
        top_ed = float(ed["edistance"].iloc[0])
    summary = {
        "mixscape_used": mixscape_ok,
        "n_edistance_targets": int(len(ed)),
        "top_edistance": top_ed,
        "n_de_contrasts": int(n_de_contrasts),
        "replicate_aware_de": bool((design or {}).get("replicate_aware")),
        "guide_confidence": (guide_flags or {}).get("confidence", "unknown"),
    }
    # Simple integrated confidence: down-weight when guides disagree or DE is exploratory-only.
    confidence = (guide_flags or {}).get("confidence", "unknown")
    if not mixscape_ok and ed.empty and n_de_contrasts == 0:
        confidence = "insufficient_evidence"
    elif confidence == "high" and not (design or {}).get("replicate_aware") and n_de_contrasts > 0:
        confidence = "moderate"
    summary["integrated_confidence"] = confidence
    return summary


def integrate_evidence(
    report: dict,
    *,
    guide_summary: dict | None = None,
    consistency_df: pd.DataFrame | None = None,
    edistances: pd.DataFrame | None = None,
    mixscape_ok: bool = False,
    n_de_contrasts: int = 0,
    design: dict | None = None,
) -> dict[str, Any]:
    """Attach confidence flags and a sensitivity snapshot onto the pipeline report."""
    guide_flags = effect_consistency_flags(guide_summary=guide_summary, consistency_df=consistency_df)
    sensitivity = sensitivity_summary(
        mixscape_ok=mixscape_ok,
        edistances=edistances,
        n_de_contrasts=n_de_contrasts,
        design=design,
        guide_flags=guide_flags,
    )
    block = {
        "guide_flags": guide_flags,
        "sensitivity": sensitivity,
        "design": design or {},
    }
    report["robustness"] = block
    report["confidence_flags"] = {
        "guide_confidence": guide_flags.get("confidence"),
        "integrated_confidence": sensitivity.get("integrated_confidence"),
        "guides_mostly_consistent": guide_flags.get("guides_mostly_consistent"),
        "replicate_aware_de": sensitivity.get("replicate_aware_de"),
    }
    report["steps"] = list(report.get("steps") or []) + ["robustness"]
    return block


__all__ = [
    "effect_consistency_flags",
    "integrate_evidence",
    "sensitivity_summary",
]
