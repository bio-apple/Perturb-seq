"""Optional tertiary guide re-assignment vs DRAGEN GMM calls.

DRAGEN secondary already GMM-calls guides into ``positive_cell_guide_assignments.csv``.
This module optionally **compares** (and optionally overrides) those calls using:

- ``max`` — per-cell argmax CRISPR UMI (simple ambient-prone baseline)
- ``gmm`` — lightweight per-guide 2-component GMM on log1p counts (sklearn)

Does **not** install CatchR / Cell Ranger Feature Barcode / pertpy GuideAssignment;
those remain documented external alternatives.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.guides import DEFAULT_CONTROL_PATTERNS, build_feature_target_map, parse_gene_target

GUIDE_REASSIGN_MODES = ("off", "compare", "apply_max", "apply_gmm")


def _counts_matrix(crispr: AnnData) -> np.ndarray:
    x = crispr.X
    if sparse.issparse(x):
        return np.asarray(x.toarray(), dtype=float)
    return np.asarray(x, dtype=float)


def assign_max_count(crispr: AnnData, *, min_umi: float = 1.0) -> pd.Series:
    """Per-cell guide = argmax CRISPR UMI; empty string if max UMI < ``min_umi``."""
    if crispr.n_vars == 0 or crispr.n_obs == 0:
        return pd.Series("", index=crispr.obs_names, dtype=object)
    mat = _counts_matrix(crispr)
    names = crispr.var_names.astype(str).to_numpy()
    max_idx = np.argmax(mat, axis=1)
    max_vals = mat[np.arange(mat.shape[0]), max_idx]
    out = np.where(max_vals >= float(min_umi), names[max_idx], "")
    return pd.Series(out, index=crispr.obs_names, dtype=object, name="guide_max")


def assign_gmm_per_guide(
    crispr: AnnData,
    *,
    min_umi: float = 1.0,
    random_state: int = 0,
) -> pd.Series:
    """Assign guides via per-feature 2-component GMM on log1p counts.

    For each guide, fit a 2-GMM; cells in the higher-mean component are candidates.
    Each cell keeps the candidate with the highest UMI (ties → first). Falls back to
    ``assign_max_count`` if sklearn is unavailable or a guide has no variance.
    """
    if crispr.n_vars == 0 or crispr.n_obs == 0:
        return pd.Series("", index=crispr.obs_names, dtype=object)
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        warnings.warn(
            "sklearn not available; guide GMM reassignment falling back to max-count",
            stacklevel=2,
        )
        return assign_max_count(crispr, min_umi=min_umi).rename("guide_gmm")

    mat = _counts_matrix(crispr)
    names = crispr.var_names.astype(str).to_numpy()
    n_cells, n_guides = mat.shape
    on = np.zeros((n_cells, n_guides), dtype=bool)
    for j in range(n_guides):
        counts = mat[:, j]
        if float(counts.max()) < float(min_umi):
            continue
        x = np.log1p(counts).reshape(-1, 1)
        if np.unique(x).size < 2:
            on[:, j] = counts >= float(min_umi)
            continue
        try:
            gm = GaussianMixture(n_components=2, random_state=int(random_state), covariance_type="full")
            gm.fit(x)
            high = int(np.argmax(gm.means_.ravel()))
            on[:, j] = gm.predict(x) == high
        except Exception:  # noqa: BLE001
            on[:, j] = counts >= float(min_umi)

    # Prefer GMM-positive guides; if none, fall back to max-count cell-wise.
    assigned = np.full(n_cells, "", dtype=object)
    for i in range(n_cells):
        cand = np.flatnonzero(on[i])
        if cand.size == 0:
            j = int(np.argmax(mat[i]))
            if mat[i, j] >= float(min_umi):
                assigned[i] = names[j]
            continue
        # Among positives, pick highest UMI.
        j = int(cand[np.argmax(mat[i, cand])])
        if mat[i, j] >= float(min_umi):
            assigned[i] = names[j]
    return pd.Series(assigned, index=crispr.obs_names, dtype=object, name="guide_gmm")


def compare_guide_assignments(
    dragen: pd.Series,
    alternatives: dict[str, pd.Series],
) -> pd.DataFrame:
    """Per-cell comparison of DRAGEN ``guide_id`` vs alternative method(s)."""
    idx = dragen.index
    frame = pd.DataFrame({"guide_dragen": dragen.astype(str).reindex(idx).fillna("")}, index=idx)
    for name, series in alternatives.items():
        col = f"guide_{name}"
        frame[col] = series.astype(str).reindex(idx).fillna("")
        frame[f"agree_dragen_{name}"] = frame["guide_dragen"] == frame[col]
    return frame


def summarize_guide_reassignment(comparison: pd.DataFrame) -> dict[str, Any]:
    """Agreement rates and simple contingency tallies for report.json."""
    n = int(len(comparison))
    summary: dict[str, Any] = {"n_cells": n, "methods": {}}
    for col in comparison.columns:
        if not col.startswith("agree_dragen_"):
            continue
        method = col.removeprefix("agree_dragen_")
        agree = comparison[col].fillna(False).astype(bool)
        alt_col = f"guide_{method}"
        n_alt = int((comparison[alt_col].astype(str) != "").sum()) if alt_col in comparison.columns else 0
        summary["methods"][method] = {
            "n_agree_with_dragen": int(agree.sum()),
            "frac_agree_with_dragen": float(agree.mean()) if n else 0.0,
            "n_assigned": n_alt,
        }
    return summary


def apply_guide_ids(
    adata: AnnData,
    guide_ids: pd.Series,
    feature_ref: pd.DataFrame | None = None,
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS,
) -> AnnData:
    """Overwrite ``guide_id`` / ``gene_target`` / ``perturbation`` from new guide IDs.

    Keeps DRAGEN originals in ``guide_id_dragen`` / ``gene_target_dragen`` when absent.
    """
    if "guide_id_dragen" not in adata.obs.columns and "guide_id" in adata.obs.columns:
        adata.obs["guide_id_dragen"] = adata.obs["guide_id"].astype(str)
    if "gene_target_dragen" not in adata.obs.columns and "gene_target" in adata.obs.columns:
        adata.obs["gene_target_dragen"] = adata.obs["gene_target"].astype(str)

    mapping: dict[str, str] = {}
    if feature_ref is not None and not feature_ref.empty:
        mapping = build_feature_target_map(feature_ref, control_patterns)

    aligned = guide_ids.astype(str).reindex(adata.obs_names).fillna("")
    adata.obs["guide_id"] = aligned.to_numpy()
    targets: list[str] = []
    for gid in aligned:
        if not gid:
            targets.append("unassigned")
            continue
        target = mapping.get(gid, parse_gene_target(gid, control_patterns))
        targets.append(target)
    adata.obs["gene_target"] = targets
    adata.obs["perturbation"] = np.where(
        adata.obs["gene_target"].astype(str).eq("unassigned"),
        "unassigned",
        np.where(adata.obs["gene_target"].astype(str).eq("NT"), "NT", "perturbed"),
    )
    # Singlet-like: one assigned guide → num_features=1; empty → 0.
    adata.obs["num_features"] = np.where(aligned.to_numpy() == "", 0, 1).astype(int)
    return adata


def plot_guide_reassignment(comparison: pd.DataFrame, figures_dir: Path) -> Path | None:
    """Bar chart of agreement fractions vs DRAGEN."""
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    methods = [c.removeprefix("agree_dragen_") for c in comparison.columns if c.startswith("agree_dragen_")]
    if not methods:
        return None
    import matplotlib.pyplot as plt

    fracs = [float(comparison[f"agree_dragen_{m}"].fillna(False).astype(bool).mean()) for m in methods]
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    ax.bar(methods, fracs, color="#4C72B0")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction agreeing with DRAGEN")
    ax.set_title("Guide reassignment vs DRAGEN")
    for i, v in enumerate(fracs):
        ax.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    out = figures_dir / "guide_reassignment_agreement.png"
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def run_guide_reassignment(
    rna: AnnData,
    crispr: AnnData,
    *,
    mode: str = "compare",
    feature_ref: pd.DataFrame | None = None,
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS,
    min_umi: float = 1.0,
    random_state: int = 0,
    tables_dir: Path | None = None,
    figures_dir: Path | None = None,
) -> dict[str, Any]:
    """Compare DRAGEN vs max/GMM; optionally apply max or GMM labels to ``rna``.

    ``mode``: ``off`` | ``compare`` | ``apply_max`` | ``apply_gmm``.
    """
    mode = str(mode or "off").strip().lower()
    if mode not in GUIDE_REASSIGN_MODES:
        raise ValueError(f"guide_reassign must be one of {GUIDE_REASSIGN_MODES}, got {mode!r}")
    if mode == "off":
        return {"skipped": True, "reason": "user_off", "mode": mode}

    # Align CRISPR to RNA cells (same barcode index).
    shared = rna.obs_names.intersection(crispr.obs_names)
    if len(shared) == 0 or crispr.n_vars == 0:
        return {
            "skipped": True,
            "reason": "no_crispr_overlap",
            "mode": mode,
            "detail": "No overlapping cell barcodes with CRISPR features for reassignment.",
        }
    crispr_sub = crispr[shared].copy()
    dragen = (
        rna.obs.loc[shared, "guide_id"].astype(str)
        if "guide_id" in rna.obs.columns
        else pd.Series("", index=shared, dtype=object)
    )
    max_ids = assign_max_count(crispr_sub, min_umi=min_umi)
    gmm_ids = assign_gmm_per_guide(crispr_sub, min_umi=min_umi, random_state=random_state)
    comparison = compare_guide_assignments(dragen, {"max": max_ids, "gmm": gmm_ids})
    summary = summarize_guide_reassignment(comparison)
    summary.update({"skipped": False, "reason": None, "mode": mode, "applied": False})

    if tables_dir is not None:
        tables_dir = Path(tables_dir)
        tables_dir.mkdir(parents=True, exist_ok=True)
        comparison.to_csv(tables_dir / "guide_reassignment_comparison.csv")
        summary_path = tables_dir / "guide_reassignment_summary.json"
        import json

        summary_path.write_text(json.dumps(summary, indent=2, default=str))
        summary["tables"] = str(tables_dir / "guide_reassignment_comparison.csv")

    if figures_dir is not None:
        plot_path = plot_guide_reassignment(comparison, figures_dir)
        if plot_path is not None:
            summary["figure"] = str(plot_path)

    if mode == "apply_max":
        apply_guide_ids(rna, max_ids.reindex(rna.obs_names).fillna(""), feature_ref, control_patterns)
        summary["applied"] = True
        summary["applied_method"] = "max"
    elif mode == "apply_gmm":
        apply_guide_ids(rna, gmm_ids.reindex(rna.obs_names).fillna(""), feature_ref, control_patterns)
        summary["applied"] = True
        summary["applied_method"] = "gmm"

    return summary


__all__ = [
    "GUIDE_REASSIGN_MODES",
    "apply_guide_ids",
    "assign_gmm_per_guide",
    "assign_max_count",
    "compare_guide_assignments",
    "plot_guide_reassignment",
    "run_guide_reassignment",
    "summarize_guide_reassignment",
]
