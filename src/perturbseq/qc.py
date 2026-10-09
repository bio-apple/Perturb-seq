from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.typing import NDArray
from scipy import sparse

# Fraction of a guide's cells removed by QC → strong cytotoxicity signal (not silent discard).
DEFAULT_CYTOTOXICITY_FRAC_REMOVED = 0.9
_MIN_NT_FOR_AWARE_QC = 3


def _to_dense_1d(matrix: Any) -> NDArray[np.floating]:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix).sum(axis=1).ravel()


def add_qc_metrics(adata: AnnData, mt_prefixes: tuple[str, ...] = ("MT-", "mt-")) -> AnnData:
    """Add ``n_counts``, ``n_genes``, ``pct_counts_mt`` to ``adata.obs`` (in place)."""
    names = pd.Index(adata.var_names.astype(str))
    ids = pd.Index(adata.var["gene_ids"].astype(str)) if "gene_ids" in adata.var else names
    mt_mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in mt_prefixes:
        mt_mask |= np.asarray(names.str.startswith(prefix), dtype=bool)
        mt_mask |= np.asarray(ids.str.startswith(prefix), dtype=bool)
    adata.var["mt"] = mt_mask
    counts = _to_dense_1d(adata.X)
    genes = _to_dense_1d(adata.X > 0)
    mito = _to_dense_1d(adata[:, mt_mask].X) if mt_mask.any() else np.zeros(adata.n_obs)
    adata.obs["n_counts"] = counts
    adata.obs["n_genes"] = genes
    adata.obs["pct_counts_mt"] = np.divide(mito, counts, out=np.zeros_like(counts, dtype=float), where=counts > 0) * 100
    return adata


def _mad(values: np.ndarray) -> float:
    med = np.median(values)
    return float(np.median(np.abs(values - med))) * 1.4826


def mad_mask(
    values: pd.Series,
    n_mads: float,
    high: bool = True,
    low: bool = True,
    reference: pd.Series | None = None,
) -> pd.Series:
    """Keep values within median ± n_mads·MAD.

    When ``reference`` is set, median/MAD are estimated from that subset only
    (perturbation-aware QC: fit on NT, apply to all cells).
    """
    array = values.to_numpy(dtype=float)
    ref = reference.to_numpy(dtype=float) if reference is not None else array
    if len(ref) == 0:
        return pd.Series(True, index=values.index)
    med = np.median(ref)
    spread = _mad(ref)
    if spread == 0:
        return pd.Series(True, index=values.index)
    keep = np.ones(len(array), dtype=bool)
    if high:
        keep &= array <= med + n_mads * spread
    if low:
        keep &= array >= med - n_mads * spread
    return pd.Series(keep, index=values.index)


def control_cell_mask(obs: pd.DataFrame, control: str = "NT") -> pd.Series | None:
    """Boolean mask for negative-control cells (``gene_target`` / ``perturbation`` == control)."""
    if "gene_target" in obs.columns:
        return obs["gene_target"].astype(str) == str(control)
    if "perturbation" in obs.columns:
        return obs["perturbation"].astype(str) == str(control)
    return None


def qc_filter_by_category(
    obs: pd.DataFrame,
    keep: pd.Series,
    *,
    control: str = "NT",
    cytotoxicity_frac: float = DEFAULT_CYTOTOXICITY_FRAC_REMOVED,
    categories: tuple[str, ...] = ("guide_id", "gene_target"),
) -> pd.DataFrame:
    """Per-guide / gene_target counts and fraction removed by the cell QC mask."""
    rows: list[dict[str, Any]] = []
    keep_aligned = keep.reindex(obs.index).fillna(False).astype(bool)
    control_tokens = {str(control), "NT", "unassigned", ""}
    for col in categories:
        if col not in obs.columns:
            continue
        labels = obs[col].astype(str)
        for value, idx in labels.groupby(labels, sort=True).groups.items():
            n_before = int(len(idx))
            n_kept = int(keep_aligned.loc[idx].sum())
            n_removed = n_before - n_kept
            frac = float(n_removed / n_before) if n_before else 0.0
            is_ctrl = str(value) in control_tokens
            row: dict[str, Any] = {
                "category": col,
                "value": str(value),
                "n_cells_before_qc": n_before,
                "n_cells_after_qc": n_kept,
                "n_cells_removed_qc": n_removed,
                "frac_removed_qc": frac,
                "cytotoxicity_qc_flag": (not is_ctrl) and frac >= float(cytotoxicity_frac),
            }
            # Carry gene_target on guide_id rows so fully QC-depleted guides can be recorded.
            if col == "guide_id" and "gene_target" in obs.columns:
                targets = obs.loc[idx, "gene_target"].astype(str)
                row["gene_target"] = str(targets.mode().iloc[0]) if len(targets) else ""
            rows.append(row)
    return pd.DataFrame(rows)


def write_qc_filter_by_guide(table: pd.DataFrame, path: Path) -> pd.DataFrame:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table


def filter_cells(
    adata: AnnData,
    n_mads: float = 5.0,
    min_cells: int = 3,
    min_genes: int = 0,
    *,
    perturbation_aware: bool = False,
    control: str = "NT",
    cytotoxicity_frac: float = DEFAULT_CYTOTOXICITY_FRAC_REMOVED,
) -> tuple[AnnData, dict[str, Any]]:
    """MAD-based cell QC. Thresholds come from the observed distributions, not fixed cutoffs.

    When ``perturbation_aware`` is True, fit MAD thresholds on NT / control cells only,
    then apply those thresholds to all cells. Prefer this for Perturb-seq so strong
    phenotypes (essential KO, apoptosis) are not used to set cutoffs that then drop them.
    """
    log: dict[str, Any] = {
        "n_cells_start": int(adata.n_obs),
        "n_genes_start": int(adata.n_vars),
        "perturbation_aware_qc": bool(perturbation_aware),
        "qc_threshold_source": "all_cells",
    }
    reference: pd.Series | None = None
    if perturbation_aware:
        ctrl = control_cell_mask(adata.obs, control=control)
        n_nt = int(ctrl.sum()) if ctrl is not None else 0
        log["n_control_cells_for_qc"] = n_nt
        if ctrl is not None and n_nt >= _MIN_NT_FOR_AWARE_QC:
            reference = adata.obs.loc[ctrl]
            log["qc_threshold_source"] = f"control:{control}"
        else:
            log["qc_threshold_source"] = "all_cells_fallback"
            log["perturbation_aware_fallback"] = (
                f"need ≥{_MIN_NT_FOR_AWARE_QC} control cells with gene_target/perturbation="
                f"{control!r}; using global MAD"
            )

    def _ref(col: str) -> pd.Series | None:
        if reference is None:
            return None
        return reference[col]

    keep = mad_mask(adata.obs["n_counts"], n_mads, reference=_ref("n_counts"))
    keep &= mad_mask(adata.obs["n_genes"], n_mads, reference=_ref("n_genes"))
    if adata.var["mt"].any():
        keep &= mad_mask(adata.obs["pct_counts_mt"], n_mads, low=False, reference=_ref("pct_counts_mt"))
    if min_genes:
        keep &= adata.obs["n_genes"] >= min_genes

    by_guide = qc_filter_by_category(
        adata.obs,
        keep,
        control=control,
        cytotoxicity_frac=cytotoxicity_frac,
    )
    log["qc_filter_by_guide"] = by_guide
    n_cyto = int(by_guide["cytotoxicity_qc_flag"].sum()) if not by_guide.empty else 0
    log["n_cytotoxicity_qc_flags"] = n_cyto

    filtered = adata[keep].copy()
    log["n_cells_after_qc"] = int(filtered.n_obs)
    gene_counts = np.asarray((filtered.X > 0).sum(axis=0)).ravel()
    filtered = filtered[:, gene_counts >= min_cells].copy()
    log["n_genes_after_qc"] = int(filtered.n_vars)
    log["n_cells_removed_qc"] = log["n_cells_start"] - log["n_cells_after_qc"]
    return filtered, log


def sample_qc_summary(adata: AnnData, sample_col: str = "sample_id") -> dict[str, Any]:
    """Per-sample cell/UMI summaries for sample-level QC reporting."""
    if sample_col not in adata.obs:
        return {"n_samples": 1, "n_cells": int(adata.n_obs)}
    rows = []
    for sample, sub in adata.obs.groupby(adata.obs[sample_col].astype(str), sort=True):
        row = {"sample_id": sample, "n_cells": int(len(sub))}
        for col in ("n_counts", "n_genes", "pct_counts_mt"):
            if col in sub.columns:
                row[f"median_{col}"] = float(sub[col].astype(float).median())
        rows.append(row)
    return {"n_samples": len(rows), "samples": rows}
