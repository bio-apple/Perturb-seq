from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.typing import NDArray
from scipy import sparse


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


def mad_mask(values: pd.Series, n_mads: float, high: bool = True, low: bool = True) -> pd.Series:
    array = values.to_numpy(dtype=float)
    med = np.median(array)
    spread = _mad(array)
    if spread == 0:
        return pd.Series(True, index=values.index)
    keep = np.ones(len(array), dtype=bool)
    if high:
        keep &= array <= med + n_mads * spread
    if low:
        keep &= array >= med - n_mads * spread
    return pd.Series(keep, index=values.index)


def filter_cells(
    adata: AnnData,
    n_mads: float = 5.0,
    min_cells: int = 3,
    min_genes: int = 0,
) -> tuple[AnnData, dict[str, int]]:
    """MAD-based cell QC. Thresholds come from the observed distributions, not fixed cutoffs."""
    log = {"n_cells_start": int(adata.n_obs), "n_genes_start": int(adata.n_vars)}
    keep = mad_mask(adata.obs["n_counts"], n_mads)
    keep &= mad_mask(adata.obs["n_genes"], n_mads)
    if adata.var["mt"].any():
        keep &= mad_mask(adata.obs["pct_counts_mt"], n_mads, low=False)
    if min_genes:
        keep &= adata.obs["n_genes"] >= min_genes
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
