"""Shared Wilcoxon DE gate used by Mixscape and Mixscale."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

SKIPPED_NO_DE = "skipped_no_de_genes"


def differential_gene_gate(
    adata: AnnData,
    *,
    control: str,
    min_de_genes: int,
    logfc_threshold: float,
    pval_cutoff: float,
    max_de_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wilcoxon each target against NT. The table stores counts, not p-values.

    A target passes only when at least ``min_de_genes`` genes clear the absolute
    log fold-change and adjusted p-value cutoffs. ``max_de_genes`` only caps the
    gene list used later for scoring.
    """
    labels = adata.obs["target_gene"].astype(str)
    targets = [str(gene) for gene in pd.unique(labels) if str(gene) != control]
    n_nt = int((labels == control).sum())
    rows: list[dict[str, Any]] = []
    gene_rows: list[dict[str, Any]] = []
    eligible: list[str] = []
    for gene in targets:
        n_cells = int((labels == gene).sum())
        if n_cells < 2 or n_nt < 2:
            rows.append(
                _gate_row(
                    gene,
                    0,
                    min_de_genes,
                    logfc_threshold,
                    pval_cutoff,
                    max_de_genes,
                    passed=False,
                )
            )
            continue
        eligible.append(gene)
    if eligible:
        work = AnnData(
            X=adata.X,
            obs=pd.DataFrame({"target_gene": labels.to_numpy()}, index=adata.obs_names),
            var=pd.DataFrame(index=adata.var_names),
        )
        if "log1p" in adata.uns:
            work.uns["log1p"] = adata.uns["log1p"]
        scored: list[dict[str, Any]] = []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                sc.tl.rank_genes_groups(
                    work,
                    groupby="target_gene",
                    groups=eligible,
                    reference=control,
                    method="wilcoxon",
                    use_raw=False,
                )
            result = work.uns["rank_genes_groups"]
            symbol_by_name = _symbol_by_var_name(adata)
            for gene in eligible:
                names = np.asarray(result["names"][gene]).astype(str)
                signed = np.asarray(result["logfoldchanges"][gene], dtype=float)
                adjusted = np.asarray(result["pvals_adj"][gene], dtype=float)
                raw_p = np.asarray(result["pvals"][gene], dtype=float)
                keep = (np.abs(signed) >= logfc_threshold) & (adjusted < pval_cutoff)
                kept_names = names[keep]
                kept_signed = signed[keep]
                kept_adjusted = adjusted[keep]
                kept_p = raw_p[keep]
                order = np.argsort(kept_p, kind="stable")
                kept_names = kept_names[order]
                kept_signed = kept_signed[order]
                kept_adjusted = kept_adjusted[order]
                n_de = int(kept_names.size)
                passed = n_de >= min_de_genes
                scored.append(
                    _gate_row(
                        gene,
                        n_de,
                        min_de_genes,
                        logfc_threshold,
                        pval_cutoff,
                        max_de_genes,
                        passed=passed,
                    )
                )
                if not passed:
                    continue
                capped = slice(0, max_de_genes)
                for rank, name, log_fc, adj_p in zip(
                    range(int(kept_names[capped].size)),
                    kept_names[capped],
                    kept_signed[capped],
                    kept_adjusted[capped],
                    strict=True,
                ):
                    gene_rows.append(
                        {
                            "target_gene": gene,
                            "gene": str(name),
                            "gene_symbol": symbol_by_name.get(str(name), str(name)),
                            "rank": rank,
                            "logfc": float(log_fc),
                            "adj_p": float(adj_p),
                        }
                    )
        except Exception:
            scored = [
                _gate_row(
                    gene,
                    0,
                    min_de_genes,
                    logfc_threshold,
                    pval_cutoff,
                    max_de_genes,
                    passed=False,
                )
                for gene in eligible
            ]
            gene_rows = []
        rows.extend(scored)
    table = pd.DataFrame(
        rows,
        columns=[
            "target_gene",
            "n_de_genes",
            "min_de_genes",
            "logfc_threshold",
            "pval_cutoff",
            "max_de_genes",
            "fine_mode",
            "passed",
            "reason",
        ],
    )
    genes = pd.DataFrame(
        gene_rows,
        columns=["target_gene", "gene", "gene_symbol", "rank", "logfc", "adj_p"],
    )
    return table, genes


def _symbol_by_var_name(adata: AnnData) -> dict[str, str]:
    """Map matrix names to gene symbols. Guide ids are not used."""
    if "gene_symbol" not in adata.var.columns:
        return {}
    mapping: dict[str, str] = {}
    for name, symbol in zip(
        adata.var_names.astype(str), adata.var["gene_symbol"].astype(str), strict=False
    ):
        mapping.setdefault(name, symbol)
    return mapping


def _gate_row(
    gene: str,
    n_de_genes: int,
    min_de_genes: int,
    logfc_threshold: float,
    pval_cutoff: float,
    max_de_genes: int,
    *,
    passed: bool,
) -> dict[str, Any]:
    return {
        "target_gene": gene,
        "n_de_genes": int(n_de_genes),
        "min_de_genes": int(min_de_genes),
        "logfc_threshold": float(logfc_threshold),
        "pval_cutoff": float(pval_cutoff),
        "max_de_genes": int(max_de_genes),
        "fine_mode": False,
        "passed": bool(passed),
        "reason": "" if passed else SKIPPED_NO_DE,
    }


def _gene_lists(frame: pd.DataFrame) -> dict[str, list[str]]:
    lists: dict[str, list[str]] = {}
    if frame.empty:
        return lists
    for gene, group in frame.groupby("target_gene", sort=False):
        ordered = group.sort_values("rank")
        lists[str(gene)] = ordered["gene"].astype(str).tolist()
    return lists


def stamp_control_label(obs: pd.DataFrame, mask: np.ndarray, column: str, label: str) -> None:
    """Write the control label without requiring it to already be a category."""
    values = obs[column].astype(str).to_numpy(copy=True)
    values[np.asarray(mask)] = label
    obs[column] = values
