"""Statistical inference: pseudobulk DE, FDR-bearing contrasts, and design checks."""

from __future__ import annotations

import pandas as pd
from anndata import AnnData
import scanpy as sc


def _require_pertpy():
    try:
        import pertpy as pt
    except ImportError as exc:
        raise ImportError(
            "pertpy is required for Mixscape/E-distance/PyDESeq2. "
            "Install with: pip install 'pertpy[de]>=1.3'"
        ) from exc
    return pt


def check_experimental_design(
    adata: AnnData,
    replicate_col: str | None = None,
    groupby: str = "gene_target",
) -> dict:
    """Summarize whether DE can use replicate-aware pseudobulk inference."""
    n_groups = int(adata.obs[groupby].astype(str).nunique()) if groupby in adata.obs else 0
    n_reps = 1
    if replicate_col and replicate_col in adata.obs:
        n_reps = int(adata.obs[replicate_col].nunique())
    replicate_aware = bool(replicate_col and n_reps >= 2)
    note = None
    if not replicate_aware:
        note = (
            "Single-sample / no biological replicates: DE is cell-level Wilcoxon and is exploratory. "
            "Do not treat cell-level p-values as biological-replicate inference."
        )
    return {
        "replicate_col": replicate_col,
        "n_replicates": n_reps,
        "n_groups": n_groups,
        "replicate_aware": replicate_aware,
        "recommended_method": "pydeseq2_pseudobulk" if replicate_aware else "wilcoxon_cell_level_exploratory",
        "note": note,
    }


def exploratory_wilcoxon(adata: AnnData, groupby: str, group: str, reference: str = "NT") -> pd.DataFrame:
    mask = adata.obs[groupby].astype(str).isin([group, reference])
    subset = adata[mask].copy()
    if subset.obs[groupby].nunique() < 2:
        return pd.DataFrame()
    sc.tl.rank_genes_groups(subset, groupby=groupby, groups=[group], reference=reference, method="wilcoxon")
    table = sc.get.rank_genes_groups_df(subset, group=group)
    table.insert(0, "contrast", f"{group}_vs_{reference}")
    table.insert(0, "method", "wilcoxon_cell_level_exploratory")
    table.insert(0, "expression_matrix", "log1p_normalized_X")
    return table


def run_deseq2_or_wilcoxon(
    counts: AnnData,
    group: str,
    reference: str = "NT",
    replicate_col: str | None = None,
    groupby: str = "mixscape_class",
) -> pd.DataFrame:
    """Pseudobulk PyDESeq2 when replicates exist; otherwise exploratory Wilcoxon."""
    mask = counts.obs[groupby].astype(str).isin([group, reference])
    subset = counts[mask].copy()
    n_reps = subset.obs[replicate_col].nunique() if replicate_col and replicate_col in subset.obs else 1
    if replicate_col and n_reps >= 2:
        pt = _require_pertpy()
        if "counts" in subset.layers:
            subset.X = subset.layers["counts"]
        pseudobulk = pt.tl.PseudobulkSpace().compute(
            subset, target_col=groupby, groups_col=replicate_col, mode="sum"
        )
        pds2 = pt.tl.PyDESeq2(adata=pseudobulk, design=f"~ {replicate_col} + {groupby}")
        pds2.fit()
        results = pds2.test_contrasts(
            pds2.contrast(column=groupby, baseline=reference, group_to_compare=group)
        )
        results = results.copy()
        results.insert(0, "contrast", f"{group}_vs_{reference}")
        results.insert(0, "method", "pydeseq2_pseudobulk")
        results.insert(0, "expression_matrix", "layers_counts_pseudobulk_sum")
        return results
    return exploratory_wilcoxon(subset, groupby=groupby, group=group, reference=reference)


__all__ = [
    "check_experimental_design",
    "exploratory_wilcoxon",
    "run_deseq2_or_wilcoxon",
]
