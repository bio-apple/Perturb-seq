"""Statistical inference: pseudobulk DE, FDR-bearing contrasts, and design checks."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pandas as pd
import scanpy as sc
from anndata import AnnData

from perturbseq._deps import require_pertpy
from perturbseq.parallel import parallel_map


def check_experimental_design(
    adata: AnnData,
    replicate_col: str | None = None,
    groupby: str = "gene_target",
) -> dict[str, Any]:
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
        pt = require_pertpy()
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


def _de_one_contrast(
    args: tuple[AnnData, str, str, str | None, str],
) -> tuple[str, pd.DataFrame | None, str | None]:
    """Worker: one group vs reference on a pre-subset AnnData (picklable top-level)."""
    subset, group, reference, replicate_col, groupby = args
    try:
        table = run_deseq2_or_wilcoxon(
            subset,
            group=group,
            reference=reference,
            replicate_col=replicate_col,
            groupby=groupby,
        )
        if table is None or table.empty:
            return str(group), None, None
        return str(group), table, None
    except Exception as exc:  # noqa: BLE001
        return str(group), None, str(exc)


def run_de_contrasts(
    counts: AnnData,
    groups: Sequence[str],
    *,
    reference: str = "NT",
    replicate_col: str | None = None,
    groupby: str = "de_group",
    min_cells: int = 10,
    n_jobs: int = 1,
) -> tuple[list[tuple[str, pd.DataFrame]], dict[str, str]]:
    """Run DE for each group vs reference; parallelizes the outer loop when n_jobs != 1.

    Each worker receives a copy of cells in {group, reference} only. CSV/plots stay
    on the caller so AnnData is never mutated across processes.
    """
    tasks: list[tuple[AnnData, str, str, str | None, str]] = []
    labels = counts.obs[groupby].astype(str)
    for group in groups:
        g = str(group)
        if int((labels == g).sum()) < min_cells:
            continue
        mask = labels.isin([g, str(reference)])
        tasks.append((counts[mask].copy(), g, str(reference), replicate_col, groupby))

    results = parallel_map(_de_one_contrast, tasks, n_jobs=n_jobs)
    tables: list[tuple[str, pd.DataFrame]] = []
    errors: dict[str, str] = {}
    for group, table, err in results:
        if err:
            errors[group] = err
        elif table is not None:
            tables.append((group, table))
    return tables, errors


__all__ = [
    "check_experimental_design",
    "exploratory_wilcoxon",
    "run_de_contrasts",
    "run_deseq2_or_wilcoxon",
]
