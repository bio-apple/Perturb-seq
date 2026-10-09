"""Domain: experimental-design checks and pseudobulk / DE contrasts (FDR-bearing).

Owns: DE contrast runners, design awareness, covariate helpers.
Does NOT own: Mixscape/E-distance (``perturbation``), evidence flags (``robustness``),
or stage wiring (``stages.statistical_inference``). Not a general “analysis” dump —
use ``perturbation`` for modeling.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from perturbseq._deps import require_pertpy
from perturbseq.parallel import parallel_map

# Machine-readable DE claim strength for report.json / downstream tools.
EvidenceLevel = Literal["exploratory", "inferential"]
EVIDENCE_LEVEL_EXPLORATORY: EvidenceLevel = "exploratory"
EVIDENCE_LEVEL_INFERENTIAL: EvidenceLevel = "inferential"

# Preferred covariate obs columns when ``de_covariates`` is enabled (available-only).
DEFAULT_DE_COVARIATE_CANDIDATES: tuple[str, ...] = (
    "phase",
    "pct_counts_mt",
    "log_n_counts",
)


def de_evidence_level(*, replicate_aware: bool) -> EvidenceLevel:
    """Map design to evidence_level: Wilcoxon/no replicates → exploratory; else inferential."""
    return EVIDENCE_LEVEL_INFERENTIAL if replicate_aware else EVIDENCE_LEVEL_EXPLORATORY


def ensure_de_covariate_columns(adata: AnnData) -> list[str]:
    """Create ``log_n_counts`` from ``n_counts`` when missing; return available candidates."""
    if "log_n_counts" not in adata.obs.columns and "n_counts" in adata.obs.columns:
        n = pd.to_numeric(adata.obs["n_counts"], errors="coerce").to_numpy(dtype=float)
        adata.obs["log_n_counts"] = np.log1p(np.nan_to_num(n, nan=0.0))
    return [c for c in DEFAULT_DE_COVARIATE_CANDIDATES if c in adata.obs.columns]


def resolve_de_covariates(
    adata: AnnData,
    covariates: Sequence[str] | bool | None,
) -> list[str]:
    """Resolve covariate column names.

    - ``None`` / ``False`` → no covariates
    - ``True`` → default candidates present on ``adata.obs``
    - sequence → intersection with available columns (after ensuring ``log_n_counts``)
    """
    ensure_de_covariate_columns(adata)
    if covariates is None or covariates is False:
        return []
    if covariates is True:
        return ensure_de_covariate_columns(adata)
    wanted = [str(c).strip() for c in covariates if str(c).strip()]
    return [c for c in wanted if c in adata.obs.columns]


def build_deseq2_design(
    groupby: str,
    replicate_col: str | None = None,
    covariates: Sequence[str] | None = None,
) -> str:
    """Build a PyDESeq2-style design formula ``~ terms + groupby``.

    Term order: optional replicate, then covariates, then the contrast factor.
    Duplicate names are dropped while preserving order.
    """
    terms: list[str] = []
    if replicate_col:
        terms.append(str(replicate_col))
    for cov in covariates or ():
        c = str(cov).strip()
        if c and c not in terms and c != groupby:
            terms.append(c)
    terms.append(str(groupby))
    # De-dupe while keeping order
    seen: set[str] = set()
    ordered: list[str] = []
    for t in terms:
        if t not in seen:
            seen.add(t)
            ordered.append(t)
    return "~ " + " + ".join(ordered)


def check_experimental_design(
    adata: AnnData,
    replicate_col: str | None = None,
    groupby: str = "gene_target",
    covariates: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Summarize whether DE can use replicate-aware pseudobulk inference."""
    n_groups = int(adata.obs[groupby].astype(str).nunique()) if groupby in adata.obs else 0
    n_reps = 1
    if replicate_col and replicate_col in adata.obs:
        n_reps = int(adata.obs[replicate_col].nunique())
    replicate_aware = bool(replicate_col and n_reps >= 2)
    covs = list(covariates or [])
    note = None
    if not replicate_aware:
        note = (
            "Single-sample / no biological replicates: DE is cell-level Wilcoxon and is exploratory. "
            "Do not treat cell-level p-values as biological-replicate inference."
        )
        if covs:
            note += (
                f" Requested covariates ({', '.join(covs)}) are recorded for provenance but "
                "are not modeled by Wilcoxon (scanpy rank_genes_groups has no covariate design)."
            )
    design_formula = None
    if replicate_aware:
        design_formula = build_deseq2_design(groupby, replicate_col=replicate_col, covariates=covs)
    level = de_evidence_level(replicate_aware=replicate_aware)
    return {
        "replicate_col": replicate_col,
        "n_replicates": n_reps,
        "n_groups": n_groups,
        "replicate_aware": replicate_aware,
        "covariates": covs,
        "design_formula": design_formula,
        "recommended_method": "pydeseq2_pseudobulk" if replicate_aware else "wilcoxon_cell_level_exploratory",
        "evidence_level": level,
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
    table.insert(0, "evidence_level", EVIDENCE_LEVEL_EXPLORATORY)
    table.insert(0, "expression_matrix", "log1p_normalized_X")
    return table


def run_deseq2_or_wilcoxon(
    counts: AnnData,
    group: str,
    reference: str = "NT",
    replicate_col: str | None = None,
    groupby: str = "mixscape_class",
    covariates: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Pseudobulk PyDESeq2 when replicates exist; otherwise exploratory Wilcoxon."""
    mask = counts.obs[groupby].astype(str).isin([group, reference])
    subset = counts[mask].copy()
    covs = resolve_de_covariates(subset, list(covariates) if covariates else None)
    n_reps = subset.obs[replicate_col].nunique() if replicate_col and replicate_col in subset.obs else 1
    if replicate_col and n_reps >= 2:
        pt = require_pertpy()
        if "counts" in subset.layers:
            subset.X = subset.layers["counts"]
        pseudobulk = pt.tl.PseudobulkSpace().compute(
            subset, target_col=groupby, groups_col=replicate_col, mode="sum"
        )
        # Attach covariates onto pseudobulk obs (mean for numeric; mode for categorical).
        usable_covs: list[str] = []
        for c in covs:
            if c in pseudobulk.obs.columns:
                usable_covs.append(c)
                continue
            if c not in subset.obs.columns:
                continue
            keys = list(
                zip(
                    pseudobulk.obs[groupby].astype(str),
                    pseudobulk.obs[replicate_col].astype(str),
                    strict=True,
                )
            )
            vals = pd.to_numeric(subset.obs[c], errors="coerce")
            tmp = subset.obs[[groupby, replicate_col]].copy()
            if vals.notna().any() and vals.isna().mean() < 0.5:
                tmp[c] = vals
                means = tmp.groupby([groupby, replicate_col], observed=False)[c].mean()
                pseudobulk.obs[c] = [float(means.get((g, r), np.nan)) for g, r in keys]
            else:
                tmp[c] = subset.obs[c].astype(str)
                mode = tmp.groupby([groupby, replicate_col], observed=False)[c].agg(
                    lambda s: s.value_counts().index[0] if len(s) else ""
                )
                pseudobulk.obs[c] = [str(mode.get((g, r), "")) for g, r in keys]
            usable_covs.append(c)

        design = build_deseq2_design(groupby, replicate_col=replicate_col, covariates=usable_covs)
        pds2 = pt.tl.PyDESeq2(adata=pseudobulk, design=design)
        pds2.fit()
        results = pds2.test_contrasts(
            pds2.contrast(column=groupby, baseline=reference, group_to_compare=group)
        )
        results = results.copy()
        results.insert(0, "contrast", f"{group}_vs_{reference}")
        results.insert(0, "method", "pydeseq2_pseudobulk")
        results.insert(0, "evidence_level", EVIDENCE_LEVEL_INFERENTIAL)
        results.insert(0, "expression_matrix", "layers_counts_pseudobulk_sum")
        results.insert(0, "design", design)
        return results
    table = exploratory_wilcoxon(subset, groupby=groupby, group=group, reference=reference)
    if not table.empty and covs:
        table.insert(0, "covariates_ignored", ",".join(covs))
    return table


def _de_one_contrast(
    args: tuple[AnnData, str, str, str | None, str, tuple[str, ...]],
) -> tuple[str, pd.DataFrame | None, str | None]:
    """Worker: one group vs reference on a pre-subset AnnData (picklable top-level)."""
    subset, group, reference, replicate_col, groupby, covariates = args
    try:
        table = run_deseq2_or_wilcoxon(
            subset,
            group=group,
            reference=reference,
            replicate_col=replicate_col,
            groupby=groupby,
            covariates=covariates,
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
    covariates: Sequence[str] | None = None,
) -> tuple[list[tuple[str, pd.DataFrame]], dict[str, str]]:
    """Run DE for each group vs reference; parallelizes the outer loop when n_jobs != 1.

    Each worker receives a copy of cells in {group, reference} only. CSV/plots stay
    on the caller so AnnData is never mutated across processes.
    """
    cov_tuple = tuple(covariates or ())
    tasks: list[tuple[AnnData, str, str, str | None, str, tuple[str, ...]]] = []
    labels = counts.obs[groupby].astype(str)
    for group in groups:
        g = str(group)
        if int((labels == g).sum()) < min_cells:
            continue
        mask = labels.isin([g, str(reference)])
        tasks.append((counts[mask].copy(), g, str(reference), replicate_col, groupby, cov_tuple))

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
    "DEFAULT_DE_COVARIATE_CANDIDATES",
    "EVIDENCE_LEVEL_EXPLORATORY",
    "EVIDENCE_LEVEL_INFERENTIAL",
    "build_deseq2_design",
    "check_experimental_design",
    "de_evidence_level",
    "ensure_de_covariate_columns",
    "exploratory_wilcoxon",
    "resolve_de_covariates",
    "run_de_contrasts",
    "run_deseq2_or_wilcoxon",
]
