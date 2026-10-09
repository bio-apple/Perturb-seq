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

METHOD_PYDESEQ2 = "pydeseq2_pseudobulk"
METHOD_PYDESEQ2_NO_BIO_REPS = "pydeseq2_pseudobulk_no_bio_reps"
METHOD_WILCOXON = "wilcoxon_cell_level_exploratory"

DEFAULT_SAMPLE_COL = "sample_id"
PSEUDO_REPLICATE_COL = "_de_pseudo_replicate"

# Preferred covariate obs columns when ``de_covariates`` is enabled (available-only).
DEFAULT_DE_COVARIATE_CANDIDATES: tuple[str, ...] = (
    "phase",
    "pct_counts_mt",
    "log_n_counts",
)


def de_evidence_level(*, replicate_aware: bool) -> EvidenceLevel:
    """Map design to evidence_level: true bio replicates → inferential; else exploratory."""
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


def _n_levels(adata: AnnData, col: str | None) -> int:
    if not col or col not in adata.obs:
        return 0
    return int(adata.obs[col].astype(str).nunique())


def assign_pseudo_replicates(
    adata: AnnData,
    *,
    groupby: str,
    n_pseudo_replicates: int,
    random_state: int = 0,
    col: str = PSEUDO_REPLICATE_COL,
) -> AnnData:
    """Randomly split cells within each ``groupby`` level into N technical pseudo-replicates.

    Labels are stored on ``obs[col]``. This does **not** create biological replicates —
    use only for exploratory shared-dispersion DE and always keep ``evidence_level=exploratory``.
    """
    if n_pseudo_replicates < 2:
        raise ValueError("n_pseudo_replicates must be >= 2")
    out = adata.copy()
    rng = np.random.default_rng(int(random_state))
    labels = np.empty(out.n_obs, dtype=object)
    groups = out.obs[groupby].astype(str)
    for g in groups.unique():
        idx = np.flatnonzero(groups.to_numpy() == g)
        if len(idx) == 0:
            continue
        assign = rng.integers(0, int(n_pseudo_replicates), size=len(idx))
        labels[idx] = [f"pseudo_{i}" for i in assign]
    out.obs[col] = pd.Categorical(labels)
    return out


def resolve_de_aggregation(
    adata: AnnData,
    *,
    replicate_col: str | None = None,
    sample_col: str = DEFAULT_SAMPLE_COL,
    n_pseudo_replicates: int | None = None,
    de_prefer_pseudobulk: bool = True,
) -> dict[str, Any]:
    """Choose aggregation column / method for DE without running the contrast.

    Priority:
    1. True biological ``replicate_col`` with ≥2 levels → inferential PyDESeq2
    2. If ``de_prefer_pseudobulk`` and ``n_pseudo_replicates`` ≥ 2 → exploratory PyDESeq2
    3. If ``de_prefer_pseudobulk`` and ``sample_col`` has ≥2 levels → exploratory PyDESeq2
    4. Else cell-level Wilcoxon (exploratory fallback)
    """
    n_reps = _n_levels(adata, replicate_col)
    replicate_aware = bool(replicate_col and n_reps >= 2)
    n_samples = _n_levels(adata, sample_col)
    n_pseudo = int(n_pseudo_replicates) if n_pseudo_replicates is not None else 0

    if replicate_aware:
        return {
            "replicate_aware": True,
            "groups_col": str(replicate_col),
            "method": METHOD_PYDESEQ2,
            "evidence_level": EVIDENCE_LEVEL_INFERENTIAL,
            "use_pseudo_replicates": False,
            "n_replicates": n_reps,
            "n_samples": n_samples,
            "n_pseudo_replicates": n_pseudo if n_pseudo >= 2 else None,
            "include_groups_in_design": True,
        }

    if de_prefer_pseudobulk and n_pseudo >= 2:
        return {
            "replicate_aware": False,
            "groups_col": PSEUDO_REPLICATE_COL,
            "method": METHOD_PYDESEQ2_NO_BIO_REPS,
            "evidence_level": EVIDENCE_LEVEL_EXPLORATORY,
            "use_pseudo_replicates": True,
            "n_replicates": n_reps if replicate_col else 1,
            "n_samples": n_samples,
            "n_pseudo_replicates": n_pseudo,
            # Artificial splits are not paired across conditions — design is ~ groupby (+ covs).
            "include_groups_in_design": False,
        }

    if de_prefer_pseudobulk and n_samples >= 2:
        return {
            "replicate_aware": False,
            "groups_col": str(sample_col),
            "method": METHOD_PYDESEQ2_NO_BIO_REPS,
            "evidence_level": EVIDENCE_LEVEL_EXPLORATORY,
            "use_pseudo_replicates": False,
            "n_replicates": n_reps if replicate_col else 1,
            "n_samples": n_samples,
            "n_pseudo_replicates": None,
            "include_groups_in_design": True,
        }

    return {
        "replicate_aware": False,
        "groups_col": None,
        "method": METHOD_WILCOXON,
        "evidence_level": EVIDENCE_LEVEL_EXPLORATORY,
        "use_pseudo_replicates": False,
        "n_replicates": n_reps if replicate_col else 1,
        "n_samples": n_samples,
        "n_pseudo_replicates": n_pseudo if n_pseudo >= 2 else None,
        "include_groups_in_design": False,
    }


def check_experimental_design(
    adata: AnnData,
    replicate_col: str | None = None,
    groupby: str = "gene_target",
    covariates: Sequence[str] | None = None,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    n_pseudo_replicates: int | None = None,
    de_prefer_pseudobulk: bool = True,
) -> dict[str, Any]:
    """Summarize DE design: prefer pseudobulk; Wilcoxon only as loud exploratory fallback."""
    n_groups = int(adata.obs[groupby].astype(str).nunique()) if groupby in adata.obs else 0
    covs = list(covariates or [])
    plan = resolve_de_aggregation(
        adata,
        replicate_col=replicate_col,
        sample_col=sample_col,
        n_pseudo_replicates=n_pseudo_replicates,
        de_prefer_pseudobulk=de_prefer_pseudobulk,
    )
    method = plan["method"]
    replicate_aware = bool(plan["replicate_aware"])
    note: str | None = None
    design_formula: str | None = None

    if method == METHOD_PYDESEQ2:
        design_formula = build_deseq2_design(
            groupby, replicate_col=plan["groups_col"], covariates=covs
        )
        note = None
    elif method == METHOD_PYDESEQ2_NO_BIO_REPS:
        design_rep = plan["groups_col"] if plan["include_groups_in_design"] else None
        design_formula = build_deseq2_design(groupby, replicate_col=design_rep, covariates=covs)
        if plan["use_pseudo_replicates"]:
            note = (
                f"No true biological replicates: aggregating cells into "
                f"{plan['n_pseudo_replicates']} technical pseudo-replicates for exploratory "
                f"PyDESeq2 (method={METHOD_PYDESEQ2_NO_BIO_REPS}). Pseudo-replicates are not a "
                "substitute for biological replicates (Squair et al. 2021)."
            )
        else:
            note = (
                f"No true biological replicate column with ≥2 levels: pseudobulk by "
                f"{plan['groups_col']!r} + {groupby!r} then exploratory PyDESeq2 "
                f"(method={METHOD_PYDESEQ2_NO_BIO_REPS}). Shared dispersion helps ranking but "
                "does not replace biological replicates (Squair et al. 2021)."
            )
    else:
        note = (
            "DE falls back to cell-level Wilcoxon (exploratory only). Treating cells as "
            "independent samples massively inflates significance — prefer biological "
            "replicates, multi-sample pseudobulk, or --pseudo-replicates N for exploratory "
            f"PyDESeq2. method={METHOD_WILCOXON}."
        )
        if covs:
            note += (
                f" Requested covariates ({', '.join(covs)}) are recorded for provenance but "
                "are not modeled by Wilcoxon (scanpy rank_genes_groups has no covariate design)."
            )
        if de_prefer_pseudobulk and (plan["n_samples"] or 0) < 2 and not plan["use_pseudo_replicates"]:
            note += (
                " Single sample_id (or missing): cannot form ≥2 pseudobulk units without "
                "n_pseudo_replicates≥2."
            )

    return {
        "replicate_col": replicate_col,
        "n_replicates": plan["n_replicates"],
        "n_samples": plan["n_samples"],
        "n_groups": n_groups,
        "replicate_aware": replicate_aware,
        "covariates": covs,
        "design_formula": design_formula,
        "recommended_method": method,
        "evidence_level": plan["evidence_level"],
        "de_prefer_pseudobulk": de_prefer_pseudobulk,
        "n_pseudo_replicates": plan["n_pseudo_replicates"],
        "pseudobulk_groups_col": plan["groups_col"],
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
    table.insert(0, "method", METHOD_WILCOXON)
    table.insert(0, "evidence_level", EVIDENCE_LEVEL_EXPLORATORY)
    table.insert(0, "expression_matrix", "log1p_normalized_X")
    return table


def _attach_covariates_to_pseudobulk(
    subset: AnnData,
    pseudobulk: AnnData,
    groupby: str,
    groups_col: str,
    covs: Sequence[str],
) -> list[str]:
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
                pseudobulk.obs[groups_col].astype(str),
                strict=True,
            )
        )
        vals = pd.to_numeric(subset.obs[c], errors="coerce")
        tmp = subset.obs[[groupby, groups_col]].copy()
        if vals.notna().any() and vals.isna().mean() < 0.5:
            tmp[c] = vals
            means = tmp.groupby([groupby, groups_col], observed=False)[c].mean()
            pseudobulk.obs[c] = [float(means.get((g, r), np.nan)) for g, r in keys]
        else:
            tmp[c] = subset.obs[c].astype(str)
            mode = tmp.groupby([groupby, groups_col], observed=False)[c].agg(
                lambda s: s.value_counts().index[0] if len(s) else ""
            )
            pseudobulk.obs[c] = [str(mode.get((g, r), "")) for g, r in keys]
        usable_covs.append(c)
    return usable_covs


def _run_pydeseq2_pseudobulk(
    subset: AnnData,
    *,
    group: str,
    reference: str,
    groupby: str,
    groups_col: str,
    covariates: Sequence[str],
    method: str,
    evidence_level: EvidenceLevel,
    include_groups_in_design: bool,
) -> pd.DataFrame:
    pt = require_pertpy()
    if "counts" in subset.layers:
        subset.X = subset.layers["counts"]
    pseudobulk = pt.tl.PseudobulkSpace().compute(
        subset, target_col=groupby, groups_col=groups_col, mode="sum"
    )
    usable_covs = _attach_covariates_to_pseudobulk(
        subset, pseudobulk, groupby, groups_col, covariates
    )
    design_rep = groups_col if include_groups_in_design else None
    design = build_deseq2_design(groupby, replicate_col=design_rep, covariates=usable_covs)
    pds2 = pt.tl.PyDESeq2(adata=pseudobulk, design=design)
    pds2.fit()
    results = pds2.test_contrasts(
        pds2.contrast(column=groupby, baseline=reference, group_to_compare=group)
    )
    results = results.copy()
    results.insert(0, "contrast", f"{group}_vs_{reference}")
    results.insert(0, "method", method)
    results.insert(0, "evidence_level", evidence_level)
    results.insert(0, "expression_matrix", "layers_counts_pseudobulk_sum")
    results.insert(0, "design", design)
    return results


def run_deseq2_or_wilcoxon(
    counts: AnnData,
    group: str,
    reference: str = "NT",
    replicate_col: str | None = None,
    groupby: str = "mixscape_class",
    covariates: Sequence[str] | None = None,
    *,
    sample_col: str = DEFAULT_SAMPLE_COL,
    n_pseudo_replicates: int | None = None,
    de_prefer_pseudobulk: bool = True,
    random_state: int = 0,
) -> pd.DataFrame:
    """Prefer pseudobulk + PyDESeq2; Wilcoxon only when ≥2 pseudobulk units are unavailable."""
    mask = counts.obs[groupby].astype(str).isin([group, reference])
    subset = counts[mask].copy()
    covs = resolve_de_covariates(subset, list(covariates) if covariates else None)
    plan = resolve_de_aggregation(
        subset,
        replicate_col=replicate_col,
        sample_col=sample_col,
        n_pseudo_replicates=n_pseudo_replicates,
        de_prefer_pseudobulk=de_prefer_pseudobulk,
    )

    if plan["method"] in {METHOD_PYDESEQ2, METHOD_PYDESEQ2_NO_BIO_REPS}:
        work = subset
        groups_col = str(plan["groups_col"])
        if plan["use_pseudo_replicates"]:
            work = assign_pseudo_replicates(
                subset,
                groupby=groupby,
                n_pseudo_replicates=int(plan["n_pseudo_replicates"]),
                random_state=random_state,
                col=PSEUDO_REPLICATE_COL,
            )
            groups_col = PSEUDO_REPLICATE_COL
        return _run_pydeseq2_pseudobulk(
            work,
            group=group,
            reference=reference,
            groupby=groupby,
            groups_col=groups_col,
            covariates=covs,
            method=str(plan["method"]),
            evidence_level=plan["evidence_level"],
            include_groups_in_design=bool(plan["include_groups_in_design"]),
        )

    table = exploratory_wilcoxon(subset, groupby=groupby, group=group, reference=reference)
    if not table.empty and covs:
        table.insert(0, "covariates_ignored", ",".join(covs))
    return table


def _de_one_contrast(
    args: tuple[
        AnnData,
        str,
        str,
        str | None,
        str,
        tuple[str, ...],
        str,
        int | None,
        bool,
        int,
    ],
) -> tuple[str, pd.DataFrame | None, str | None]:
    """Worker: one group vs reference on a pre-subset AnnData (picklable top-level)."""
    (
        subset,
        group,
        reference,
        replicate_col,
        groupby,
        covariates,
        sample_col,
        n_pseudo_replicates,
        de_prefer_pseudobulk,
        random_state,
    ) = args
    try:
        table = run_deseq2_or_wilcoxon(
            subset,
            group=group,
            reference=reference,
            replicate_col=replicate_col,
            groupby=groupby,
            covariates=covariates,
            sample_col=sample_col,
            n_pseudo_replicates=n_pseudo_replicates,
            de_prefer_pseudobulk=de_prefer_pseudobulk,
            random_state=random_state,
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
    sample_col: str = DEFAULT_SAMPLE_COL,
    n_pseudo_replicates: int | None = None,
    de_prefer_pseudobulk: bool = True,
    random_state: int = 0,
) -> tuple[list[tuple[str, pd.DataFrame]], dict[str, str]]:
    """Run DE for each group vs reference; parallelizes the outer loop when n_jobs != 1.

    Each worker receives a copy of cells in {group, reference} only. CSV/plots stay
    on the caller so AnnData is never mutated across processes.
    """
    cov_tuple = tuple(covariates or ())
    tasks: list[
        tuple[AnnData, str, str, str | None, str, tuple[str, ...], str, int | None, bool, int]
    ] = []
    labels = counts.obs[groupby].astype(str)
    for group in groups:
        g = str(group)
        if int((labels == g).sum()) < min_cells:
            continue
        mask = labels.isin([g, str(reference)])
        tasks.append(
            (
                counts[mask].copy(),
                g,
                str(reference),
                replicate_col,
                groupby,
                cov_tuple,
                sample_col,
                n_pseudo_replicates,
                de_prefer_pseudobulk,
                int(random_state),
            )
        )

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
    "DEFAULT_SAMPLE_COL",
    "EVIDENCE_LEVEL_EXPLORATORY",
    "EVIDENCE_LEVEL_INFERENTIAL",
    "METHOD_PYDESEQ2",
    "METHOD_PYDESEQ2_NO_BIO_REPS",
    "METHOD_WILCOXON",
    "PSEUDO_REPLICATE_COL",
    "assign_pseudo_replicates",
    "build_deseq2_design",
    "check_experimental_design",
    "de_evidence_level",
    "ensure_de_covariate_columns",
    "exploratory_wilcoxon",
    "resolve_de_aggregation",
    "resolve_de_covariates",
    "run_de_contrasts",
    "run_deseq2_or_wilcoxon",
]
