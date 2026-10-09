"""Domain: Mixscape, E-distance, secondary distances, perturbation clustering.

Owns: perturbation-effect modeling algorithms and related helpers.
Does NOT own: pseudobulk DE / FDR (``statistics``), confidence flags (``robustness``),
or stage wiring (``stages.perturbation_modeling``). Guide QC math lives in ``guide_qc``
(re-exported hooks only).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from perturbseq._deps import require_pertpy
from perturbseq.guide_qc import compute_gene_guide_consistency, run_guide_qc
from perturbseq.parallel import parallel_map, resolve_n_jobs

# Secondary metrics tried after E-distance. Wasserstein needs OTT-JAX; skip if absent.
DEFAULT_SECONDARY_DISTANCE_METRICS: tuple[str, ...] = ("mmd", "wasserstein")
# Peidli et al.: E-test is underpowered at small n; ~50–100 safer, ~200 more stable.
DEFAULT_ETEST_POWER_MIN_CELLS = 50
# Bootstrap replicates for E-distance CI over cells (0 = skip). Modest default for speed.
DEFAULT_N_BOOTSTRAP = 100
DEFAULT_EDISTANCE_CI_LEVEL = 0.95


def _pairwise_mean_norm(a: np.ndarray, b: np.ndarray) -> float:
    """Mean Euclidean distance between all pairs of rows in a and b."""
    a2 = np.sum(a * a, axis=1)[:, None]
    b2 = np.sum(b * b, axis=1)[None, :]
    d2 = np.maximum(a2 + b2 - 2.0 * (a @ b.T), 0.0)
    return float(np.sqrt(d2).mean())


def energy_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Scalar energy distance between two point clouds (embedding rows)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[0] == 0 or y.shape[0] == 0:
        return float("nan")
    return (
        2.0 * _pairwise_mean_norm(x, y)
        - _pairwise_mean_norm(x, x)
        - _pairwise_mean_norm(y, y)
    )


def _bootstrap_edistance_one_group(
    args: tuple[str, np.ndarray, np.ndarray, int, float, int],
) -> tuple[str, float, float]:
    """Worker: percentile CI for energy distance by resampling cells in both groups."""
    group, group_emb, ctrl_emb, n_bootstrap, ci_level, seed = args
    rng = np.random.default_rng(int(seed))
    n_g = int(group_emb.shape[0])
    n_c = int(ctrl_emb.shape[0])
    if n_g < 2 or n_c < 2 or n_bootstrap <= 0:
        return str(group), float("nan"), float("nan")
    stats = np.empty(int(n_bootstrap), dtype=float)
    for i in range(int(n_bootstrap)):
        g_idx = rng.integers(0, n_g, size=n_g)
        c_idx = rng.integers(0, n_c, size=n_c)
        stats[i] = energy_distance(group_emb[g_idx], ctrl_emb[c_idx])
    alpha = 1.0 - float(ci_level)
    lo, hi = np.quantile(stats, [alpha / 2.0, 1.0 - alpha / 2.0])
    return str(group), float(lo), float(hi)


def annotate_edistance_bootstrap_ci(
    edistances: pd.DataFrame,
    adata: AnnData,
    *,
    groupby: str = "gene_target",
    contrast: str = "NT",
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    ci_level: float = DEFAULT_EDISTANCE_CI_LEVEL,
    random_state: int = 0,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Add ``edistance_ci_low`` / ``edistance_ci_high`` via cell bootstrap (percentile CI).

    ``n_bootstrap <= 0`` leaves the frame unchanged (no CI columns). Uses the same
    energy-distance definition as the unit tests / scPerturb spirit; does not require
    pertpy. Parallelizes over groups when ``n_jobs != 1``.
    """
    out = edistances.copy()
    n_boot = int(n_bootstrap)
    if n_boot <= 0 or out.empty:
        return out
    if "X_pca" not in adata.obsm:
        out["edistance_ci_low"] = np.nan
        out["edistance_ci_high"] = np.nan
        out["n_bootstrap"] = 0
        out["ci_level"] = float(ci_level)
        return out

    labels = adata.obs[groupby].astype(str)
    embedding = np.asarray(adata.obsm["X_pca"], dtype=float)
    ctrl_mask = labels == str(contrast)
    ctrl_emb = embedding[ctrl_mask.to_numpy()]
    tasks: list[tuple[str, np.ndarray, np.ndarray, int, float, int]] = []
    for i, group in enumerate(out.index.astype(str)):
        g_emb = embedding[(labels == group).to_numpy()]
        # Distinct seed per group for reproducibility across n_jobs.
        seed = int(random_state) + 1_000_003 * (i + 1)
        tasks.append((str(group), g_emb, ctrl_emb, n_boot, float(ci_level), seed))

    results = parallel_map(_bootstrap_edistance_one_group, tasks, n_jobs=n_jobs)
    lo_map = {g: lo for g, lo, _ in results}
    hi_map = {g: hi for g, _, hi in results}
    out["edistance_ci_low"] = [lo_map.get(str(g), float("nan")) for g in out.index.astype(str)]
    out["edistance_ci_high"] = [hi_map.get(str(g), float("nan")) for g in out.index.astype(str)]
    out["n_bootstrap"] = n_boot
    out["ci_level"] = float(ci_level)
    return out


def _attach_edistance_ci_to_etest(
    etest: pd.DataFrame,
    edistances: pd.DataFrame,
) -> pd.DataFrame:
    """Copy bootstrap CI columns from edistance table onto matching E-test rows."""
    out = etest.copy()
    for col in ("edistance_ci_low", "edistance_ci_high", "n_bootstrap", "ci_level"):
        if col in edistances.columns:
            out[col] = edistances[col].reindex(out.index).to_numpy()
    return out


def estimate_mixscape_cost(
    n_cells: int,
    n_targets: int,
    mixscape_max_targets: int,
) -> dict[str, Any]:
    """Rough O(n_cells × n_targets) cost hint for auto-skip transparency (not a benchmark)."""
    n_cells = int(max(0, n_cells))
    n_targets = int(max(0, n_targets))
    max_t = int(max(1, mixscape_max_targets))
    work = n_cells * n_targets
    work_at_max = n_cells * max_t
    # Heuristic: ~50 float features × 8 bytes × cells × targets, floored at 0.5 GB.
    approx_memory_gb = max(0.5, round(work * 50 * 8 / 1e9, 2))
    return {
        "note": (
            "Heuristic cost ∝ O(n_cells × n_targets); not a wall-clock benchmark. "
            "Genome-scale Mixscape is typically memory- and time-heavy versus the "
            f"default threshold (mixscape_max_targets={max_t})."
        ),
        "n_cells": n_cells,
        "n_targets": n_targets,
        "mixscape_max_targets": max_t,
        "approx_work_units": int(work),
        "approx_work_units_at_threshold": int(work_at_max),
        "approx_relative_to_threshold": round(n_targets / max_t, 2),
        "approx_memory_hint_gb": approx_memory_gb,
    }


def filter_cells_for_mixscape_targets(
    adata: AnnData,
    control: str,
    targets: list[str] | tuple[str, ...],
    *,
    gene_col: str = "gene_target",
) -> AnnData:
    """Keep control + selected gene targets only (subset Mixscape path)."""
    wanted = {str(control), *[str(t) for t in targets]}
    labels = adata.obs[gene_col].astype(str)
    keep = labels.isin(wanted)
    if not bool(keep.any()):
        raise ValueError(
            f"No cells for Mixscape subset: control={control!r}, targets={list(targets)[:10]}"
        )
    return adata[keep].copy()


def select_mixscape_targets_from_edistance(
    edistances: pd.DataFrame,
    *,
    top_n: int,
    control: str,
    exclude: set[str] | None = None,
) -> list[str]:
    """Pick top-N gene targets by E-distance rank (cheap pre-Mixscape ranking)."""
    if edistances is None or edistances.empty or top_n <= 0:
        return []
    skip = {str(control), "unassigned", ""} | {str(x) for x in (exclude or set())}
    ranked: list[str] = []
    for name in edistances.index.astype(str):
        if name in skip:
            continue
        ranked.append(name)
        if len(ranked) >= int(top_n):
            break
    return ranked


def merge_mixscape_annotations(full: AnnData, subset: AnnData) -> list[str]:
    """Copy mixscape_* obs columns from a subset run onto the full object."""
    copied: list[str] = []
    for col in subset.obs.columns:
        if not str(col).startswith("mixscape"):
            continue
        if col not in full.obs.columns:
            full.obs[col] = pd.Series(pd.NA, index=full.obs_names, dtype=object)
        full.obs[col] = full.obs[col].astype(object)
        full.obs.loc[subset.obs_names, col] = subset.obs[col].astype(object).to_numpy()
        copied.append(str(col))
    return copied


def run_mixscape(
    adata: AnnData,
    control: str = "NT",
    split_by: str | None = None,
    n_neighbors: int = 20,
    perturbation_type: str = "KO",
) -> AnnData:
    """Run pertpy Mixscape (signature → KO/NP). Requires ``pertpy``."""
    pt = require_pertpy()
    if control not in set(adata.obs["perturbation"].astype(str)) and control not in set(
        adata.obs["gene_target"].astype(str)
    ):
        raise ValueError(f"Control '{control}' not found in perturbation/gene_target labels")
    n_control = int((adata.obs["gene_target"].astype(str) == control).sum())
    if n_control < 10:
        raise ValueError(f"Only {n_control} control cells; Mixscape needs a larger NT pool")
    mixscape = pt.tl.Mixscape()
    signature_kwargs = {
        "pert_key": "perturbation",
        "control": control,
        "n_neighbors": min(n_neighbors, max(2, n_control)),
    }
    if split_by and split_by in adata.obs:
        signature_kwargs["split_by"] = split_by
    mixscape.perturbation_signature(adata, **signature_kwargs)
    mixscape.mixscape(
        adata,
        pert_key="gene_target",
        control=control,
        layer="X_pert",
        perturbation_type=perturbation_type,
    )
    try:
        mixscape.lda(adata, pert_key="gene_target", control=control, layer="X_pert")
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"Mixscape LDA skipped: {exc}", stacklevel=2)
    try:
        mixscale = pt.tl.Mixscale()
        mixscale.mixscale(adata, "gene_target", control, layer="X_pert")
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"Mixscale skipped: {exc}", stacklevel=2)
    return adata


def _etest_one_group(
    args: tuple[np.ndarray, np.ndarray, str, str, int],
) -> pd.DataFrame:
    """Worker: DistanceTest for one group vs contrast on a PCA embedding slice."""
    embedding, labels, group, contrast, n_perms = args
    pt = require_pertpy()
    from anndata import AnnData as _AnnData

    obs = pd.DataFrame({"_group": labels})
    ad = _AnnData(obs=obs)
    ad.obsm["X_pca"] = np.asarray(embedding)
    etest = pt.tl.DistanceTest("edistance", n_perms=n_perms, obsm_key="X_pca")
    return etest(ad, groupby="_group", contrast=contrast, show_progressbar=False)


def _reapply_etest_padj(table: pd.DataFrame, method: str = "holm-sidak") -> pd.DataFrame:
    """Recompute multiple-testing correction across concatenated per-group tests."""
    out = table.copy()
    pcol = None
    for cand in ("pvalue", "p_value", "pval"):
        if cand in out.columns:
            pcol = cand
            break
    if pcol is None:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        return out
    pvals = pd.to_numeric(out[pcol], errors="coerce").to_numpy()
    ok = np.isfinite(pvals)
    padj = np.full(len(pvals), np.nan, dtype=float)
    if ok.any():
        _, adj, _, _ = multipletests(pvals[ok], method=method)
        padj[ok] = adj
    for col in ("pvalue_adj", "padj", "p_adj"):
        if col in out.columns:
            out[col] = padj
            break
    else:
        out["pvalue_adj"] = padj
    return out


def annotate_etest_power(
    etest: pd.DataFrame,
    n_cells: pd.Series | dict[str, Any],
    power_min_cells: int = DEFAULT_ETEST_POWER_MIN_CELLS,
) -> pd.DataFrame:
    """Tag E-test rows with sample-size power flags.

    ``significant_adj`` remains the raw adjusted call. ``significant_adj_reported``
    is True only when significant *and* not ``low_power`` (n_cells < threshold).
    """
    out = etest.copy()
    if isinstance(n_cells, dict):
        cell_map = pd.Series(n_cells, dtype=float)
    else:
        cell_map = pd.to_numeric(n_cells, errors="coerce")
    out["n_cells"] = cell_map.reindex(out.index).to_numpy()
    n = pd.to_numeric(out["n_cells"], errors="coerce")
    out["low_power"] = (n < int(power_min_cells)) | ~np.isfinite(n)
    if "significant_adj" not in out.columns:
        padj = None
        for col in ("pvalue_adj", "padj", "p_adj"):
            if col in out.columns:
                padj = pd.to_numeric(out[col], errors="coerce")
                break
        if padj is not None:
            out["significant_adj"] = padj < 0.05
        else:
            out["significant_adj"] = False
    sig = out["significant_adj"].fillna(False).astype(bool)
    out["significant_adj_reported"] = sig & ~out["low_power"].astype(bool)
    return out


def _onesided_distance_frame(
    adata: AnnData,
    metric: str,
    *,
    groupby: str,
    contrast: str,
    min_cells: int,
    n_jobs: int,
    pca_source: str,
) -> pd.DataFrame:
    """Compute onesided distances for one pertpy metric on the same group filter as E-distance."""
    pt = require_pertpy()
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)
    counts = adata.obs[groupby].astype(str).value_counts()
    keep_groups = set(counts[counts >= min_cells].index)
    keep_groups.add(contrast)
    subset = adata[adata.obs[groupby].astype(str).isin(keep_groups)].copy()
    distance = pt.tl.Distance(metric, obsm_key="X_pca")
    onesided_kwargs: dict[str, Any] = {
        "groupby": groupby,
        "selected_group": contrast,
        "show_progressbar": False,
    }
    resolved = resolve_n_jobs(n_jobs)
    onesided_kwargs["n_jobs"] = resolved
    try:
        series = distance.onesided_distances(subset, **onesided_kwargs)
    except TypeError:
        onesided_kwargs.pop("n_jobs", None)
        series = distance.onesided_distances(subset, **onesided_kwargs)
    frame = (
        series.drop(contrast, errors="ignore")
        .sort_values(ascending=False)
        .rename(metric)
        .to_frame()
    )
    frame["n_cells"] = counts.reindex(frame.index)
    frame["embedding"] = "X_pca"
    frame["pca_source"] = pca_source
    frame["metric"] = metric
    return frame


def run_secondary_distances(
    adata: AnnData,
    *,
    groupby: str = "gene_target",
    contrast: str = "NT",
    min_cells: int = 10,
    metrics: tuple[str, ...] | list[str] = DEFAULT_SECONDARY_DISTANCE_METRICS,
    pca_source: str = "log1p_hvg",
    n_jobs: int = 1,
) -> tuple[dict[str, pd.DataFrame], list[dict[str, Any]]]:
    """Compute secondary distance metrics (MMD, Wasserstein, …) with graceful skips.

    Returns ``(tables_by_metric, statuses)`` where each status mirrors
    ``skipped`` + ``reason`` (+ optional ``detail``).
    """
    tables: dict[str, pd.DataFrame] = {}
    statuses: list[dict[str, Any]] = []
    for metric in metrics:
        metric = str(metric).strip()
        if not metric or metric == "edistance":
            continue
        try:
            tables[metric] = _onesided_distance_frame(
                adata,
                metric,
                groupby=groupby,
                contrast=contrast,
                min_cells=min_cells,
                n_jobs=n_jobs,
                pca_source=pca_source,
            )
            statuses.append({"metric": metric, "skipped": False, "reason": None})
        except ImportError as exc:
            reason = "missing_jax" if metric == "wasserstein" else "missing_dependency"
            detail = str(exc) or f"{metric} requires an optional dependency that is not installed"
            if metric == "wasserstein" and "jax" not in detail.lower() and "ott" not in detail.lower():
                # Still likely OTT-JAX for wasserstein.
                reason = "missing_jax"
            statuses.append(
                {"metric": metric, "skipped": True, "reason": reason, "detail": detail}
            )
            warnings.warn(f"Secondary distance '{metric}' skipped: {detail}", stacklevel=2)
        except Exception as exc:  # noqa: BLE001
            detail = str(exc)
            reason = "missing_jax" if (
                metric == "wasserstein"
                and any(tok in detail.lower() for tok in ("jax", "ott", "module"))
            ) else "failed"
            statuses.append(
                {"metric": metric, "skipped": True, "reason": reason, "detail": detail}
            )
            warnings.warn(f"Secondary distance '{metric}' skipped: {detail}", stacklevel=2)
    return tables, statuses


def combine_distance_tables(
    edistances: pd.DataFrame,
    secondary: dict[str, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    """Wide table: edistance + secondary metric columns (+ shared n_cells / provenance)."""
    out = edistances.copy()
    if "edistance" not in out.columns and out.shape[1] >= 1:
        # Defensive: already named.
        pass
    secondary = secondary or {}
    for metric, frame in secondary.items():
        if frame is None or frame.empty or metric not in frame.columns:
            continue
        out[metric] = frame[metric].reindex(out.index)
    return out


def run_edistance(
    adata: AnnData,
    groupby: str = "gene_target",
    contrast: str = "NT",
    n_perms: int = 1000,
    min_cells: int = 10,
    etest_top_n: int = 10,
    etest_random_n: int = 10,
    max_control_cells: int = 1000,
    random_state: int = 0,
    pca_source: str = "log1p_hvg",
    n_jobs: int = 1,
    power_min_cells: int = DEFAULT_ETEST_POWER_MIN_CELLS,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    ci_level: float = DEFAULT_EDISTANCE_CI_LEVEL,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """E-distance on ``obsm['X_pca']``. ``pca_source`` is recorded for report provenance only.

    When ``n_jobs != 1``, E-test runs one DistanceTest per group in parallel and
    re-applies holm-sidak across groups. ``n_jobs=1`` keeps a single joint
    DistanceTest call (original behavior). pertpy DistanceTest has no internal
    n_jobs hook; only the outer group loop is parallelized here.

    By default, bootstrap percentile CIs over cells are attached
    (``edistance_ci_low`` / ``edistance_ci_high``; ``n_bootstrap=0`` disables).
    E-test rows are annotated with ``n_cells``, ``low_power``, and
    ``significant_adj_reported`` (significance only trusted when not low-power).
    """
    edistances = _onesided_distance_frame(
        adata,
        "edistance",
        groupby=groupby,
        contrast=contrast,
        min_cells=min_cells,
        n_jobs=n_jobs,
        pca_source=pca_source,
    )
    # Drop helper column used by secondary metrics; keep edistance-focused schema.
    if "metric" in edistances.columns:
        edistances = edistances.drop(columns=["metric"])
    edistances = annotate_edistance_bootstrap_ci(
        edistances,
        adata,
        groupby=groupby,
        contrast=contrast,
        n_bootstrap=n_bootstrap,
        ci_level=ci_level,
        random_state=random_state,
        n_jobs=n_jobs,
    )
    pt = require_pertpy()
    counts = adata.obs[groupby].astype(str).value_counts()
    keep_groups = set(counts[counts >= min_cells].index)
    keep_groups.add(contrast)
    subset = adata[adata.obs[groupby].astype(str).isin(keep_groups)].copy()
    resolved = resolve_n_jobs(n_jobs)
    etest_results = None
    try:
        rng = np.random.default_rng(random_state)
        ranked = list(edistances.index)
        tested = ranked[:etest_top_n]
        rest = ranked[etest_top_n:]
        if rest and etest_random_n:
            tested += list(rng.choice(rest, size=min(etest_random_n, len(rest)), replace=False))
        control_names = subset.obs_names[subset.obs[groupby].astype(str) == contrast].to_numpy()
        if len(control_names) > max_control_cells:
            control_names = rng.choice(control_names, size=max_control_cells, replace=False)
        labels_all = subset.obs[groupby].astype(str)
        embedding = np.asarray(subset.obsm["X_pca"])
        name_to_idx = {n: i for i, n in enumerate(subset.obs_names)}
        ctrl_idx = np.array([name_to_idx[n] for n in control_names], dtype=int)

        if resolved == 1:
            tested_names = subset.obs_names[labels_all.isin(tested)].to_numpy()
            etest_cells = np.concatenate([tested_names, control_names])
            etest = pt.tl.DistanceTest("edistance", n_perms=n_perms, obsm_key="X_pca")
            etest_results = etest(
                subset[etest_cells].copy(),
                groupby=groupby,
                contrast=contrast,
                show_progressbar=False,
            )
        else:
            tasks = []
            for group in tested:
                g_idx = np.array(
                    [name_to_idx[n] for n in subset.obs_names[labels_all == group]],
                    dtype=int,
                )
                if g_idx.size == 0:
                    continue
                idx = np.concatenate([g_idx, ctrl_idx])
                emb = embedding[idx]
                labs = np.concatenate(
                    [np.full(g_idx.size, str(group)), np.full(ctrl_idx.size, str(contrast))]
                )
                tasks.append((emb, labs, str(group), str(contrast), int(n_perms)))
            pieces = parallel_map(_etest_one_group, tasks, n_jobs=n_jobs)
            if pieces:
                etest_results = _reapply_etest_padj(pd.concat(pieces, axis=0))
                # Preserve deterministic row order matching `tested`.
                order = [g for g in tested if g in etest_results.index]
                etest_results = etest_results.reindex(order)
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"E-test skipped: {exc}", stacklevel=2)
    if etest_results is not None:
        etest_results = annotate_etest_power(
            etest_results,
            edistances["n_cells"] if "n_cells" in edistances.columns else counts,
            power_min_cells=power_min_cells,
        )
        etest_results = _attach_edistance_ci_to_etest(etest_results, edistances)
    return edistances, etest_results


# Configurable perturbation-space construction (CLI --perturbation-space).
# pca_silhouette = historical default (mean PCA + graph clustering).
PERTURBATION_SPACE_METHODS: tuple[str, ...] = ("pca_silhouette", "kmeans", "lr_classifier")


@dataclass
class PerturbationSpaceResult:
    """Embeddings + cluster labels from ``build_perturbation_space``."""

    clusters: pd.DataFrame
    embeddings: pd.DataFrame | None
    method: str
    backend: str
    metadata: dict[str, Any] = field(default_factory=dict)


def _ensure_pca(adata: AnnData) -> None:
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)


def _mean_embedding_profiles(
    adata: AnnData,
    groupby: str,
    *,
    embedding_key: str = "X_pca",
) -> AnnData:
    """One observation per group = mean of ``obsm[embedding_key]`` rows."""
    labels = adata.obs[groupby].astype(str)
    emb = np.asarray(adata.obsm[embedding_key], dtype=float)
    groups = list(pd.unique(labels))
    rows = []
    for g in groups:
        mask = (labels == g).to_numpy()
        rows.append(emb[mask].mean(axis=0))
    profiles = AnnData(X=np.asarray(rows, dtype=float))
    profiles.obs_names = pd.Index([str(g) for g in groups], name=groupby)
    profiles.obs[groupby] = profiles.obs_names.astype(str)
    return profiles


def _cluster_profiles_leiden(profiles: AnnData, *, n_neighbors: int = 5) -> AnnData:
    """Graph clustering on a per-perturbation profile AnnData (``use_rep='X'``)."""
    if profiles.n_obs < 2:
        profiles.obs["pert_cluster"] = "0"
        return profiles
    k = min(int(n_neighbors), max(2, profiles.n_obs - 1))
    sc.pp.neighbors(profiles, use_rep="X", n_neighbors=k)
    try:
        sc.tl.leiden(profiles, resolution=1.0, flavor="igraph", n_iterations=2, key_added="pert_cluster")
    except Exception:
        profiles.obs["pert_cluster"] = "0"
    return profiles


def _cluster_profiles_kmeans(
    profiles: AnnData,
    *,
    n_clusters: int | None = None,
    random_state: int = 0,
) -> tuple[AnnData, dict[str, Any]]:
    """KMeans on profile rows; optional silhouette pick of k in [2, max]."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    meta: dict[str, Any] = {}
    n = int(profiles.n_obs)
    if n < 2:
        profiles.obs["pert_cluster"] = "0"
        meta["n_clusters"] = 1
        return profiles, meta
    x = np.asarray(profiles.X, dtype=float)
    max_k = min(8, n - 1)
    if n_clusters is not None:
        k = max(2, min(int(n_clusters), n))
        if k >= n:
            k = max(1, n - 1) if n > 1 else 1
        chosen = int(k)
        meta["k_selection"] = "user"
    else:
        best_k, best_score = 2, -1.0
        for k in range(2, max_k + 1):
            labels = KMeans(n_clusters=k, n_init=10, random_state=random_state).fit_predict(x)
            if len(set(labels)) < 2:
                continue
            score = float(silhouette_score(x, labels))
            if score > best_score:
                best_k, best_score = k, score
        chosen = int(best_k)
        meta["k_selection"] = "silhouette"
        meta["silhouette_score"] = best_score
    if chosen <= 1:
        profiles.obs["pert_cluster"] = "0"
        meta["n_clusters"] = 1
        return profiles, meta
    labels = KMeans(n_clusters=chosen, n_init=10, random_state=random_state).fit_predict(x)
    profiles.obs["pert_cluster"] = pd.Series(labels, index=profiles.obs_names).astype(str)
    meta["n_clusters"] = int(chosen)
    return profiles, meta


def _embeddings_frame(profiles: AnnData, groupby: str) -> pd.DataFrame:
    x = np.asarray(profiles.X, dtype=float)
    cols = [f"dim_{i}" for i in range(x.shape[1])]
    frame = pd.DataFrame(x, index=profiles.obs_names.astype(str), columns=cols)
    frame.index.name = groupby
    return frame


def _clusters_frame(profiles: AnnData, groupby: str) -> pd.DataFrame:
    table = profiles.obs[["pert_cluster"]].copy()
    table.index = profiles.obs_names.astype(str)
    table.index.name = groupby
    return table


def _space_pca_silhouette(
    adata: AnnData,
    *,
    groupby: str,
    n_neighbors: int,
) -> PerturbationSpaceResult:
    """Mean PCA per perturbation + neighbors/Leiden (historical default).

    Prefers pertpy ``PseudobulkSpace`` when importable; otherwise local mean of
    ``obsm['X_pca']`` (same biological intent). Method name keeps the proposal /
    docs label ``pca_silhouette``; clustering uses Leiden on the profile graph.
    """
    _ensure_pca(adata)
    backend = "local"
    try:
        pt = require_pertpy()
        profiles = pt.tl.PseudobulkSpace().compute(
            adata, target_col=groupby, embedding_key="X_pca", mode="mean"
        )
        backend = "pertpy"
    except ImportError:
        profiles = _mean_embedding_profiles(adata, groupby, embedding_key="X_pca")
    profiles = _cluster_profiles_leiden(profiles, n_neighbors=n_neighbors)
    return PerturbationSpaceResult(
        clusters=_clusters_frame(profiles, groupby),
        embeddings=_embeddings_frame(profiles, groupby),
        method="pca_silhouette",
        backend=backend,
        metadata={"embedding": "X_pca", "clusterer": "leiden", "n_neighbors": int(n_neighbors)},
    )


def _space_kmeans(
    adata: AnnData,
    *,
    groupby: str,
    n_clusters: int | None,
    random_state: int,
) -> PerturbationSpaceResult:
    """KMeans on per-perturbation mean PCA profiles (pertpy KMeansSpace-like intent).

    Local implementation: mean ``X_pca`` → sklearn KMeans (silhouette-selected k
    unless ``n_clusters`` is set). Does not call cell-level ``pt.tl.KMeansSpace``,
    which labels cells rather than producing one embedding per perturbation.
    """
    _ensure_pca(adata)
    profiles = _mean_embedding_profiles(adata, groupby, embedding_key="X_pca")
    profiles, km_meta = _cluster_profiles_kmeans(
        profiles, n_clusters=n_clusters, random_state=random_state
    )
    return PerturbationSpaceResult(
        clusters=_clusters_frame(profiles, groupby),
        embeddings=_embeddings_frame(profiles, groupby),
        method="kmeans",
        backend="local",
        metadata={
            "embedding": "X_pca_mean",
            "clusterer": "kmeans",
            "note": (
                "Local KMeans on perturbation-mean PCA profiles "
                "(pertpy KMeansSpace-like response clustering)."
            ),
            **km_meta,
        },
    )


def _lr_coefficients_local(
    adata: AnnData,
    *,
    groupby: str,
    embedding_key: str = "X_pca",
    random_state: int = 0,
    max_iter: int = 500,
) -> AnnData:
    """One-vs-rest logistic regression; coefficient rows = perturbation embeddings."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import LabelEncoder

    labels = adata.obs[groupby].astype(str)
    x = np.asarray(adata.obsm[embedding_key], dtype=float)
    classes = sorted(labels.unique())
    if len(classes) < 2:
        profiles = AnnData(X=np.zeros((len(classes), x.shape[1]), dtype=float))
        profiles.obs_names = pd.Index(classes, name=groupby)
        return profiles
    enc = LabelEncoder()
    y = enc.fit_transform(labels.to_numpy())
    # Multinomial / OVR coefficients: shape (n_classes, n_features) for multinomial.
    clf = LogisticRegression(
        multi_class="multinomial",
        solver="lbfgs",
        max_iter=int(max_iter),
        random_state=int(random_state),
    )
    clf.fit(x, y)
    coef = np.asarray(clf.coef_, dtype=float)
    # sklearn may drop a class column in some binary cases; align to enc.classes_.
    order = list(enc.classes_)
    if coef.shape[0] == 1 and len(order) == 2:
        # Binary: single coef row for class 1; class 0 ≈ -coef.
        coef = np.vstack([-coef[0], coef[0]])
    profiles = AnnData(X=coef)
    profiles.obs_names = pd.Index([str(c) for c in order], name=groupby)
    profiles.obs[groupby] = profiles.obs_names.astype(str)
    return profiles


def _space_lr_classifier(
    adata: AnnData,
    *,
    groupby: str,
    n_clusters: int | None,
    random_state: int,
) -> PerturbationSpaceResult:
    """LR coefficients as perturbation embedding (pertpy LRClassifierSpace-like).

    Prefers ``pt.tl.LRClassifierSpace`` when available; otherwise local multinomial
    logistic regression on ``X_pca``. Clusters coefficient rows with KMeans.
    """
    _ensure_pca(adata)
    try:
        pt = require_pertpy()
        profiles = pt.tl.LRClassifierSpace().compute(
            adata,
            target_col=groupby,
            embedding_key="X_pca",
            random_state=int(random_state),
        )
        backend = "pertpy"
    except ImportError:
        profiles = _lr_coefficients_local(
            adata, groupby=groupby, embedding_key="X_pca", random_state=random_state
        )
        backend = "local"
    except Exception as exc:  # noqa: BLE001 — prefer local over failing the stage
        warnings.warn(
            f"pertpy LRClassifierSpace failed ({exc}); using local logistic regression",
            stacklevel=2,
        )
        profiles = _lr_coefficients_local(
            adata, groupby=groupby, embedding_key="X_pca", random_state=random_state
        )
        backend = "local"
    profiles, km_meta = _cluster_profiles_kmeans(
        profiles, n_clusters=n_clusters, random_state=random_state
    )
    return PerturbationSpaceResult(
        clusters=_clusters_frame(profiles, groupby),
        embeddings=_embeddings_frame(profiles, groupby),
        method="lr_classifier",
        backend=backend,
        metadata={
            "embedding": "lr_coefficients",
            "clusterer": "kmeans",
            "note": (
                "Logistic-regression coefficients as perturbation embedding "
                "(pertpy LRClassifierSpace when available, else local sklearn)."
            ),
            **km_meta,
        },
    )


def build_perturbation_space(
    adata: AnnData,
    method: str = "pca_silhouette",
    *,
    groupby: str = "gene_target",
    n_neighbors: int = 5,
    n_clusters: int | None = None,
    random_state: int = 0,
) -> PerturbationSpaceResult:
    """Build a perturbation-level embedding and cluster labels.

    Parameters
    ----------
    method
        ``pca_silhouette`` (default) — mean PCA + Leiden (backward-compatible).
        ``kmeans`` — KMeans on mean PCA profiles (response-pattern clustering).
        ``lr_classifier`` — LR coefficient embedding + KMeans (marker-like space).
    """
    method = str(method).strip().lower()
    if method not in PERTURBATION_SPACE_METHODS:
        raise ValueError(
            f"perturbation_space must be one of {PERTURBATION_SPACE_METHODS}, got {method!r}"
        )
    if method == "pca_silhouette":
        return _space_pca_silhouette(adata, groupby=groupby, n_neighbors=n_neighbors)
    if method == "kmeans":
        return _space_kmeans(
            adata, groupby=groupby, n_clusters=n_clusters, random_state=random_state
        )
    return _space_lr_classifier(
        adata, groupby=groupby, n_clusters=n_clusters, random_state=random_state
    )


def cluster_perturbations(
    adata: AnnData,
    groupby: str = "gene_target",
    n_neighbors: int = 5,
    *,
    method: str = "pca_silhouette",
    n_clusters: int | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    """Cluster perturbations; thin wrapper over ``build_perturbation_space``."""
    result = build_perturbation_space(
        adata,
        method=method,
        groupby=groupby,
        n_neighbors=n_neighbors,
        n_clusters=n_clusters,
        random_state=random_state,
    )
    return result.clusters


__all__ = [
    "DEFAULT_EDISTANCE_CI_LEVEL",
    "DEFAULT_ETEST_POWER_MIN_CELLS",
    "DEFAULT_N_BOOTSTRAP",
    "DEFAULT_SECONDARY_DISTANCE_METRICS",
    "PERTURBATION_SPACE_METHODS",
    "PerturbationSpaceResult",
    "annotate_edistance_bootstrap_ci",
    "annotate_etest_power",
    "build_perturbation_space",
    "cluster_perturbations",
    "combine_distance_tables",
    "compute_gene_guide_consistency",
    "energy_distance",
    "estimate_mixscape_cost",
    "filter_cells_for_mixscape_targets",
    "merge_mixscape_annotations",
    "run_edistance",
    "run_guide_qc",
    "run_mixscape",
    "run_secondary_distances",
    "select_mixscape_targets_from_edistance",
]
