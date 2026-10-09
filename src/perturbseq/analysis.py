from __future__ import annotations

import warnings

import numpy as np
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


def run_mixscape(
    adata: AnnData,
    control: str = "NT",
    split_by: str | None = None,
    n_neighbors: int = 20,
    perturbation_type: str = "KO",
) -> AnnData:
    pt = _require_pertpy()
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
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    pt = _require_pertpy()
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)
    counts = adata.obs[groupby].astype(str).value_counts()
    keep_groups = set(counts[counts >= min_cells].index)
    keep_groups.add(contrast)
    subset = adata[adata.obs[groupby].astype(str).isin(keep_groups)].copy()
    distance = pt.tl.Distance("edistance", obsm_key="X_pca")
    edistances = (
        distance.onesided_distances(subset, groupby=groupby, selected_group=contrast, show_progressbar=False)
        .drop(contrast, errors="ignore")
        .sort_values(ascending=False)
        .rename("edistance")
        .to_frame()
    )
    edistances["n_cells"] = counts.reindex(edistances.index)
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
        tested_names = subset.obs_names[subset.obs[groupby].astype(str).isin(tested)].to_numpy()
        etest_cells = np.concatenate([tested_names, control_names])
        etest = pt.tl.DistanceTest("edistance", n_perms=n_perms, obsm_key="X_pca")
        etest_results = etest(
            subset[etest_cells].copy(),
            groupby=groupby,
            contrast=contrast,
            show_progressbar=False,
        )
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"E-test skipped: {exc}", stacklevel=2)
    return edistances, etest_results


def exploratory_wilcoxon(adata: AnnData, groupby: str, group: str, reference: str = "NT") -> pd.DataFrame:
    mask = adata.obs[groupby].astype(str).isin([group, reference])
    subset = adata[mask].copy()
    if subset.obs[groupby].nunique() < 2:
        return pd.DataFrame()
    sc.tl.rank_genes_groups(subset, groupby=groupby, groups=[group], reference=reference, method="wilcoxon")
    table = sc.get.rank_genes_groups_df(subset, group=group)
    table.insert(0, "contrast", f"{group}_vs_{reference}")
    table.insert(0, "method", "wilcoxon_cell_level_exploratory")
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
        return results
    return exploratory_wilcoxon(subset, groupby=groupby, group=group, reference=reference)


def cluster_perturbations(adata: AnnData, groupby: str = "gene_target", n_neighbors: int = 5) -> pd.DataFrame:
    pt = _require_pertpy()
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)
    profiles = pt.tl.PseudobulkSpace().compute(adata, target_col=groupby, embedding_key="X_pca", mode="mean")
    k = min(n_neighbors, max(2, profiles.n_obs - 1))
    sc.pp.neighbors(profiles, use_rep="X", n_neighbors=k)
    try:
        sc.tl.leiden(profiles, resolution=1.0, flavor="igraph", n_iterations=2, key_added="pert_cluster")
    except Exception:
        profiles.obs["pert_cluster"] = "0"
    table = profiles.obs[["pert_cluster"]].copy()
    table.index.name = groupby
    return table
