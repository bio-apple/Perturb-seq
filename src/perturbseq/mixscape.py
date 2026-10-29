"""Mixscape classification and control-relative perturbation signature."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from anndata import AnnData
from scipy import sparse

from perturbseq.assign import STATUS_NT, control_cell_mask, perturbation_analysis_mask
from perturbseq.config import DEFAULT_SIGNATURE_NDIMS, ZERO_GUIDE_LABEL, control_spec

PERTURBED = "perturbed"
NON_PERTURBED = "non_perturbed"
AMBIGUOUS = "ambiguous"
NOT_EVALUATED = "not_evaluated"

# Seurat Mixscape vignette used CalcPerturbSig ndims=40; this pipeline keeps
# preprocess PCA and signature width aligned at DEFAULT_SIGNATURE_NDIMS (30).
# num.neighbors=20. logfc_threshold passed to pertpy is 0.1 (API default is 0.25).
# https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscape.html
# https://satijalab.org/seurat/articles/mixscape_vignette
VIGNETTE_CONTROL_NEIGHBORS = 20


def perturbation_type(config: dict[str, Any]) -> str:
    mode = config["experiment"]["crispr_mode"]
    return {"CRISPR_KO": "KO", "CRISPRi": "KD", "CRISPRa": "OE"}[mode]


def mixscape_call_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Arguments for ``pertpy.tools.Mixscape.mixscape``.

    Reads ``mixscape.*`` only. Mixscale thresholds are not used. pertpy's
    ``mixscape`` has no ``max_de_genes`` argument.
    """
    mix = config.get("mixscape", {})
    perturb = config.get("perturb", {})
    if not isinstance(mix, dict):
        mix = {}
    if not isinstance(perturb, dict):
        perturb = {}
    return {
        "min_de_genes": int(mix.get("min_de_genes", perturb.get("mixscape_min_de_genes", 5))),
        "logfc_threshold": float(mix.get("logfc_threshold", perturb.get("mixscape_logfc", 0.1))),
        "pval_cutoff": float(mix.get("pval_cutoff", perturb.get("mixscape_pval", 0.05))),
    }


def _stamp_control_label(obs, mask: np.ndarray, column: str, label: str) -> None:
    """Write the control category onto ``pert_key`` before the pertpy call."""
    values = obs[column].astype(str).to_numpy(copy=True)
    values[np.asarray(mask)] = label
    obs[column] = values


def signature_n_dims(n_obs: int, n_vars: int, requested: int = DEFAULT_SIGNATURE_NDIMS) -> int:
    """Cap ``n_dims`` by the representation width and the number of cells."""
    return max(1, min(int(requested), max(int(n_obs) - 1, 1), int(n_vars)))


def neighbor_n_dims(work: AnnData, requested: int) -> int:
    """``n_dims`` for ``perturbation_signature``.

    Does not scale ``X`` and does not replace ``X_pca``. pertpy uses an existing
    ``X_pca`` for neighbors (``use_rep=None``) and otherwise computes PCA with
    its default parameters. The signature itself is unscaled ``X``.
    """
    if "X_pca" in work.obsm:
        width = int(work.obsm["X_pca"].shape[1])
    else:
        width = int(work.n_vars)
    return signature_n_dims(work.n_obs, width, requested)


def require_unscaled_log1p(adata: AnnData, step: str) -> None:
    """pertpy subtracts ``expm1(X)``; scaled or raw counts are the wrong input."""
    if "log1p" not in adata.uns:
        raise ValueError(f"{step} 需要未 scale 的 log1p 矩阵（adata.X）。")
    matrix = adata.X
    if sparse.issparse(matrix):
        values = matrix.data
    else:
        values = np.asarray(matrix).ravel()
    if values.size and float(np.min(values)) < 0:
        raise ValueError(f"{step} 的 adata.X 含负值，不能把 scale 后的矩阵交给 perturbation_signature。")


def run_mixscape(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """``perturbation_signature`` then ``mixscape`` on analysis-eligible cells ∪ control.

    Input follows pertpy: unscaled log1p ``.X`` for the signature and for DE
    (``de_layer=None``), ``layer="X_pert"`` for the GMM. Every target in the
    mask is passed through. pertpy labels a target NP when it has fewer than
    ``min_de_genes`` DE genes. Posterior probability above 0.5 is the perturbed
    class. The signature UMAP is stored separately from the expression UMAP.
    """
    from pertpy.tools import Mixscape

    spec = control_spec(config)
    control = spec["control_label"]
    control_cells = control_cell_mask(adata.obs, spec)
    adata.uns["control_source"] = spec["control_source"]
    adata.uns["control_label"] = control
    adata.uns["control_temporary"] = bool(spec["temporary"])
    adata.uns["control_phrase"] = spec["control_phrase"]
    adata.uns["control_note"] = spec["note"]
    adata.uns["control_fallback_from"] = str(spec.get("fallback_from", ""))
    if spec["temporary"] and int(control_cells.sum()) == 0:
        adata.obs["mixscape_class"] = NOT_EVALUATED
        adata.obs["mixscape_class_global"] = NOT_EVALUATED
        adata.obs["mixscape_class_p"] = np.nan
        adata.obs["mixscape_status"] = NOT_EVALUATED
        adata.uns["mixscape_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        _stamp_control_label(adata.obs, control_cells, "target_gene", control)
    # Embedding may include all QC-pass cells; scoring stays on analysis_eligible ∪ control.
    pert_mask = perturbation_analysis_mask(adata.obs, spec)
    mixscape_settings = config.get("mixscape", {})
    seed = int(config["runtime"]["seed"])
    require_unscaled_log1p(adata, "Mixscape")
    params = mixscape_call_parameters(config)
    labels = adata.obs["target_gene"].astype(str)
    adata.obs["mixscape_class"] = NOT_EVALUATED
    adata.obs["mixscape_class_global"] = NOT_EVALUATED
    adata.obs["mixscape_class_p"] = np.nan
    adata.obs["mixscape_status"] = NOT_EVALUATED
    label = perturbation_type(config)
    if not bool(np.any(pert_mask)):
        adata.uns["mixscape_scope"] = "unstratified"
        adata.uns["mixscape_umap_status"] = "no_analysis_cells"
        adata.uns["mixscape_parameters"] = {
            "pert_key": "target_gene",
            "control": control,
            "split_by": "",
            "layer": "X_pert",
            "de_layer": None,
            "posterior_cutoff": 0.5,
            "perturbation_type": label,
            "random_state": seed,
            **params,
        }
        return adata
    targets = {
        str(gene)
        for gene in labels[pert_mask].unique()
        if str(gene) != control
    }
    if not targets:
        nt = (labels == control).to_numpy() & pert_mask
        adata.obs.loc[nt, "mixscape_class"] = control
        adata.obs.loc[nt, "mixscape_status"] = (
            ZERO_GUIDE_LABEL if control == ZERO_GUIDE_LABEL else STATUS_NT
        )
        adata.uns["mixscape_scope"] = "unstratified"
        adata.uns["mixscape_umap_status"] = "no_target_cells"
        adata.uns["mixscape_parameters"] = {
            "pert_key": "target_gene",
            "control": control,
            "split_by": "",
            "layer": "X_pert",
            "de_layer": None,
            "posterior_cutoff": 0.5,
            "perturbation_type": label,
            "random_state": seed,
            **params,
        }
        return adata
    keep = pert_mask
    # Full unscaled log1p matrix. pertpy reads X for the signature and for DE,
    # and X_pert for the GMM. Do not subset to HVGs or scale X first.
    work = adata[keep].copy()
    requested_dims = int(mixscape_settings.get("n_dims", DEFAULT_SIGNATURE_NDIMS))
    n_dims = neighbor_n_dims(work, requested_dims)
    n_nt = int((work.obs["target_gene"].astype(str) == control).sum())
    if n_nt < 1:
        raise ValueError(f"没有对照细胞（{control}），无法计算 perturbation signature。")
    requested_neighbors = int(
        mixscape_settings.get("n_control_neighbors", VIGNETTE_CONTROL_NEIGHBORS)
    )
    n_neighbors = min(requested_neighbors, n_nt, max(work.n_obs - 1, 1))
    if n_neighbors < requested_neighbors:
        adata.uns["mixscape_neighbor_note"] = (
            f"n_control_neighbors reduced from {requested_neighbors} to {n_neighbors}"
        )
    split_by = real_cell_type_split(work, config)
    signature_args: dict[str, Any] = {
        "pert_key": "target_gene",
        "control": control,
        "n_neighbors": n_neighbors,
        "n_dims": n_dims,
    }
    if split_by:
        signature_args["split_by"] = split_by
    model = Mixscape()
    model.perturbation_signature(work, **signature_args)
    embedding, umap_status = _signature_umap(work, seed, n_dims)
    # layer="X_pert" scores the signature; de_layer=None runs DE on unscaled X.
    # https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscape.html
    mix_args: dict[str, Any] = {
        "pert_key": "target_gene",
        "control": control,
        "layer": "X_pert",
        "de_layer": None,
        "min_de_genes": params["min_de_genes"],
        "logfc_threshold": params["logfc_threshold"],
        "pval_cutoff": params["pval_cutoff"],
        "iter_num": int(mixscape_settings.get("iter_num", 10)),
        "perturbation_type": label,
        "random_state": seed,
        "scale": True,
    }
    if split_by:
        mix_args["split_by"] = split_by
    model.mixscape(work, **mix_args)
    _write_signature_rows(adata, work)
    full_embedding = np.full((adata.n_obs, 2), np.nan)
    if embedding is not None:
        positions = adata.obs_names.get_indexer(work.obs_names)
        full_embedding[positions] = embedding
        adata.obsm["X_umap_mixscape"] = full_embedding
    adata.obs.loc[work.obs_names, "mixscape_class"] = (
        work.obs["mixscape_class"].astype(str).to_numpy()
    )
    adata.obs.loc[work.obs_names, "mixscape_class_global"] = (
        work.obs["mixscape_class_global"].astype(str).to_numpy()
    )
    probability_key = f"mixscape_class_p_{label.lower()}"
    if probability_key in work.obs:
        adata.obs.loc[work.obs_names, "mixscape_class_p"] = work.obs[probability_key].to_numpy()
    adata.obs.loc[work.obs_names, "mixscape_status"] = [
        _normalize_mixscape(value, control, label)
        for value in work.obs["mixscape_class"].astype(str)
    ]
    if "mixscape" in work.uns:
        adata.uns["mixscape"] = work.uns["mixscape"]
    adata.uns["mixscape_scope"] = "unstratified" if not split_by else str(split_by)
    adata.uns["mixscape_umap_status"] = umap_status if embedding is not None else umap_status
    adata.uns["mixscape_parameters"] = {
        "pert_key": "target_gene",
        "control": control,
        "n_neighbors": n_neighbors,
        "n_dims": n_dims,
        "split_by": split_by or "",
        "layer": "X_pert",
        "posterior_cutoff": 0.5,
        "perturbation_type": label,
        "random_state": seed,
        "umap_neighbors_metric": "cosine",
        "umap_key": "X_umap_mixscape",
        "iter_num": int(mixscape_settings.get("iter_num", 10)),
        # Recorded for provenance; call omits test_method so pertpy default applies.
        "test_method": "wilcoxon",
        "scale": True,
        "control_relative": (
            "perturbation_signature on unscaled log1p X (all genes); "
            "neighbors from existing X_pca when present (use_rep=None); "
            "mixscape layer=X_pert, de_layer=None (DE on X), scale=True"
        ),
        **params,
    }
    return adata


def real_cell_type_split(adata: AnnData, config: dict[str, Any]) -> str | None:
    """Return a real cell-type column. Never invent a replicate split."""
    key = None
    for section in ("preprocess", "perturb", "mixscape"):
        candidate = config.get(section, {}).get("cell_type_key")
        if candidate:
            key = str(candidate)
            break
    if key is None and "cell_type" in adata.obs:
        labels = set(adata.obs["cell_type"].astype(str).unique())
        if labels - {"ALL", "nan", "None", "", "<NA>"}:
            key = "cell_type"
    if key is None:
        return None
    if key not in adata.obs:
        raise ValueError(f"cell_type_key={key} 不在 obs 中。")
    labels = set(adata.obs[key].astype(str).unique())
    if len(labels - {"ALL", "nan", "None", "", "<NA>"}) < 2:
        return None
    return key


def _signature_umap(work: AnnData, seed: int, n_dims: int) -> tuple[np.ndarray | None, str]:
    """PCA and UMAP on the perturbation signature.

    Seurat Mixscape vignette runs ``ScaleData(do.scale=FALSE, do.center=TRUE)`` on
    PRTB then ``RunPCA`` / ``RunUMAP`` on dims 1:40; this pipeline uses the shared
    ``n_dims`` (default 30). Scanpy PCA zero-centers genes (equivalent to
    center-only). Do not gene-scale here: that would add a global Z-score on top
    of the control-neighbor signature and diverge from the vignette.
    """
    import scanpy as sc

    if "X_pert" not in work.layers or work.n_obs < 4:
        return None, "not_run_too_few_cells"
    vis = work.copy()
    vis.X = vis.layers["X_pert"].copy()
    try:
        n_comps = max(2, min(int(n_dims), vis.n_obs - 1, vis.n_vars - 1))
        sc.pp.pca(vis, n_comps=n_comps, zero_center=True, random_state=seed)
        n_neighbors = min(15, vis.n_obs - 1)
        sc.pp.neighbors(vis, n_neighbors=n_neighbors, metric="cosine", random_state=seed)
        sc.tl.umap(vis, random_state=seed)
    except Exception as exc:
        return None, f"signature_umap_failed:{exc}"
    return np.asarray(vis.obsm["X_umap"], dtype=float), "signature_space"


def _write_signature_rows(adata: AnnData, work: AnnData) -> None:
    parent_genes = adata.var_names.get_indexer(work.var_names)
    parent_cells = adata.obs_names.get_indexer(work.obs_names)
    if np.any(parent_genes < 0) or np.any(parent_cells < 0):
        raise ValueError("X_pert 无法对齐回原矩阵。")
    pert = work.layers["X_pert"]
    matrix = pert.tocoo() if sparse.issparse(pert) else sparse.coo_matrix(np.asarray(pert))
    rows = parent_cells[np.asarray(matrix.row)]
    columns = parent_genes[np.asarray(matrix.col)]
    full = sparse.coo_matrix((matrix.data, (rows, columns)), shape=(adata.n_obs, adata.n_vars))
    adata.layers["X_pert"] = full.tocsr()


def _normalize_mixscape(value: str, control: str, label: str) -> str:
    if value == control or value.endswith(f" {control}"):
        if control == ZERO_GUIDE_LABEL:
            return ZERO_GUIDE_LABEL
        return STATUS_NT
    if value.endswith(f" {label}") or value == label:
        return PERTURBED
    if value.endswith(" NP") or value == "NP":
        return NON_PERTURBED
    return AMBIGUOUS


def write_mixscape_figures(
    output_dir: Path,
    config: dict[str, Any],
    adata: AnnData,
) -> None:
    """Write Mixscape / Mixscape-UMAP PNGs at the end of the mixscape step."""
    from perturbseq.config import control_spec
    from perturbseq.plotting import (
        draw_mixscape_class_counts,
        draw_mixscape_class_global_umap,
        draw_mixscape_target_panels,
        draw_target_de_counts,
        draw_target_in_de_gene_counts,
        draw_umap,
        open_figure_store,
        retire_removed_mixscape_fraction,
    )

    output_dir = Path(output_dir)
    retire_removed_mixscape_fraction(output_dir)
    store = open_figure_store(output_dir, config, adata)
    label = perturbation_type(config)
    contrast = control_spec(config)["control_label"]
    if "X_umap" in adata.obsm and "mixscape_class" in adata.obs:
        draw_umap(store, adata, "mixscape_class", "umap_mixscape_class")
    else:
        store.skip("umap_mixscape_class", "missing_fields")
    if "mixscape_class" in adata.obs and "X_umap_mixscape" in adata.obsm:
        draw_umap(store, adata, "mixscape_class", "umap_mixscape", embedding="X_umap_mixscape")
    else:
        reason = (
            "mixscape_not_run" if "mixscape_class" not in adata.obs else "signature_umap_not_run"
        )
        store.skip("umap_mixscape", reason)
    draw_mixscape_class_global_umap(store, adata, label)
    draw_mixscape_class_counts(store, adata, label)
    plotting = config.get("plotting", {})
    draw_mixscape_target_panels(
        store,
        adata,
        label,
        control=contrast,
        seed=int(config["runtime"]["seed"]),
        max_genes=int(plotting.get("mixscape_heatmap_max_genes", 20)),
        max_cells=int(plotting.get("mixscape_heatmap_max_cells", 300)),
        pval_cutoff=float(config.get("mixscape", {}).get("pval_cutoff", 0.05)),
        lfc_line=float(config.get("mixscape", {}).get("logfc_threshold", 0.1)),
    )
    draw_target_de_counts(
        store,
        adata,
        control=contrast,
        min_de_genes=int(config.get("mixscape", {}).get("min_de_genes", 5)),
    )
    draw_target_in_de_gene_counts(
        store,
        adata,
        control=contrast,
        min_de_genes=int(config.get("mixscape", {}).get("min_de_genes", 5)),
    )
