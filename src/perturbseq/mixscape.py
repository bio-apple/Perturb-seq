"""Mixscape classification and control-relative perturbation signature."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.assign import STATUS_NT, STATUS_TARGET, control_cell_mask
from perturbseq.config import DEFAULT_SIGNATURE_NDIMS, ZERO_GUIDE_LABEL, control_spec

PERTURBED = "perturbed"
NON_PERTURBED = "non_perturbed"
AMBIGUOUS = "ambiguous"
NOT_EVALUATED = "not_evaluated"

# Seurat Mixscape vignette used CalcPerturbSig ndims=40; this pipeline keeps
# preprocess PCA and signature width aligned at DEFAULT_SIGNATURE_NDIMS (30).
# num.neighbors=20. RunMixscape omits logfc.threshold; function default is 0.25.
# https://satijalab.org/seurat/articles/mixscape_vignette
# https://satijalab.org/seurat/reference/runmixscape
# DEFAULT_SIGNATURE_NDIMS is imported from config (re-exported for default checks).
VIGNETTE_CONTROL_NEIGHBORS = 20


def perturbation_type(config: dict[str, Any]) -> str:
    mode = config["experiment"]["crispr_mode"]
    return {"CRISPR_KO": "KO", "CRISPRi": "KD", "CRISPRa": "OE"}.get(mode, "PT")


def mixscape_gate_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Mixscape classification thresholds from the Seurat vignette.

    ``logfc_threshold`` defaults to RunMixscape / MixscapeHeatmap's 0.25.
    Mixscale's 0.2 cutoff is not read here. ``fine_mode`` stays off.
    """
    mix = config.get("mixscape", {})
    perturb = config.get("perturb", {})
    if not isinstance(mix, dict):
        mix = {}
    if not isinstance(perturb, dict):
        perturb = {}
    return {
        "min_de_genes": int(mix.get("min_de_genes", perturb.get("mixscape_min_de_genes", 5))),
        "logfc_threshold": float(mix.get("logfc_threshold", perturb.get("mixscape_logfc", 0.25))),
        "pval_cutoff": float(mix.get("pval_cutoff", perturb.get("mixscape_pval", 0.05))),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
        "fine_mode": False,
    }


def signature_n_dims(n_obs: int, n_vars: int, requested: int = DEFAULT_SIGNATURE_NDIMS) -> int:
    """CalcPerturbSig ndims. Cap only when the matrix is smaller than the request."""
    return max(1, min(int(requested), int(n_obs) - 1, int(n_vars)))


def run_mixscape(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Gate on DE genes, then perturbation_signature and Mixscape for passing targets.

    No replicate column is invented. A real cell-type column with more than one
    label is the only optional split. Targets below ``min_de_genes`` stay
    ``skipped_no_de_genes`` and are not forced to NP. Posterior probability
    above 0.5 is the software's perturbed class. The signature UMAP is stored
    separately from the expression UMAP.
    """
    from pertpy.tools import Mixscape

    from perturbseq.de_gate import SKIPPED_NO_DE, differential_gene_gate, stamp_control_label

    spec = control_spec(config)
    control = spec["control_label"]
    control_cells = control_cell_mask(adata.obs, spec)
    if spec["temporary"] and int(control_cells.sum()) == 0:
        adata.obs["mixscape_class"] = NOT_EVALUATED
        adata.obs["mixscape_class_global"] = NOT_EVALUATED
        adata.obs["mixscape_class_p"] = np.nan
        adata.obs["mixscape_status"] = NOT_EVALUATED
        adata.uns["mixscape_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        stamp_control_label(adata.obs, control_cells, "target_gene", control)
        adata.uns["control_source"] = spec["control_source"]
        adata.uns["control_label"] = control
        adata.uns["control_note"] = spec["note"]
    mixscape_settings = config.get("mixscape", {})
    seed = int(config["runtime"]["seed"])
    if "highly_variable" not in adata.var:
        raise ValueError("Mixscape 需要先完成高变基因选择。")
    hvg = adata.var["highly_variable"].to_numpy()
    if int(np.sum(hvg)) < 2:
        raise ValueError("高变基因不足，无法计算 perturbation signature。")
    params = mixscape_gate_parameters(config)
    gate, gene_table = differential_gene_gate(
        adata,
        control=control,
        min_de_genes=int(params["min_de_genes"]),
        logfc_threshold=float(params["logfc_threshold"]),
        pval_cutoff=float(params["pval_cutoff"]),
        max_de_genes=int(params["max_de_genes"]),
    )
    adata.uns["mixscape_de_gate"] = gate
    adata.uns["mixscape_de_genes"] = gene_table
    labels = adata.obs["target_gene"].astype(str)
    adata.obs["mixscape_class"] = NOT_EVALUATED
    adata.obs["mixscape_class_global"] = NOT_EVALUATED
    adata.obs["mixscape_class_p"] = np.nan
    adata.obs["mixscape_status"] = NOT_EVALUATED
    skipped_genes = set(gate.loc[~gate["passed"].astype(bool), "target_gene"].astype(str))
    skipped = labels.isin(skipped_genes)
    adata.obs.loc[skipped, "mixscape_class"] = SKIPPED_NO_DE
    adata.obs.loc[skipped, "mixscape_status"] = SKIPPED_NO_DE
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    label = perturbation_type(config)
    if not passed:
        nt = labels == control
        adata.obs.loc[nt, "mixscape_class"] = control
        adata.obs.loc[nt, "mixscape_status"] = (
            ZERO_GUIDE_LABEL if control == ZERO_GUIDE_LABEL else STATUS_NT
        )
        adata.uns["mixscape_scope"] = "unstratified"
        adata.uns["mixscape_umap_status"] = "no_targets_passed_de_gate"
        adata.uns["mixscape_parameters"] = {
            "pert_key": "target_gene",
            "control": control,
            "split_by": "",
            "layer": "X_pert",
            "posterior_cutoff": 0.5,
            "perturbation_type": label,
            "random_state": seed,
            "fine_mode": False,
            **params,
        }
        return adata
    keep = labels.isin(passed | {control}).to_numpy()
    work = adata[keep, hvg].copy()
    requested_dims = int(mixscape_settings.get("n_dims", DEFAULT_SIGNATURE_NDIMS))
    n_dims = _prepare_signature_dims(work, requested_dims, seed)
    n_nt = int((work.obs["target_gene"].astype(str) == control).sum())
    if n_nt < 1:
        raise ValueError("没有 NT 细胞，无法计算 perturbation signature。")
    requested_neighbors = int(
        mixscape_settings.get("n_control_neighbors", VIGNETTE_CONTROL_NEIGHBORS)
    )
    n_neighbors = min(requested_neighbors, n_nt, max(work.n_obs - 1, 1))
    if n_neighbors < requested_neighbors:
        adata.uns["mixscape_neighbor_note"] = (
            f"n_control_neighbors reduced from {requested_neighbors} to {n_neighbors}"
        )
    split_by = _real_cell_type_split(work, config)
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
    # scale=True matches Seurat RunMixscape's scale.data on the PRTB assay for GMM.
    # Do not add a second gene-level Z-score versus pooled NT/zero_guide means on X.
    mix_args: dict[str, Any] = {
        "pert_key": "target_gene",
        "control": control,
        "layer": "X_pert",
        "de_layer": None,
        "min_de_genes": params["min_de_genes"],
        "logfc_threshold": params["logfc_threshold"],
        "pval_cutoff": params["pval_cutoff"],
        "iter_num": int(mixscape_settings.get("iter_num", 10)),
        "test_method": "wilcoxon",
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
        "test_method": "wilcoxon",
        "fine_mode": False,
        "scale": True,
        "control_relative": (
            "X_pert = log1p(mean of k control neighbors in linear space) - cell log1p; "
            "GMM runs on gene-scaled X_pert DE genes (pertpy scale=True); "
            "signature UMAP PCA zero-centers X_pert without a pooled-NT Z-score layer"
        ),
        **params,
    }
    return adata


def _real_cell_type_split(adata: AnnData, config: dict[str, Any]) -> str | None:
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


def _prepare_signature_dims(work: AnnData, requested: int, seed: int) -> int:
    """Width of the PCA used to find NT neighbors.

    Pertpy uses ``X_pca`` when the gene count exceeds Scanpy's PC threshold.
    A shorter preprocess PCA is rebuilt on a z-scored copy so the search can use
    the requested width (aligned with ``preprocess.n_pcs``, default 30).
    ``work.X`` stays log-normalized: the signature itself is the difference of
    unscaled expression. The cap applies only when the matrix itself is smaller.
    """
    import scanpy as sc

    n_dims = signature_n_dims(work.n_obs, work.n_vars, requested)
    if work.n_vars <= int(sc.settings.N_PCS):
        return n_dims
    pca_limit = min(work.n_obs - 1, work.n_vars - 1)
    target = min(n_dims, pca_limit)
    if target < 2:
        return n_dims
    have = int(work.obsm["X_pca"].shape[1]) if "X_pca" in work.obsm else 0
    if have < target:
        scaled = work.copy()
        sc.pp.scale(scaled)
        sc.pp.pca(scaled, n_comps=target, random_state=seed)
        work.obsm["X_pca"] = np.asarray(scaled.obsm["X_pca"])
        have = int(work.obsm["X_pca"].shape[1])
    return min(n_dims, have)


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
    if value == control:
        return STATUS_NT
    return AMBIGUOUS


def mixscape_fraction(adata: AnnData) -> pd.DataFrame:
    targeting = adata.obs[adata.obs["perturbation_status"] == STATUS_TARGET]
    rows = []
    for target, group in targeting.groupby("target_gene", sort=False):
        counts = group["mixscape_status"].value_counts()
        total = int(len(group))
        rows.append(
            {
                "target_gene": target,
                "n_cells": total,
                "n_perturbed": int(counts.get(PERTURBED, 0)),
                "n_non_perturbed": int(counts.get(NON_PERTURBED, 0)),
                "n_ambiguous": int(counts.get(AMBIGUOUS, 0)),
                "fraction_perturbed": float(counts.get(PERTURBED, 0) / total) if total else 0.0,
                "fraction_non_perturbed": float(counts.get(NON_PERTURBED, 0) / total)
                if total
                else 0.0,
                "fraction_ambiguous": float(counts.get(AMBIGUOUS, 0) / total) if total else 0.0,
            }
        )
    return pd.DataFrame(rows)
