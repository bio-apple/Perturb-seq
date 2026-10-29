"""Mixscale perturbation-strength scoring.

Independent of Mixscape. Input follows pertpy Mixscale: unscaled log1p ``.X``,
``perturbation_signature`` into ``layers['X_pert']``, then ``mixscale`` with
``layer='X_pert'`` and ``de_layer=None``. The signature is stored as
``layers['X_pert_mixscale']`` so it does not collide with Mixscape.

https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscale.html
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.assign import control_cell_mask, perturbation_analysis_mask
from perturbseq.config import DEFAULT_SIGNATURE_NDIMS, control_spec

# Mixscale tutorial CalcPerturbSig uses ndims=40 / num.neighbors=20; this pipeline
# keeps signature width aligned with preprocess PCA (DEFAULT_SIGNATURE_NDIMS=30).
MIXSCALE_LAYER = "X_pert_mixscale"
VIGNETTE_CONTROL_NEIGHBORS = 20


def mixscale_call_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Arguments for ``pertpy.tools.Mixscale.mixscale``.

    Reads ``mixscale.*`` only. Mixscape thresholds are not used.
    """
    mix = config.get("mixscale", {})
    if not isinstance(mix, dict):
        mix = {}
    return {
        "min_de_genes": int(mix.get("min_de_genes", 5)),
        "logfc_threshold": float(mix.get("logfc_threshold", 0.1)),
        "pval_cutoff": float(mix.get("pval_cutoff", 0.05)),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
    }


def _signature_n_dims(n_obs: int, n_vars: int, requested: int) -> int:
    """Cap ``n_dims`` by the representation width and the number of cells."""
    return max(1, min(int(requested), max(int(n_obs) - 1, 1), int(n_vars)))


def _neighbor_n_dims(work: AnnData, requested: int) -> int:
    """``n_dims`` for ``perturbation_signature``.

    Does not scale ``X`` and does not replace ``X_pca``. pertpy uses an existing
    ``X_pca`` for neighbors (``use_rep=None``) and otherwise computes PCA with
    its default parameters.
    """
    if "X_pca" in work.obsm:
        width = int(work.obsm["X_pca"].shape[1])
    else:
        width = int(work.n_vars)
    return _signature_n_dims(work.n_obs, width, requested)


def _require_unscaled_log1p(adata: AnnData) -> None:
    """pertpy subtracts ``expm1(X)``; scaled or raw counts are the wrong input."""
    if "log1p" not in adata.uns:
        raise ValueError("Mixscale 需要未 scale 的 log1p 矩阵（adata.X）。")
    matrix = adata.X
    if sparse.issparse(matrix):
        values = matrix.data
    else:
        values = np.asarray(matrix).ravel()
    if values.size and float(np.min(values)) < 0:
        raise ValueError("Mixscale 的 adata.X 含负值，不能把 scale 后的矩阵交给 perturbation_signature。")


def _stamp_control_label(obs: pd.DataFrame, mask: np.ndarray, column: str, label: str) -> None:
    """Write the control category onto ``pert_key`` before the pertpy call."""
    values = obs[column].astype(str).to_numpy(copy=True)
    values[np.asarray(mask)] = label
    obs[column] = values


def _cell_type_split(adata: AnnData, config: dict[str, Any]) -> str | None:
    """Return a real cell-type column. Does not read Mixscape settings."""
    key = None
    for section in ("preprocess", "perturb", "mixscale"):
        block = config.get(section, {})
        if not isinstance(block, dict):
            continue
        candidate = block.get("cell_type_key")
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


def _perturbation_label(config: dict[str, Any]) -> str:
    mode = config["experiment"]["crispr_mode"]
    return {"CRISPR_KO": "KO", "CRISPRi": "KD", "CRISPRa": "OE"}[mode]


def run_mixscale(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Score analysis-eligible targets on a Mixscale-owned signature.

    Does not read Mixscape classes, ``layers['X_pert']``, or a shared DE table.
    pertpy writes 0 for control cells and 1 for targets that do not have enough
    DE genes. Those official scores are kept.
    """
    from pertpy.tools import Mixscale

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
        adata.uns["mixscale_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        _stamp_control_label(adata.obs, control_cells, "target_gene", control)
    pert_mask = perturbation_analysis_mask(adata.obs, spec)
    _require_unscaled_log1p(adata)

    params = mixscale_call_parameters(config)
    scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    labels = adata.obs["target_gene"].astype(str)
    seed = int(config["runtime"]["seed"])
    mixscale_settings = config.get("mixscale", {})
    if not isinstance(mixscale_settings, dict):
        mixscale_settings = {}
    split_by_used = ""
    skipped: list[str] = []
    gene_lists: dict[str, list[str]] = {}

    targets = {
        str(gene) for gene in labels[pert_mask].unique() if str(gene) != control
    }
    if targets and bool(np.any(pert_mask)):
        work = adata[pert_mask].copy()
        requested_dims = int(mixscale_settings.get("n_dims", DEFAULT_SIGNATURE_NDIMS))
        n_dims = _neighbor_n_dims(work, requested_dims)
        n_nt = int((work.obs["target_gene"].astype(str) == control).sum())
        if n_nt < 1:
            raise ValueError("没有对照细胞，无法计算 Mixscale perturbation signature。")
        requested_neighbors = int(
            mixscale_settings.get("n_control_neighbors", VIGNETTE_CONTROL_NEIGHBORS)
        )
        n_neighbors = min(requested_neighbors, n_nt, max(work.n_obs - 1, 1))
        split_by = _cell_type_split(work, config)
        split_by_used = split_by or ""
        signature_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "n_neighbors": n_neighbors,
            "n_dims": n_dims,
        }
        if split_by:
            signature_args["split_by"] = split_by
        # https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscale.html
        Mixscale().perturbation_signature(work, **signature_args)
        mix_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "layer": "X_pert",
            "de_layer": None,
            "min_de_genes": int(params["min_de_genes"]),
            "max_de_genes": int(params["max_de_genes"]),
            "logfc_threshold": float(params["logfc_threshold"]),
            "pval_cutoff": float(params["pval_cutoff"]),
            "fine_mode": False,
            "random_state": seed,
            "scale": True,
        }
        if split_by:
            mix_args["split_by"] = split_by
        Mixscale().mixscale(work, **mix_args)
        if "X_pert" not in work.layers:
            raise ValueError("Mixscale perturbation_signature 没有写入 X_pert。")
        work.layers[MIXSCALE_LAYER] = work.layers.pop("X_pert")
        _write_mixscale_signature_rows(adata, work)
        written = pd.Series(work.obs["mixscale_score"].to_numpy(dtype=float), index=work.obs_names)
        scores.loc[written.index] = written.to_numpy()
        raw_genes = work.uns.get("mixscale_de_genes", {})
        if isinstance(raw_genes, dict):
            gene_lists = {
                str(gene): [str(name) for name in list(markers)]
                for gene, markers in raw_genes.items()
            }
        skipped = sorted(gene for gene, markers in gene_lists.items() if len(markers) == 0)
        if "mixscale" in work.uns:
            adata.uns["mixscale"] = work.uns["mixscale"]

    adata.obs["mixscale_score"] = scores.to_numpy(dtype=float)
    adata.uns["mixscale_de_genes"] = gene_lists
    adata.uns["mixscale_unscored_targets"] = skipped
    adata.uns["mixscale_parameters"] = {
        "layer": MIXSCALE_LAYER,
        "control": control,
        "split_by": split_by_used,
        "meaning": "perturbation_strength_not_editing_success",
        "scale": True,
        "fine_mode": False,
        "test_method": "wilcoxon",
        "independent_of_mixscape": True,
        "reference": "https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscale.html",
        "control_relative": (
            "perturbation_signature on unscaled log1p X (all genes) writes X_pert; "
            "mixscale(layer='X_pert', de_layer=None) scores that signature and "
            f"detects DE on X; stored as {MIXSCALE_LAYER}. "
            "Unscored target cells keep pertpy's fallback score 1; controls are 0."
        ),
        **params,
    }
    return adata


def _write_mixscale_signature_rows(adata: AnnData, work: AnnData) -> None:
    parent_genes = adata.var_names.get_indexer(work.var_names)
    parent_cells = adata.obs_names.get_indexer(work.obs_names)
    if np.any(parent_genes < 0) or np.any(parent_cells < 0):
        raise ValueError(f"{MIXSCALE_LAYER} 无法对齐回原矩阵。")
    pert = work.layers[MIXSCALE_LAYER]
    matrix = pert.tocoo() if sparse.issparse(pert) else sparse.coo_matrix(np.asarray(pert))
    rows = parent_cells[np.asarray(matrix.row)]
    columns = parent_genes[np.asarray(matrix.col)]
    full = sparse.coo_matrix((matrix.data, (rows, columns)), shape=(adata.n_obs, adata.n_vars))
    adata.layers[MIXSCALE_LAYER] = full.tocsr()


def write_mixscale_figures(
    output_dir: Path, config: dict[str, Any], adata: AnnData
) -> None:
    """Write Mixscale ridge / scatter / DE-count / DE-heatmap PNGs."""
    from perturbseq.plotting import (
        draw_mixscale_de_heatmaps,
        draw_mixscale_ridge,
        draw_mixscale_scatter,
        draw_target_de_counts,
        draw_target_in_de_gene_counts,
        open_figure_store,
    )

    output_dir = Path(output_dir)
    store = open_figure_store(output_dir, config, adata)
    contrast = control_spec(config)["control_label"]
    plotting = config.get("plotting", {})
    params = mixscale_call_parameters(config)
    prtb = _perturbation_label(config)
    draw_mixscale_ridge(store, adata, contrast, perturbation_label=prtb)
    draw_mixscale_scatter(store, adata, contrast, perturbation_label=prtb)
    draw_target_de_counts(
        store,
        adata,
        control=contrast,
        min_de_genes=int(params["min_de_genes"]),
        source="mixscale",
    )
    draw_target_in_de_gene_counts(
        store,
        adata,
        control=contrast,
        min_de_genes=int(params["min_de_genes"]),
        source="mixscale",
    )
    draw_mixscale_de_heatmaps(
        store,
        adata,
        control=contrast,
        seed=int(config["runtime"]["seed"]),
        max_genes=int(plotting.get("mixscape_heatmap_max_genes", 20)),
        max_cells=int(plotting.get("mixscape_heatmap_max_cells", 300)),
        perturbation_label=prtb,
    )
