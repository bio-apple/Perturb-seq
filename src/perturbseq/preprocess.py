"""Normalize one sample and build a visualization embedding."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from perturbseq.config import DEFAULT_EMBEDDING_NDIMS, mito_prefix


def preprocess(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Normalize all genes, then PCA / neighbors / Leiden / UMAP on top HVGs.

    Library-size normalization and ``log1p`` run on the full gene matrix so
    later effect tables still see every gene. Highly variable genes are then
    selected (default ``n_hvg=2000``, Seurat flavor). Mitochondrial genes are
    dropped from that set and replaced by the next non-mitochondrial genes.
    PCA is fit on the z-scored HVG slice only; neighbors, Leiden clustering and
    UMAP use that PCA. ``X`` stays sparse and log-normalized. The embedding is
    descriptive. Perturbation statistics do not use the UMAP coordinates.
    Harmony runs only when ``batch_key`` names a column with more than one
    level. A missing Harmony install stops that run instead of skipping the
    correction.
    """
    settings = config["preprocess"]
    seed = int(config["runtime"]["seed"])
    target_sum = float(settings["normalize_target_sum"])
    n_hvg = min(int(settings["n_hvg"]), adata.n_vars)
    _normalize_from_counts(adata, target_sum)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    _drop_mitochondrial_hvg(adata, config, n_hvg)
    usable = int(adata.var["highly_variable"].sum())
    if usable < 2:
        raise ValueError("高变基因数不足以计算 PCA。")
    n_comps = min(
        int(settings.get("n_pcs", DEFAULT_EMBEDDING_NDIMS)),
        usable - 1,
        adata.n_obs - 1,
    )
    if n_comps < 2:
        raise ValueError("细胞数或高变基因数不足以计算 PCA。")
    _pca_on_scaled_hvg(adata, n_comps, seed)
    batch_key = settings.get("batch_key")
    representation = "X_pca"
    if batch_key:
        if batch_key not in adata.obs:
            raise ValueError(f"preprocess.batch_key={batch_key} 不在 obs 中。")
        labels = adata.obs[batch_key].astype(str)
        n_levels = int(labels[~labels.isin(["nan", "None", "<NA>", ""])].nunique())
        if n_levels < 2:
            adata.uns["batch_correction"] = "not_run_fewer_than_two_batch_levels"
        else:
            try:
                sc.external.pp.harmony_integrate(adata, key=batch_key, random_state=seed)
            except ImportError as exc:
                raise RuntimeError(
                    "配置了批次校正，但当前环境没有 Harmony。安装 harmonypy 后重跑。"
                ) from exc
            representation = "X_pca_harmony"
            adata.uns["batch_correction"] = "harmony"
    n_neighbors = min(int(settings["n_neighbors"]), adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=representation, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=float(settings["leiden_resolution"]),
        flavor="igraph",
        random_state=seed,
        key_added="leiden",
    )
    adata.uns["preprocess"] = {
        "normalize_target_sum": target_sum,
        "n_hvg_requested": int(settings["n_hvg"]),
        "n_hvg_used": usable,
        "hvg_flavor": "seurat",
        "pca_on": "scaled_highly_variable_genes",
        "neighbors_umap_leiden_on": representation,
        "n_pcs": n_comps,
        "n_neighbors": n_neighbors,
        "leiden_resolution": float(settings["leiden_resolution"]),
        "note": (
            "Library-size normalization and log1p use all genes; "
            "PCA, neighbors, Leiden and UMAP use the top highly variable genes. "
            "HVG scale is gene-wise across all analysis cells, not relative to "
            "NT/zero_guide. Perturbation assessment uses Mixscape X_pert instead."
        ),
    }
    assign_cell_type(adata, config)
    return adata


def assign_cell_type(adata: AnnData, config: dict[str, Any]) -> str | None:
    """Use a real cell-type column, otherwise mark every cell ALL.

    Returns the obs column used for cell type, or None when unstratified.
    """
    annotation = config.get("input", {}).get("cell_annotations")
    if annotation:
        path = Path(str(annotation))
        if not path.exists():
            raise ValueError(f"细胞注释表不存在: {path}")
        table = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",")
        if not {"cell_barcode", "cell_type"}.issubset(table.columns):
            raise ValueError("细胞注释表需要 cell_barcode 和 cell_type。")
        mapped = table.drop_duplicates("cell_barcode").set_index("cell_barcode")["cell_type"]
        labels = adata.obs_names.astype(str).map(mapped)
        if labels.notna().any():
            adata.obs["cell_type"] = labels.astype("string")
            adata.uns["cell_type_scope"] = "annotated"
            adata.uns["cell_type_col"] = "cell_type"
            adata.uns["cell_type_note"] = "annotation table; unmapped cells stay missing"
            return "cell_type"
    key = None
    for section in ("preprocess", "perturb", "mixscape"):
        candidate = config.get(section, {}).get("cell_type_key")
        if candidate:
            key = str(candidate)
            break
    if key:
        if key not in adata.obs.columns:
            raise ValueError(f"cell_type_key={key} 不在 obs 中。")
        adata.obs["cell_type"] = adata.obs[key].astype(str)
        adata.uns["cell_type_scope"] = "annotated"
        adata.uns["cell_type_col"] = key
        adata.uns["cell_type_note"] = "annotated"
        return key
    if "cell_type" in adata.obs.columns:
        labels = set(adata.obs["cell_type"].astype(str).unique())
        if labels - {"ALL", "nan", "None", "", "<NA>"}:
            adata.uns["cell_type_scope"] = "annotated"
            adata.uns["cell_type_col"] = "cell_type"
            adata.uns["cell_type_note"] = "annotated"
            return "cell_type"
    adata.obs["cell_type"] = "ALL"
    adata.uns["cell_type_scope"] = "unstratified"
    adata.uns["cell_type_col"] = "cell_type"
    adata.uns["cell_type_note"] = "unstratified; no cell-type annotation"
    return None


def _normalize_from_counts(adata: AnnData, target_sum: float) -> None:
    if "counts" not in adata.layers:
        raise ValueError("预处理需要 layers['counts']。")
    if adata.layers["counts"] is adata.X:
        adata.layers["counts"] = adata.X.copy()
    else:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def _drop_mitochondrial_hvg(adata: AnnData, config: dict[str, Any], n_hvg: int) -> None:
    prefix = mito_prefix(config)
    if "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    mito = symbols.str.startswith(prefix).to_numpy()
    selected = adata.var["highly_variable"].to_numpy().copy()
    if not mito.any() or not np.any(selected & mito):
        return
    selected[mito] = False
    need = n_hvg - int(selected.sum())
    if need > 0 and "dispersions_norm" in adata.var.columns:
        scores = np.array(adata.var["dispersions_norm"].to_numpy(dtype=np.float64), copy=True)
        scores[~np.isfinite(scores)] = -np.inf
        scores[selected | mito] = -np.inf
        order = np.argsort(-scores, kind="mergesort")
        for index in order:
            if need <= 0 or not np.isfinite(scores[index]):
                break
            selected[index] = True
            need -= 1
    adata.var["highly_variable"] = selected


def _pca_on_scaled_hvg(adata: AnnData, n_comps: int, seed: int) -> None:
    """Scale only the HVG slice, then store PCA on the parent object.

    Scaling zero-centers and densifies. Doing that on the HVG columns avoids
    turning the full sparse log-normalized matrix into a dense array.
    """
    mask = adata.var["highly_variable"].to_numpy()
    work = adata[:, mask].copy()
    sc.pp.scale(work)
    sc.pp.pca(work, n_comps=n_comps, random_state=seed)
    adata.obsm["X_pca"] = work.obsm["X_pca"]
    components = np.asarray(work.varm["PCs"])
    loadings = np.zeros((adata.n_vars, n_comps), dtype=components.dtype)
    loadings[mask] = components
    adata.varm["PCs"] = loadings
    params = dict(work.uns["pca"]["params"])
    params["mask_var"] = "highly_variable"
    params["use_highly_variable"] = True
    params["scaled"] = True
    adata.uns["pca"] = {
        "params": params,
        "variance": work.uns["pca"]["variance"],
        "variance_ratio": work.uns["pca"]["variance_ratio"],
    }


def write_preprocess_figures(
    output_dir: Path, config: dict[str, Any], adata: AnnData
) -> None:
    """Write HVG / PCA / cluster / cell-cycle PNGs after preprocess."""
    from perturbseq.config import control_spec, default_markers
    from perturbseq.plotting import (
        draw_annotation_cell_cycle_panel,
        draw_dotplot_markers,
        draw_marker_detection_heatmap,
        draw_pca_elbow,
        draw_pca_scatter,
        draw_umap,
        draw_variable_features,
        open_figure_store,
        retire_removed_umaps,
    )

    output_dir = Path(output_dir)
    store = open_figure_store(output_dir, config, adata)
    retire_removed_umaps(output_dir)
    draw_variable_features(store, adata)
    draw_pca_elbow(store, adata)
    draw_pca_scatter(store, adata)
    spec = control_spec(config)
    control_label = str(spec["control_label"])
    if spec["temporary"]:
        control_values = (control_label,)
        control_title = "zero-guide"
    else:
        control_values = tuple(
            dict.fromkeys(
                value for value in ("nt", control_label) if value and value not in {"nan", "None"}
            )
        )
        control_title = "NT"
    draw_annotation_cell_cycle_panel(
        store, adata, control_values=control_values, control_title=control_title
    )
    markers = default_markers(config)
    draw_marker_detection_heatmap(store, adata, markers)
    draw_dotplot_markers(store, adata, markers)
    scope = str(adata.uns.get("cell_type_scope", "unstratified"))
    labels = set(adata.obs["cell_type"].astype(str).unique()) if "cell_type" in adata.obs else set()
    if "X_umap" in adata.obsm and scope == "annotated" and labels - {"ALL", "nan", "None", "<NA>"}:
        draw_umap(store, adata, "cell_type", "umap_cell_type")
    else:
        store.skip("umap_cell_type", "unstratified_cell_type")
