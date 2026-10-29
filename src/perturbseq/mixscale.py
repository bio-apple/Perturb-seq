"""Mixscale perturbation-strength scoring.

Independent of the Mixscape classification step. Follows the Mixscale tutorial:
CalcPerturbSig-equivalent signature, then RunMixscale on that assay.

https://satijalab.github.io/Mixscale/
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
from perturbseq.de_gate import differential_gene_gate, stamp_control_label

# Mixscale tutorial CalcPerturbSig uses ndims=40 / num.neighbors=20; this pipeline
# keeps signature width aligned with preprocess PCA (DEFAULT_SIGNATURE_NDIMS=30).
MIXSCALE_LAYER = "X_pert_mixscale"
VIGNETTE_CONTROL_NEIGHBORS = 20


def de_gate_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Thresholds passed to Mixscale DE gate / ``RunMixscale``.

    Project default ``logfc_threshold`` is 0.1 (|logFC|≥0.1); the Mixscale
    tutorial originally used 0.2. Also ``min.de.genes = 5``,
    ``max.de.genes = 100``, ``fine.mode = F``, ``split.by = NULL``.
    ``pval.cutoff`` default is 0.05. Mixscape's log fold-change is not read.
    """
    mix = config.get("mixscale", {})
    return {
        "min_de_genes": int(mix.get("min_de_genes", 5)),
        "logfc_threshold": float(mix.get("logfc_threshold", 0.1)),
        "pval_cutoff": float(mix.get("pval_cutoff", 0.05)),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
        "fine_mode": False,
    }


def run_mixscale(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Score targets that passed the Mixscale DE gate on a Mixscale-owned signature.

    Computes ``layers['X_pert_mixscale']`` here (tutorial CalcPerturbSig), then
    RunMixscale. Does not read Mixscape classes or ``layers['X_pert']``.
    """
    from pertpy.tools import Mixscale

    from perturbseq.mixscape import neighbor_n_dims, real_cell_type_split, require_unscaled_log1p

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
        stamp_control_label(adata.obs, control_cells, "target_gene", control)
    # Embedding may include all QC-pass cells; scoring stays on analysis_eligible ∪ control.
    pert_mask = perturbation_analysis_mask(adata.obs, spec)
    require_unscaled_log1p(adata, "Mixscale")

    params = de_gate_parameters(config)
    gate_adata = adata[pert_mask] if not bool(np.all(pert_mask)) else adata
    gate, gene_table = differential_gene_gate(
        gate_adata,
        control=control,
        min_de_genes=int(params["min_de_genes"]),
        logfc_threshold=float(params["logfc_threshold"]),
        pval_cutoff=float(params["pval_cutoff"]),
        max_de_genes=int(params["max_de_genes"]),
    )
    adata.uns["mixscale_de_gate"] = gate
    adata.uns["mixscale_de_genes"] = gene_table
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    passed.discard(control)
    ready = set(passed)
    skipped = sorted(set(gate["target_gene"].astype(str)) - ready - {control})
    scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    labels = adata.obs["target_gene"].astype(str)
    min_de_genes = int(params["min_de_genes"])
    seed = int(config["runtime"]["seed"])
    mixscale_settings = config.get("mixscale", {})
    split_by_used = ""

    if ready:
        keep = labels.isin(ready | {control}).to_numpy() & pert_mask
        # Same input contract as Mixscape: unscaled log1p X, all genes.
        work = adata[keep].copy()
        requested_dims = int(
            mixscale_settings.get(
                "n_dims",
                config.get("mixscape", {}).get("n_dims", DEFAULT_SIGNATURE_NDIMS),
            )
        )
        n_dims = neighbor_n_dims(work, requested_dims)
        n_nt = int((work.obs["target_gene"].astype(str) == control).sum())
        if n_nt < 1:
            raise ValueError("没有对照细胞，无法计算 Mixscale perturbation signature。")
        requested_neighbors = int(
            mixscale_settings.get(
                "n_control_neighbors",
                config.get("mixscape", {}).get("n_control_neighbors", VIGNETTE_CONTROL_NEIGHBORS),
            )
        )
        n_neighbors = min(requested_neighbors, n_nt, max(work.n_obs - 1, 1))
        split_by = real_cell_type_split(work, config)
        split_by_used = split_by or ""
        signature_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "n_neighbors": n_neighbors,
            "n_dims": n_dims,
        }
        if split_by:
            signature_args["split_by"] = split_by
        # perturbation_signature writes X_pert from unscaled X. mixscale scores
        # that layer and finds DE genes on X (de_layer=None). de_genes_by_target
        # is omitted so the DE test is not skipped.
        # https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscale.html
        Mixscale().perturbation_signature(work, **signature_args)
        mix_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "layer": "X_pert",
            "de_layer": None,
            "min_de_genes": min_de_genes,
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
        scored_targets = {str(gene) for gene in work.uns.get("mixscale", {})}
        used = work.uns.get("mixscale_de_genes", {})
        work_labels = labels.reindex(written.index).astype(str)
        for gene in ready:
            markers = list(used.get(gene, ()))
            if gene not in scored_targets or len(markers) < min_de_genes:
                written.loc[work_labels == gene] = np.nan
        scores.loc[written.index] = written.to_numpy()
        if "mixscale" in work.uns:
            adata.uns["mixscale"] = work.uns["mixscale"]

    scores.loc[~labels.isin(ready | {control})] = np.nan
    adata.obs["mixscale_score"] = scores.to_numpy(dtype=float)
    adata.uns["mixscale_unscored_targets"] = skipped
    adata.uns["mixscale_parameters"] = {
        "layer": MIXSCALE_LAYER,
        "control": control,
        "split_by": split_by_used,
        "meaning": "perturbation_strength_not_editing_success",
        "scale": True,
        "test_method": "wilcoxon",
        "independent_of_mixscape": True,
        "reference": "https://pertpy.scverse.org/en/latest/api/tools/pertpy.tools.Mixscale.html",
        "control_relative": (
            "perturbation_signature on unscaled log1p X (all genes) writes X_pert; "
            "mixscale(layer='X_pert', de_layer=None) scores that signature and "
            f"detects DE on X; stored as {MIXSCALE_LAYER}"
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
    from perturbseq.mixscape import perturbation_type
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
    params = de_gate_parameters(config)
    prtb = perturbation_type(config)
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
