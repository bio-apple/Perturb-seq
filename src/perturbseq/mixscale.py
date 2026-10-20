"""Mixscale perturbation-strength scoring on X_pert."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from pertpy.tools import Mixscale

from perturbseq.assign import control_cell_mask
from perturbseq.config import control_spec
from perturbseq.de_gate import _gene_lists, differential_gene_gate, stamp_control_label


def de_gate_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Thresholds the Mixscale tutorial passes to ``RunMixscale``.

    The tutorial calls ``min.de.genes = 5``, ``logfc.threshold = 0.2``,
    ``max.de.genes = 100``, and ``fine.mode = F``, with ``split.by = NULL``.
    ``pval.cutoff`` is omitted there; the function default is 0.05.
    Mixscape's own log fold-change is not read here.
    """
    mix = config.get("mixscale", {})
    return {
        "min_de_genes": int(mix.get("min_de_genes", 5)),
        "logfc_threshold": float(mix.get("logfc_threshold", 0.2)),
        "pval_cutoff": float(mix.get("pval_cutoff", 0.05)),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
        "fine_mode": False,
    }


def run_mixscale(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Score only targets that passed the Mixscale DE gate.

    Control cells are 0. Targets with fewer than ``min_de_genes`` keep a missing
    score rather than the software fallback of 1. No replicate split is invented.
    The Mixscape classification gate is left as Mixscape wrote it.
    These scores are not written into the descriptive effect tables.
    """
    spec = control_spec(config)
    control = spec["control_label"]
    control_cells = control_cell_mask(adata.obs, spec)
    if spec["temporary"] and int(control_cells.sum()) == 0:
        adata.uns["mixscale_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        stamp_control_label(adata.obs, control_cells, "target_gene", control)
    if "X_pert" not in adata.layers:
        raise ValueError("Mixscale 需要 layers['X_pert']。请先跑 mixscape 步骤。")
    params = de_gate_parameters(config)
    gate, gene_table = differential_gene_gate(
        adata,
        control=control,
        min_de_genes=int(params["min_de_genes"]),
        logfc_threshold=float(params["logfc_threshold"]),
        pval_cutoff=float(params["pval_cutoff"]),
        max_de_genes=int(params["max_de_genes"]),
    )
    adata.uns["mixscale_de_gate"] = gate
    adata.uns["mixscale_de_genes"] = gene_table
    gene_lists = _gene_lists(gene_table)
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    passed.discard(control)
    ready = set(passed)
    skipped = sorted(set(gate["target_gene"].astype(str)) - ready - {control})
    scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    labels = adata.obs["target_gene"].astype(str)
    min_de_genes = int(params["min_de_genes"])
    if ready:
        keep = labels.isin(ready | {control}).to_numpy()
        if "highly_variable" in adata.var and int(adata.var["highly_variable"].sum()) >= 2:
            work = adata[keep, adata.var["highly_variable"].to_numpy()].copy()
        else:
            work = adata[keep].copy()
        available = set(work.var_names.astype(str))
        supplied = {
            gene: [name for name in gene_lists.get(gene, []) if name in available][
                : int(params["max_de_genes"])
            ]
            for gene in ready
        }
        # scale=True: gene-scale X_pert, then z-score the projection versus control cells.
        mix_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "layer": "X_pert",
            "min_de_genes": min_de_genes,
            "max_de_genes": int(params["max_de_genes"]),
            "logfc_threshold": float(params["logfc_threshold"]),
            "pval_cutoff": float(params["pval_cutoff"]),
            "fine_mode": False,
            "random_state": int(config["runtime"]["seed"]),
            "scale": True,
        }
        if supplied and all(len(genes) >= min_de_genes for genes in supplied.values()):
            mix_args["de_genes_by_target"] = supplied
        Mixscale().mixscale(work, **mix_args)
        written = pd.Series(work.obs["mixscale_score"].to_numpy(dtype=float), index=work.obs_names)
        scored_targets = {str(gene) for gene in work.uns.get("mixscale", {})}
        used = work.uns.get("mixscale_de_genes", {})
        work_labels = labels.reindex(written.index).astype(str)
        for gene in ready:
            markers = list(used.get(gene, ()))
            if gene not in scored_targets or len(markers) < min_de_genes:
                written.loc[work_labels == gene] = np.nan
        scores.loc[written.index] = written.to_numpy()
    scores.loc[~labels.isin(ready | {control})] = np.nan
    adata.obs["mixscale_score"] = scores.to_numpy(dtype=float)
    adata.uns["mixscale_unscored_targets"] = skipped
    adata.uns["mixscale_parameters"] = {
        "layer": "X_pert",
        "control": control,
        "split_by": "",
        "meaning": "perturbation_strength_not_editing_success",
        "scale": True,
        "control_relative": (
            "scores on X_pert; projection z-scored to the control distribution "
            "(NT or temporary zero_guide); control cells get score 0"
        ),
        **params,
    }
    return adata
