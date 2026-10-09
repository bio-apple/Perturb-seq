"""Shared helpers for pipeline stage modules.

Owns: Mixscape mode resolution, embedding provenance, E-distance I/O helpers,
stage bookkeeping markers used inside stage bodies.

Does NOT own: orchestration, checkpoints, resume, or domain algorithms
(``perturbation``, ``statistics``, ``preprocessing``, …).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from anndata import AnnData

from perturbseq.perturbation import combine_distance_tables, run_secondary_distances
from perturbseq.plots import plot_edistance

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def _mixscape_subset(adata: AnnData, control: str) -> AnnData:
    if "mixscape_class_global" not in adata.obs:
        return adata
    keep = adata.obs["mixscape_class_global"].astype(str).isin([control, "KO"])
    return adata[keep].copy()


def _effective_mixscape_mode(config: PipelineConfig) -> str:
    """Resolve mixscape_mode with legacy --skip-mixscape / --force-mixscape flags."""
    mode = (config.mixscape_mode or "auto").strip().lower()
    if config.skip_mixscape or mode == "skip":
        return "skip"
    if config.force_mixscape or mode == "force":
        return "force"
    if mode == "subset":
        return "subset"
    # Explicit target list / top-n implies subset even under auto.
    if config.mixscape_targets or config.mixscape_top_n:
        return "subset"
    return "auto"


def _ko_label(gene: str, perturbation_type: str) -> str:
    return f"{gene} {perturbation_type}"


def _embedding_provenance(adata: AnnData, pca_source: str, config: PipelineConfig) -> dict[str, Any]:
    """Always-present embedding provenance keys for report.json → matrix_provenance."""
    if "highly_variable" in adata.var.columns:
        n_hvg = int(adata.var["highly_variable"].sum())
    else:
        n_hvg = int(config.n_top_genes)
    if "X_pca" in adata.obsm:
        n_pcs = int(adata.obsm["X_pca"].shape[1])
    else:
        n_pcs = int(config.n_pcs)
    return {"pca_source": str(pca_source), "n_hvg": n_hvg, "n_pcs": n_pcs}


def _write_edistance_outputs(
    edistances: pd.DataFrame,
    etest: pd.DataFrame | None,
    tables_dir: Path,
    figures_dir: Path,
    secondary: dict[str, pd.DataFrame] | None = None,
    secondary_status: list[dict[str, Any]] | None = None,
) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    edistances.to_csv(tables_dir / "edistance.csv")
    plot_edistance(edistances, figures_dir)
    if etest is not None:
        etest.to_csv(tables_dir / "etest.csv")
    secondary = secondary or {}
    for metric, frame in secondary.items():
        if frame is not None and not frame.empty:
            frame.to_csv(tables_dir / f"distance_{metric}.csv")
    combine_distance_tables(edistances, secondary).to_csv(tables_dir / "distances.csv")
    if secondary_status is not None:
        (tables_dir / "distance_metrics_status.json").write_text(
            json.dumps(secondary_status, indent=2, default=str)
        )


def _compute_secondary_distances(
    adata: AnnData,
    config: PipelineConfig,
    pca_source: str,
) -> tuple[dict[str, pd.DataFrame], list[dict[str, Any]]]:
    metrics = tuple(config.secondary_distance_metrics or ())
    if not metrics:
        return {}, []
    return run_secondary_distances(
        adata,
        groupby="gene_target",
        contrast=config.control,
        min_cells=config.min_cells_per_pert,
        metrics=metrics,
        pca_source=pca_source,
        n_jobs=config.n_jobs,
    )


def _pca_from_layer(adata: AnnData, layer: str, n_pcs: int) -> None:
    import scanpy as sc

    tmp = adata.copy()
    tmp.X = tmp.layers[layer]
    n_comps = min(n_pcs, max(2, tmp.n_obs - 1), max(2, tmp.n_vars - 1))
    sc.pp.pca(tmp, n_comps=n_comps)
    adata.obsm["X_pca"] = tmp.obsm["X_pca"]


def _mark_stage(report: dict, stage: str) -> None:
    stages = list(report.get("stages") or [])
    if stage not in stages:
        stages.append(stage)
    report["stages"] = stages
