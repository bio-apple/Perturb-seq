"""Figures for one sample. UMAP is a display of the PCA neighborhood."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "PingFang SC"
plt.rcParams["axes.unicode_minus"] = False
import pandas as pd
from anndata import AnnData



def _save(fig: plt.Figure, path: Path) -> None:
    _save_many(fig, [path])


def _save_many(fig: plt.Figure, paths: list[Path]) -> None:
    fig.tight_layout()
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=120)
    plt.close(fig)


def qc_histograms(adata: AnnData, singlet: np.ndarray, figures: Path) -> None:
    obs = adata.obs.loc[singlet]
    panels = [
        ("n_genes", "nFeature_RNA"),
        ("total_counts", "nCount_RNA"),
        ("pct_mito", "percent.mt"),
        ("pct_ribo", "percent.ribo"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.8))
    for ax, (column, title) in zip(axes, panels):
        ax.hist(obs[column].to_numpy(), bins=80, color="#4C78A8")
        ax.set_title(title)
        ax.set_xlabel(title)
    axes[0].set_ylabel("Single-guide cells")
    _save(fig, figures / "qc_metrics.png")


def guide_charts(n_features: pd.DataFrame, frequency: pd.DataFrame, figures: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.2))
    bars = ax.bar(n_features["num_features"].astype(str), n_features["n_cells"], color="#4C78A8")
    ax.bar_label(bars, fmt="%d", padding=2, fontsize=8)
    ax.set_ylim(0, float(n_features["n_cells"].max()) * 1.14)
    ax.set_xlabel("Guides per cell")
    ax.set_ylabel("Passing cells")
    ax.set_title("Guides per cell")
    _save_many(fig, [figures / "guides_per_cell.png", figures / "guide_n_features.png"])

    ordered = frequency.sort_values(["n_cells", "guide"], ascending=[False, True])
    fig, ax = plt.subplots(figsize=(10, 4.2))
    if not ordered.empty:
        ax.bar(range(len(ordered)), ordered["n_cells"].to_numpy(), width=1.0, color="#4C78A8", linewidth=0)
        peak = ordered.iloc[0]
        if len(ordered) <= 40:
            ax.set_xticks(range(len(ordered)), ordered["guide"], rotation=90)
        else:
            ax.set_xticks([])
            ax.annotate(
                f"{peak['guide']}  {int(peak['n_cells']):,}",
                xy=(0, float(peak["n_cells"])),
                xytext=(10, 2),
                textcoords="offset points",
                va="bottom",
                ha="left",
                fontsize=8,
            )
        ax.set_yscale("log")
        ax.set_ylim(0.8, float(peak["n_cells"]) * 2.2)
    ax.set_xlabel("Guides ordered by cell count")
    ax.set_ylabel("Passing cells")
    ax.set_title("Cells per guide")
    _save(fig, figures / "guide_cell_counts.png")

    top = frequency.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bars = ax.barh(top["guide"], top["transcript_sum"], color="#F58518")
    ax.bar_label(bars, labels=[f"{int(value):,}" for value in top["transcript_sum"]], padding=3, fontsize=8)
    ax.set_xlim(0, float(top["transcript_sum"].max()) * 1.28)
    ax.set_xlabel("Transcript sum")
    ax.set_title("Top 10 guides")
    _save(fig, figures / "guide_top10.png")


def _umap_scatter(adata: AnnData, column: str, path: Path, legend: bool) -> None:
    plotted = adata.obs[column].astype(str)
    plotted = plotted[plotted != ""]
    if plotted.empty or "X_umap" not in adata.obsm:
        return
    xy = np.asarray(adata.obsm["X_umap"])
    positions = adata.obs_names.get_indexer(plotted.index)
    fig, ax = plt.subplots(figsize=(7, 6))
    labels = pd.Categorical(plotted)
    ax.scatter(
        xy[positions, 0],
        xy[positions, 1],
        c=labels.codes,
        s=2,
        cmap="tab20",
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(column)
    if legend and labels.categories.size <= 12:
        handles = [
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=plt.cm.tab20(i / max(labels.categories.size, 1)), label=str(cat))
            for i, cat in enumerate(labels.categories)
        ]
        ax.legend(handles=handles, frameon=False, markerscale=3, fontsize=8)
    _save(fig, path)


def umaps(adata: AnnData, figures: Path) -> None:
    _umap_scatter(adata, "cluster", figures / "umap_cluster.png", legend=True)
    _umap_scatter(adata, "perturbation", figures / "umap_perturbation.png", legend=False)
    _umap_scatter(adata, "icm_class", figures / "umap_icm_class.png", legend=True)
    _umap_scatter(adata, "mixscape_class", figures / "umap_mixscape.png", legend=True)


def composition_bars(composition: pd.DataFrame, figures: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(composition["cluster"].astype(str), composition["proportion"], color="#54A24B")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Proportion of single-guide cells")
    ax.set_title("Composition")
    _save(fig, figures / "composition.png")


def perturbation_composition_heatmap(table: pd.DataFrame, figures: Path) -> None:
    if table.empty:
        return
    values = table.drop(columns=["perturbation"]).to_numpy(dtype=float)
    fig_h = max(4.0, min(40.0, 0.18 * len(table)))
    fig, ax = plt.subplots(figsize=(8, fig_h))
    image = ax.imshow(values, aspect="auto", interpolation="none", cmap="viridis")
    ax.set_xticks(range(values.shape[1]), [column for column in table.columns if column != "perturbation"], rotation=90)
    if len(table) <= 40:
        ax.set_yticks(range(len(table)), table["perturbation"])
    else:
        ax.set_yticks([])
    ax.set_title("Perturbation composition")
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    _save(fig, figures / "perturbation_composition.png")

