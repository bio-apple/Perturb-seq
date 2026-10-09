from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

sc.settings.figdir = Path(".")


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_qc(adata: AnnData, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    for ax, key in zip(axes, ("n_counts", "n_genes", "pct_counts_mt"), strict=True):
        ax.hist(adata.obs[key], bins=40, color="#4C72B0")
        ax.set_title(key)
        ax.set_xlabel(key)
        ax.set_ylabel("cells")
    _save(fig, out_dir / "qc_histograms.png")


def plot_guide_composition(adata: AnnData, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    adata.obs["num_features"].value_counts().sort_index().plot(kind="bar", ax=axes[0], color="#4C72B0")
    axes[0].set_title("guides per cell")
    axes[0].set_xlabel("num_features")
    adata.obs["gene_target"].value_counts().head(20).plot(kind="bar", ax=axes[1], color="#DD8452")
    axes[1].set_title("top gene targets")
    fig.tight_layout()
    _save(fig, out_dir / "guide_composition.png")


def plot_umap(adata: AnnData, out_dir: Path, color: list[str] | None = None) -> None:
    keys = [key for key in (color or ["leiden", "gene_target", "perturbation"]) if key in adata.obs]
    if "X_umap" not in adata.obsm or not keys:
        return
    sc.pl.umap(adata, color=keys, show=False, wspace=0.4)
    fig = plt.gcf()
    _save(fig, out_dir / "umap.png")


def plot_edistance(edistances: pd.DataFrame, out_dir: Path) -> None:
    if edistances.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    top = edistances["edistance"].head(30).iloc[::-1]
    top.plot(kind="barh", ax=ax, color="#4C72B0")
    ax.set_xlabel("E-distance to NT")
    ax.set_title("Perturbation effect size")
    fig.tight_layout()
    _save(fig, out_dir / "edistance.png")


def plot_volcano(table: pd.DataFrame, out_dir: Path, name: str) -> None:
    if table.empty:
        return
    lfc_col = next((c for c in ("log2FoldChange", "logfoldchanges") if c in table.columns), None)
    p_col = next((c for c in ("padj", "pvals_adj") if c in table.columns), None)
    if lfc_col is None or p_col is None:
        return
    plot_df = table.dropna(subset=[lfc_col, p_col]).copy()
    pvals = pd.to_numeric(plot_df[p_col], errors="coerce").clip(lower=1e-300)
    plot_df["neglog10p"] = -np.log10(pvals)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(plot_df[lfc_col], plot_df["neglog10p"], s=8, c="#4C72B0", alpha=0.6)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel(lfc_col)
    ax.set_ylabel("-log10(padj)")
    ax.set_title(name)
    _save(fig, out_dir / f"volcano_{name}.png")
