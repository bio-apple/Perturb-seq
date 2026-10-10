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


def _label_bar_counts(ax) -> None:
    """Annotate each bar with its integer count; rotate when many bars crowd the axis."""
    n_bars = 0
    for container in ax.containers:
        values = getattr(container, "datavalues", None)
        if values is None:
            continue
        n_bars = max(n_bars, len(values))
        labels = [f"{int(round(float(v)))}" for v in values]
        crowded = n_bars > 8
        ax.bar_label(
            container,
            labels=labels,
            fontsize=6.5 if crowded else 8,
            rotation=90 if crowded else 0,
            padding=2,
        )
    if n_bars and ax.containers:
        ymin, ymax = ax.get_ylim()
        if ymax > ymin:
            ax.set_ylim(ymin, ymax * (1.22 if n_bars > 8 else 1.12))


def plot_qc(adata: AnnData, out_dir: Path) -> None:
    """QC histograms, violins, and metric scatters (必做 cell QC figures)."""
    keys = ("n_counts", "n_genes", "pct_counts_mt")
    present = [k for k in keys if k in adata.obs]
    if not present:
        return

    fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 3.5))
    if len(present) == 1:
        axes = [axes]
    for ax, key in zip(axes, present, strict=True):
        ax.hist(adata.obs[key], bins=40, color="#4C72B0")
        ax.set_title(key)
        ax.set_xlabel(key)
        ax.set_ylabel("cells")
    _save(fig, out_dir / "qc_histograms.png")

    fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 3.8))
    if len(present) == 1:
        axes = [axes]
    for ax, key in zip(axes, present, strict=True):
        ax.violinplot(adata.obs[key].astype(float).to_numpy(), showmeans=True, showextrema=True)
        ax.set_title(key)
        ax.set_xticks([1])
        ax.set_xticklabels([key])
        ax.set_ylabel(key)
    fig.suptitle("QC violin", fontsize=11)
    fig.tight_layout()
    _save(fig, out_dir / "qc_violin.png")

    scatter_pairs = [
        ("n_counts", "n_genes", "qc_scatter_counts_genes.png"),
        ("n_counts", "pct_counts_mt", "qc_scatter_counts_mt.png"),
    ]
    for x_key, y_key, filename in scatter_pairs:
        if x_key not in adata.obs or y_key not in adata.obs:
            continue
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(
            adata.obs[x_key].astype(float),
            adata.obs[y_key].astype(float),
            s=6,
            alpha=0.35,
            c="#4C72B0",
            linewidths=0,
        )
        ax.set_xlabel(x_key)
        ax.set_ylabel(y_key)
        ax.set_title(f"{y_key} vs {x_key}")
        _save(fig, out_dir / filename)


def plot_qc_cell_counts(counts: pd.DataFrame | dict, out_dir: Path) -> None:
    """Before/after cell counts across QC / singlet (and optional later) stages."""
    if isinstance(counts, dict):
        table = pd.DataFrame(
            [{"stage": k, "n_cells": int(v)} for k, v in counts.items() if v is not None]
        )
    else:
        table = counts.copy()
    if table.empty or "stage" not in table.columns or "n_cells" not in table.columns:
        return
    table = table.sort_values("stage", kind="stable")
    fig, ax = plt.subplots(figsize=(max(5, 0.9 * len(table)), 3.8))
    table.set_index("stage")["n_cells"].plot(kind="bar", ax=ax, color="#4C72B0")
    ax.set_ylabel("n_cells")
    ax.set_title("Cell counts before / after filters")
    _label_bar_counts(ax)
    fig.tight_layout()
    _save(fig, out_dir / "qc_cell_counts.png")


def plot_guide_composition(adata: AnnData, out_dir: Path) -> None:
    """Guides-per-cell distribution + top gene targets (gRNA QC)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    adata.obs["num_features"].value_counts().sort_index().plot(kind="bar", ax=axes[0], color="#4C72B0")
    axes[0].set_title("guides per cell")
    axes[0].set_xlabel("num_features")
    _label_bar_counts(axes[0])
    adata.obs["gene_target"].value_counts().head(20).plot(kind="bar", ax=axes[1], color="#DD8452")
    axes[1].set_title("top gene targets")
    _label_bar_counts(axes[1])
    fig.tight_layout()
    _save(fig, out_dir / "guide_composition.png")


def plot_guide_cell_counts(
    adata: AnnData,
    out_dir: Path,
    *,
    top_n: int = 40,
    filename: str = "guide_cell_counts.png",
) -> None:
    """Barplot of cells per guide_id (必做 gRNA QC)."""
    if "guide_id" not in adata.obs:
        return
    counts = adata.obs["guide_id"].astype(str).value_counts().head(top_n)
    if counts.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.28 * len(counts)), 4))
    counts.plot(kind="bar", ax=ax, color="#4C72B0")
    ax.set_title(f"Cells per guide (top {len(counts)})")
    ax.set_xlabel("guide_id")
    ax.set_ylabel("n_cells")
    _label_bar_counts(ax)
    fig.tight_layout()
    _save(fig, out_dir / filename)


def plot_umap(adata: AnnData, out_dir: Path, color: list[str] | None = None, filename: str = "umap.png") -> None:
    keys = [key for key in (color or ["leiden", "gene_target", "perturbation"]) if key in adata.obs]
    if "X_umap" not in adata.obsm or not keys:
        return
    sc.pl.umap(adata, color=keys, show=False, wspace=0.4)
    fig = plt.gcf()
    fig.suptitle("UMAP (visualization only; not perturbation-effect evidence)", fontsize=10, y=1.02)
    _save(fig, out_dir / filename)


def plot_perturbation_umap(
    adata: AnnData,
    out_dir: Path,
    *,
    control: str = "NT",
    filename: str = "umap_perturbation.png",
) -> None:
    """UMAP colored by guide, target gene, and NTC vs targeting (必做 perturbation map)."""
    if "X_umap" not in adata.obsm:
        return
    obs = adata.obs
    if "is_ntc" not in obs.columns:
        if "gene_target" in obs.columns:
            adata.obs["is_ntc"] = obs["gene_target"].astype(str).eq(str(control))
        elif "perturbation" in obs.columns:
            adata.obs["is_ntc"] = obs["perturbation"].astype(str).eq(str(control))
    # Cap guide categories for readable legends on large libraries.
    if "guide_id" in adata.obs and adata.obs["guide_id"].nunique() > 30:
        top = set(adata.obs["guide_id"].astype(str).value_counts().head(29).index)
        adata.obs["guide_id_plot"] = adata.obs["guide_id"].astype(str).where(
            adata.obs["guide_id"].astype(str).isin(top), other="other"
        )
        guide_key = "guide_id_plot"
    else:
        guide_key = "guide_id"
    keys = [k for k in (guide_key, "gene_target", "is_ntc") if k in adata.obs]
    if not keys:
        return
    plot_umap(adata, out_dir, color=keys, filename=filename)


def plot_composition_by_cluster(
    table: pd.DataFrame,
    out_dir: Path,
    *,
    filename: str = "composition_by_perturbation.png",
    top_perturbations: int = 20,
) -> None:
    """Stacked proportions of perturbations within Leiden clusters."""
    if table is None or table.empty:
        return
    needed = {"cluster", "perturbation", "n_cells"}
    if not needed.issubset(table.columns):
        return
    work = table.copy()
    work["cluster"] = work["cluster"].astype(str)
    work["perturbation"] = work["perturbation"].astype(str)
    top = (
        work.groupby("perturbation", observed=True)["n_cells"].sum().sort_values(ascending=False).head(top_perturbations).index
    )
    work["perturbation_plot"] = work["perturbation"].where(work["perturbation"].isin(top), other="other")
    pivot = (
        work.groupby(["cluster", "perturbation_plot"], observed=True)["n_cells"]
        .sum()
        .unstack(fill_value=0)
    )
    if pivot.empty:
        return
    frac = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(frac)), 4.2))
    frac.plot(kind="bar", stacked=True, ax=ax, width=0.85)
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel("cluster")
    ax.set_title("Perturbation composition within clusters (descriptive)")
    ax.legend(title="perturbation", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    fig.tight_layout()
    _save(fig, out_dir / filename)


def plot_de_heatmap(
    tables: list[tuple[str, pd.DataFrame]],
    out_dir: Path,
    *,
    top_n: int = 20,
    filename: str = "de_heatmap.png",
) -> bool:
    """Heatmap of top |LFC| genes across DE contrasts (when suitable control / DE exists)."""
    if not tables:
        return False
    gene_scores: dict[str, float] = {}
    lfc_maps: dict[str, dict[str, float]] = {}
    for name, table in tables:
        if table is None or table.empty:
            continue
        lfc_col = next((c for c in ("log2FoldChange", "logfoldchanges") if c in table.columns), None)
        gene_col = next((c for c in ("gene", "names", "gene_symbol") if c in table.columns), None)
        if lfc_col is None:
            continue
        work = table.copy()
        if gene_col is None:
            work = work.reset_index()
            gene_col = work.columns[0]
        work[lfc_col] = pd.to_numeric(work[lfc_col], errors="coerce")
        work = work.dropna(subset=[lfc_col, gene_col])
        if work.empty:
            continue
        work["gene"] = work[gene_col].astype(str)
        best = work.reindex(work[lfc_col].abs().sort_values(ascending=False).index).drop_duplicates("gene")
        lfc_maps[str(name)] = dict(zip(best["gene"], best[lfc_col], strict=False))
        for gene, lfc in zip(best["gene"].head(top_n), best[lfc_col].head(top_n), strict=False):
            gene_scores[gene] = max(gene_scores.get(gene, 0.0), abs(float(lfc)))
    if not lfc_maps or not gene_scores:
        return False
    genes = [g for g, _ in sorted(gene_scores.items(), key=lambda kv: -kv[1])[:top_n]]
    mat = pd.DataFrame(
        {pert: [lfc_maps[pert].get(g, np.nan) for g in genes] for pert in lfc_maps},
        index=genes,
    )
    if mat.isna().all().all():
        return False
    fig, ax = plt.subplots(figsize=(max(4.5, 0.55 * mat.shape[1] + 2), max(4, 0.28 * mat.shape[0] + 1.5)))
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(list(mat.columns), rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(list(mat.index), fontsize=7)
    ax.set_title("DE log2FC heatmap (top |LFC| genes)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log2FC")
    fig.tight_layout()
    _save(fig, out_dir / filename)
    return True


def plot_mixscape_lda(
    adata: AnnData,
    out_dir: Path,
    *,
    control: str = "NT",
    perturbation_type: str = "KO",
    filename: str = "lda_umap.png",
) -> bool:
    """Best-effort MixscapeLDA viz (Seurat vignette analogue). Requires pertpy + ``uns['mixscape_lda']``."""
    if "mixscape_lda" not in adata.uns or "mixscape_class" not in adata.obs:
        return False
    try:
        from perturbseq._deps import require_pertpy

        pt = require_pertpy()
        fig = pt.tl.Mixscape().plot_lda(
            adata,
            control=control,
            perturbation_type=perturbation_type,
            return_fig=True,
        )
    except Exception:  # noqa: BLE001
        return False
    if fig is None:
        fig = plt.gcf()
    _save(fig, out_dir / filename)
    return True


def plot_cell_annotation(adata: AnnData, out_dir: Path) -> None:
    keys = [k for k in ("leiden", "phase", "cell_state", "cluster_annotation") if k in adata.obs]
    if keys:
        plot_umap(adata, out_dir, color=keys, filename="umap_cell_annotation.png")
    if "phase" in adata.obs and "cell_state" in adata.obs:
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
        adata.obs["phase"].astype(str).value_counts().plot(kind="bar", ax=axes[0], color="#4C72B0")
        axes[0].set_title("cell cycle phase")
        _label_bar_counts(axes[0])
        adata.obs["cell_state"].astype(str).value_counts().plot(kind="bar", ax=axes[1], color="#DD8452")
        axes[1].set_title("cell state")
        _label_bar_counts(axes[1])
        fig.tight_layout()
        _save(fig, out_dir / "cell_annotation_counts.png")


def plot_edistance(edistances: pd.DataFrame, out_dir: Path, filename: str = "edistance.png") -> None:
    if edistances.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    top = edistances["edistance"].head(30).iloc[::-1]
    top.plot(kind="barh", ax=ax, color="#4C72B0")
    ax.set_xlabel("E-distance to NT")
    source = ""
    if "pca_source" in edistances.columns and len(edistances):
        source = f" on X_pca ({edistances['pca_source'].iloc[0]})"
    ax.set_title(f"E-distance effect size{source}")
    fig.tight_layout()
    _save(fig, out_dir / filename)


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


def plot_guide_qc(guide_df: pd.DataFrame, consistency_df: pd.DataFrame, out_dir: Path) -> None:
    """Guide-level QC: cells/UMI vs target effect, and per-gene guide consistency."""
    if guide_df is None or guide_df.empty:
        return
    non_ctrl = guide_df.loc[~guide_df["is_control"]].copy() if "is_control" in guide_df.columns else guide_df
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    if not non_ctrl.empty and "n_cells" in non_ctrl.columns:
        axes[0].hist(non_ctrl["n_cells"], bins=40, color="#4C72B0")
        axes[0].set_title("cells per guide")
        axes[0].set_xlabel("n_cells")
    if (
        not non_ctrl.empty
        and "median_guide_umi" in non_ctrl.columns
        and "target_expr_log2fc_vs_control" in non_ctrl.columns
    ):
        axes[1].scatter(
            non_ctrl["median_guide_umi"],
            non_ctrl["target_expr_log2fc_vs_control"],
            s=10,
            alpha=0.5,
            c="#DD8452",
        )
        axes[1].axhline(0, color="grey", lw=0.8)
        axes[1].set_xlabel("median guide UMI")
        axes[1].set_ylabel("target log2FC vs control")
        axes[1].set_title("assignment vs target effect")
    if consistency_df is not None and not consistency_df.empty and "guides_consistent" in consistency_df.columns:
        counts = consistency_df["guides_consistent"].map({True: "consistent", False: "inconsistent"}).fillna("unknown")
        counts.value_counts().reindex(["consistent", "inconsistent", "unknown"]).fillna(0).plot(
            kind="bar", ax=axes[2], color="#55A868"
        )
        axes[2].set_title("guide consistency per gene")
        _label_bar_counts(axes[2])
    fig.tight_layout()
    _save(fig, out_dir / "guide_qc.png")
