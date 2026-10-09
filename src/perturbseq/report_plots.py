"""Biologist-facing report figures: effect summary, guide consistency, target validation.

PNG files under ``figures/`` (not base64). Skip with explicit notes when inputs missing.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy import sparse

DEFAULT_REPORT_PLOT_TOP_N = 15
DEFAULT_CORR_N_GENES = 200
DEFAULT_MIN_GUIDE_CELLS = 10

EFFECT_SUMMARY_NAME = "perturbation_effect_summary.png"
GUIDE_CONSISTENCY_PREFIX = "guide_consistency_"
TARGET_VALIDATION_PREFIX = "target_validation_"


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _safe_stem(name: str) -> str:
    cleaned = re.sub(r"[^\w.\-]+", "_", str(name).strip())
    return cleaned or "unnamed"


def _to_dense_block(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray(), dtype=np.float64)
    return np.asarray(matrix, dtype=np.float64)


def _guide_col(adata: AnnData) -> str | None:
    for col in ("guide_id", "feature_call", "guide"):
        if col in adata.obs.columns:
            return col
    return None


def _select_genes_for_corr(adata: AnnData, n_genes: int) -> list[str]:
    if "highly_variable" in adata.var.columns:
        hv = adata.var_names[adata.var["highly_variable"].astype(bool)].tolist()
        if hv:
            return hv[:n_genes]
    # Variance on a cell subsample for speed
    n_obs = adata.n_obs
    rng = np.random.default_rng(0)
    idx = rng.choice(n_obs, size=min(n_obs, 500), replace=False) if n_obs else np.array([], dtype=int)
    if len(idx) == 0:
        return list(adata.var_names[:n_genes])
    block = _to_dense_block(adata.X[idx, :])
    var = np.nanvar(block, axis=0)
    order = np.argsort(-var)
    return [str(adata.var_names[i]) for i in order[:n_genes]]


def _mean_profile(adata: AnnData, cell_mask: np.ndarray, gene_idx: np.ndarray) -> np.ndarray:
    if not cell_mask.any() or len(gene_idx) == 0:
        return np.full(len(gene_idx), np.nan, dtype=np.float64)
    block = _to_dense_block(adata.X[cell_mask][:, gene_idx])
    return block.mean(axis=0)


def plot_perturbation_effect_summary(
    edistance: pd.DataFrame | None,
    etest: pd.DataFrame | None,
    out_dir: Path,
    *,
    filename: str = EFFECT_SUMMARY_NAME,
) -> Path | None:
    """Volcano-like: x=E-distance, y=-log10(p), size=n_cells."""
    if edistance is None or edistance.empty:
        return None
    ed = edistance.copy()
    if "gene_target" not in ed.columns:
        ed = ed.rename(columns={ed.columns[0]: "gene_target"})
    if "edistance" not in ed.columns:
        return None

    et = etest.copy() if etest is not None and not etest.empty else None
    if et is not None:
        et = et.reset_index()
        # index name may be gene_target / Unnamed / first col
        if "gene_target" not in et.columns:
            et = et.rename(columns={et.columns[0]: "gene_target"})
        keep = [c for c in ("gene_target", "pvalue", "pvalue_adj", "significant_adj", "significant_adj_reported") if c in et.columns]
        ed = ed.merge(et[keep], on="gene_target", how="left")

    plot_df = ed.dropna(subset=["edistance"]).copy()
    if plot_df.empty:
        return None

    pcol = next((c for c in ("pvalue_adj", "pvalue") if c in plot_df.columns), None)
    if pcol is None:
        plot_df["neglog10p"] = 0.0
        ylabel = "-log10(p) (unavailable)"
    else:
        pvals = pd.to_numeric(plot_df[pcol], errors="coerce").clip(lower=1e-300)
        plot_df["neglog10p"] = -np.log10(pvals.fillna(1.0))
        ylabel = f"-log10({pcol})"

    n_cells = (
        pd.to_numeric(plot_df["n_cells"], errors="coerce")
        if "n_cells" in plot_df.columns
        else pd.Series(np.full(len(plot_df), 10.0), index=plot_df.index)
    )
    n_fill = float(n_cells.median()) if n_cells.notna().any() else 10.0
    n_max = max(float(n_cells.max()) if n_cells.notna().any() else 1.0, 1.0)
    sizes = (20 + 80 * (n_cells.fillna(n_fill) / n_max)).to_numpy()

    sig = None
    for col in ("significant_adj_reported", "significant_adj"):
        if col in plot_df.columns:
            sig = plot_df[col].map(lambda v: bool(v) if pd.notna(v) else False)
            break

    x = plot_df["edistance"].to_numpy(dtype=float)
    y = plot_df["neglog10p"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    if sig is not None:
        for flag, color, label in ((True, "#C44E52", "sig_adj_reported"), (False, "#4C72B0", "not reported")):
            m = sig.to_numpy() if flag else (~sig).to_numpy()
            if not m.any():
                continue
            ax.scatter(
                x[m],
                y[m],
                s=sizes[m],
                c=color,
                alpha=0.75,
                edgecolors="white",
                linewidths=0.4,
                label=label,
            )
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.scatter(
            x,
            y,
            s=sizes,
            c="#4C72B0",
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
        )

    # Label top points by E-distance
    top = plot_df.nlargest(min(12, len(plot_df)), "edistance")
    for _, row in top.iterrows():
        ax.annotate(
            str(row["gene_target"]),
            (row["edistance"], row["neglog10p"]),
            textcoords="offset points",
            xytext=(4, 3),
            fontsize=7,
            alpha=0.85,
        )

    ax.axhline(-np.log10(0.05), color="grey", lw=0.8, ls="--", alpha=0.7)
    ax.set_xlabel("E-distance vs control")
    ax.set_ylabel(ylabel)
    ax.set_title("Perturbation effect summary")
    fig.tight_layout()
    path = Path(out_dir) / filename
    _save(fig, path)
    return path


def plot_guide_consistency_matrix(
    adata: AnnData,
    gene: str,
    out_dir: Path,
    *,
    control: str = "NT",
    min_cells: int = DEFAULT_MIN_GUIDE_CELLS,
    n_genes: int = DEFAULT_CORR_N_GENES,
    guide_col: str | None = None,
) -> Path | None:
    """Heatmap of Pearson correlation of per-guide log2FC vectors for one gene."""
    guide_col = guide_col or _guide_col(adata)
    if guide_col is None or "gene_target" not in adata.obs.columns:
        return None

    obs = adata.obs
    gene_mask = obs["gene_target"].astype(str) == gene
    ctrl_mask = obs["gene_target"].astype(str) == control
    if not gene_mask.any() or not ctrl_mask.any():
        return None

    guides = (
        obs.loc[gene_mask, guide_col]
        .astype(str)
        .value_counts()
    )
    guides = guides[guides >= min_cells]
    if len(guides) < 2:
        return None

    gene_names = _select_genes_for_corr(adata, n_genes)
    gene_idx = np.array([adata.var_names.get_loc(g) for g in gene_names if g in adata.var_names], dtype=int)
    if len(gene_idx) < 5:
        return None

    ctrl_mean = _mean_profile(adata, ctrl_mask.to_numpy(), gene_idx)
    ctrl_mean = np.maximum(ctrl_mean, 0.0) + 1e-6

    vectors: dict[str, np.ndarray] = {}
    for gid, _n in guides.items():
        mask = gene_mask & (obs[guide_col].astype(str) == gid)
        mean = _mean_profile(adata, mask.to_numpy(), gene_idx)
        if not np.isfinite(mean).any():
            continue
        lfc = np.log2((np.maximum(mean, 0.0) + 1e-6) / ctrl_mean)
        if np.nanstd(lfc) < 1e-12:
            continue
        vectors[str(gid)] = lfc

    if len(vectors) < 2:
        return None

    labels = list(vectors.keys())
    mat = np.vstack([vectors[k] for k in labels])
    # Pairwise Pearson; pandas corr handles NaNs
    corr = pd.DataFrame(mat.T, columns=labels).corr(method="pearson")

    fig, ax = plt.subplots(figsize=(max(3.5, 0.7 * len(labels) + 1.5), max(3.0, 0.7 * len(labels) + 1.2)))
    sns.heatmap(
        corr,
        ax=ax,
        vmin=-1,
        vmax=1,
        cmap="RdBu_r",
        annot=len(labels) <= 8,
        fmt=".2f",
        square=True,
        cbar_kws={"shrink": 0.8, "label": "Pearson r"},
    )
    ax.set_title(f"Guide consistency — {gene}\n(log2FC vs {control})")
    fig.tight_layout()
    path = Path(out_dir) / f"{GUIDE_CONSISTENCY_PREFIX}{_safe_stem(gene)}.png"
    _save(fig, path)
    return path


def plot_target_gene_validation(
    adata: AnnData,
    gene: str,
    out_dir: Path,
    *,
    control: str = "NT",
) -> Path | None:
    """Violin: control cells vs targeting cells for the target gene's own expression."""
    if "gene_target" not in adata.obs.columns or gene not in adata.var_names:
        return None
    gidx = adata.var_names.get_loc(gene)
    if not isinstance(gidx, (int, np.integer)):
        gidx = int(gidx.start)

    obs = adata.obs
    ctrl_mask = obs["gene_target"].astype(str) == control
    pert_mask = obs["gene_target"].astype(str) == gene
    if not ctrl_mask.any() or not pert_mask.any():
        return None

    ctrl_vals = _to_dense_block(adata.X[ctrl_mask.to_numpy(), gidx]).ravel()
    pert_vals = _to_dense_block(adata.X[pert_mask.to_numpy(), gidx]).ravel()
    plot_df = pd.DataFrame(
        {
            "expression": np.concatenate([ctrl_vals, pert_vals]),
            "group": ([f"{control} (n={len(ctrl_vals)})"] * len(ctrl_vals))
            + ([f"{gene} (n={len(pert_vals)})"] * len(pert_vals)),
        }
    )

    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    sns.violinplot(
        data=plot_df,
        x="group",
        y="expression",
        hue="group",
        ax=ax,
        inner="box",
        cut=0,
        palette=["#4C72B0", "#DD8452"],
        legend=False,
    )
    ax.set_xlabel("")
    ax.set_ylabel(f"{gene} expression")
    ax.set_title(f"Target validation — {gene}")
    fig.tight_layout()
    path = Path(out_dir) / f"{TARGET_VALIDATION_PREFIX}{_safe_stem(gene)}.png"
    _save(fig, path)
    return path


def _top_genes_for_plots(
    edistance: pd.DataFrame | None,
    guide_qc: pd.DataFrame | None,
    *,
    top_n: int,
    min_guides: int = 2,
    min_cells: int = DEFAULT_MIN_GUIDE_CELLS,
    multi_guide_only: bool = False,
) -> list[str]:
    """Rank genes by E-distance (fallback: total cells in guide_qc)."""
    scores: dict[str, float] = {}
    if edistance is not None and not edistance.empty:
        ed = edistance.copy()
        if "gene_target" not in ed.columns:
            ed = ed.rename(columns={ed.columns[0]: "gene_target"})
        if "edistance" in ed.columns:
            for _, row in ed.iterrows():
                scores[str(row["gene_target"])] = float(row["edistance"]) if pd.notna(row["edistance"]) else 0.0

    multi_guide: set[str] = set()
    if guide_qc is not None and not guide_qc.empty and "gene_target" in guide_qc.columns:
        gq = guide_qc.loc[~guide_qc["is_control"].astype(bool)] if "is_control" in guide_qc.columns else guide_qc
        for gene, sub in gq.groupby(gq["gene_target"].astype(str)):
            n_ok = int((pd.to_numeric(sub["n_cells"], errors="coerce").fillna(0) >= min_cells).sum()) if "n_cells" in sub.columns else len(sub)
            if n_ok >= min_guides:
                multi_guide.add(str(gene))
            if str(gene) not in scores:
                scores[str(gene)] = float(pd.to_numeric(sub["n_cells"], errors="coerce").fillna(0).sum()) if "n_cells" in sub.columns else float(len(sub))

    genes = sorted(scores, key=lambda g: (-scores[g], g))
    if multi_guide_only:
        genes = [g for g in genes if g in multi_guide]
    return genes[:top_n]


def generate_report_plots(
    output_dir: Path | str,
    *,
    tables: dict[str, Any] | None = None,
    adata: AnnData | None = None,
    h5ad: Path | str | None = None,
    control: str = "NT",
    top_n: int = DEFAULT_REPORT_PLOT_TOP_N,
    min_guide_cells: int = DEFAULT_MIN_GUIDE_CELLS,
    n_genes_corr: int = DEFAULT_CORR_N_GENES,
) -> dict[str, Any]:
    """Create report PNGs; return status with paths and skip reasons."""
    from perturbseq.report import load_result_tables

    output_dir = Path(output_dir)
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    if tables is None:
        tables = load_result_tables(output_dir)

    status: dict[str, Any] = {
        "skipped": False,
        "reason": None,
        "generated": [],
        "skipped_plots": [],
        "effect_summary": None,
        "guide_consistency": {},
        "target_validation": {},
    }

    # 1) Effect summary (tables only)
    effect_path = plot_perturbation_effect_summary(
        tables.get("edistance"),
        tables.get("etest"),
        figures,
    )
    if effect_path is not None:
        rel = f"figures/{effect_path.name}"
        status["effect_summary"] = rel
        status["generated"].append(rel)
    else:
        status["skipped_plots"].append(
            {"plot": "perturbation_effect_summary", "reason": "missing_edistance", "detail": "edistance.csv absent or empty"}
        )

    # Resolve expression data
    loaded_adata = adata
    if loaded_adata is None and h5ad is not None:
        h5ad_path = Path(h5ad)
        if h5ad_path.is_file():
            import scanpy as sc

            loaded_adata = sc.read_h5ad(h5ad_path)
        else:
            status["skipped_plots"].append(
                {"plot": "h5ad_load", "reason": "h5ad_not_found", "detail": str(h5ad_path)}
            )
    if loaded_adata is None:
        # Auto-discover tertiary h5ad next to report
        candidates = sorted(output_dir.glob("*.tertiary.h5ad"))
        if candidates:
            import scanpy as sc

            loaded_adata = sc.read_h5ad(candidates[0])

    if loaded_adata is None:
        note = "h5ad absent — guide consistency matrices and target validation violins skipped"
        status["skipped_plots"].append(
            {"plot": "guide_consistency", "reason": "h5ad_missing", "detail": note}
        )
        status["skipped_plots"].append(
            {"plot": "target_validation", "reason": "h5ad_missing", "detail": note}
        )
        return status

    # 2) Guide consistency (multi-guide genes)
    consistency_genes = _top_genes_for_plots(
        tables.get("edistance"),
        tables.get("guide_qc"),
        top_n=top_n,
        min_cells=min_guide_cells,
        multi_guide_only=True,
    )
    if not consistency_genes:
        status["skipped_plots"].append(
            {
                "plot": "guide_consistency",
                "reason": "no_multi_guide_genes",
                "detail": f"Need ≥2 guides with ≥{min_guide_cells} cells each",
            }
        )
    for gene in consistency_genes:
        path = plot_guide_consistency_matrix(
            loaded_adata,
            gene,
            figures,
            control=control,
            min_cells=min_guide_cells,
            n_genes=n_genes_corr,
        )
        if path is not None:
            rel = f"figures/{path.name}"
            status["guide_consistency"][gene] = rel
            status["generated"].append(rel)
        else:
            status["skipped_plots"].append(
                {
                    "plot": f"guide_consistency_{gene}",
                    "reason": "insufficient_guides_or_genes",
                    "detail": "Could not build log2FC vectors for ≥2 guides",
                }
            )

    # 3) Target validation (top-N by E-distance / cells)
    validation_genes = _top_genes_for_plots(
        tables.get("edistance"),
        tables.get("guide_qc"),
        top_n=top_n,
        multi_guide_only=False,
    )
    if not validation_genes:
        status["skipped_plots"].append(
            {"plot": "target_validation", "reason": "no_genes", "detail": "No ranked perturbations for violin plots"}
        )
    for gene in validation_genes:
        path = plot_target_gene_validation(loaded_adata, gene, figures, control=control)
        if path is not None:
            rel = f"figures/{path.name}"
            status["target_validation"][gene] = rel
            status["generated"].append(rel)
        else:
            status["skipped_plots"].append(
                {
                    "plot": f"target_validation_{gene}",
                    "reason": "gene_not_in_matrix_or_empty_groups",
                    "detail": f"{gene} missing from var_names or no NT/targeting cells",
                }
            )

    return status


__all__ = [
    "DEFAULT_REPORT_PLOT_TOP_N",
    "EFFECT_SUMMARY_NAME",
    "generate_report_plots",
    "plot_guide_consistency_matrix",
    "plot_perturbation_effect_summary",
    "plot_target_gene_validation",
]
