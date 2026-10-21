"""E-test from the pertpy Distance Tests tutorial.

The notebook calls ``DistanceTest("edistance")`` on ``X_pca`` and plots each
group's permutation p-value against its E-distance to the contrast. It does
not draw a pairwise heatmap.

https://pertpy.scverse.org/en/latest/tutorials/notebooks/distance_tests.html
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from pertpy.tools import DistanceTest

from perturbseq.assign import control_cell_mask
from perturbseq.config import ZERO_GUIDE_LABEL, control_spec

INFERENCE = "within_sample_distance_test"
INFERENCE_NOTE = (
    "单样本内距离检验（within-sample distance test）："
    "Monte Carlo 置换打乱的是这一个生物学样本中的细胞，不是生物学重复。"
    "对这些置换 p 值做多重校正，不会得到生物学重复水平上的 FDR。"
    "功效和样本间可重复性无法估计。"
    "p 值、校正后 p 值和显著性不写入 effect_main.tsv 或 effect_sensitivity.tsv。"
)
_TUTORIAL_N_PERMS = 1000
_TUTORIAL_PADJ = 0.0015
_OBSM_KEY = "X_pca"


def run_edistance_test(
    adata: AnnData,
    *,
    contrast: str,
    groupby: str = "target_gene",
    groups: list[str] | None = None,
    n_perms: int = _TUTORIAL_N_PERMS,
    padj_threshold: float = _TUTORIAL_PADJ,
    obsm_key: str = _OBSM_KEY,
    selection_rule: str = "",
    show_progressbar: bool = False,
) -> pd.DataFrame:
    """Test each perturbation group against ``contrast`` with pertpy's E-test.

    Cells are the permutation units, matching the tutorial. The result is a
    within-sample distance test for this one biological sample.
    """
    if obsm_key not in adata.obsm:
        raise ValueError(f"缺少 {obsm_key}。教程要求先在 log-normalized 表达上做 PCA。")
    view, counts = _pca_view(
        adata, groupby=groupby, contrast=contrast, groups=groups, obsm_key=obsm_key
    )
    etest = DistanceTest(
        "edistance",
        n_perms=n_perms,
        obsm_key=obsm_key,
        padj_threshold=padj_threshold,
    )
    tab = etest(view, groupby=groupby, contrast=contrast, show_progressbar=show_progressbar)
    tab = tab.copy()
    tab.index = tab.index.astype(str)
    n_contrast = int(counts.get(contrast, 0))
    n_pcs = int(view.obsm[obsm_key].shape[1])
    rows = []
    for group in tab.index.astype(str):
        record = tab.loc[group]
        rows.append(
            {
                "group": group,
                "contrast": contrast,
                "pair": f"{group} vs {contrast}",
                "distance": float(record["distance"]),
                "pvalue": float(record["pvalue"]),
                "significant": bool(record["significant"]),
                "pvalue_adj": float(record["pvalue_adj"]),
                "significant_adj": bool(record["significant_adj"]),
                "n_cells": int(counts.get(group, 0)),
                "n_cells_contrast": n_contrast,
                "n_perms": int(n_perms),
                "metric": "edistance",
                "obsm_key": obsm_key,
                "n_pcs": n_pcs,
                "padj_threshold": float(padj_threshold),
                "correction": str(etest.correction),
                "inference": INFERENCE,
                "inference_note": INFERENCE_NOTE,
                "selection_rule": selection_rule,
            }
        )
    return pd.DataFrame(rows)


def write_edistance_results(
    output_dir: Path,
    adata: AnnData,
    config: dict[str, Any],
    *,
    n_perms: int = _TUTORIAL_N_PERMS,
    padj_threshold: float = _TUTORIAL_PADJ,
    show_progressbar: bool = False,
) -> pd.DataFrame:
    """Write the E-test table and scatter under ``E-distance/``."""
    figures, tables = _output_dirs(output_dir)
    spec = control_spec(config)
    contrast = spec["control_label"]
    n_control = int(control_cell_mask(adata.obs, spec).sum()) if adata.n_obs else 0
    if n_control == 0:
        reason = "no_zero_guide_cells" if spec["temporary"] else "no_designed_nt_cells"
        print(f"  E-distance skipped: {reason}", flush=True)
        _clear_edistance_skip_markers(figures)
        return pd.DataFrame()
    if not contrast:
        print("  E-distance skipped: missing_control_label", flush=True)
        _clear_edistance_skip_markers(figures)
        return pd.DataFrame()
    if _OBSM_KEY not in adata.obsm:
        print("  E-distance skipped: missing_X_pca", flush=True)
        _clear_edistance_skip_markers(figures)
        return pd.DataFrame()
    gate, source = _load_de_gate(adata, output_dir)
    groups, rule = _select_groups(gate, contrast, source)
    print(
        f"  E-distance: {len(groups) if groups is not None else 'all'} groups vs {contrast}, "
        f"n_perms={n_perms}",
        flush=True,
    )
    table = run_edistance_test(
        adata,
        contrast=contrast,
        groups=groups,
        n_perms=n_perms,
        padj_threshold=padj_threshold,
        selection_rule=rule,
        show_progressbar=show_progressbar,
    )
    # Lean outputs: PNG only; table stays in the perturb checkpoint via caller if needed.
    if not table.empty:
        plot_edistance_test(table, figures / "edistance_etest.png")
    _clear_edistance_skip_markers(figures)
    stale_table = tables / "edistance_etest.tsv"
    if stale_table.exists():
        stale_table.unlink()
    return table


def _clear_edistance_skip_markers(figures: Path) -> None:
    skipped = figures / "edistance_etest.skipped.txt"
    if skipped.exists():
        skipped.unlink()


def plot_edistance_test(table: pd.DataFrame, path: Path) -> None:
    """Scatter of E-test p-value versus E-distance, as in the tutorial notebook."""
    if table.empty:
        raise ValueError("E-test 表是空的，不能作图。")
    contrast = str(table["contrast"].iloc[0])
    plotted = table[table["group"].astype(str) != contrast].copy()
    # Color by adjusted significance when present: matches padj_threshold in the table.
    hue_col = "significant_adj" if "significant_adj" in plotted.columns else "significant"
    plotted[hue_col] = _as_bool(plotted[hue_col])
    path.parent.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams["font.sans-serif"] = ["Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    figure, axis = plt.subplots(figsize=(7.2, 5.6))
    with sns.axes_style("darkgrid"):
        sns.scatterplot(
            data=plotted,
            x="pvalue",
            y="distance",
            hue=hue_col,
            palette={True: "green", False: "red"},
            ax=axis,
        )
    shown = "zero-guide cells" if contrast == ZERO_GUIDE_LABEL else contrast
    axis.set_title(f"E-test on expression X_pca vs {shown} (within-sample; not replicates)")
    axis.set_xlabel("within-sample p-value (cell permutation; not biological replicates)")
    axis.set_ylabel(f"E-distance to {shown}")
    legend = axis.get_legend()
    if legend is not None:
        legend.set_title(
            "within-sample sig. (adj)" if hue_col == "significant_adj" else "within-sample sig."
        )
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    if path.stat().st_size <= 0:
        raise RuntimeError(f"图文件为空: {path}")
    for stale in (
        path.with_suffix(".pdf"),
        path.with_suffix(".txt"),
        path.parent / f"{path.stem}.skipped.txt",
    ):
        if stale.exists():
            stale.unlink()


def _pca_view(
    adata: AnnData,
    *,
    groupby: str,
    contrast: str,
    groups: list[str] | None,
    obsm_key: str,
) -> tuple[AnnData, dict[str, int]]:
    labels = adata.obs[groupby].astype(str)
    mask = np.ones(adata.n_obs, dtype=bool)
    if "n_targets" in adata.obs.columns:
        n_targets = pd.to_numeric(adata.obs["n_targets"], errors="coerce")
        keep_cells = n_targets.eq(1)
        if contrast == ZERO_GUIDE_LABEL:
            keep_cells = keep_cells | labels.eq(contrast)
        mask &= keep_cells.to_numpy()
    elif "n_guides" in adata.obs.columns:
        n_guides = pd.to_numeric(adata.obs["n_guides"], errors="coerce")
        keep_cells = n_guides.eq(1)
        if contrast == ZERO_GUIDE_LABEL:
            keep_cells = keep_cells | labels.eq(contrast)
        mask &= keep_cells.to_numpy()
    if groups is not None:
        mask &= labels.isin(set(groups) | {contrast}).to_numpy()
    indices = np.flatnonzero(mask)
    kept = labels.iloc[indices]
    if contrast not in set(kept):
        raise ValueError(f"对照 {contrast} 在分析细胞里没有细胞。")
    if kept.nunique() < 2:
        raise ValueError("没有可与对照比较的扰动组。")
    embedding = np.asarray(adata.obsm[obsm_key][indices], dtype=np.float64)
    if embedding.ndim != 2 or embedding.shape[1] < 1:
        raise ValueError(f"{obsm_key} 不是可用的嵌入。")
    obs = pd.DataFrame({groupby: kept.to_numpy()})
    view = AnnData(X=np.zeros((obs.shape[0], 1), dtype=np.float32), obs=obs)
    view.obsm[obsm_key] = embedding
    counts = {str(name): int(count) for name, count in kept.value_counts().items()}
    return view, counts


def _load_de_gate(adata: AnnData, output_dir: Path) -> tuple[pd.DataFrame | None, str]:
    gate = adata.uns.get("mixscape_de_gate") if adata.uns is not None else None
    if isinstance(gate, pd.DataFrame) and not gate.empty:
        return gate, "adata.uns['mixscape_de_gate']"
    for relative in (
        Path("Mixscape") / "tables" / "mixscape_de_gate.tsv",
        Path("DGE") / "Mixscape" / "tables" / "mixscape_de_gate.tsv",
        Path("DGE") / "mixscape" / "tables" / "mixscape_de_gate.tsv",
        Path("DGE") / "tables" / "mixscape_de_gate.tsv",
        Path("tables") / "mixscape_de_gate.tsv",
    ):
        path = output_dir / relative
        if path.is_file():
            return pd.read_csv(path, sep="\t"), str(relative)
    return None, ""


def _select_groups(
    gate: pd.DataFrame | None, contrast: str, source: str
) -> tuple[list[str] | None, str]:
    base = "教程图是每个组相对对照的 E-test 散点，不是两两距离热图。"
    if gate is None or not {"passed", "target_gene"}.issubset(gate.columns):
        return None, (
            f"{base}未使用差异基因门槛，对象里每一个有细胞的 target_gene 都相对 {contrast} 检验。"
        )
    passed = gate.loc[_as_bool(gate["passed"]), "target_gene"].astype(str)
    chosen = sorted({gene for gene in passed if gene and gene != contrast})
    if not chosen:
        return None, (f"{base}差异基因门槛里没有 passed 的靶基因，改为检验对象中的全部分组。")
    return chosen, (
        f"{base}全量靶基因的两两矩阵不是教程的主结果，"
        f"因此使用 {source} 中 passed 为 true 的 {len(chosen)} 个靶基因，加上对照 {contrast}。"
        "距离用对象里已有的 X_pca，即 log-normalized 高变基因的 PCA。"
    )


def _output_dirs(output_dir: Path) -> tuple[Path, Path]:
    figures = output_dir / "E-distance" / "figures"
    tables = output_dir / "E-distance" / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    return figures, tables


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.strip().str.lower().isin(["true", "1"])
