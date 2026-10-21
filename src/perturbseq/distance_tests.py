"""E-test on expression ``X_pca`` (pertpy tutorial metric, memory-safe runner).

The pertpy tutorial calls ``DistanceTest("edistance")`` on ``X_pca``. That path
precomputes a dense pairwise distance matrix for every group∪contrast (~N²
floats). With a large temporary zero-guide contrast this OOMs.

This module keeps the same statistic (energy distance + cell-label Monte Carlo)
but:

* never materializes an N×N matrix (uses pertpy's euclidean kernel-sum edistance);
* tests one group at a time;
* optionally downsamples a huge contrast to ``max_contrast_cells`` (seeded);
* reads ``n_perms`` / thresholds from config.

https://pertpy.scverse.org/en/latest/tutorials/notebooks/distance_tests.html
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from pertpy.tools import Distance
from statsmodels.stats.multitest import multipletests

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
_DEFAULT_MAX_CONTRAST = 4000
_OBSM_KEY = "X_pca"
_CORRECTION = "holm-sidak"


def run_edistance_test(
    adata: AnnData,
    *,
    contrast: str,
    groupby: str = "target_gene",
    groups: list[str] | None = None,
    n_perms: int = _TUTORIAL_N_PERMS,
    padj_threshold: float = _TUTORIAL_PADJ,
    obsm_key: str = _OBSM_KEY,
    max_contrast_cells: int = _DEFAULT_MAX_CONTRAST,
    seed: int = 0,
    selection_rule: str = "",
    show_progressbar: bool = False,
) -> pd.DataFrame:
    """Test each perturbation group against ``contrast`` with energy distance.

    Cells are the permutation units, matching the tutorial. The result is a
    within-sample distance test for this one biological sample.
    """
    del show_progressbar  # kept for call-site compatibility; streaming has no bar
    if obsm_key not in adata.obsm:
        raise ValueError(f"缺少 {obsm_key}。教程要求先在 log-normalized 表达上做 PCA。")
    view, counts = _pca_view(
        adata, groupby=groupby, contrast=contrast, groups=groups, obsm_key=obsm_key
    )
    view, counts, subsample_note = _subsample_contrast(
        view,
        groupby=groupby,
        contrast=contrast,
        max_contrast_cells=max_contrast_cells,
        seed=seed,
        counts=counts,
    )
    tab = _estream_etest(
        view,
        groupby=groupby,
        contrast=contrast,
        n_perms=n_perms,
        padj_threshold=padj_threshold,
        obsm_key=obsm_key,
        seed=seed,
    )
    n_contrast = int(counts.get(contrast, 0))
    n_pcs = int(view.obsm[obsm_key].shape[1])
    rule = selection_rule
    if subsample_note:
        rule = f"{rule} {subsample_note}".strip()
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
                "correction": _CORRECTION,
                "inference": INFERENCE,
                "inference_note": INFERENCE_NOTE,
                "selection_rule": rule,
            }
        )
    return pd.DataFrame(rows)


def write_edistance_results(
    output_dir: Path,
    adata: AnnData,
    config: dict[str, Any],
    *,
    n_perms: int | None = None,
    padj_threshold: float | None = None,
    show_progressbar: bool = False,
) -> pd.DataFrame:
    """Write the E-test table and scatter under ``E-distance/``."""
    figures, tables = _output_dirs(output_dir)
    spec = control_spec(config)
    contrast = spec["control_label"]
    settings = _edistance_settings(config, n_perms=n_perms, padj_threshold=padj_threshold)
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
        f"n_perms={settings['n_perms']}, max_contrast_cells={settings['max_contrast_cells']}",
        flush=True,
    )
    table = run_edistance_test(
        adata,
        contrast=contrast,
        groups=groups,
        n_perms=int(settings["n_perms"]),
        padj_threshold=float(settings["padj_threshold"]),
        max_contrast_cells=int(settings["max_contrast_cells"]),
        seed=int(settings["seed"]),
        selection_rule=rule,
        show_progressbar=show_progressbar,
    )
    if not table.empty:
        plot_edistance_test(table, figures / "edistance_etest.png")
    _clear_edistance_skip_markers(figures)
    legacy_png = figures / "figures" / "edistance_etest.png"
    if legacy_png.is_file():
        legacy_png.unlink()
    stale_table = tables / "edistance_etest.tsv"
    if stale_table.exists():
        stale_table.unlink()
    return table


def _edistance_settings(
    config: dict[str, Any],
    *,
    n_perms: int | None,
    padj_threshold: float | None,
) -> dict[str, Any]:
    block = config.get("edistance") if isinstance(config.get("edistance"), dict) else {}
    runtime = config.get("runtime") if isinstance(config.get("runtime"), dict) else {}
    run = config.get("run") if isinstance(config.get("run"), dict) else {}
    resolved_perms = int(block.get("n_perms", _TUTORIAL_N_PERMS)) if n_perms is None else int(n_perms)
    resolved_padj = (
        float(block.get("padj_threshold", _TUTORIAL_PADJ))
        if padj_threshold is None
        else float(padj_threshold)
    )
    return {
        "n_perms": resolved_perms,
        "padj_threshold": resolved_padj,
        "max_contrast_cells": int(block.get("max_contrast_cells", _DEFAULT_MAX_CONTRAST)),
        "seed": int(runtime.get("seed", run.get("random_seed", 0))),
    }


def _estream_etest(
    adata: AnnData,
    *,
    groupby: str,
    contrast: str,
    n_perms: int,
    padj_threshold: float,
    obsm_key: str,
    seed: int,
) -> pd.DataFrame:
    """One-group-at-a-time E-test without pairwise N×N matrices."""
    metric = Distance("edistance", obsm_key=obsm_key).metric_fct
    embedding = np.asarray(adata.obsm[obsm_key], dtype=np.float64, order="C")
    labels = adata.obs[groupby].astype(str).to_numpy()
    contrast_idx = np.flatnonzero(labels == contrast)
    if contrast_idx.size == 0:
        raise ValueError(f"对照 {contrast} 在分析细胞里没有细胞。")
    contrast_xy = embedding[contrast_idx]
    groups = [str(name) for name in pd.unique(labels) if str(name) != contrast]
    if not groups:
        raise ValueError("没有可与对照比较的扰动组。")

    distances = np.full(len(groups), np.nan, dtype=float)
    pvalues = np.ones(len(groups), dtype=float)
    for index, group in enumerate(groups):
        group_idx = np.flatnonzero(labels == group)
        if group_idx.size == 0:
            continue
        group_xy = embedding[group_idx]
        observed = float(metric(group_xy, contrast_xy))
        distances[index] = observed
        combined = np.vstack((group_xy, contrast_xy))
        n_group = int(group_xy.shape[0])
        n_total = int(combined.shape[0])
        rng = np.random.default_rng(_group_seed(seed, group))
        n_extreme = 0
        for _ in range(int(n_perms)):
            order = rng.permutation(n_total)
            permuted = float(metric(combined[order[:n_group]], combined[order[n_group:]]))
            if permuted >= observed:
                n_extreme += 1
        pvalues[index] = max(n_extreme, 1) / float(n_perms)
        del combined, group_xy

    significant_adj, pvalue_adj, _, _ = multipletests(
        pvalues, alpha=padj_threshold, method=_CORRECTION
    )
    tab = pd.DataFrame(
        {
            "distance": distances,
            "pvalue": pvalues,
            "significant": pvalues < padj_threshold,
            "pvalue_adj": pvalue_adj,
            "significant_adj": significant_adj,
        },
        index=pd.Index(groups, name=groupby),
    )
    tab.loc[contrast, "distance"] = 0.0
    tab.loc[contrast, "pvalue"] = 1.0
    tab.loc[contrast, "significant"] = False
    tab.loc[contrast, "pvalue_adj"] = 1.0
    tab.loc[contrast, "significant_adj"] = False
    return tab


def _group_seed(seed: int, group: str) -> int:
    payload = f"{int(seed)}\0{group}".encode()
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2**31 - 1)


def _subsample_contrast(
    adata: AnnData,
    *,
    groupby: str,
    contrast: str,
    max_contrast_cells: int,
    seed: int,
    counts: dict[str, int],
) -> tuple[AnnData, dict[str, int], str]:
    """Cap huge temporary-control sets so E-test stays in memory/time budget."""
    if max_contrast_cells <= 0:
        return adata, counts, ""
    labels = adata.obs[groupby].astype(str)
    contrast_mask = labels.eq(contrast).to_numpy()
    n_contrast = int(contrast_mask.sum())
    if n_contrast <= int(max_contrast_cells):
        return adata, counts, ""
    rng = np.random.default_rng(int(seed))
    keep_contrast = rng.choice(
        np.flatnonzero(contrast_mask), size=int(max_contrast_cells), replace=False
    )
    keep = np.flatnonzero(~contrast_mask)
    indices = np.sort(np.concatenate([keep, keep_contrast]))
    view = adata[indices].copy()
    kept = view.obs[groupby].astype(str)
    updated = {str(name): int(count) for name, count in kept.value_counts().items()}
    note = (
        f"对照 {contrast} 从 {n_contrast} 细胞随机下采样到 {max_contrast_cells} "
        f"（edistance.max_contrast_cells；seed={seed}），以控制 pairwise 置换内存与时间。"
    )
    print(f"  E-distance: {note}", flush=True)
    return view, updated, note


def _clear_edistance_skip_markers(figures: Path) -> None:
    for skipped in (
        figures / "edistance_etest.skipped.txt",
        figures / "figures" / "edistance_etest.skipped.txt",
    ):
        if skipped.exists():
            skipped.unlink()


def plot_edistance_test(table: pd.DataFrame, path: Path) -> None:
    """Scatter of E-test p-value versus E-distance, as in the tutorial notebook."""
    if table.empty:
        raise ValueError("E-test 表是空的，不能作图。")
    contrast = str(table["contrast"].iloc[0])
    plotted = table[table["group"].astype(str) != contrast].copy()
    hue_col = "significant_adj" if "significant_adj" in plotted.columns else "significant"
    plotted[hue_col] = _as_bool(plotted[hue_col])
    path.parent.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
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
    shown = "zero-guide" if contrast == ZERO_GUIDE_LABEL else contrast
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
    # Slim AnnData: tiny X + labels + PCA (avoid copying the expression matrix).
    obs = pd.DataFrame({groupby: kept.to_numpy()})
    view = AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs,
    )
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
    figures = output_dir / "E-distance"
    tables = output_dir / "E-distance" / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    return figures, tables


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.strip().str.lower().isin(["true", "1"])
