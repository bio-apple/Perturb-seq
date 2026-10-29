"""Save PNG figures and the table behind each figure."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, NamedTuple

import matplotlib

matplotlib.use("Agg", force=False)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy import sparse

from perturbseq.assign import STATUS_TARGET, control_cell_mask
from perturbseq.config import ZERO_GUIDE_LABEL

TEMPORARY_CONTROL_QC_NAME = "temporary_control_qc"
ZERO_GUIDE_QC_GROUP = "zero-guide"
TARGETING_QC_GROUP = "targeting"
_TEMPORARY_SPEC = {"temporary": True, "control_label": ZERO_GUIDE_LABEL}

# Assign-step QC: full ledger / pre gene-mito filter (or descriptive temporary control).
_QC_FIGURES = frozenset(
    {
        "qc_metrics_panel",
        "qc_scatter_panel",
        TEMPORARY_CONTROL_QC_NAME,
        "guides_per_cell",
        "cells_per_target_bar",
        "cells_per_guide_name_bar",
        "target_cell_counts",
    }
)
# Preprocess figures on QC-pass cells → UMAP/ (not qc/).
_UMAP_QC_PASS_FIGURES = frozenset(
    {
        "variable_features",
        "pca_elbow",
        "marker_detection_heatmap",
        "dotplot_markers",
    }
)
# Library-D / ggplot2-like defaults (figures under xihu/Library-D/figures).
_HIST_FILL = "#305080"
_THRESHOLD_LINE = "#D62728"
_GRID_MAJOR = "#E0E0E0"
_GRID_MINOR = "#F0F0F0"
# Seurat pbmc3k VlnPlot default (ggplot2 hue 1 / salmon).
_SEURAT_VIOLIN_FILL = "#F8766D"
_QC_PANEL_NAME = "qc_metrics_panel"
_QC_SCATTER_PANEL_NAME = "qc_scatter_panel"
_LEGACY_QC_HIST_NAMES = ("hist_ncount", "hist_nfeature", "hist_percent_mt_gex")
_LEGACY_QC_SCATTER_NAMES = ("qc_scatter", "scatter_ncount_nfeature_gex")
# Seurat DoHeatmap-like Expression: purple → black → yellow.
_MARKER_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "seurat_doheatmap",
    ["#FF00FF", "#000000", "#FFFF00"],
)
# Seurat DotPlot-like Average Expression: lightgrey → blue.
_DOTPLOT_CMAP = LinearSegmentedColormap.from_list(
    "seurat_dotplot",
    ["#D3D3D3", "#9ECAE1", "#3182BD", "#08519C", "#0000FF"],
)
# Mixscape DE gate / class / 2x2 panels. Signature UMAPs stay under UMAP/ by name.
_MIXSCAPE_FIGURES = frozenset(
    {
        "target_de_counts",
        "target_in_de_gene_counts",
        "target_in_de",
        "mixscape_target_panels",
        "mixscape_class_global_bar",
    }
)
_MIXSCALE_FIGURES = frozenset({"mixscale_score_ridge", "mixscale_scatter"})


def _output_bucket(name: str, section: str) -> str:
    # Signature UMAPs may pass section "mixscape"; the image still goes to UMAP/.
    if name.startswith("umap_") or name in _UMAP_QC_PASS_FIGURES:
        return "umap"
    # Mixscale section / heatmap stems win before shared DE-count names.
    if (
        section in {"mixscale", "mixscale_target_in_de"}
        or name in _MIXSCALE_FIGURES
        or name.endswith("_mixscale_heatmap")
    ):
        return "mixscale"
    if (
        section in {"mixscape", "target_in_de"}
        or name in _MIXSCAPE_FIGURES
        or name.endswith("_mixscape_2x2")
    ):
        return "mixscape"
    if name in _QC_FIGURES:
        return "qc"
    return "main"


def _analysis_output_root(output_dir: Path, bucket: str) -> Path:
    if bucket == "mixscape":
        return output_dir / "Mixscape"
    if bucket == "mixscale":
        return output_dir / "mixscale"
    raise ValueError(f"不是分析输出桶: {bucket}")


class FigureStore:
    def __init__(self, output_dir: Path, config: dict[str, Any], context: dict[str, str]) -> None:
        self.output_dir = output_dir
        self.plotting = config["plotting"]
        self.singler_enabled = config.get("singler", {}).get("enabled", True) is not False
        self.context = context
        self.plot_data = output_dir / "plot_data"
        self.plot_data.mkdir(parents=True, exist_ok=True)
        theme = str(self.plotting.get("theme", "white"))
        sns.set_theme(style=theme if theme in {"white", "whitegrid", "ticks"} else "white")
        matplotlib.rcParams.update(
            {
                "font.family": "sans-serif",
                "font.sans-serif": [
                    "Arial",
                    "Helvetica",
                    "DejaVu Sans",
                    "Heiti SC",
                    "Arial Unicode MS",
                ],
                "axes.unicode_minus": False,
                "axes.facecolor": "white",
                "figure.facecolor": "white",
                "axes.edgecolor": "black",
                "axes.linewidth": 0.8,
                "axes.grid": True,
                "grid.color": _GRID_MAJOR,
                "grid.linewidth": 0.7,
                "xtick.direction": "out",
                "ytick.direction": "out",
                "legend.frameon": False,
            }
        )

    def figure_dir(self, name: str, section: str) -> Path:
        bucket = _output_bucket(name, section)
        if bucket == "umap":
            return self.output_dir / "UMAP"
        if bucket == "qc":
            return self.output_dir / "qc"
        if bucket in {"mixscape", "mixscale"}:
            root = _analysis_output_root(self.output_dir, bucket)
            if bucket == "mixscape" and section == "target_in_de":
                return root / "target_in_de"
            if bucket == "mixscale" and section == "mixscale_target_in_de":
                return root / "target_in_de"
            return root
        return self.output_dir / section

    def data_dir(self, name: str, section: str = "") -> Path:
        bucket = _output_bucket(name, section)
        if bucket == "umap":
            return self.output_dir / "UMAP" / "tables"
        if bucket == "qc":
            return self.output_dir / "qc" / "tables"
        if bucket in {"mixscape", "mixscale"}:
            return _analysis_output_root(self.output_dir, bucket) / "plot_data"
        return self.plot_data

    def save(
        self, figure: Figure, name: str, section: str, data: pd.DataFrame, method: str
    ) -> None:
        folder = self.figure_dir(name, section)
        folder.mkdir(parents=True, exist_ok=True)
        if getattr(figure, "perturbseq_manual_layout", False):
            figure.subplots_adjust(left=0.22, right=0.9, bottom=0.14, top=0.96)
        elif not getattr(figure, "perturbseq_skip_tight_layout", False):
            figure.tight_layout()
        path = folder / f"{name}.png"
        figure.savefig(path, dpi=int(self.plotting["dpi"]), bbox_inches="tight")
        if path.stat().st_size <= 0:
            raise RuntimeError(f"图文件为空: {path}")
        for stale in (folder / f"{name}.pdf", folder / f"{name}.txt"):
            if stale.exists():
                stale.unlink()
        plt.close(figure)
        data_folder = self.data_dir(name, section)
        data_folder.mkdir(parents=True, exist_ok=True)
        stale = data_folder / f"{name}.skipped.txt"
        if stale.exists():
            stale.unlink()
        legacy = self.plot_data / f"{name}.skipped.txt"
        if legacy != stale and legacy.exists():
            legacy.unlink()
        if self.plotting.get("save_plot_data", False):
            data_folder.mkdir(parents=True, exist_ok=True)
            data.to_csv(data_folder / f"{name}.tsv", sep="\t", index=False)
        else:
            stale_tsv = data_folder / f"{name}.tsv"
            if stale_tsv.exists():
                stale_tsv.unlink()
        bucket = _output_bucket(name, section)
        self._retire_legacy_plot(name, section, bucket)

    def skip(self, name: str, reason: str, section: str = "") -> None:
        # Lean outputs: do not leave *.skipped.txt on disk (cleanup would remove them anyway).
        data_folder = self.data_dir(name, section)
        for stale in (
            data_folder / f"{name}.skipped.txt",
            data_folder / f"{name}.tsv",
            self.plot_data / f"{name}.skipped.txt",
        ):
            if stale.exists():
                stale.unlink()
        bucket = _output_bucket(name, section)
        self._retire_legacy_plot(name, section, bucket)

    def _retire_legacy_plot(self, name: str, section: str, bucket: str) -> None:
        """Remove products left under pre-bucket figure or table paths."""
        relatives = [
            Path("figures") / "plot" / f"{name}.png",
            Path("figures") / "plot" / f"{name}.txt",
            Path("figures") / "plot" / f"{name}.pdf",
            Path("figures") / "embedding" / f"{name}.png",
            Path("figures") / "embedding" / f"{name}.txt",
            Path("figures") / "embedding" / f"{name}.pdf",
            Path("figures") / "mixscape" / f"{name}.png",
            Path("figures") / "mixscape" / f"{name}.txt",
            Path("figures") / "mixscape" / f"{name}.pdf",
            Path("figures") / "mixscale" / f"{name}.png",
            Path("figures") / "de" / f"{name}.png",
            Path("figures") / section / f"{name}.png" if section else Path("figures") / f"{name}.png",
            Path("plot_data") / f"{name}.tsv",
            Path("plot_data") / f"{name}.skipped.txt",
            # Prior layout nested PNGs under */figures/.
            Path("UMAP") / "figures" / f"{name}.png",
            Path("qc") / "figures" / f"{name}.png",
            Path("Mixscape") / "figures" / f"{name}.png",
            Path("Mixscape") / "figures" / "target_in_de" / f"{name}.png",
            Path("mixscale") / "figures" / f"{name}.png",
            Path("wilcoxon") / "figures" / f"{name}.png",
            # Retired output roots (old descriptive / E-test steps).
            Path("DGE") / f"{name}.png",
            Path("DGE") / "figures" / f"{name}.png",
            Path("DGE") / "descriptive" / f"{name}.png",
            Path("E-distance") / f"{name}.png",
            Path("E-distance") / "figures" / f"{name}.png",
        ]
        # QC-pass preprocess plots moved from qc/ → UMAP/.
        if bucket == "umap" and name in _UMAP_QC_PASS_FIGURES:
            relatives.extend(
                (
                    Path("qc") / f"{name}.png",
                    Path("qc") / f"{name}.pdf",
                    Path("qc") / f"{name}.txt",
                    Path("qc") / "tables" / f"{name}.tsv",
                    Path("qc") / "tables" / f"{name}.skipped.txt",
                )
            )
        keep = {
            self.figure_dir(name, section) / f"{name}.png",
            self.data_dir(name, section) / f"{name}.tsv",
            self.data_dir(name, section) / f"{name}.skipped.txt",
        }
        for relative in relatives:
            stale = self.output_dir / relative
            if stale.is_file() and stale not in keep:
                stale.unlink()

    def _figure(
        self, width: float | None = None, height: float | None = None
    ) -> tuple[Figure, Axes]:
        figure, axis = plt.subplots(
            figsize=(
                width or float(self.plotting["width"]),
                height or float(self.plotting["height"]),
            )
        )
        return figure, axis

    def palette(self, key: str, fallback: str = "#4C78A8") -> str:
        return str(self.plotting.get("palette", {}).get(key, fallback))

def figure_context(
    output_dir: Path, config: dict[str, Any], adata: AnnData | None
) -> dict[str, str]:
    """Shared figure caption context for ``FigureStore``."""
    import json

    from perturbseq.config import control_spec, expected_expression, mito_prefix

    _direction, sentence = expected_expression(config)
    spec = control_spec(config)
    thresholds: dict[str, Any] = {}
    if adata is not None:
        stored = adata.uns.get("qc_thresholds")
        if isinstance(stored, dict):
            thresholds = dict(stored)
    if not thresholds:
        threshold_path = Path(output_dir) / "tables" / "qc_thresholds.json"
        if threshold_path.exists():
            thresholds = json.loads(threshold_path.read_text())
    n_cells = int(adata.n_obs) if adata is not None else 0
    return {
        "sample": str(config["experiment"]["sample_id"]),
        "organism": str(config["experiment"]["organism"]),
        "crispr_mode": str(config["experiment"]["crispr_mode"]),
        "n_cells": str(n_cells),
        "mito_prefix": mito_prefix(config),
        "qc": json.dumps(thresholds, ensure_ascii=False, default=_json_default),
        "direction": sentence,
        "seed": str(config["runtime"]["seed"]),
        "control_source": spec["control_source"],
        "control_label": spec["control_label"],
        "control_phrase": spec["control_phrase"],
    }


def _json_default(value: Any) -> Any:
    """Serialize DataFrames / numpy scalars for figure_context qc metadata."""
    if isinstance(value, pd.DataFrame):
        return value.replace({np.nan: None}).to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.replace({np.nan: None}).tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def open_figure_store(
    output_dir: Path, config: dict[str, Any], adata: AnnData | None = None
) -> FigureStore:
    return FigureStore(Path(output_dir), config, figure_context(output_dir, config, adata))


def retire_removed_umaps(output_dir: Path) -> None:
    """Remove non-current UMAP PNG/TSV names if present under the output root."""
    names = (
        "umap_leiden",
        "umap_assignment",
        "umap_target",
        "umap_status",
        "pca_scatter",
    )
    root = Path(output_dir)
    relatives: list[Path] = []
    for name in names:
        relatives.extend(
            (
                Path("UMAP") / f"{name}.png",
                Path("UMAP") / "figures" / f"{name}.png",
                Path("UMAP") / "tables" / f"{name}.tsv",
                Path("UMAP") / "tables" / f"{name}.skipped.txt",
                Path("figures") / "plot" / f"{name}.png",
                Path("figures") / "embedding" / f"{name}.png",
                Path("plot_data") / f"{name}.tsv",
                Path("plot_data") / f"{name}.skipped.txt",
            )
        )
    for relative in relatives:
        stale = root / relative
        if stale.is_file():
            stale.unlink()


def retire_removed_mixscape_fraction(output_dir: Path) -> None:
    """Remove retired mixscape_fraction_bar products if present under the output root."""
    name = "mixscape_fraction_bar"
    root = Path(output_dir)
    for relative in (
        Path("Mixscape") / f"{name}.png",
        Path("Mixscape") / "figures" / f"{name}.png",
        Path("Mixscape") / "plot_data" / f"{name}.tsv",
        Path("Mixscape") / "plot_data" / f"{name}.skipped.txt",
        Path("figures") / "mixscape" / f"{name}.png",
        Path("figures") / "mixscape" / f"{name}.txt",
        Path("figures") / "mixscape" / f"{name}.pdf",
        Path("plot_data") / f"{name}.tsv",
        Path("plot_data") / f"{name}.skipped.txt",
    ):
        stale = root / relative
        if stale.is_file():
            stale.unlink()


def _label_count_bars(axis: Axes, counts, *, horizontal: bool = False) -> None:
    """Put integer counts on bars when there are few enough to read."""
    values = [int(value) for value in counts]
    if not values or len(values) > 40:
        return
    peak = max(max(values), 1)
    if horizontal:
        axis.set_xlim(0, peak * 1.22)
        for patch, count in zip(axis.patches, values, strict=True):
            axis.text(
                patch.get_width(),
                patch.get_y() + patch.get_height() / 2,
                str(count),
                va="center",
                ha="left",
                fontsize=8,
                clip_on=False,
            )
        return
    axis.set_ylim(0, peak * 1.28)
    for patch, count in zip(axis.patches, values, strict=True):
        axis.text(
            patch.get_x() + patch.get_width() / 2,
            patch.get_height(),
            str(count),
            ha="center",
            va="bottom",
            fontsize=8,
            clip_on=False,
        )


def _contrast_name(control: str) -> str:
    if control == "zero_guide":
        return "zero-guide cells"
    return control or "NT"


def _control_positions(status: np.ndarray) -> np.ndarray:
    return np.flatnonzero(np.isin(np.asarray(status).astype(str), ["nt", "zero_guide"]))


def guides_per_cell_bars(n_guides: pd.Series) -> tuple[Figure, pd.DataFrame]:
    """Bar chart of guide counts. Each bar is labeled with the count and the percent."""
    values = n_guides.dropna().astype(int)
    total = int(values.shape[0])
    table = values.value_counts().sort_index().rename_axis("n_guides").reset_index(name="n_cells")
    table["percent"] = (100.0 * table["n_cells"] / total) if total else 0.0
    figure, axis = plt.subplots(figsize=(max(6, 0.7 * max(len(table), 1)), 5))
    axis.bar(table["n_guides"].astype(str), table["n_cells"], color="#4C78A8")
    ceiling = float(table["n_cells"].max()) if len(table) else 1.0
    axis.set_ylim(0, ceiling * 1.45)
    for bar, count, percent in zip(axis.patches, table["n_cells"], table["percent"], strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(count)}\n{percent:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    axis.set_title("Guides per cell")
    axis.set_xlabel("Guides per cell")
    axis.set_ylabel("Cells")
    return figure, table


def draw_guides_per_cell(store: FigureStore, n_guides: pd.Series) -> None:
    figure, table = guides_per_cell_bars(n_guides)
    store.save(
        figure,
        "guides_per_cell",
        "qc",
        table,
        "guides per cell; bar labels are cell count and percent of the assignment table",
    )


def draw_cells_per_target(store: FigureStore, targets: pd.DataFrame, min_cells: int) -> None:
    """Bars by reference ``target_gene_name`` (ledger ``target_gene``); name under each bar."""
    _draw_cells_per_label_bar(
        store,
        targets,
        min_cells=min_cells,
        figure_name="cells_per_target_bar",
        label_column="target_gene",
        unit_word="targets",
        threshold_key="min_cells_per_target",
        threshold_note=f"eligibility keeps n_cells ≥ {min_cells}",
        method_prefix=(
            "cells per target_gene_name for every reference/assigned target "
            "(incl. zeros; no analysis_eligible / min-cell filter)"
        ),
    )


def draw_cells_per_guide_name(
    store: FigureStore, guides: pd.DataFrame, min_cells: int
) -> None:
    """Bars by reference ``name`` (ledger ``guide_name``); name under each bar."""
    _draw_cells_per_label_bar(
        store,
        guides,
        min_cells=min_cells,
        figure_name="cells_per_guide_name_bar",
        label_column="guide_name",
        unit_word="guides",
        threshold_key="min_cells_per_guide",
        threshold_note="guide-level count reference only; not a cell exclusion gate",
        method_prefix=(
            "cells per feature_barcode_reference name for every reference/assigned "
            "guide (incl. zeros; no analysis_eligible filter)"
        ),
    )


def _draw_cells_per_label_bar(
    store: FigureStore,
    table: pd.DataFrame,
    *,
    min_cells: int,
    figure_name: str,
    label_column: str,
    unit_word: str,
    threshold_key: str,
    threshold_note: str,
    method_prefix: str,
) -> None:
    if table.empty or label_column not in table.columns or "n_cells" not in table.columns:
        store.skip(figure_name, "no_targets" if unit_word == "targets" else "no_guides")
        return
    shown = table.sort_values("n_cells", ascending=False, kind="mergesort")
    names = shown[label_column].astype(str).to_numpy()
    n_cells = shown["n_cells"].astype(float)
    values = n_cells.to_numpy()
    min_obs = float(np.min(values))
    median_cells = float(np.median(values))
    max_obs = float(np.max(values))
    ge_mask = values >= float(min_cells)
    n_ge_threshold = int(np.sum(ge_mask))
    n_cells_ge = int(np.sum(values[ge_mask]))
    total_cells = int(np.sum(values))
    n_items = int(len(shown))
    width = max(8.0, min(48.0, 0.16 * n_items + 2.0))
    height = max(6.0, min(10.0, 4.5 + 0.012 * n_items))
    figure, axis = store._figure(width=width, height=height)
    figure.perturbseq_skip_tight_layout = True
    axis.bar(np.arange(n_items), n_cells, color=store.palette("targeting"))
    axis.axhline(
        min_cells,
        color=store.palette("low_power"),
        linestyle="--",
        linewidth=1.0,
        label=f"{threshold_key} = {min_cells} ({threshold_note})",
    )
    axis.axhline(
        min_obs,
        color=store.palette("nt", "#7F7F7F"),
        linestyle=":",
        linewidth=1.0,
        label=f"min = {min_obs:g}",
    )
    axis.axhline(
        median_cells,
        color="black",
        linestyle="-",
        linewidth=1.0,
        label=f"median = {median_cells:g}",
    )
    axis.axhline(
        max_obs,
        color=store.palette("qc_failed", "#D62728"),
        linestyle="-.",
        linewidth=1.0,
        label=f"max = {max_obs:g}",
    )
    axis.set_ylabel("cells")
    tick_fs = 7 if n_items <= 60 else (5 if n_items <= 150 else 4)
    axis.set_xticks(
        np.arange(n_items), names, rotation=45, ha="right", rotation_mode="anchor", fontsize=tick_fs
    )
    if n_items <= 40:
        _label_count_bars(axis, shown["n_cells"])
    handles, labels = axis.get_legend_handles_labels()
    handles.append(Line2D([], [], linestyle="none", marker="", color="none"))
    labels.append(f"total cells = {total_cells}")
    handles.append(Line2D([], [], linestyle="none", marker="", color="none"))
    labels.append(f"{unit_word} with n_cells ≥ {min_cells}: {n_ge_threshold} / {n_items}")
    handles.append(Line2D([], [], linestyle="none", marker="", color="none"))
    labels.append(f"sum n_cells (those {unit_word}): {n_cells_ge}")
    axis.legend(handles, labels, frameon=False, fontsize=8, loc="upper right")
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.28, top=0.92)
    store.save(
        figure,
        figure_name,
        "qc",
        table,
        f"{method_prefix}; dashed = {threshold_key} "
        f"({threshold_note}); "
        f"dotted/solid/dashdot = min/median/max ({min_obs:g}/{median_cells:g}/{max_obs:g}); "
        f"total cells = {total_cells}; "
        f"n_cells ≥ {min_cells}: {n_ge_threshold}/{n_items}; "
        f"sum n_cells on those {unit_word}: {n_cells_ge}; "
        f"x tick = {label_column}",
    )


def _mito_column(metrics: pd.DataFrame) -> str | None:
    if "mito_frac" in metrics.columns:
        return "mito_frac"
    if "pct_counts_mito" in metrics.columns:
        return "pct_counts_mito"
    return None


def _as_percent(values: pd.Series) -> tuple[pd.Series, str]:
    """Map a 0–1 fraction onto 0–100. Values already above 1 stay as they are."""
    numeric = values.astype(float)
    finite = numeric[np.isfinite(numeric)]
    if finite.empty or float(finite.max()) <= 1:
        return numeric * 100, "fraction_times_100"
    return numeric, "already_percent"


def percent_mt_values(metrics: pd.DataFrame) -> tuple[pd.Series, str]:
    """Return mitochondrial percentage on a 0–100 scale."""
    column = _mito_column(metrics)
    if column is None:
        raise KeyError("mito_frac")
    return _as_percent(metrics[column])


def _style_library_d(axis: Axes) -> None:
    """Match Library-D histograms: white face, light grid, full black spines."""
    axis.set_facecolor("white")
    axis.grid(True, which="major", color=_GRID_MAJOR, linewidth=0.7, zorder=0)
    axis.grid(True, which="minor", color=_GRID_MINOR, linewidth=0.45, zorder=0)
    axis.minorticks_on()
    axis.set_axisbelow(True)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    axis.tick_params(colors="black", labelsize=9)


def _hist_bins(store: FigureStore) -> int:
    raw = store.plotting.get("hist_bins", 50)
    try:
        bins = int(raw)
    except (TypeError, ValueError):
        bins = 50
    return max(bins, 10)


def _draw_qc_histogram(
    store: FigureStore,
    *,
    name: str,
    values: np.ndarray,
    xlabel: str,
    title: str,
    method: str,
    threshold: float | None = None,
    threshold_label: str = "",
) -> None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        store.skip(name, "missing_fields")
        return
    figure, axis = store._figure(width=8.4, height=5.6)
    bins = _hist_bins(store)
    axis.hist(
        finite,
        bins=bins,
        color=store.palette("hist_fill", _HIST_FILL),
        edgecolor="white",
        linewidth=0.2,
        zorder=2,
    )
    if threshold is not None and np.isfinite(threshold):
        axis.axvline(
            threshold,
            color=_THRESHOLD_LINE,
            linestyle="--",
            linewidth=1.4,
            zorder=3,
            label=threshold_label or f"threshold={threshold:g}",
        )
        axis.legend(frameon=False, fontsize=8, loc="upper right")
    axis.set_title(title, loc="left", fontsize=12)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Cells")
    _style_library_d(axis)
    table = pd.DataFrame({xlabel: finite})
    store.save(figure, name, "qc", table, method)


def _style_seurat_violin_axis(axis: Axes) -> None:
    """Seurat/ggplot2 classic: white face, black spines, no grid on violin panels."""
    axis.set_facecolor("white")
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    axis.tick_params(colors="black", labelsize=9)


def _draw_seurat_violin(
    axis: Axes,
    values: np.ndarray,
    *,
    identity: str,
    feature: str,
    seed: int,
    max_points: int = 8000,
) -> None:
    """Salmon violin + black jitter, matching Seurat pbmc3k VlnPlot."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        axis.set_title(feature)
        axis.set_xlabel("Identity")
        axis.set_ylabel(feature)
        _style_seurat_violin_axis(axis)
        return
    frame = pd.DataFrame({"Identity": identity, feature: finite})
    sns.violinplot(
        data=frame,
        x="Identity",
        y=feature,
        ax=axis,
        color=_SEURAT_VIOLIN_FILL,
        inner=None,
        cut=0,
        linewidth=0.8,
        density_norm="width",
        order=[identity],
    )
    for collection in axis.collections:
        collection.set_edgecolor("black")
        collection.set_linewidth(0.8)
    shown = frame
    if len(frame) > max_points:
        shown = frame.sample(max_points, random_state=seed)
    sns.stripplot(
        data=shown,
        x="Identity",
        y=feature,
        ax=axis,
        color="black",
        size=1.15,
        alpha=0.35,
        jitter=0.28,
        order=[identity],
        legend=False,
    )
    axis.set_title(feature, fontsize=11)
    axis.set_xlabel("Identity")
    axis.set_ylabel(feature)
    plt.setp(axis.get_xticklabels(), rotation=45, ha="right")
    _style_seurat_violin_axis(axis)


def _draw_panel_histogram(
    axis: Axes,
    *,
    values: np.ndarray,
    xlabel: str,
    title: str,
    bins: int,
    fill: str,
    threshold: float | None = None,
    threshold_label: str = "",
) -> None:
    finite = values[np.isfinite(values)]
    if finite.size:
        axis.hist(
            finite,
            bins=bins,
            color=fill,
            edgecolor="white",
            linewidth=0.2,
            zorder=2,
        )
    if threshold is not None and np.isfinite(threshold):
        axis.axvline(
            threshold,
            color=_THRESHOLD_LINE,
            linestyle="--",
            linewidth=1.4,
            zorder=3,
            label=threshold_label or f"threshold={threshold:g}",
        )
        axis.legend(frameon=False, fontsize=8, loc="upper right")
    axis.set_title(title, loc="left", fontsize=11)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Cells")
    _style_library_d(axis)


def draw_qc_metrics_panel(
    store: FigureStore,
    metrics: pd.DataFrame,
    config: dict[str, Any],
    *,
    population: str = "qc metrics",
) -> None:
    """2×3 QC panel: Seurat-style violins (top) + Library-D histograms (bottom).

    Violin order matches Seurat pbmc3k: nFeature_RNA, nCount_RNA, percent.mt.
    Histogram row: nCount_RNA, nFeature_RNA, percent.mt.
    """
    name = _QC_PANEL_NAME
    _remove_legacy_qc_hists(store.output_dir)
    mito = _mito_column(metrics)
    needed = {"total_counts", "n_genes_by_counts"}
    if metrics.empty or mito is None or not needed.issubset(metrics.columns):
        store.skip(name, "missing_fields")
        return
    scaled, note = percent_mt_values(metrics)
    ncount = pd.to_numeric(metrics["total_counts"], errors="coerce").to_numpy(dtype=float)
    nfeature = pd.to_numeric(metrics["n_genes_by_counts"], errors="coerce").to_numpy(
        dtype=float
    )
    percent_mt = scaled.to_numpy(dtype=float)
    if not (
        np.isfinite(ncount).any() and np.isfinite(nfeature).any() and np.isfinite(percent_mt).any()
    ):
        store.skip(name, "missing_fields")
        return
    sample = str(config["experiment"]["sample_id"])
    seed = int(config["runtime"]["seed"])
    max_mito = config.get("qc", {}).get("max_mito_frac")
    mito_line = float(max_mito) * 100.0 if max_mito is not None else None
    bins = _hist_bins(store)
    fill = store.palette("hist_fill", _HIST_FILL)
    figure, axes = plt.subplots(2, 3, figsize=(12.6, 8.8))
    _draw_seurat_violin(
        axes[0, 0], nfeature, identity=sample, feature="nFeature_RNA", seed=seed
    )
    _draw_seurat_violin(
        axes[0, 1], ncount, identity=sample, feature="nCount_RNA", seed=seed
    )
    _draw_seurat_violin(
        axes[0, 2], percent_mt, identity=sample, feature="percent.mt", seed=seed
    )
    _draw_panel_histogram(
        axes[1, 0],
        values=ncount,
        xlabel="nCount_RNA",
        title="nCount_RNA",
        bins=bins,
        fill=fill,
    )
    _draw_panel_histogram(
        axes[1, 1],
        values=nfeature,
        xlabel="nFeature_RNA",
        title="nFeature_RNA",
        bins=bins,
        fill=fill,
    )
    _draw_panel_histogram(
        axes[1, 2],
        values=percent_mt,
        xlabel="percent.mt",
        title="percent.mt",
        bins=bins,
        fill=fill,
        threshold=mito_line,
        threshold_label=f"max_mito_frac={mito_line:g}%" if mito_line is not None else "",
    )
    table = pd.DataFrame(
        {
            "nCount_RNA": ncount,
            "nFeature_RNA": nfeature,
            "percent.mt": percent_mt,
        }
    )
    store.save(
        figure,
        name,
        "qc",
        table,
        (
            "2x3 QC panel; top=Seurat-style VlnPlot (nFeature, nCount, percent.mt); "
            f"bottom=Library-D histograms; mito scale={note}; bins={bins}; {population}"
        ),
    )


def draw_qc_histograms(
    store: FigureStore,
    metrics: pd.DataFrame,
    config: dict[str, Any],
    *,
    population: str = "qc metrics",
) -> None:
    """Backward-compatible alias: write the combined ``qc_metrics_panel`` PNG."""
    draw_qc_metrics_panel(store, metrics, config, population=population)


def _remove_legacy_qc_hists(output_dir: Path) -> None:
    """Drop standalone hist_*.png replaced by qc_metrics_panel."""
    root = Path(output_dir)
    for name in _LEGACY_QC_HIST_NAMES:
        for path in (
            root / "qc" / f"{name}.png",
            root / "qc" / "figures" / f"{name}.png",
            root / "qc" / "tables" / f"{name}.tsv",
            root / "qc" / "tables" / f"{name}.skipped.txt",
            root / "figures" / "plot" / f"{name}.png",
            root / "plot_data" / f"{name}.tsv",
            root / "plot_data" / f"{name}.skipped.txt",
        ):
            if path.is_file():
                path.unlink()


def _remove_legacy_qc_scatters(output_dir: Path) -> None:
    """Drop standalone qc_scatter / scatter_ncount_nfeature_gex replaced by panel."""
    root = Path(output_dir)
    for name in _LEGACY_QC_SCATTER_NAMES:
        for path in (
            root / "qc" / f"{name}.png",
            root / "qc" / "figures" / f"{name}.png",
            root / "qc" / "tables" / f"{name}.tsv",
            root / "qc" / "tables" / f"{name}.skipped.txt",
            root / "figures" / "plot" / f"{name}.png",
            root / "plot_data" / f"{name}.tsv",
            root / "plot_data" / f"{name}.skipped.txt",
        ):
            if path.is_file():
                path.unlink()


def remove_qc_violin_products(output_dir: Path) -> None:
    """Delete retired qc_violin / assignment_status_bar / standalone hist/scatter products."""
    root = Path(output_dir)
    for name in (
        "qc_violin",
        "assignment_status_bar",
        *_LEGACY_QC_HIST_NAMES,
        *_LEGACY_QC_SCATTER_NAMES,
    ):
        for path in (
            root / "qc" / f"{name}.png",
            root / "qc" / "figures" / f"{name}.png",
            root / "qc" / "tables" / f"{name}.tsv",
            root / "qc" / "tables" / f"{name}.skipped.txt",
            root / "figures" / "plot" / f"{name}.png",
            root / "plot_data" / f"{name}.tsv",
            root / "plot_data" / f"{name}.skipped.txt",
        ):
            if path.is_file():
                path.unlink()


def remove_temporary_control_qc(output_dir: Path) -> None:
    """Delete a temporary-control comparison left by an earlier zero-guide run."""
    name = TEMPORARY_CONTROL_QC_NAME
    for path in (
        output_dir / "qc" / f"{name}.png",
        output_dir / "qc" / "figures" / f"{name}.png",
        output_dir / "qc" / "tables" / f"{name}.tsv",
        output_dir / "qc" / "tables" / f"{name}.skipped.txt",
    ):
        if path.is_file():
            path.unlink()


def _flag(values: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).to_numpy(dtype=bool)
    return values.astype(str).str.strip().str.lower().isin(["true", "1"]).to_numpy()


def _temporary_control_masks(ledger: pd.DataFrame) -> tuple[np.ndarray, np.ndarray] | str:
    """QC-pass zero-guide cells versus analysis-eligible targeting cells."""
    qc_columns = {"total_counts", "n_genes_by_counts"}
    if not qc_columns.issubset(ledger.columns) or _mito_column(ledger) is None:
        return "missing_qc_columns"
    identity = {
        "assignment_status",
        "n_guides",
        "qc_pass",
        "perturbation_status",
        "analysis_eligible",
    }
    if not identity.issubset(ledger.columns):
        return "missing_fields"
    control = control_cell_mask(ledger, _TEMPORARY_SPEC) & _flag(ledger["qc_pass"])
    status = ledger["perturbation_status"].astype(str).eq(STATUS_TARGET).to_numpy()
    targeting = status & _flag(ledger["analysis_eligible"])
    reasons = []
    if not bool(np.any(control)):
        reasons.append("no_zero_guide_cells")
    if not bool(np.any(targeting)):
        reasons.append("no_targeting_cells")
    if reasons:
        return ";".join(reasons)
    return control, targeting


def _distribution_row(group: str, metric: str, values: np.ndarray, n_cells: int) -> dict[str, Any]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        median = mad = q1 = q3 = float("nan")
    else:
        median = float(np.median(finite))
        mad = float(np.median(np.abs(finite - median)))
        q1 = float(np.percentile(finite, 25))
        q3 = float(np.percentile(finite, 75))
    return {
        "group": group,
        "metric": metric,
        "n_cells": n_cells,
        "median": median,
        "mad": mad,
        "q1": q1,
        "q3": q3,
    }


def _qc_distribution_summary(shown: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group in (ZERO_GUIDE_QC_GROUP, TARGETING_QC_GROUP):
        block = shown[shown["group"].astype(str) == group]
        n_cells = int(len(block))
        for metric in ("total_counts", "n_genes_by_counts", "percent.mt"):
            values = pd.to_numeric(block[metric], errors="coerce").to_numpy(dtype=float)
            rows.append(_distribution_row(group, metric, values, n_cells))
    return pd.DataFrame(rows)


def _draw_grouped_metric(
    axis: Axes,
    panel: pd.DataFrame,
    column: str,
    order: list[str],
    palette: dict[str, str],
) -> None:
    kwargs: dict[str, Any] = {
        "data": panel,
        "x": "group",
        "y": column,
        "order": order,
        "hue": "group",
        "hue_order": order,
        "palette": palette,
        "ax": axis,
    }
    try:
        sns.boxplot(**kwargs, legend=False, fliersize=2, linewidth=0.8)
    except TypeError:
        sns.boxplot(**kwargs, fliersize=2, linewidth=0.8)
        legend = axis.get_legend()
        if legend is not None:
            legend.remove()
    _style_library_d(axis)


def _temporary_control_qc_figure(store: FigureStore, shown: pd.DataFrame) -> Figure:
    order = [ZERO_GUIDE_QC_GROUP, TARGETING_QC_GROUP]
    palette = {
        ZERO_GUIDE_QC_GROUP: store.palette("no_guide", "#BDBDBD"),
        TARGETING_QC_GROUP: store.palette("targeting", "#1F77B4"),
    }
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    panels = (
        ("total_counts", "UMI (total_counts)"),
        ("n_genes_by_counts", "genes (n_genes_by_counts)"),
        ("percent.mt", "percent.mt"),
    )
    for axis, (column, title) in zip(axes, panels, strict=True):
        numeric = pd.to_numeric(shown[column], errors="coerce")
        panel = shown.loc[np.isfinite(numeric), ["group", column]].copy()
        panel[column] = pd.to_numeric(panel[column], errors="coerce")
        axis.set_title(title)
        axis.set_ylabel(column)
        if not panel.empty:
            _draw_grouped_metric(axis, panel, column, order, palette)
        axis.set_xlabel("")
        axis.tick_params(axis="x", labelrotation=18)
        for label in axis.get_xticklabels():
            label.set_ha("right")
    return figure


def draw_temporary_control_qc(store: FigureStore, ledger: pd.DataFrame) -> None:
    """Compare temporary zero-guide cells with targeting cells. Descriptive only."""
    name = TEMPORARY_CONTROL_QC_NAME
    try:
        masks = _temporary_control_masks(ledger)
        if isinstance(masks, str):
            remove_temporary_control_qc(store.output_dir)
            store.skip(name, masks)
            return
        control, targeting = masks
        mito = _mito_column(ledger)
        if mito is None:
            remove_temporary_control_qc(store.output_dir)
            store.skip(name, "missing_qc_columns")
            return
        keep = control | targeting
        shown = ledger.loc[keep, ["total_counts", "n_genes_by_counts", mito]].copy()
        for column in ("total_counts", "n_genes_by_counts", mito):
            shown[column] = pd.to_numeric(shown[column], errors="coerce")
        shown.insert(0, "group", np.where(control[keep], ZERO_GUIDE_QC_GROUP, TARGETING_QC_GROUP))
        scaled, note = percent_mt_values(shown)
        shown["percent.mt"] = scaled.to_numpy()
        figure = _temporary_control_qc_figure(store, shown)
        store.save(
            figure,
            name,
            "qc",
            _qc_distribution_summary(shown),
            (
                "descriptive QC; zero-guide cells (temporary control) vs targeting; "
                f"percent.mt {note}; no pass/fail threshold"
            ),
        )
    except Exception as exc:
        print(f"  temporary control QC comparison failed: {exc}", flush=True)
        remove_temporary_control_qc(store.output_dir)
        store.skip(name, f"comparison_failed: {exc}")


_SINGLER_REFERENCES = {
    "mouse": "MouseRNAseqData",
    "human": "HumanPrimaryCellAtlasData",
}
SINGLER_IMPLEMENTATION = "R SingleR"

_SINGLER_RSCRIPT = r"""
options(timeout = 600)
# Inherit OMP/OpenBLAS limits from the Python parent (runtime.n_threads).
n_threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")))
if (is.na(n_threads) || n_threads < 1) n_threads <- 1L
Sys.setenv(
  OMP_NUM_THREADS = as.character(n_threads),
  OPENBLAS_NUM_THREADS = as.character(n_threads),
  MKL_NUM_THREADS = as.character(n_threads)
)
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(n_threads)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("参数不成对")
opt <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[[i]]))
  opt[[key]] <- args[[i + 1]]
}
organism <- opt$organism
reference_name <- opt$reference_name
if (is.null(organism) || is.null(reference_name) || is.null(opt$test) || is.null(opt$out)) {
  stop("缺少 SingleR 参数")
}
known <- c(mouse = "MouseRNAseqData", human = "HumanPrimaryCellAtlasData")
if (!(organism %in% names(known))) stop("experiment.organism 只能是 human 或 mouse")

read_matrix <- function(path) {
  frame <- utils::read.delim(path, check.names = FALSE, sep = "\t")
  genes <- as.character(frame[[1]])
  if (anyDuplicated(genes)) stop("表达矩阵基因名重复")
  mat <- as.matrix(frame[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

reference_matrix <- if (is.null(opt$reference_matrix)) "" else opt$reference_matrix
if (nzchar(reference_matrix)) {
  reference_organism <- if (is.null(opt$reference_organism)) "" else opt$reference_organism
  if (reference_organism != organism) {
    stop(sprintf("参考物种 %s 与输入物种 %s 不一致", reference_organism, organism))
  }
  if (reference_name %in% as.character(known) && known[[organism]] != reference_name) {
    stop(sprintf("参考 %s 与输入物种 %s 不一致", reference_name, organism))
  }
  ref <- read_matrix(reference_matrix)
  label_frame <- utils::read.delim(opt$reference_labels, sep = "\t", check.names = FALSE)
  labels <- as.character(label_frame[[1]])
  if (length(labels) != ncol(ref)) stop("参考标签数量和样本数不一致")
} else {
  if (reference_name != known[[organism]]) {
    stop(sprintf("参考 %s 与输入物种 %s 不一致", reference_name, organism))
  }
  ref <- if (organism == "mouse") {
    celldex::MouseRNAseqData()
  } else {
    celldex::HumanPrimaryCellAtlasData()
  }
  labels <- as.character(SummarizedExperiment::colData(ref)$label.main)
  if (anyNA(labels) || any(!nzchar(labels))) stop("参考 label.main 有缺失")
}

test <- read_matrix(opt$test)
pred <- SingleR::SingleR(test = test, ref = ref, labels = labels)
score_mat <- as.matrix(pred$scores)
pruned <- as.character(pred$pruned.labels)
pruned[is.na(pruned)] <- ""
out <- data.frame(
  cluster = rownames(pred),
  label = as.character(pred$labels),
  delta_next = as.numeric(pred$delta.next),
  pruned_label = pruned,
  max_score = apply(score_mat, 1, max),
  singler_version = as.character(packageVersion("SingleR")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
"""


class SingleRReference(NamedTuple):
    """One SingleR reference. ``matrix`` is genes by samples; empty means load celldex."""

    name: str
    organism: str
    matrix: np.ndarray | None = None
    labels: tuple[str, ...] = ()
    features: tuple[str, ...] = ()


def singler_reference_name(organism: str) -> str:
    """celldex reference for this organism. Mouse never maps to the human atlas."""
    key = str(organism).strip().lower()
    try:
        return _SINGLER_REFERENCES[key]
    except KeyError:
        raise ValueError(f"experiment.organism 只能是 human 或 mouse，收到 {organism}") from None


def cluster_annotations(
    adata: AnnData,
    organism: str,
    *,
    reference: SingleRReference | None = None,
) -> tuple[pd.Series, pd.DataFrame]:
    """Annotate each Leiden cluster with SingleR and store ``obs['singler_label']``.

    Expression is the per-cluster mean of ``adata.X``. The reference has to be the
    same species as ``organism``. This does not write ``obs['cell_type']``.
    """
    species = str(organism).strip().lower()
    expected = singler_reference_name(species)
    if reference is None:
        reference = SingleRReference(name=expected, organism=species)
    else:
        ref_species = str(reference.organism).strip().lower()
        if ref_species != species:
            raise ValueError(f"参考物种 {reference.organism} 与输入物种 {organism} 不一致")
        if reference.name in _SINGLER_REFERENCES.values() and reference.name != expected:
            raise ValueError(f"参考 {reference.name} 与输入物种 {organism} 不一致")
    if "leiden" not in adata.obs:
        raise ValueError("缺少 leiden。")
    profiles = _cluster_mean_expression(adata)
    label_field = "label.main" if reference.name in _SINGLER_REFERENCES.values() else "provided"
    predicted = _run_singler(profiles, reference, species)
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    order = _cluster_order(clusters)
    by_cluster = predicted.set_index(predicted["cluster"].astype(str))
    missing = [cluster for cluster in order if cluster not in by_cluster.index]
    if missing:
        raise ValueError(f"SingleR 没有返回簇 {missing[0]}")
    mapping = {cluster: str(by_cluster.loc[cluster, "label"]) for cluster in order}
    adata.obs["singler_label"] = [mapping[cluster] for cluster in clusters]
    rows = []
    for cluster in order:
        record = by_cluster.loc[cluster]
        pruned = str(record["pruned_label"])
        rows.append(
            {
                "leiden": cluster,
                "singler_label": mapping[cluster],
                "pruned_label": pruned if pruned else pd.NA,
                "delta_next": float(record["delta_next"]),
                "max_score": float(record["max_score"]),
                "n_cells": int((clusters == cluster).sum()),
                "reference": reference.name,
                "organism": species,
                "label_field": label_field,
                "implementation": SINGLER_IMPLEMENTATION,
                "singler_version": str(record["singler_version"]),
            }
        )
    labels = pd.Series(adata.obs["singler_label"].to_numpy(), index=adata.obs_names)
    return labels, pd.DataFrame(rows)


def _cluster_order(clusters: np.ndarray) -> list[str]:
    return sorted(set(clusters), key=lambda value: int(value) if str(value).isdigit() else value)


def _cluster_mean_expression(adata: AnnData) -> pd.DataFrame:
    """Genes by Leiden cluster. Values are the mean of ``adata.X`` inside the cluster."""
    from scipy import sparse

    if "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    else:
        symbols = adata.var_names.astype(str).to_numpy()
    keep: list[int] = []
    names: list[str] = []
    seen: set[str] = set()
    for index, gene in enumerate(symbols):
        if gene in seen or gene in {"", "nan", "None"}:
            continue
        seen.add(gene)
        keep.append(index)
        names.append(gene)
    if not keep:
        raise ValueError("矩阵里没有基因符号。")
    if len(keep) == adata.n_vars and keep[-1] == len(keep) - 1:
        chosen = adata.X
    else:
        chosen = adata.X[:, keep]
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    rows = []
    for cluster in _cluster_order(clusters):
        block = chosen[clusters == cluster]
        if sparse.issparse(block):
            values = np.asarray(block.mean(axis=0)).ravel()
        else:
            values = np.asarray(block, dtype=float).mean(axis=0)
        rows.append(values)
    return pd.DataFrame(np.vstack(rows), index=_cluster_order(clusters), columns=names).T


def _run_singler(
    profiles: pd.DataFrame, reference: SingleRReference, organism: str
) -> pd.DataFrame:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise ValueError("没有找到 Rscript，无法运行 SingleR。")
    with tempfile.TemporaryDirectory(prefix="perturbseq-singler-") as folder:
        root = Path(folder)
        test_path = root / "test.tsv"
        out_path = root / "labels.tsv"
        script_path = root / "singler_cluster.R"
        table = profiles.reset_index(names="gene")
        table.to_csv(test_path, sep="\t", index=False)
        script_path.write_text(_SINGLER_RSCRIPT)
        command = [
            rscript,
            str(script_path),
            "--test",
            str(test_path),
            "--out",
            str(out_path),
            "--organism",
            organism,
            "--reference-name",
            reference.name,
        ]
        if reference.matrix is not None:
            matrix = np.asarray(reference.matrix, dtype=float)
            if matrix.ndim != 2:
                raise ValueError("参考矩阵必须是基因×样本。")
            if (
                len(reference.features) != matrix.shape[0]
                or len(reference.labels) != matrix.shape[1]
            ):
                raise ValueError("参考矩阵和基因或标签长度不一致。")
            ref_path = root / "reference.tsv"
            pd.DataFrame(matrix, index=list(reference.features)).reset_index(names="gene").to_csv(
                ref_path, sep="\t", index=False
            )
            label_path = root / "reference_labels.tsv"
            pd.DataFrame({"label": list(reference.labels)}).to_csv(
                label_path, sep="\t", index=False
            )
            command.extend(
                [
                    "--reference-matrix",
                    str(ref_path),
                    "--reference-labels",
                    str(label_path),
                    "--reference-organism",
                    str(reference.organism).strip().lower(),
                ]
            )
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout).strip()
            message = detail.splitlines()[-1] if detail else "SingleR 失败"
            raise ValueError(message)
        predicted = pd.read_csv(out_path, sep="\t")
    if predicted.empty or predicted["cluster"].astype(str).duplicated().any():
        raise ValueError("SingleR 没有返回每个簇一个标签。")
    return predicted


def leiden_resolution(adata: AnnData) -> float | None:
    """Resolution stored by scanpy.tl.leiden, if the object still has it."""
    block = adata.uns.get("leiden")
    if not isinstance(block, dict):
        return None
    params = block.get("params")
    if not isinstance(params, dict) or "resolution" not in params:
        return None
    return float(params["resolution"])


_ANNOTATION_CELL_CYCLE_PANEL = "umap_annotation_cell_cycle_panel"


def _remove_legacy_annotation_umaps(output_dir: Path) -> None:
    """Drop standalone Leiden/SingleR and cell-cycle UMAPs replaced by the panel."""
    umap_roots = (Path(output_dir) / "UMAP", Path(output_dir) / "UMAP" / "figures")
    for name in ("umap_cluster_annotation", "umap_nt_perturbed_cell_cycle"):
        for root in umap_roots:
            for path in (
                root / f"{name}.png",
                root / f"{name}.pdf",
                root / f"{name}.txt",
            ):
                if path.is_file():
                    path.unlink()
        for path in (
            Path(output_dir) / "UMAP" / "tables" / f"{name}.tsv",
            Path(output_dir) / "UMAP" / "tables" / f"{name}.skipped.txt",
            Path(output_dir) / "plot_data" / f"{name}.tsv",
            Path(output_dir) / "plot_data" / f"{name}.skipped.txt",
        ):
            if path.is_file():
                path.unlink()


def draw_cluster_annotation(store: FigureStore, adata: AnnData) -> None:
    """Deprecated standalone name: write the combined annotation + cell-cycle panel."""
    draw_annotation_cell_cycle_panel(store, adata)


def draw_annotation_cell_cycle_panel(
    store: FigureStore,
    adata: AnnData,
    *,
    control_values: tuple[str, ...] | None = None,
    control_title: str | None = None,
) -> None:
    """Two-row UMAP panel: Leiden+SingleR (top) and NT/perturbed/cell-cycle (bottom)."""
    name = _ANNOTATION_CELL_CYCLE_PANEL
    _remove_legacy_annotation_umaps(store.output_dir)
    if not store.singler_enabled:
        store.skip(name, "singler_disabled")
        return
    if "leiden" not in adata.obs or "X_umap" not in adata.obsm:
        store.skip(name, "missing_fields")
        return
    if "phase" not in adata.obs.columns:
        store.skip(name, "missing_cell_cycle_phase")
        return
    organism = str(store.context.get("organism", ""))
    try:
        labels, summary = cluster_annotations(adata, organism)
    except ValueError as exc:
        store.skip(name, str(exc))
        return
    if control_values is None:
        control_values = ("nt",)
    if not control_values:
        store.skip(name, "missing_control_values")
        return
    if control_title is None:
        if set(control_values) == {ZERO_GUIDE_LABEL}:
            control_title = "zero-guide"
        elif set(control_values) <= {"nt", ZERO_GUIDE_LABEL}:
            control_title = "NT / zero-guide"
        else:
            control_title = "NT" if set(control_values) == {"nt"} else " / ".join(control_values)

    resolution = leiden_resolution(adata)
    resolution_text = "not stored" if resolution is None else f"{resolution:g}"
    reference_name = str(summary["reference"].iloc[0])
    label_field = str(summary["label_field"].iloc[0])
    version = str(summary["singler_version"].iloc[0])
    summary = summary.copy()
    summary["leiden_resolution"] = resolution if resolution is not None else pd.NA
    summary["clustering"] = "Leiden igraph"
    summary["annotation_method"] = (
        f"{SINGLER_IMPLEMENTATION} {version} on cluster-mean log-normalized expression; "
        f"reference={reference_name} {label_field}"
    )
    summary["control_title"] = control_title
    summary["control_values"] = "|".join(control_values)

    coordinates = np.asarray(adata.obsm["X_umap"], dtype=float)
    table = pd.DataFrame(
        {
            "UMAP1": coordinates[:, 0],
            "UMAP2": coordinates[:, 1],
            "leiden": adata.obs["leiden"].astype(str).to_numpy(),
            "singler_label": labels.to_numpy(),
            "perturbation_status": adata.obs["perturbation_status"].astype(str).to_numpy()
            if "perturbation_status" in adata.obs
            else np.array([""] * adata.n_obs),
            "phase": adata.obs["phase"].astype(str).to_numpy(),
        }
    )
    x_span = float(np.max(coordinates[:, 0]) - np.min(coordinates[:, 0])) or 1.0
    y_span = float(np.max(coordinates[:, 1]) - np.min(coordinates[:, 1])) or 1.0
    xlim = (
        float(np.min(coordinates[:, 0]) - 0.05 * x_span),
        float(np.max(coordinates[:, 0]) + 0.05 * x_span),
    )
    ylim = (
        float(np.min(coordinates[:, 1]) - 0.05 * y_span),
        float(np.max(coordinates[:, 1]) + 0.05 * y_span),
    )

    figure = plt.figure(figsize=(16.8, 10.8))
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.28)
    top = grid[0].subgridspec(1, 2, wspace=0.28)
    bottom = grid[1].subgridspec(1, 3, wspace=0.35)
    ax_leiden = figure.add_subplot(top[0, 0])
    ax_singler = figure.add_subplot(top[0, 1], sharex=ax_leiden, sharey=ax_leiden)
    ax_pert = figure.add_subplot(bottom[0, 0])
    ax_control = figure.add_subplot(bottom[0, 1])
    ax_cycle = figure.add_subplot(bottom[0, 2])

    sns.scatterplot(
        data=table,
        x="UMAP1",
        y="UMAP2",
        hue="leiden",
        s=4,
        linewidth=0,
        ax=ax_leiden,
        legend=False,
    )
    sns.scatterplot(
        data=table,
        x="UMAP1",
        y="UMAP2",
        hue="singler_label",
        s=4,
        linewidth=0,
        ax=ax_singler,
    )
    ax_leiden.set_title(f"Leiden (igraph)\nresolution = {resolution_text}")
    ax_singler.set_title(f"SingleR ({reference_name})\n{label_field}")
    for axis in (ax_leiden, ax_singler):
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
    for cluster, group in table.groupby("leiden"):
        ax_singler.text(
            group["UMAP1"].median(),
            group["UMAP2"].median(),
            str(cluster),
            fontsize=8,
            ha="center",
            va="center",
        )
    ax_singler.legend(title="singler_label", frameon=False, fontsize=7, loc="best")

    status = table["perturbation_status"]
    phase = table["phase"]
    control_mask = status.isin(control_values).to_numpy()
    perturbed_mask = status.isin(("targeting",)).to_numpy()
    phase_known = phase.isin(["G1", "G2M", "S"]).to_numpy()
    for axis in (ax_pert, ax_control, ax_cycle):
        _prepare_umap_axis(axis, xlim, ylim)
    ax_control.scatter(
        coordinates[control_mask, 0],
        coordinates[control_mask, 1],
        c=_NT_COLOR,
        s=3.0,
        linewidths=0,
        rasterized=True,
    )
    ax_control.set_title(f"Perturbation Status\n{control_title}", fontsize=14, pad=8)
    ax_pert.scatter(
        coordinates[perturbed_mask, 0],
        coordinates[perturbed_mask, 1],
        c=_PERTURBED_COLOR,
        s=3.0,
        linewidths=0,
        rasterized=True,
    )
    ax_pert.set_title("Perturbed", fontsize=14, pad=8)
    for label in ("G1", "S", "G2M"):
        mask = phase_known & phase.eq(label).to_numpy()
        ax_cycle.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            c=_PHASE_COLORS[label],
            s=3.0,
            linewidths=0,
            rasterized=True,
        )
    ax_cycle.set_title("Cell Cycle Phase", fontsize=14, pad=8)
    ax_cycle.legend(
        handles=[_proxy(_PHASE_COLORS[label], label) for label in ("G1", "G2M", "S")],
        frameon=False,
        loc="lower left",
        fontsize=8,
    )
    figure.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.93, hspace=0.30)
    figure.perturbseq_skip_tight_layout = True
    n_cells = int(adata.n_obs)
    scope = str(adata.uns.get("embedding_cell_scope", "preprocess_entry"))
    method = (
        f"2-row UMAP panel on {n_cells} {scope} cells (QC-pass preprocess entry; "
        f"not analysis_eligible-only); top=Leiden+SingleR (resolution={resolution_text}; "
        f"{summary['annotation_method'].iloc[0]}); "
        f"bottom=perturbed / {control_title} / cell-cycle phase; "
        "obs['cell_type'] unchanged"
    )
    store.save(figure, name, "plot", summary, method)


def draw_qc_scatter_panel(
    store: FigureStore,
    metrics: pd.DataFrame,
    seed: int,
    *,
    population: str = "qc metrics",
) -> None:
    """1×2 QC scatter panel: log10(nCount) vs percent.mt and log10(nCount) vs log10(nFeature)."""
    name = _QC_SCATTER_PANEL_NAME
    _remove_legacy_qc_scatters(store.output_dir)
    mito = _mito_column(metrics)
    needed = {"total_counts", "n_genes_by_counts"}
    if metrics.empty or mito is None or not needed.issubset(metrics.columns):
        store.skip(name, "missing_fields")
        return
    limit = int(store.plotting.get("scatter_max_points", 50000))
    shown = metrics
    if len(metrics) > limit:
        shown = metrics.sample(limit, random_state=seed)
    scaled_all, note = percent_mt_values(metrics)
    scaled = scaled_all.loc[shown.index]
    ncount = pd.to_numeric(shown["total_counts"], errors="coerce").to_numpy(dtype=float)
    nfeature = pd.to_numeric(shown["n_genes_by_counts"], errors="coerce").to_numpy(dtype=float)
    percent_mt = scaled.to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ncount = np.where(np.isfinite(ncount) & (ncount > 0), np.log10(ncount), np.nan)
        log_nfeature = np.where(
            np.isfinite(nfeature) & (nfeature > 0), np.log10(nfeature), np.nan
        )
    fill = store.palette("hist_fill", _HIST_FILL)
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 5.0))
    axes[0].scatter(
        log_ncount,
        percent_mt,
        s=8,
        linewidths=0,
        c=fill,
        alpha=0.55,
        rasterized=True,
        zorder=2,
    )
    axes[0].set_xlabel("log10(nCount_RNA)")
    axes[0].set_ylabel("percent.mt")
    axes[0].set_title("log10(nCount) vs percent.mt", loc="left", fontsize=11)
    axes[1].scatter(
        log_ncount,
        log_nfeature,
        s=8,
        linewidths=0,
        c=fill,
        alpha=0.55,
        rasterized=True,
        zorder=2,
    )
    axes[1].set_xlabel("log10(nCount_RNA)")
    axes[1].set_ylabel("log10(nFeature_RNA)")
    axes[1].set_title("log10(nCount) vs log10(nFeature)", loc="left", fontsize=11)
    for axis in axes:
        _style_library_d(axis)
    table = pd.DataFrame(
        {
            "nCount_RNA": ncount,
            "nFeature_RNA": nfeature,
            "log10_nCount_RNA": log_ncount,
            "log10_nFeature_RNA": log_nfeature,
            "percent.mt": percent_mt,
        }
    )
    store.save(
        figure,
        name,
        "qc",
        table,
        (
            "1x2 QC scatter panel; log10(nCount) vs percent.mt / log10(nFeature) "
            f"({note}); {population}; seed={seed}"
        ),
    )


def _gene_symbol_series(adata: AnnData) -> pd.Series:
    if "gene_symbol" in adata.var.columns:
        return adata.var["gene_symbol"].astype(str)
    return pd.Series(adata.var_names.astype(str), index=adata.var_names)


def _resolve_marker_indices(
    adata: AnnData, markers: list[str]
) -> tuple[list[str], list[int], list[str]]:
    """Map requested markers to var indices (case-insensitive symbol match)."""
    symbols = _gene_symbol_series(adata)
    folded = {str(symbol).casefold(): index for index, symbol in enumerate(symbols)}
    found: list[str] = []
    indices: list[int] = []
    missing: list[str] = []
    for marker in markers:
        key = str(marker).casefold()
        if key in folded:
            found.append(str(symbols.iloc[folded[key]]))
            indices.append(int(folded[key]))
        else:
            missing.append(str(marker))
    return found, indices, missing


def _cluster_annotation_column(adata: AnnData) -> str | None:
    """Prefer SingleR, then ScType, then cell_type, for ``cluster: annotation`` rows."""
    skip = {"ALL", "nan", "None", "<NA>", "", "Unknown"}
    for column in ("singler_label", "sctype_classification", "cell_type"):
        if column not in adata.obs.columns:
            continue
        values = adata.obs[column].astype(str)
        if set(values.unique()) - skip:
            return column
    return None


def _cluster_display_labels(adata: AnnData) -> pd.Series:
    """Leiden id, optionally ``{leiden}: {annotation}`` (SingleR / ScType / cell_type)."""
    clusters = adata.obs["leiden"].astype(str)
    column = _cluster_annotation_column(adata)
    if column is None:
        return clusters
    types = adata.obs[column].astype(str)
    skip = {"ALL", "nan", "None", "<NA>", ""}
    mapping: dict[str, str] = {}
    frame = adata.obs.assign(_cluster=clusters, _annot=types)
    for cluster, block in frame.groupby("_cluster", sort=False):
        majority = block["_annot"].value_counts()
        majority = majority[~majority.index.isin(skip)]
        if majority.empty:
            mapping[str(cluster)] = str(cluster)
        else:
            mapping[str(cluster)] = f"{cluster}: {majority.index[0]}"
    return clusters.map(mapping).astype(str)


def _marker_expression_block(
    adata: AnnData, indices: list[int]
) -> tuple[np.ndarray, np.ndarray]:
    """Return (detection 0/1, expression) matrices shaped (n_cells, n_genes)."""
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    block = counts[:, indices]
    if sparse.issparse(block):
        detected = (block > 0).astype(np.float64).toarray()
    else:
        detected = (np.asarray(block) > 0).astype(np.float64)
    expr = adata.X[:, indices]
    if sparse.issparse(expr):
        expression = expr.toarray().astype(np.float64)
    else:
        expression = np.asarray(expr, dtype=np.float64)
    return detected, expression


def _marker_by_cluster_table(
    adata: AnnData, markers: list[str], indices: list[int]
) -> pd.DataFrame:
    labels = _cluster_display_labels(adata).to_numpy()
    detected, expression = _marker_expression_block(adata, indices)
    rows: list[dict[str, Any]] = []
    try:
        order = sorted(
            set(labels),
            key=lambda value: (int(str(value).split(":", 1)[0]), str(value)),
        )
    except ValueError:
        order = sorted(set(labels))
    for cluster in order:
        mask = labels == cluster
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        for gene_i, gene in enumerate(markers):
            frac = float(detected[mask, gene_i].mean())
            mean_expr = float(expression[mask, gene_i].mean())
            rows.append(
                {
                    "cluster": cluster,
                    "gene": gene,
                    "n_cells": n_cells,
                    "fraction_detected": frac,
                    "mean_expression": mean_expr,
                }
            )
    return pd.DataFrame(rows)


def draw_marker_detection_heatmap(store: FigureStore, adata: AnnData, markers: list[str]) -> None:
    """Fraction detected by cluster; Seurat DoHeatmap colors on a white panel."""
    name = "marker_detection_heatmap"
    if "leiden" not in adata.obs:
        store.skip(name, "missing_leiden")
        return
    found, indices, missing = _resolve_marker_indices(adata, markers)
    if not found:
        store.skip(name, "markers_not_in_adata")
        return
    table = _marker_by_cluster_table(adata, found, indices)
    if table.empty:
        store.skip(name, "missing_fields")
        return
    matrix = table.pivot(index="cluster", columns="gene", values="fraction_detected")
    matrix = matrix.reindex(columns=found)
    sample = str(store.context.get("sample", "sample"))
    annot_col = _cluster_annotation_column(adata) or "leiden"
    # Long annotated labels need left margin; widen with gene count.
    label_width = float(max((len(str(label)) for label in matrix.index), default=8))
    height = max(4.5, 0.48 * len(matrix) + 2.2)
    width = max(9.0, 0.55 * len(found) + 3.8 + min(label_width * 0.08, 4.0))
    figure, axis = store._figure(width=width, height=height)
    figure.perturbseq_skip_tight_layout = True
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    image = axis.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap=_MARKER_HEATMAP_CMAP,
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
    )
    axis.set_xticks(np.arange(len(found)))
    axis.set_xticklabels(found, rotation=45, ha="right", fontsize=9, color="black")
    axis.set_yticks(np.arange(len(matrix.index)))
    axis.set_yticklabels(matrix.index.astype(str), fontsize=9, color="black")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.set_title(f"{sample} marker detection by cluster", loc="left", fontsize=12, color="black")
    axis.set_xticks(np.arange(-0.5, len(found), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    axis.grid(which="minor", color="#E0E0E0", linestyle="-", linewidth=1.0, zorder=3)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.tick_params(colors="black", length=0)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.9)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colorbar.set_label("Fraction detected", color="black")
    colorbar.ax.yaxis.set_tick_params(color="black")
    colorbar.ax.set_facecolor("white")
    colorbar.outline.set_edgecolor("black")
    plt.setp(colorbar.ax.yaxis.get_ticklabels(), color="black")
    left = min(0.42, 0.12 + label_width * 0.012)
    figure.subplots_adjust(left=left, right=0.88, bottom=0.18, top=0.90)
    method = (
        f"fraction detected per Leiden cluster; row labels from {annot_col} "
        "as 'leiden: annotation' when available; "
        "cmap=Seurat DoHeatmap (purple→black→yellow); white panel; "
        f"markers={','.join(found)}"
    )
    if missing:
        method += f"; missing={','.join(missing)}"
    store.save(figure, name, "qc", table, method)


def draw_dotplot_markers(store: FigureStore, adata: AnnData, markers: list[str]) -> None:
    """Seurat DotPlot-like markers on a white panel; size=% expressed; color=lightgrey→blue.

    Cluster rows use the same ``leiden: annotation`` labels as
    ``marker_detection_heatmap`` (SingleR / ScType / cell_type when available).
    """
    name = "dotplot_markers"
    if "leiden" not in adata.obs:
        store.skip(name, "missing_leiden")
        return
    found, indices, missing = _resolve_marker_indices(adata, markers)
    if not found:
        store.skip(name, "markers_not_in_adata")
        return
    table = _marker_by_cluster_table(adata, found, indices)
    if table.empty:
        store.skip(name, "missing_fields")
        return
    clusters = list(dict.fromkeys(table["cluster"].astype(str)))
    try:
        clusters = sorted(
            clusters, key=lambda value: (int(str(value).split(":", 1)[0]), str(value))
        )
    except ValueError:
        clusters = sorted(clusters)
    cluster_pos = {cluster: index for index, cluster in enumerate(clusters)}
    gene_pos = {gene: index for index, gene in enumerate(found)}
    annot_col = _cluster_annotation_column(adata) or "leiden"
    label_width = float(max((len(str(label)) for label in clusters), default=8))
    height = max(4.5, 0.42 * len(clusters) + 2.0)
    width = max(9.5, 0.55 * len(found) + 4.5 + min(label_width * 0.08, 4.0))
    figure, axis = plt.subplots(figsize=(width, height))
    figure.perturbseq_skip_tight_layout = True
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    x = table["gene"].map(gene_pos).to_numpy(dtype=float)
    y = table["cluster"].map(cluster_pos).to_numpy(dtype=float)
    frac = table["fraction_detected"].to_numpy(dtype=float)
    # Seurat-like size scale: empty → full circle.
    sizes = 8.0 + 220.0 * np.clip(frac, 0.0, 1.0)
    points = axis.scatter(
        x,
        y,
        s=sizes,
        c=table["mean_expression"].to_numpy(dtype=float),
        cmap=_DOTPLOT_CMAP,
        edgecolors="none",
        linewidths=0.0,
        zorder=3,
    )
    axis.set_xticks(np.arange(len(found)))
    axis.set_xticklabels(found, rotation=45, ha="right", color="black", fontsize=8)
    axis.set_yticks(np.arange(len(clusters)))
    axis.set_yticklabels(clusters, color="black", fontsize=8)
    axis.set_xlim(-0.6, len(found) - 0.4)
    axis.set_ylim(-0.6, len(clusters) - 0.4)
    axis.invert_yaxis()
    axis.tick_params(colors="black", length=0, pad=3)
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    sample = str(store.context.get("sample", "sample"))
    axis.set_title(f"{sample} marker expression by cluster", loc="left", fontsize=12, color="black")
    left = min(0.08 + label_width * 0.006, 0.38)
    figure.subplots_adjust(left=left, right=0.72, bottom=0.18, top=0.90)
    # Seurat order on the right: Percent Expressed (size) above Average Expression (color).
    size_fracs = (0.0, 0.25, 0.5, 0.75, 1.0)
    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#9ECAE1",
            markeredgecolor="none",
            markersize=max(3.0, np.sqrt(8.0 + 220.0 * value) / 2.2),
            label=str(int(value * 100)),
        )
        for value in size_fracs
    ]
    legend = axis.legend(
        handles=size_handles,
        title="Percent Expressed",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        labelcolor="black",
        fontsize=8,
        title_fontsize=9,
        frameon=False,
        borderaxespad=0.0,
        handletextpad=0.6,
    )
    plt.setp(legend.get_title(), color="black")
    axis.add_artist(legend)
    colorbar = figure.colorbar(points, ax=axis, fraction=0.04, pad=0.12, shrink=0.55, anchor=(0.0, 0.15))
    colorbar.set_label("Average Expression", color="black", fontsize=9)
    colorbar.ax.yaxis.set_tick_params(color="black", labelsize=8)
    colorbar.ax.set_facecolor("white")
    colorbar.outline.set_edgecolor("black")
    plt.setp(colorbar.ax.yaxis.get_ticklabels(), color="black")
    method = (
        "Seurat DotPlot-like; white panel; size=Percent Expressed; "
        "color=Average Expression (lightgrey→blue); "
        f"row labels from {annot_col} as 'leiden: annotation' when available; "
        f"markers={','.join(found)}"
    )
    if missing:
        method += f"; missing={','.join(missing)}"
    store.save(figure, name, "qc", table, method)


def draw_umap(
    store: FigureStore,
    adata: AnnData,
    column: str,
    name: str,
    *,
    embedding: str = "X_umap",
) -> None:
    if embedding not in adata.obsm or column not in adata.obs:
        reason = "signature_umap_not_run" if embedding == "X_umap_mixscape" else "missing_fields"
        store.skip(name, reason)
        return
    coordinates = np.asarray(adata.obsm[embedding], dtype=float)
    if not np.isfinite(coordinates).all():
        store.skip(name, "non_finite_coordinates")
        return
    labels = adata.obs[column].astype(str)
    limit = int(store.plotting["max_legend_categories"])
    if labels.nunique() > limit:
        keep = set(labels.value_counts().head(limit - 1).index)
        labels = labels.where(labels.isin(keep), "Other")
    table = pd.DataFrame(
        {"umap_1": coordinates[:, 0], "umap_2": coordinates[:, 1], column: labels.to_numpy()}
    )
    figure, axis = store._figure()
    sns.scatterplot(data=table, x="umap_1", y="umap_2", hue=column, s=6, linewidth=0, ax=axis)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")
    if embedding == "X_umap_mixscape":
        axis.set_title(f"{column} on Mixscape signature UMAP (X_pert); not expression space")
        space_note = "Mixscape signature space (X_umap_mixscape from X_pert)"
    else:
        axis.set_title(f"{column} on expression UMAP; descriptive, not a distance test")
        space_note = f"descriptive expression UMAP from {embedding}"
    store.save(
        figure,
        name,
        "embedding",
        table,
        f"{space_note}; not used as a statistical distance",
    )


def draw_variable_features(store: FigureStore, adata: AnnData) -> None:
    """Seurat VariableFeaturePlot: mean-dispersion with the top variable genes labeled."""
    if "highly_variable" not in adata.var or "means" not in adata.var:
        store.skip("variable_features", "missing_fields")
        return
    dispersion_key = (
        "dispersions_norm" if "dispersions_norm" in adata.var.columns else "dispersions"
    )
    if dispersion_key not in adata.var.columns:
        store.skip("variable_features", "missing_fields")
        return
    genes = (
        adata.var["gene_symbol"].astype(str)
        if "gene_symbol" in adata.var.columns
        else adata.var_names.astype(str)
    )
    table = pd.DataFrame(
        {
            "gene": genes.to_numpy(),
            "mean": adata.var["means"].to_numpy(dtype=float),
            "dispersion": adata.var[dispersion_key].to_numpy(dtype=float),
            "highly_variable": adata.var["highly_variable"].to_numpy(),
        }
    )
    figure, axis = store._figure()
    stable = table[~table["highly_variable"].astype(bool)]
    variable = table[table["highly_variable"].astype(bool)]
    axis.scatter(stable["mean"], stable["dispersion"], s=8, c="#BDBDBD", linewidths=0)
    axis.scatter(variable["mean"], variable["dispersion"], s=8, c="#D62728", linewidths=0)
    labeled = variable.nlargest(10, "dispersion")
    for record in labeled.itertuples(index=False):
        axis.text(record.mean, record.dispersion, str(record.gene), fontsize=7)
    axis.set_xlabel("mean expression")
    axis.set_ylabel(dispersion_key)
    store.save(
        figure,
        "variable_features",
        "embedding",
        table,
        f"VariableFeaturePlot; top {len(labeled)} variable genes labeled",
    )


def draw_pca_elbow(store: FigureStore, adata: AnnData) -> None:
    """Seurat ElbowPlot: percent variance explained by each principal component."""
    variance = adata.uns.get("pca", {}).get("variance_ratio") if "pca" in adata.uns else None
    if variance is None:
        store.skip("pca_elbow", "missing_fields")
        return
    ratio = np.asarray(variance, dtype=float)
    if ratio.size == 0 or not np.isfinite(ratio).all():
        store.skip("pca_elbow", "missing_fields")
        return
    table = pd.DataFrame({"pc": np.arange(1, ratio.size + 1), "variance_explained": ratio})
    figure, axis = store._figure()
    axis.plot(table["pc"], table["variance_explained"] * 100, marker="o")
    axis.set_xlabel("PC")
    axis.set_ylabel("variance explained (%)")
    store.save(figure, "pca_elbow", "embedding", table, "ElbowPlot of PCA variance ratio")


_GLOBAL_COLORS = {
    "KO": "#1F77B4",
    "KD": "#1F77B4",
    "OE": "#1F77B4",
    "PT": "#1F77B4",
    "NP": "#FF7F0E",
    "control": "#2CA02C",
    "ambiguous": "#7F7F7F",
    "not_evaluated": "#9467BD",
}


_UNCLASSIFIED_MIXSCAPE = {"not_evaluated", "skipped_no_de_genes", "nan", "None", "<NA>", ""}


def _global_mixscape_labels(status: pd.Series, perturbation_label: str) -> pd.Series:
    mapping = {
        "nt": "control",
        "non_perturbed": "NP",
        "perturbed": perturbation_label,
        "ambiguous": "ambiguous",
        "not_evaluated": "not_evaluated",
        "skipped_no_de_genes": "skipped_no_de_genes",
    }
    return status.astype(str).map(lambda value: mapping.get(value, value))


def _classified_mixscape_labels(status: pd.Series, perturbation_label: str) -> pd.Series:
    labels = _global_mixscape_labels(status, perturbation_label)
    return labels[~labels.isin(_UNCLASSIFIED_MIXSCAPE)]


def _global_order(labels: pd.Series, perturbation_label: str) -> list[str]:
    preferred = [perturbation_label, "NP", "control"]
    present = list(dict.fromkeys(labels.astype(str)))
    ordered = [label for label in preferred if label in present]
    ordered.extend(label for label in present if label not in ordered)
    return ordered


def draw_mixscape_class_global_umap(
    store: FigureStore, adata: AnnData, perturbation_label: str
) -> None:
    """Signature-space UMAP colored KO / NP / control, matching the reference figure."""
    name = "umap_mixscape_class_global"
    if "mixscape_status" not in adata.obs:
        store.skip(name, "mixscape_not_run")
        return
    if "X_umap_mixscape" not in adata.obsm:
        reason = str(adata.uns.get("mixscape_umap_status", "signature_umap_not_run"))
        if reason == "signature_space":
            reason = "signature_umap_not_run"
        store.skip(name, reason)
        return
    coordinates = np.asarray(adata.obsm["X_umap_mixscape"], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape != (adata.n_obs, 2):
        store.skip(name, "missing_fields")
        return
    labels = _global_mixscape_labels(adata.obs["mixscape_status"], perturbation_label)
    finite = np.isfinite(coordinates).all(axis=1)
    classified = ~labels.isin(_UNCLASSIFIED_MIXSCAPE).to_numpy()
    keep = finite & classified
    if not keep.any():
        store.skip(name, "non_finite_coordinates" if not finite.any() else "no_classified_cells")
        return
    coordinates = coordinates[keep]
    labels = labels.iloc[np.flatnonzero(keep)].reset_index(drop=True)
    order = _global_order(labels, perturbation_label)
    table = pd.DataFrame(
        {
            "UMAP1": coordinates[:, 0],
            "UMAP2": coordinates[:, 1],
            "mixscape_class_global": labels.to_numpy(),
        }
    )
    figure, axis = store._figure()
    palette = {label: _GLOBAL_COLORS.get(label, "#BDBDBD") for label in order}
    sns.scatterplot(
        data=table,
        x="UMAP1",
        y="UMAP2",
        hue="mixscape_class_global",
        hue_order=order,
        palette=palette,
        s=6,
        linewidth=0,
        ax=axis,
    )
    axis.set_title("mixscape_class_global on signature UMAP (X_pert); not expression space")
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    counts = labels.value_counts()
    summary = ", ".join(f"{label}={int(counts.get(label, 0))}" for label in order)
    store.save(
        figure,
        name,
        "mixscape",
        table,
        f"X_umap_mixscape signature space (control-relative X_pert); classes kept: {summary}",
    )


def draw_mixscape_class_counts(store: FigureStore, adata: AnnData, perturbation_label: str) -> None:
    """Cell counts for the global Mixscape classes KO, NP and control."""
    name = "mixscape_class_global_bar"
    if adata is None or "mixscape_status" not in adata.obs:
        store.skip(name, "mixscape_not_run")
        return
    labels = _classified_mixscape_labels(adata.obs["mixscape_status"], perturbation_label)
    if labels.empty:
        store.skip(name, "no_classified_cells")
        return
    order = _global_order(labels, perturbation_label)
    counts = labels.value_counts()
    table = pd.DataFrame(
        {"mixscape_class_global": order, "n_cells": [int(counts.get(label, 0)) for label in order]}
    )
    figure, axis = store._figure()
    colors = [_GLOBAL_COLORS.get(label, "#BDBDBD") for label in order]
    axis.bar(table["mixscape_class_global"], table["n_cells"], color=colors)
    _label_count_bars(axis, table["n_cells"])
    axis.set_ylabel("cells")
    axis.set_xlabel("mixscape_class_global")
    summary = ", ".join(f"{label}={int(counts.get(label, 0))}" for label in order)
    store.save(figure, name, "mixscape", table, f"global class counts; {summary}")


def _row_matrix(adata: AnnData, rows: np.ndarray, genes: np.ndarray, layer: str | None):
    from scipy import sparse

    source = adata.layers[layer] if layer else adata.X
    block = source[rows][:, genes]
    if sparse.issparse(block):
        return np.asarray(block.toarray(), dtype=float)
    return np.asarray(block, dtype=float)


def _group_index(index: np.ndarray, cap: int, seed: int) -> np.ndarray:
    if index.size <= cap:
        return index
    rng = np.random.default_rng(seed)
    chosen = np.sort(rng.choice(index.size, size=cap, replace=False))
    return index[chosen]


def mixscape_heatmap_cell_order(
    np_index: np.ndarray, nt_index: np.ndarray, ko_index: np.ndarray
) -> tuple[np.ndarray, list[str]]:
    """Stack heatmap columns as contiguous NP, NT, then KO blocks."""
    np_index = np.asarray(np_index)
    nt_index = np.asarray(nt_index)
    ko_index = np.asarray(ko_index)
    columns = np.concatenate([np_index, nt_index, ko_index])
    labels = ["NP"] * int(np_index.size) + ["NT"] * int(nt_index.size) + ["KO"] * int(ko_index.size)
    return columns, labels


def gapped_heatmap_matrix(
    values: np.ndarray, labels: list[str], gap: int = 1
) -> tuple[np.ndarray, list[tuple[str, int, int]]]:
    """Insert blank columns between NP, NT and KO. ``values`` is genes by cells."""
    pieces: list[np.ndarray] = []
    spans: list[tuple[str, int, int]] = []
    cursor = 0
    index = 0
    while index < len(labels):
        label = labels[index]
        end = index + 1
        while end < len(labels) and labels[end] == label:
            end += 1
        width = end - index
        spans.append((label, cursor, cursor + width))
        pieces.append(values[:, index:end])
        cursor += width
        index = end
        if index < len(labels):
            pieces.append(np.full((values.shape[0], gap), np.nan))
            cursor += gap
    return np.concatenate(pieces, axis=1), spans


def _heatmap_cmap():
    cmap = plt.get_cmap("plasma_r").copy()
    cmap.set_bad("#FFFFFF")
    return cmap


def _annotate_heatmap_groups(
    axis: Axes,
    spans: list[tuple[str, int, int]],
    n_genes: int,
    display_names: dict[str, str],
) -> None:
    """Label NP / NT / KO (or CRISPR-mode label) inside the axes."""
    axis.set_ylim(n_genes - 0.5, -2.0)
    for key, start, end in spans:
        width = end - start
        axis.barh(
            -0.55,
            width,
            left=start - 0.5,
            height=0.45,
            color="#E6E6E6",
            edgecolor="#333333",
            align="center",
            linewidth=0.4,
        )
        label = str(display_names.get(key, key))
        axis.text(
            (start + end - 1) / 2.0,
            -1.25,
            label,
            ha="center",
            va="center",
            rotation=0,
            color="black",
            fontsize=8,
            clip_on=True,
        )


def _status_labels(status: pd.Series) -> pd.Series:
    return status.astype(str).str.replace("non-perturbed", "non_perturbed", regex=False)


def _score_frame(adata: AnnData, target: str) -> pd.DataFrame | None:
    block = adata.uns.get("mixscape", {}).get(target)
    if not isinstance(block, dict) or not block:
        return None
    frames = [
        frame for frame in block.values() if isinstance(frame, pd.DataFrame) and "pvec" in frame
    ]
    if not frames:
        return None
    scores = pd.concat(frames)
    scores = scores[~scores.index.duplicated(keep="first")]
    return scores


def _cell_level_volcano(
    adata: AnnData, target: str, control: str, lfc_line: float, pval_cutoff: float
) -> pd.DataFrame | None:
    import scanpy as sc

    labels = adata.obs["target_gene"].astype(str)
    mask = labels.isin([target, control]).to_numpy()
    if int(((labels == target) & mask).sum()) < 2 or int((labels == control).sum()) < 2:
        return None
    subset = adata[mask].copy()
    sc.tl.rank_genes_groups(
        subset,
        groupby="target_gene",
        groups=[target],
        reference=control,
        method="wilcoxon",
        use_raw=False,
    )
    result = subset.uns["rank_genes_groups"]
    genes = np.asarray(result["names"][target]).astype(str)
    log_fc = np.asarray(result["logfoldchanges"][target], dtype=float)
    p_value = np.asarray(result["pvals"][target], dtype=float)
    adj_p = np.asarray(result["pvals_adj"][target], dtype=float)
    table = pd.DataFrame(
        {"gene": genes, "log_fc": log_fc, "p_value": p_value, "adj_p_value": adj_p}
    )
    up = (table["log_fc"] >= lfc_line) & (table["adj_p_value"] < pval_cutoff)
    down = (table["log_fc"] <= -lfc_line) & (table["adj_p_value"] < pval_cutoff)
    table["direction"] = np.where(up, "Up", np.where(down, "Down", "not DE"))
    positive = table.loc[table["adj_p_value"] > 0, "adj_p_value"]
    floor = float(positive.min() / 10) if len(positive) else 1e-300
    plotted = table["adj_p_value"].where(table["adj_p_value"] > 0, floor)
    table["neg_log10_adj_p"] = -np.log10(plotted)
    return table


def _order_by_mixscape_prob(
    index: np.ndarray, probability: np.ndarray, cap: int, seed: int
) -> np.ndarray:
    """Downsample then order by Mixscape KO posterior (high → low), like order.by.prob."""
    chosen = _group_index(np.asarray(index, dtype=int), cap, seed)
    if chosen.size == 0:
        return chosen
    finite = np.isfinite(probability[chosen])
    ordered = chosen[finite][np.argsort(-probability[chosen[finite]], kind="mergesort")]
    missing = chosen[~finite]
    return np.concatenate([ordered, missing]) if missing.size else ordered


def _draw_posterior(axis: Axes, adata: AnnData, target: str, perturbation_label: str) -> None:
    status = _status_labels(adata.obs["mixscape_status"])
    targets = adata.obs["target_gene"].astype(str)
    probability = adata.obs["mixscape_class_p"].astype(float)
    prtb = str(perturbation_label or "KO")
    # Short NT / NP / KO ticks avoid overlap under the shared 2x2 title.
    control_name = "zero-guide" if (status == "zero_guide").any() else "NT"
    groups = {
        control_name: probability[status.isin(["nt", "zero_guide"])],
        "NP": probability[(targets == target) & (status == "non_perturbed")],
        prtb: probability[(targets == target) & (status == "perturbed")],
    }
    # Drop empty classes (e.g. no KO) so stripplot still draws NP+NT.
    order = [name for name, values in groups.items() if len(values) > 0]
    long = pd.DataFrame(
        {
            "Identity": np.concatenate([[name] * len(groups[name]) for name in order]),
            "probability": np.concatenate([groups[name].to_numpy() for name in order]),
        }
    )
    sns.stripplot(
        data=long,
        x="Identity",
        y="probability",
        order=order,
        color="black",
        size=2,
        ax=axis,
        jitter=0.25,
    )
    axis.set_ylim(-0.05, 1.05)
    missing = [name for name in groups if name not in order]
    title = "Posterior probability"
    if missing:
        title = f"Posterior probability (no {'/'.join(missing)})"
    axis.set_title(title, fontsize=11, pad=8)
    axis.set_xlabel("")
    axis.set_ylabel("P(perturbed)", fontsize=9, labelpad=10)
    axis.tick_params(axis="x", labelrotation=0, labelsize=9)
    axis.tick_params(axis="y", labelsize=8)
    axis.yaxis.set_major_locator(plt.MaxNLocator(6))


def _draw_perturb_score(
    axis: Axes, adata: AnnData, target: str, perturbation_label: str, seed: int
) -> None:
    scores = _score_frame(adata, target)
    if scores is None:
        raise KeyError(target)
    status = _status_labels(adata.obs["mixscape_status"])
    present = status.reindex(scores.index)
    label = np.where(
        present.eq("perturbed"),
        f"{target} {perturbation_label}",
        np.where(present.eq("non_perturbed"), f"{target} NP", "NT"),
    )
    frame = pd.DataFrame({"perturbation score": scores["pvec"].to_numpy(dtype=float), "mix": label})
    palette = {
        "NT": "#7d7d7d",
        f"{target} NP": "#c9c9c9",
        f"{target} {perturbation_label}": "#EE6A50",
    }
    order = [name for name in palette if (frame["mix"] == name).any()]
    sns.kdeplot(
        data=frame,
        x="perturbation score",
        hue="mix",
        hue_order=order,
        palette=palette,
        common_norm=False,
        fill=True,
        ax=axis,
        legend=True,
        warn_singular=False,
    )
    rng = np.random.default_rng(seed)
    ymax = axis.get_ylim()[1]
    jitter = rng.uniform(-0.04 * ymax, 0.04 * ymax, size=len(frame))
    axis.scatter(
        frame["perturbation score"],
        jitter,
        s=4,
        c=[palette.get(name, "#333333") for name in frame["mix"]],
        linewidths=0,
        alpha=0.45,
    )
    axis.set_xlabel("perturbation score", fontsize=9)
    axis.set_ylabel("Cell density", fontsize=9)
    axis.set_title("Perturbation score", fontsize=11, pad=8)
    axis.tick_params(labelsize=8)
    legend = axis.get_legend()
    if legend is not None:
        legend.set_title("")
        legend.set_frame_on(False)
        for text in legend.get_texts():
            text.set_fontsize(7)


def _draw_volcano_axis(
    axis: Axes,
    table: pd.DataFrame,
    lfc_line: float,
    pval_cutoff: float,
    *,
    target: str = "",
    control: str = "",
    max_labels: int = 3,
) -> None:
    colors = {"not DE": "#BDBDBD", "Up": "#E41A1C", "Down": "#377EB8"}
    for direction, color in colors.items():
        part = table[table["direction"] == direction]
        axis.scatter(
            part["log_fc"], part["neg_log10_adj_p"], s=8, c=color, linewidths=0, label=direction
        )
    axis.axvline(lfc_line, color="black", linestyle="--", linewidth=0.6)
    axis.axvline(-lfc_line, color="black", linestyle="--", linewidth=0.6)
    if 0 < pval_cutoff < 1:
        axis.axhline(-np.log10(pval_cutoff), color="black", linestyle="--", linewidth=0.6)
    de = table[table["direction"] != "not DE"]
    labeled = de.nsmallest(int(max_labels), "adj_p_value") if not de.empty else de
    y_max = float(np.nanmax(table["neg_log10_adj_p"])) if len(table) else 1.0
    if not np.isfinite(y_max) or y_max <= 0:
        y_max = 1.0
    for offset, record in enumerate(labeled.itertuples(index=False)):
        # Alternate left/right so labels do not stack on the same side.
        x_off = 6 if offset % 2 == 0 else -6
        y_off = 8 + 10 * (offset // 2)
        axis.annotate(
            str(record.gene),
            (record.log_fc, record.neg_log10_adj_p),
            textcoords="offset points",
            xytext=(x_off, y_off),
            fontsize=6,
            color="#333333",
            ha="left" if x_off > 0 else "right",
            clip_on=False,
        )
    axis.set_ylim(-0.05 * y_max, y_max * 1.18)
    if target and control:
        axis.set_xlabel(f"logFC ({target} vs {_contrast_name(control)})", fontsize=8)
        axis.set_title("Cell-level Wilcoxon", fontsize=11, pad=8)
    else:
        axis.set_xlabel("log_fc", fontsize=9)
        axis.set_title("Volcano", fontsize=11, pad=8)
    axis.set_ylabel("-log10 adj_p", fontsize=9)
    axis.tick_params(labelsize=8)
    axis.legend(frameon=False, fontsize=7, loc="lower right")


def _resolve_gene_indices(
    adata: AnnData, gene_names: list[str]
) -> tuple[np.ndarray, list[str]]:
    """Map gene symbols / var names to var indices and display labels."""
    var_names = adata.var_names.astype(str).to_numpy()
    symbols = (
        adata.var["gene_symbol"].astype(str).to_numpy()
        if "gene_symbol" in adata.var.columns
        else var_names
    )
    name_to_index = {name: index for index, name in enumerate(var_names)}
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    ordered: list[int] = []
    labels: list[str] = []
    seen: set[int] = set()
    for gene_name in gene_names:
        index = name_to_index.get(gene_name)
        if index is None:
            index = symbol_to_index.get(gene_name)
        if index is None or index in seen:
            continue
        seen.add(index)
        ordered.append(int(index))
        labels.append(str(symbols[index]))
    return np.asarray(ordered, dtype=int), labels


def _mixscape_heatmap_gene_indices(
    adata: AnnData,
    volcano: pd.DataFrame,
    *,
    max_genes: int,
    lfc_line: float,
    pval_cutoff: float,
    balanced: bool = False,
) -> tuple[np.ndarray, list[str]]:
    """Gene pick aligned to Seurat ``MixscapeHeatmap``.

    Seurat calls ``FindMarkers(ident.1=NT, ident.2=target)`` on RNA, so
    ``avg_log2FC`` is NT vs target (= ``-volcano.log_fc``). Vignette uses
    ``balanced=FALSE`` → NT-high genes only (target down), ``p_val < cutoff``,
    top ``max.genes`` by raw p-value. ``balanced=TRUE`` takes up to
    ``max.genes`` from each side.
    """
    if not isinstance(volcano, pd.DataFrame) or volcano.empty:
        return np.asarray([], dtype=int), []
    if "log_fc" not in volcano.columns or "gene" not in volcano.columns:
        return np.asarray([], dtype=int), []
    frame = volcano.copy()
    # Seurat MixscapeHeatmap uses raw p_val; fall back to adj_p if absent.
    p_col = "p_value" if "p_value" in frame.columns else (
        "adj_p_value" if "adj_p_value" in frame.columns else None
    )
    if p_col is None:
        return np.asarray([], dtype=int), []
    frame["_nt_fc"] = -frame["log_fc"].astype(float)
    frame["_p"] = frame[p_col].astype(float)
    frame = frame.sort_values("_p", ascending=True, kind="mergesort")
    pos = frame.loc[
        (frame["_nt_fc"] > float(lfc_line)) & (frame["_p"] < float(pval_cutoff)), "gene"
    ].astype(str)
    neg = frame.loc[
        (frame["_nt_fc"] < -float(lfc_line)) & (frame["_p"] < float(pval_cutoff)), "gene"
    ].astype(str)
    names: list[str] = []
    if balanced:
        names.extend(pos.head(int(max_genes)).tolist())
        names.extend(neg.head(int(max_genes)).tolist())
    else:
        # Vignette: balanced=FALSE → only NT-high (KO/KD depleted) genes.
        names.extend(pos.head(int(max_genes)).tolist())
    return _resolve_gene_indices(adata, names)


def _heatmap_panel(
    axis: Axes,
    adata: AnnData,
    target: str,
    perturbation_label: str,
    *,
    volcano: pd.DataFrame,
    seed: int,
    max_genes: int,
    max_cells: int,
    lfc_line: float = 0.25,
    pval_cutoff: float = 0.05,
) -> int:
    """Heatmap panel aligned to Seurat MixscapeHeatmap (RNA, order.by.prob)."""
    status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
    targets = adata.obs["target_gene"].astype(str).to_numpy()
    probability = adata.obs["mixscape_class_p"].to_numpy(dtype=float)
    nt_index = _control_positions(status)
    control_block = "zero-guide cells" if np.any(status == "zero_guide") else "NT"
    ko = np.flatnonzero((targets == target) & (status == "perturbed"))
    np_cells = np.flatnonzero((targets == target) & (status == "non_perturbed"))
    if nt_index.size == 0 or (ko.size == 0 and np_cells.size == 0):
        axis.set_axis_off()
        axis.set_title(f"{target}: missing NP/{control_block}/KO")
        return 0
    genes, display = _mixscape_heatmap_gene_indices(
        adata,
        volcano,
        max_genes=max_genes,
        lfc_line=lfc_line,
        pval_cutoff=pval_cutoff,
        balanced=False,
    )
    if genes.size == 0:
        axis.set_axis_off()
        axis.set_title(f"{target}: no MixscapeHeatmap genes")
        return 0
    # Seurat MixscapeHeatmap: assay=RNA (not PRTB / X_pert); order.by.prob=TRUE.
    columns, group_labels = mixscape_heatmap_cell_order(
        _order_by_mixscape_prob(np_cells, probability, max_cells, seed + 1),
        _order_by_mixscape_prob(nt_index, probability, max_cells, seed + 2),
        _order_by_mixscape_prob(ko, probability, max_cells, seed),
    )
    if control_block != "NT":
        group_labels = [control_block if label == "NT" else label for label in group_labels]
    values = _row_matrix(adata, columns, genes, None).T
    center = values.mean(axis=1, keepdims=True)
    scale = values.std(axis=1, keepdims=True)
    scale[scale == 0] = 1
    scaled = (values - center) / scale
    counts = [group_labels.count(label) for label in ("NP", "NT", "KO") if group_labels.count(label)]
    largest = max(counts) if counts else 1
    gap = max(6, int(round(0.08 * largest)))
    gapped, spans = gapped_heatmap_matrix(scaled, group_labels, gap=gap)
    axis.imshow(
        np.ma.masked_invalid(gapped),
        aspect="auto",
        cmap=_heatmap_cmap(),
        vmin=-2,
        vmax=2,
        interpolation="none",
    )
    axis.set_yticks(np.arange(genes.size))
    axis.set_yticklabels(display, fontsize=6)
    axis.set_xticks([])
    _annotate_heatmap_groups(
        axis,
        spans,
        int(genes.size),
        {"NP": "NP", "NT": "NT", "KO": str(perturbation_label or "KO")},
    )
    title = f"Top {genes.size} DE (RNA, NT-high)"
    if ko.size == 0:
        title = f"Top {genes.size} DE (RNA, NT-high; no KO)"
    axis.set_title(title, fontsize=11, pad=10)
    return int(genes.size)


_DE_GENE_COLUMNS = ["target_gene", "gene", "gene_symbol", "rank", "logfc", "adj_p"]


def _de_gene_frame(adata: AnnData, source: str) -> pd.DataFrame | None:
    """DE genes stored by one method. ``None`` means that method did not store a list.

    Mixscale keeps pertpy's ``uns['mixscale_de_genes']`` mapping (var names, raw-p
    order, already capped). Mixscape's pertpy call does not return gene names.
    """
    if source == "mixscale":
        raw = adata.uns.get("mixscale_de_genes")
    elif source == "mixscape":
        raw = adata.uns.get("mixscape_de_genes")
    else:
        raise ValueError(f"unknown DE source: {source}")
    if raw is None:
        return None
    if isinstance(raw, pd.DataFrame):
        return raw
    if not isinstance(raw, dict):
        return None
    symbols: dict[str, str] = {}
    if "gene_symbol" in adata.var.columns:
        for name, symbol in zip(
            adata.var_names.astype(str), adata.var["gene_symbol"].astype(str), strict=False
        ):
            symbols.setdefault(name, symbol)
    rows = []
    for target, genes in raw.items():
        for rank, gene in enumerate(list(genes)):
            name = str(gene)
            rows.append(
                {
                    "target_gene": str(target),
                    "gene": name,
                    "gene_symbol": symbols.get(name, name),
                    "rank": rank,
                    "logfc": pd.NA,
                    "adj_p": pd.NA,
                }
            )
    return pd.DataFrame(rows, columns=_DE_GENE_COLUMNS)


def target_in_de_table(
    adata: AnnData, *, control: str, source: str = "mixscape"
) -> pd.DataFrame:
    """Targets with a stored DE list, and whether each symbol is exactly in that list."""
    columns = ["target_gene", "target_in_de", "n_de_genes", "logfc", "adj_p"]
    genes = _de_gene_frame(adata, source)
    if not isinstance(genes, pd.DataFrame) or genes.empty or "target_gene" not in genes.columns:
        return pd.DataFrame(columns=columns)
    symbol_col = "gene_symbol" if "gene_symbol" in genes.columns else "gene"
    rows = []
    for target, subset in genes.groupby(genes["target_gene"].astype(str), sort=False):
        target = str(target)
        if target == control or target == "":
            continue
        names = set(subset[symbol_col].astype(str))
        if "gene" in subset.columns:
            names.update(subset["gene"].astype(str))
        matched = subset[
            subset[symbol_col].astype(str).eq(target)
            | subset.get("gene", pd.Series(index=subset.index, dtype=str)).astype(str).eq(target)
        ]
        hit = target in names
        logfc = pd.NA
        adj_p = pd.NA
        if hit and not matched.empty and "logfc" in matched.columns:
            numeric = pd.to_numeric(matched["logfc"], errors="coerce")
            if numeric.notna().any():
                logfc = float(numeric.dropna().iloc[0])
        if hit and not matched.empty and "adj_p" in matched.columns:
            numeric = pd.to_numeric(matched["adj_p"], errors="coerce")
            if numeric.notna().any():
                adj_p = float(numeric.dropna().iloc[0])
        rows.append(
            {
                "target_gene": target,
                "target_in_de": hit,
                "n_de_genes": int(len(subset)),
                "logfc": logfc,
                "adj_p": adj_p,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def export_target_in_de(store: FigureStore, adata: AnnData, *, control: str) -> pd.DataFrame:
    """Copy Mixscape 2x2 figures whose target symbol is in that target's DE list.

    Every passing target stays in ``Mixscape/``. This folder is the subset.
    """
    table = target_in_de_table(adata, control=control, source="mixscape")
    if store.plotting.get("save_plot_data", False):
        data = store.data_dir("target_in_de", "target_in_de")
        data.mkdir(parents=True, exist_ok=True)
        table.to_csv(data / "target_in_de.tsv", sep="\t", index=False)
    destination = store.figure_dir("target_in_de", "target_in_de")
    destination.mkdir(parents=True, exist_ok=True)
    source = store.figure_dir("panel_mixscape_2x2", "mixscape")
    hits = set(table.loc[table["target_in_de"].astype(bool), "target_gene"].astype(str))
    keep = {f"{_safe(target)}_mixscape_2x2" for target in hits}
    for stem in keep:
        src = source / f"{stem}.png"
        if src.exists():
            shutil.copy2(src, destination / src.name)
    for path in destination.glob("*_mixscape_2x2.*"):
        if path.suffix != ".png" or path.stem not in keep:
            path.unlink()
    return table


def export_mixscale_target_in_de(
    store: FigureStore, adata: AnnData, *, control: str
) -> pd.DataFrame:
    """Copy Mixscale DE heatmaps whose target symbol is in that target's DE list.

    Passing-target heatmaps stay in ``mixscale/``. ``mixscale/target_in_de/`` is the subset.
    """
    table = target_in_de_table(adata, control=control, source="mixscale")
    if store.plotting.get("save_plot_data", False):
        data = store.data_dir("target_in_de", "mixscale_target_in_de")
        data.mkdir(parents=True, exist_ok=True)
        table.to_csv(data / "target_in_de.tsv", sep="\t", index=False)
    destination = store.figure_dir("target_in_de", "mixscale_target_in_de")
    destination.mkdir(parents=True, exist_ok=True)
    source = store.figure_dir("panel_mixscale_heatmap", "mixscale")
    hits = set(table.loc[table["target_in_de"].astype(bool), "target_gene"].astype(str))
    keep = {f"{_safe(target)}_mixscale_heatmap" for target in hits}
    for stem in keep:
        src = source / f"{stem}.png"
        if src.exists():
            shutil.copy2(src, destination / src.name)
    for path in destination.glob("*_mixscale_heatmap.*"):
        if path.suffix != ".png" or path.stem not in keep:
            path.unlink()
    return table


def analysis_target_genes(adata: AnnData, control: str) -> list[str]:
    """Unique targeting genes in the analysis object. The control label is excluded."""
    if "target_gene" not in adata.obs:
        return []
    labels = adata.obs["target_gene"].astype(str)
    if "perturbation_status" in adata.obs:
        labels = labels[adata.obs["perturbation_status"].astype(str).eq("targeting")]
    excluded = {str(control), "nt", "zero_guide", "no_guide", "nan", "None", "<NA>", "NA", ""}
    names = {gene for gene in labels if gene not in excluded}
    return sorted(names)


def target_de_count_table(
    adata: AnnData, *, control: str, min_de_genes: int = 5, source: str = "mixscape"
) -> pd.DataFrame:
    """Counts of all targets and targets whose own method kept enough DE genes.

    Mixscape does not return the DE gene list, so its table has two rows: every
    target, and targets present in ``uns['mixscape']`` (pertpy ran the GMM).
    Mixscale uses ``uns['mixscale_de_genes']``.
    """
    targets = analysis_target_genes(adata, control)
    target_set = set(targets)
    label = "Mixscape" if source == "mixscape" else "Mixscale"
    genes = _de_gene_frame(adata, source)
    if source == "mixscape" and genes is None:
        tested: list[str] = []
        block = adata.uns.get("mixscape")
        if isinstance(block, dict):
            tested = [str(gene) for gene in block if str(gene) in target_set]
        return pd.DataFrame(
            {
                "key": ["n_targets", "n_targets_with_de"],
                "label": ["all targets", f"{label} tested (enough DE genes)"],
                "n": [len(targets), len(tested)],
            }
        )
    with_de: list[str] = []
    in_own_list: list[str] = []
    if isinstance(genes, pd.DataFrame) and not genes.empty and "target_gene" in genes.columns:
        counts = genes.groupby(genes["target_gene"].astype(str)).size()
        with_de = [
            str(gene)
            for gene, count in counts.items()
            if int(count) >= int(min_de_genes) and str(gene) in target_set
        ]
        symbol_col = "gene_symbol" if "gene_symbol" in genes.columns else "gene"
        for gene in with_de:
            subset = genes.loc[genes["target_gene"].astype(str).eq(gene)]
            names = set(subset[symbol_col].astype(str))
            if "gene" in subset.columns:
                names.update(subset["gene"].astype(str))
            if gene in names:
                in_own_list.append(gene)
    return pd.DataFrame(
        {
            "key": ["n_targets", "n_targets_with_de", "n_targets_in_own_de"],
            "label": [
                "all targets",
                f"{label} scored DE genes",
                "target in own DE list",
            ],
            "n": [len(targets), len(with_de), len(in_own_list)],
        }
    )


def target_de_count_bars(
    counts: pd.DataFrame, *, source: str = "mixscape"
) -> tuple[Figure, pd.DataFrame]:
    """Three blue bars with integer labels. Order is total, with DE, target-in-DE."""
    label = "Mixscape" if source == "mixscape" else "Mixscale"
    frame = counts.copy().reset_index(drop=True)
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.bar(frame["label"].astype(str), frame["n"].astype(int), color="#4C78A8")
    ceiling = float(frame["n"].max()) if len(frame) else 1.0
    axis.set_ylim(0, max(ceiling, 1.0) * 1.22)
    for bar, count in zip(axis.patches, frame["n"], strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    axis.set_ylabel("n targets")
    axis.set_title(
        f"{label} DE-gene counts (cell-level; not a biological-replicate test)"
    )
    return figure, frame


def draw_target_de_counts(
    store: FigureStore,
    adata: AnnData,
    *,
    control: str,
    min_de_genes: int = 5,
    source: str = "mixscape",
) -> None:
    name = "target_de_counts"
    section = "mixscale" if source == "mixscale" else "plot"
    label = "Mixscape" if source == "mixscape" else "Mixscale"
    if adata is None or "target_gene" not in adata.obs:
        store.skip(name, "missing_fields", section=section)
        return
    if source == "mixscale" and "mixscale_score" not in adata.obs:
        store.skip(name, "mixscale_not_run", section=section)
        return
    if source == "mixscape" and "mixscape_class" not in adata.obs:
        store.skip(name, "mixscape_not_run", section=section)
        return
    frame = target_de_count_table(
        adata, control=control, min_de_genes=min_de_genes, source=source
    )
    figure, frame = target_de_count_bars(frame, source=source)
    store.save(
        figure,
        name,
        section,
        frame,
        f"{label} target counts from that method's own DE result; "
        "not a biological-replicate test; integer labels; blue bars",
    )


def target_in_de_gene_count_table(
    adata: AnnData, *, control: str, min_de_genes: int = 5, source: str = "mixscape"
) -> pd.DataFrame:
    """Targets whose stored DE list contains their own symbol.

    Mixscale lists come from pertpy and are already capped at ``max_de_genes``.
    Mixscape does not return gene names, so this table is empty for that source.
    """
    columns = ["target_gene", "n_de_genes"]
    table = target_in_de_table(adata, control=control, source=source)
    if table.empty:
        return pd.DataFrame(columns=columns)
    excluded = {"", "nt", "zero_guide", "no_guide", "nan", "None", "<NA>", "NA", str(control)}
    keep = table["target_in_de"].astype(bool)
    keep &= table["n_de_genes"].astype(int) >= int(min_de_genes)
    keep &= ~table["target_gene"].astype(str).isin(excluded)
    frame = table.loc[keep, columns].copy()
    frame["n_de_genes"] = frame["n_de_genes"].astype(int)
    ordered = frame.sort_values(
        ["n_de_genes", "target_gene"], ascending=[False, True], kind="mergesort"
    )
    return ordered.reset_index(drop=True)


def target_in_de_gene_count_bars(
    counts: pd.DataFrame, *, source: str = "mixscape"
) -> tuple[Figure, pd.DataFrame]:
    """One bar per qualifying target. Height and the integer label are full n_de_genes."""
    label = "Mixscape" if source == "mixscape" else "Mixscale"
    frame = counts.copy().reset_index(drop=True)
    figure, axis = plt.subplots(figsize=(max(8.0, 0.7 * max(len(frame), 1)), 5.6))
    axis.bar(
        frame["target_gene"].astype(str),
        frame["n_de_genes"].astype(int),
        color="#4C78A8",
    )
    ceiling = float(frame["n_de_genes"].max()) if len(frame) else 1.0
    axis.set_ylim(0, max(ceiling, 1.0) * 1.28)
    for bar, count in zip(axis.patches, frame["n_de_genes"], strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    axis.set_ylabel(f"n DE genes ({label} list)")
    axis.set_title(
        f"Targets with own gene in {label} DE list (not a between-group inference)"
    )
    if len(frame) > 40:
        axis.set_xticks([])
        axis.set_xlabel(f"{len(frame)} targets; names are in the table")
    else:
        axis.set_xlabel("target gene")
        plt.setp(axis.get_xticklabels(), rotation=40, ha="right", fontsize=10)
    legend = axis.get_legend()
    if legend is not None:
        legend.remove()
    return figure, frame


def draw_target_in_de_gene_counts(
    store: FigureStore,
    adata: AnnData,
    *,
    control: str,
    min_de_genes: int = 5,
    source: str = "mixscape",
) -> None:
    name = "target_in_de_gene_counts"
    section = "mixscale" if source == "mixscale" else "plot"
    label = "Mixscape" if source == "mixscape" else "Mixscale"
    if adata is None or "target_gene" not in adata.obs:
        store.skip(name, "missing_fields", section=section)
        return
    if _de_gene_frame(adata, source) is None:
        reason = "de_genes_not_returned" if source == "mixscape" else "mixscale_not_run"
        store.skip(name, reason, section=section)
        return
    frame = target_in_de_gene_count_table(
        adata, control=control, min_de_genes=min_de_genes, source=source
    )
    if frame.empty:
        store.skip(name, "no_targets_in_own_de", section=section)
        return
    figure, frame = target_in_de_gene_count_bars(frame, source=source)
    store.save(
        figure,
        name,
        section,
        frame,
        f"{label} targets whose stored DE list contains their own symbol; "
        "bar height is the stored list length; not a biological-replicate test",
    )


def draw_mixscape_target_panels(
    store: FigureStore,
    adata: AnnData,
    perturbation_label: str,
    *,
    control: str,
    seed: int,
    max_genes: int = 20,
    max_cells: int = 300,
    pval_cutoff: float = 0.05,
    lfc_line: float = 0.1,
) -> None:
    """One PNG per classified target.

    With Mixscape scores: 2×2 (posterior, score, volcano, heatmap).
    Without scores (common when every cell is NP): 1×3 (posterior, volcano, heatmap).
    Heatmap follows Seurat MixscapeHeatmap: RNA expression, NT-high DE genes
    (balanced=FALSE), raw p cutoff, cells ordered by KO posterior within NP|NT|KO.
    """
    if "mixscape_status" not in adata.obs or "mixscape_class_p" not in adata.obs:
        store.skip("mixscape_target_panels", "mixscape_not_run")
        return
    status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
    targets = adata.obs["target_gene"].astype(str).to_numpy()
    n_nt = int(np.isin(status, ["nt", "zero_guide"]).sum())
    lfc_line = float(lfc_line)
    rows = []
    classified = np.isin(status, ["perturbed", "non_perturbed"])
    for target in pd.unique(targets[classified]):
        n_ko = int(((targets == target) & (status == "perturbed")).sum())
        n_np = int(((targets == target) & (status == "non_perturbed")).sum())
        has_score = _score_frame(adata, str(target)) is not None
        record = {
            "target_gene": target,
            "n_ko": n_ko,
            "n_np": n_np,
            "n_nt": n_nt,
            "has_score": has_score,
        }
        # KO and/or score table may be absent (all-NP targets); still plot with NT+(KO|NP).
        if n_nt < 1 or (n_ko < 1 and n_np < 1):
            record.update(status="skipped", reason="missing_nt_or_class")
            rows.append(record)
            continue
        try:
            volcano = _cell_level_volcano(adata, str(target), control, lfc_line, pval_cutoff)
        except Exception as exc:
            record.update(status="skipped", reason=f"volcano_failed:{exc}")
            rows.append(record)
            continue
        if volcano is None:
            record.update(status="skipped", reason="volcano_not_estimable")
            rows.append(record)
            continue
        if has_score:
            figure, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=False)
            figure.perturbseq_skip_tight_layout = True
            posterior_ax, score_ax = axes[0, 0], axes[0, 1]
            volcano_ax, heat_ax = axes[1, 0], axes[1, 1]
            _draw_posterior(posterior_ax, adata, str(target), perturbation_label)
            _draw_perturb_score(score_ax, adata, str(target), perturbation_label, seed)
            layout = "2x2"
            adjust = dict(left=0.10, right=0.96, bottom=0.08, top=0.90, wspace=0.42, hspace=0.52)
        else:
            # No score table: drop the score panel → 1×3.
            figure, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=False)
            figure.perturbseq_skip_tight_layout = True
            posterior_ax, volcano_ax, heat_ax = axes[0], axes[1], axes[2]
            _draw_posterior(posterior_ax, adata, str(target), perturbation_label)
            layout = "1x3"
            adjust = dict(left=0.06, right=0.98, bottom=0.14, top=0.86, wspace=0.35)
        _draw_volcano_axis(
            volcano_ax,
            volcano,
            lfc_line,
            pval_cutoff,
            target=str(target),
            control=str(control),
            max_labels=3,
        )
        n_genes = _heatmap_panel(
            heat_ax,
            adata,
            str(target),
            perturbation_label,
            volcano=volcano,
            seed=seed,
            max_genes=max_genes,
            max_cells=max_cells,
            lfc_line=lfc_line,
            pval_cutoff=pval_cutoff,
        )
        contrast = _contrast_name(str(control))
        title = f"{target} vs {contrast}"
        extras = []
        if n_ko < 1:
            extras.append("no KO cells")
        if not has_score:
            extras.append("no score → 1×3")
        if extras:
            title = f"{title} ({'; '.join(extras)})"
        figure.suptitle(title, fontsize=14, y=0.98)
        figure.subplots_adjust(**adjust)
        note_parts = []
        if n_ko < 1:
            note_parts.append("no KO cells (NP+NT only)")
        if not has_score:
            note_parts.append("no Mixscape score table; layout 1x3")
        note = f"; {'; '.join(note_parts)}" if note_parts else ""
        store.save(
            figure,
            f"{_safe(str(target))}_mixscape_2x2",
            "mixscape",
            volcano,
            (
                f"cell-level Wilcoxon {target} vs {control}; not a biological-replicate test; "
                f"layout={layout}; MixscapeHeatmap-aligned RNA heatmap genes={n_genes} "
                f"(NT-high, raw p, order.by.prob); "
                f"lfc={lfc_line}; pval_cutoff={pval_cutoff}{note}"
            ),
        )
        record.update(
            status="plotted", reason="", n_heatmap_genes=n_genes, layout=layout
        )
        rows.append(record)
    selection = pd.DataFrame(rows)
    if store.plotting.get("save_plot_data", False):
        panels = store.data_dir("mixscape_target_panels", "mixscape")
        panels.mkdir(parents=True, exist_ok=True)
        selection.to_csv(panels / "mixscape_target_panels.tsv", sep="\t", index=False)
    export_target_in_de(store, adata, control=control)
    if selection.empty or (selection["status"] != "plotted").all():
        if selection.empty or not (selection["status"] == "plotted").any():
            store.skip("mixscape_target_panels", "no_nt_or_classified_cells")


def _passed_de_targets(adata: AnnData, control: str) -> set[str] | None:
    """Targets pertpy Mixscale kept a DE list for. None means Mixscale did not run."""
    if "mixscale_score" not in adata.obs and adata.uns.get("mixscale_de_genes") is None:
        return None
    frame = _de_gene_frame(adata, "mixscale")
    if frame is None:
        return None
    if frame.empty or "target_gene" not in frame.columns:
        return set()
    passed = set(frame["target_gene"].astype(str))
    passed.discard(str(control))
    passed.discard("")
    return passed


def _skipped_de_targets(adata: AnnData) -> list[str]:
    stored = adata.uns.get("mixscale_unscored_targets", [])
    return [str(gene) for gene in stored]


def _mixscale_de_gene_indices(
    adata: AnnData, target: str, *, max_genes: int
) -> tuple[np.ndarray, list[str]]:
    """Indices / labels aligned to Mixscale_DoHeatmap feature selection.

    Tutorial: top ``max_genes`` DE by p-value, then reorder those by log2FC
    descending; ``features = c(target, top_DEG)`` so the target row is first
    and does not consume the DE budget.
    """
    genes = _de_gene_frame(adata, "mixscale")
    if not isinstance(genes, pd.DataFrame) or genes.empty:
        return np.asarray([], dtype=int), []
    subset = genes[genes["target_gene"].astype(str) == str(target)].copy()
    if subset.empty:
        return np.asarray([], dtype=int), []
    # pertpy orders Mixscale DE genes by raw p. Reorder by logFC only when stored.
    if "rank" in subset.columns:
        subset = subset.sort_values("rank", ascending=True, kind="mergesort")
    elif "adj_p" in subset.columns and pd.to_numeric(subset["adj_p"], errors="coerce").notna().any():
        subset = subset.sort_values("adj_p", ascending=True, kind="mergesort")
    top = subset.head(int(max_genes)).copy()
    if "logfc" in top.columns and pd.to_numeric(top["logfc"], errors="coerce").notna().any():
        top = top.sort_values("logfc", ascending=False, kind="mergesort")
    var_names = adata.var_names.astype(str).to_numpy()
    symbols = (
        adata.var["gene_symbol"].astype(str).to_numpy()
        if "gene_symbol" in adata.var.columns
        else var_names
    )
    name_to_index = {name: index for index, name in enumerate(var_names)}
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    ordered: list[int] = []
    labels: list[str] = []
    seen: set[int] = set()

    def _push(index: int, label: str) -> None:
        if index in seen:
            return
        seen.add(index)
        ordered.append(index)
        labels.append(label)

    target_index = symbol_to_index.get(str(target))
    if target_index is None:
        target_index = name_to_index.get(str(target))
    if target_index is not None:
        _push(int(target_index), str(target))
    gene_col = "gene" if "gene" in top.columns else "gene_symbol"
    symbol_col = "gene_symbol" if "gene_symbol" in top.columns else gene_col
    for record in top.itertuples(index=False):
        gene_name = str(getattr(record, gene_col))
        symbol = str(getattr(record, symbol_col))
        index = name_to_index.get(gene_name)
        if index is None:
            index = symbol_to_index.get(symbol)
        if index is None:
            continue
        _push(int(index), symbol if symbol else gene_name)
    return np.asarray(ordered, dtype=int), labels


def _mixscale_classified_masks(
    adata: AnnData, target: str, control: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    """NT / NP / KO index arrays for ridge/scatter; ``True`` when Mixscape classes exist."""
    labels = adata.obs["target_gene"].astype(str).to_numpy()
    if "mixscape_status" not in adata.obs:
        nt = np.flatnonzero(labels == str(control))
        ko = np.flatnonzero(labels == str(target))
        return nt, np.asarray([], dtype=int), ko, False
    status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
    nt = _control_positions(status)
    if nt.size == 0:
        nt = np.flatnonzero(labels == str(control))
    np_cells = np.flatnonzero((labels == str(target)) & (status == "non_perturbed"))
    ko = np.flatnonzero((labels == str(target)) & (status == "perturbed"))
    if np_cells.size == 0 and ko.size == 0:
        ko = np.flatnonzero(labels == str(target))
        return nt, np.asarray([], dtype=int), ko, False
    return nt, np_cells, ko, True


def _order_by_mixscale_score(
    index: np.ndarray, scores: np.ndarray, cap: int, seed: int
) -> np.ndarray:
    chosen = _group_index(np.asarray(index, dtype=int), cap, seed)
    if chosen.size == 0:
        return chosen
    finite = np.isfinite(scores[chosen])
    ordered = chosen[finite][np.argsort(scores[chosen[finite]], kind="mergesort")]
    missing = chosen[~finite]
    return np.concatenate([ordered, missing]) if missing.size else ordered


def _mixscale_heatmap_blocks(
    adata: AnnData,
    target: str,
    control: str,
    *,
    seed: int,
    max_cells: int,
) -> tuple[np.ndarray, list[str], str]:
    """NT | all target cells (Mixscale_DoHeatmap); within-block by ``mixscale_score``."""
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").to_numpy(dtype=float)
    labels = adata.obs["target_gene"].astype(str).to_numpy()
    if "mixscape_status" in adata.obs:
        status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
        nt = _control_positions(status)
        if nt.size == 0:
            nt = np.flatnonzero(labels == str(control))
    else:
        nt = np.flatnonzero(labels == str(control))
    # Official Mixscale_DoHeatmap: gene identity, not Mixscape NP/KO.
    prtb = np.flatnonzero(labels == str(target))
    nt_ord = _order_by_mixscale_score(nt, scores, max_cells, seed + 2)
    prtb_ord = _order_by_mixscale_score(prtb, scores, max_cells, seed)
    if nt_ord.size == 0 or prtb_ord.size == 0:
        return np.asarray([], dtype=int), [], "missing_nt_or_perturbed"
    columns = np.concatenate([nt_ord, prtb_ord])
    control_block = "zero-guide" if str(control) == ZERO_GUIDE_LABEL else "NT"
    group_labels = [control_block] * int(nt_ord.size) + [str(target)] * int(prtb_ord.size)
    return columns, group_labels, "mixscale_nt_perturbed"


def _draw_mixscale_de_heatmap(
    axis: Axes,
    adata: AnnData,
    target: str,
    *,
    control: str,
    seed: int,
    max_genes: int,
    max_cells: int,
    perturbation_label: str = "KO",
) -> int:
    """Mixscale_DoHeatmap-aligned panel: RNA, NT|target, score-ordered, target+top DE."""
    del perturbation_label  # grouping is gene identity, not CRISPR KO/NP label
    if "mixscale_score" not in adata.obs or "target_gene" not in adata.obs:
        axis.set_axis_off()
        axis.set_title(f"{target}: mixscale_not_run")
        return 0
    columns, group_labels, mode = _mixscale_heatmap_blocks(
        adata, target, control, seed=seed, max_cells=max_cells
    )
    if columns.size == 0:
        axis.set_axis_off()
        axis.set_title(f"{target}: missing NT or PRTB cells")
        return 0
    genes, gene_labels = _mixscale_de_gene_indices(adata, target, max_genes=max_genes)
    if genes.size == 0:
        axis.set_axis_off()
        axis.set_title(f"{target}: no Mixscale DE genes")
        return 0
    values = _row_matrix(adata, columns, genes, None).T
    center = values.mean(axis=1, keepdims=True)
    scale = values.std(axis=1, keepdims=True)
    scale[scale == 0] = 1
    scaled = (values - center) / scale
    unique_labels = list(dict.fromkeys(group_labels))
    counts = [group_labels.count(key) for key in unique_labels]
    largest = max(counts) if counts else int(columns.size)
    gap = max(6, int(round(0.08 * largest)))
    gapped, spans = gapped_heatmap_matrix(scaled, group_labels, gap=gap)
    axis.imshow(
        np.ma.masked_invalid(gapped),
        aspect="auto",
        cmap=_heatmap_cmap(),
        vmin=-2,
        vmax=2,
        interpolation="none",
    )
    axis.set_yticks(np.arange(genes.size))
    axis.set_yticklabels(gene_labels, fontsize=6)
    axis.set_xticks([])
    _annotate_heatmap_groups(
        axis,
        spans,
        int(genes.size),
        {label: label for label in unique_labels},
    )
    axis.set_title(f"{target} Mixscale_DoHeatmap (n={genes.size})")
    return int(genes.size)


def draw_mixscale_de_heatmaps(
    store: FigureStore,
    adata: AnnData,
    *,
    control: str,
    seed: int,
    max_genes: int = 20,
    max_cells: int = 300,
    perturbation_label: str = "KO",
) -> None:
    """One Mixscale DE heatmap per target pertpy scored; copy target-in-own-DE.

    Aligned to Mixscale_DoHeatmap: NT | all target cells ordered by
    ``mixscale_score``; features = target gene + top ``max_genes`` DE (by p, then
    logFC desc). Ridge/scatter may still use Mixscape NT/NP/KO when available.
    """
    if adata is None or "mixscale_score" not in adata.obs or "target_gene" not in adata.obs:
        store.skip("mixscale_de_heatmaps", "mixscale_not_run", section="mixscale")
        return
    if adata.uns.get("mixscale_de_genes") is None:
        store.skip("mixscale_de_heatmaps", "mixscale_not_run", section="mixscale")
        return
    passed = _passed_de_targets(adata, control)
    if passed is None:
        store.skip("mixscale_de_heatmaps", "mixscale_not_run", section="mixscale")
        return
    if not passed:
        store.skip("mixscale_de_heatmaps", "no_scored_targets", section="mixscale")
        return
    prtb = str(perturbation_label or "KO")
    rows = []
    for target in sorted(passed):
        record = {"target_gene": target, "status": "plotted", "reason": "", "n_heatmap_genes": 0}
        figure, axis = plt.subplots(figsize=(10, 6))
        figure.perturbseq_manual_layout = True
        try:
            n_genes = _draw_mixscale_de_heatmap(
                axis,
                adata,
                str(target),
                control=control,
                seed=seed,
                max_genes=max_genes,
                max_cells=max_cells,
                perturbation_label=prtb,
            )
        except Exception as exc:
            plt.close(figure)
            record.update(status="skipped", reason=f"heatmap_failed:{exc}")
            rows.append(record)
            continue
        if n_genes == 0:
            plt.close(figure)
            record.update(status="skipped", reason="no_genes_or_cells")
            rows.append(record)
            continue
        _, _, mode = _mixscale_heatmap_blocks(
            adata, str(target), control, seed=seed, max_cells=max_cells
        )
        store.save(
            figure,
            f"{_safe(str(target))}_mixscale_heatmap",
            "mixscale",
            pd.DataFrame(
                {
                    "target_gene": [target],
                    "n_heatmap_genes": [n_genes],
                    "control": [control],
                    "grouping": [mode],
                    "perturbation_label": [prtb],
                }
            ),
            (
                f"Mixscale_DoHeatmap-aligned {target} vs {_contrast_name(control)}; "
                f"grouping={mode} (NT|target by gene identity); "
                f"cells ordered by mixscale_score within block; "
                f"features=target+top{max_genes} DE in pertpy raw-p order; "
                f"heatmap genes={n_genes}; not a biological-replicate test"
            ),
        )
        record["n_heatmap_genes"] = n_genes
        record["grouping"] = mode
        rows.append(record)
    selection = pd.DataFrame(rows)
    if store.plotting.get("save_plot_data", False):
        panels = store.data_dir("mixscale_de_heatmaps", "mixscale")
        panels.mkdir(parents=True, exist_ok=True)
        selection.to_csv(panels / "mixscale_de_heatmaps.tsv", sep="\t", index=False)
    export_mixscale_target_in_de(store, adata, control=control)
    if selection.empty or not (selection["status"] == "plotted").any():
        store.skip("mixscale_de_heatmaps", "no_heatmaps_plotted", section="mixscale")


def _mixscale_ridge_identities(
    adata: AnnData, control: str, passed: set[str], perturbation_label: str
) -> tuple[list[str], dict[str, np.ndarray], str]:
    """Ridge row labels and boolean masks; prefer NT / NP / KO when Mixscape exists."""
    labels = adata.obs["target_gene"].astype(str).to_numpy()
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").to_numpy(dtype=float)
    prtb = str(perturbation_label or "KO")
    identities: list[str] = []
    masks: dict[str, np.ndarray] = {}

    def _add(name: str, mask: np.ndarray) -> None:
        values = scores[mask]
        if np.isfinite(values).any():
            identities.append(name)
            masks[name] = mask

    used_mixscape = False
    if "mixscape_status" in adata.obs:
        status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
        control_mask = np.isin(status, ["nt", "zero_guide"])
        if not control_mask.any():
            control_mask = labels == str(control)
        _add("NT" if str(control) != ZERO_GUIDE_LABEL else "zero-guide", control_mask)
        for gene in sorted(passed):
            np_mask = (labels == gene) & (status == "non_perturbed")
            ko_mask = (labels == gene) & (status == "perturbed")
            if np_mask.any() or ko_mask.any():
                used_mixscape = True
                _add(f"{gene} NP", np_mask)
                _add(f"{gene} {prtb}", ko_mask)
            else:
                _add(gene, labels == gene)
    else:
        _add(_contrast_name(control) if control == ZERO_GUIDE_LABEL else str(control), labels == str(control))
        for gene in sorted(passed):
            _add(gene, labels == gene)
    mode = "mixscape_nt_np_ko" if used_mixscape else "mixscale_by_target"
    return identities, masks, mode


def draw_mixscale_ridge(
    store: FigureStore,
    adata: AnnData,
    control: str,
    perturbation_label: str = "KO",
) -> None:
    """Ridge of mixscale_score; NT/NP/KO rows when Mixscape labels exist.

    Missing PNGs are expected when pertpy scores no target
    (common with a temporary ``zero_guide`` control and weak/noisy effects).
    """
    name = "mixscale_score_ridge"
    if adata is None or "mixscale_score" not in adata.obs or "target_gene" not in adata.obs:
        store.skip(name, "mixscale_not_run")
        return
    from scipy.stats import gaussian_kde

    passed = _passed_de_targets(adata, control)
    if passed is None:
        store.skip(name, "mixscale_not_run")
        return
    if not passed:
        store.skip(name, "no_scored_targets")
        return
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").to_numpy(dtype=float)
    identities, masks, mode = _mixscale_ridge_identities(
        adata, control, passed, perturbation_label
    )
    if len(identities) < 2:
        store.skip(name, "no_targets_with_finite_scores")
        return
    finite_values = np.concatenate(
        [scores[masks[name]][np.isfinite(scores[masks[name]])] for name in identities]
    )
    if finite_values.size == 0:
        store.skip(name, "no_targets_with_finite_scores")
        return
    grid = np.linspace(float(finite_values.min()) - 0.5, float(finite_values.max()) + 0.5, 200)
    figure, axis = plt.subplots(figsize=(8, max(4, 0.32 * len(identities))))
    figure.perturbseq_manual_layout = True
    cmap = plt.get_cmap("turbo")
    rows = []
    for index, identity in enumerate(identities):
        values = scores[masks[identity]]
        values = values[np.isfinite(values)]
        color = cmap(index / max(len(identities) - 1, 1))
        if values.size >= 2 and float(np.std(values)) > 0:
            density = gaussian_kde(values)(grid)
            density = density / density.max()
            axis.fill_between(grid, index, index + density * 0.85, color=color, linewidth=0)
            axis.plot(grid, index + density * 0.85, color="black", linewidth=0.4)
        elif values.size:
            axis.plot([values[0], values[0]], [index, index + 0.8], color=color, linewidth=2)
        rows.append(
            {
                "identity": identity,
                "n_cells": int(values.size),
                "median_score": float(np.median(values)) if values.size else np.nan,
            }
        )
    axis.set_yticks(range(len(identities)))
    axis.set_yticklabels(identities, fontsize=7)
    axis.set_xlabel("mixscale_score")
    axis.set_ylabel("Identity")
    axis.set_title(f"mixscale_score; grouping={mode}; control = {_contrast_name(control)}")
    axis.set_ylim(-0.2, len(identities))
    table = pd.DataFrame(rows)
    skipped = _skipped_de_targets(adata)
    skipped_text = ",".join(skipped) if skipped else "none"
    store.save(
        figure,
        name,
        "mixscale",
        table,
        "Mixscale ridge; prefer NT/NP/"
        f"{perturbation_label} when Mixscape labels exist else one row per target; "
        "score is perturbation strength, not editing success; no legend; "
        f"grouping={mode}; skipped_targets={skipped_text}",
    )


def _mixscale_expression_bins(
    scores: np.ndarray, expression: np.ndarray, nbin: int = 10
) -> pd.DataFrame:
    """Bin means from Mixscale_ScatterPlot.

    Cutoffs follow R ``quantile(pvec, probs = seq(0, 1, 1/nbin))`` (type 7, the
    same rule as ``numpy.quantile(..., method="linear")``). A bin is the closed
    interval between adjacent cutoffs, so a score on a cutoff is in both bins.
    """
    values = np.asarray(scores, dtype=float)
    observed = np.asarray(expression, dtype=float)
    if values.shape != observed.shape:
        raise ValueError("mixscale 分数和表达长度不一致。")
    probs = np.arange(nbin + 1, dtype=float) / nbin
    cutoffs = np.quantile(values, probs, method="linear")
    rows = []
    for index in range(nbin):
        low = float(cutoffs[index])
        high = float(cutoffs[index + 1])
        chosen = (values >= low) & (values <= high)
        mean = float(np.mean(observed[chosen])) if chosen.any() else np.nan
        rows.append(
            {
                "bin": index + 1,
                "mean_expression": mean,
                "n_cells": int(chosen.sum()),
                "score_low": low,
                "score_high": high,
            }
        )
    return pd.DataFrame(rows)


def _log_expression(adata: AnnData, gene: str) -> tuple[np.ndarray | None, str]:
    """One gene from ``adata.X``, the log-normalized matrix (Seurat RNA ``data``).

    Exact ``gene_symbol`` or var name wins. If both miss, one case-insensitive
    ``gene_symbol`` is used. Two symbols that differ only by case stay missing.
    This lookup is only the scatter y-axis. DE-list membership stays exact.
    """
    column, match = _expression_column(adata, gene)
    if column is None:
        return None, "missing"
    block = adata.X[:, column]
    if hasattr(block, "toarray"):
        block = block.toarray()
    return np.ravel(np.asarray(block, dtype=float)), match


def _expression_column(adata: AnnData, gene: str) -> tuple[int | None, str]:
    symbols = None
    if "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str).to_numpy()
        exact = np.flatnonzero(symbols == gene)
        if exact.size:
            return int(exact[0]), "exact"
    names = adata.var_names.astype(str).to_numpy()
    exact_names = np.flatnonzero(names == gene)
    if exact_names.size:
        return int(exact_names[0]), "exact"
    if symbols is None:
        return None, "missing"
    folded = np.array([value.casefold() for value in symbols], dtype=object)
    hits = np.flatnonzero(folded == gene.casefold())
    if hits.size == 1:
        return int(hits[0]), "case_insensitive_symbol"
    return None, "missing"


def _fallback_only(scores: np.ndarray) -> bool:
    finite = scores[np.isfinite(scores)]
    return bool(finite.size) and bool(np.all(np.abs(finite - 1.0) <= 1e-8))


def draw_mixscale_scatter(
    store: FigureStore,
    adata: AnnData,
    control: str,
    nbin: int = 10,
    perturbation_label: str = "KO",
) -> None:
    """Connected scatter of target-gene expression across mixscale score bins.

    One facet per target pertpy Mixscale scored. When Mixscape NP/KO labels exist, each facet
    draws separate NP and KO curves (NT mean as a dashed reference). Without
    Mixscape, one curve per target (Mixscale_ScatterPlot). Fallback score ``1``
    targets are omitted. ``nbin`` defaults to 10.
    """
    name = "mixscale_scatter"
    if adata is None or "mixscale_score" not in adata.obs or "target_gene" not in adata.obs:
        store.skip(name, "mixscale_not_run")
        return
    nbin = int(nbin)
    if nbin < 1:
        store.skip(name, "invalid_nbin")
        return
    passed = _passed_de_targets(adata, control)
    if passed is None:
        store.skip(name, "mixscale_not_run")
        return
    if not passed:
        store.skip(name, "no_scored_targets")
        return
    labels = adata.obs["target_gene"].astype(str).to_numpy()
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").to_numpy(dtype=float)
    control = str(control)
    prtb = str(perturbation_label or "KO")
    status = (
        _status_labels(adata.obs["mixscape_status"]).to_numpy()
        if "mixscape_status" in adata.obs
        else None
    )
    fallback: list[str] = []
    unscored: list[str] = []
    missing: list[str] = []
    case_insensitive: list[str] = []
    frames: list[pd.DataFrame] = []
    nt_refs: dict[str, float] = {}
    for gene in sorted(passed):
        if gene == control:
            continue
        chosen = labels == gene
        gene_scores = scores[chosen]
        if not np.isfinite(gene_scores).any():
            unscored.append(gene)
            continue
        if _fallback_only(gene_scores):
            fallback.append(gene)
            continue
        expression, match = _log_expression(adata, gene)
        if expression is None:
            missing.append(gene)
            continue
        if match == "case_insensitive_symbol":
            case_insensitive.append(gene)
        control_mask = labels == control
        if status is not None:
            control_mask = np.isin(status, ["nt", "zero_guide"]) | control_mask
        if control_mask.any() and np.isfinite(expression[control_mask]).any():
            nt_refs[gene] = float(np.nanmean(expression[control_mask]))
        class_masks: list[tuple[str, np.ndarray]] = []
        if status is not None:
            np_mask = chosen & (status == "non_perturbed") & np.isfinite(scores)
            ko_mask = chosen & (status == "perturbed") & np.isfinite(scores)
            if np_mask.any():
                class_masks.append(("NP", np_mask))
            if ko_mask.any():
                class_masks.append((prtb, ko_mask))
        if not class_masks:
            class_masks.append(("all", chosen & np.isfinite(scores)))
        for class_name, mask in class_masks:
            if int(mask.sum()) < 2:
                continue
            binned = _mixscale_expression_bins(scores[mask], expression[mask], nbin)
            binned.insert(0, "target_gene", gene)
            binned.insert(1, "class", class_name)
            frames.append(binned)
    if not frames:
        store.skip(name, "no_scored_de_targets")
        return
    table = pd.concat(frames, ignore_index=True)
    genes = table["target_gene"].drop_duplicates().tolist()
    ncol = int(np.ceil(np.sqrt(len(genes))))
    nrow = int(np.ceil(len(genes) / ncol))
    figure, axes = plt.subplots(
        nrow, ncol, figsize=(3.4 * ncol, 2.7 * nrow), squeeze=False
    )
    figure.perturbseq_skip_tight_layout = True
    class_colors = {"NP": "#c9c9c9", prtb: "#EE6A50", "all": "#4C78A8"}
    for index, gene in enumerate(genes):
        axis = axes[index // ncol, index % ncol]
        piece = table[table["target_gene"] == gene]
        for class_name, sub in piece.groupby("class", sort=False):
            axis.plot(
                sub["bin"],
                sub["mean_expression"],
                color=class_colors.get(str(class_name), "#4C78A8"),
                marker="o",
                linewidth=1.4,
                markersize=4.5,
                label=str(class_name),
            )
        if gene in nt_refs:
            axis.axhline(nt_refs[gene], color="#7d7d7d", linestyle="--", linewidth=0.8, label="NT")
        axis.set_title(gene, fontsize=10, pad=6)
        axis.set_xlim(0.5, nbin + 0.5)
        axis.tick_params(labelsize=8)
        if piece["class"].nunique() > 1 or gene in nt_refs:
            axis.legend(frameon=False, fontsize=6, loc="best")
    for index in range(len(genes), nrow * ncol):
        axes[index // ncol, index % ncol].axis("off")
    figure.supxlabel("bins for perturbation scores", y=0.02)
    figure.supylabel("expression levels of target", x=0.02)
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.10, top=0.92, wspace=0.35, hspace=0.40)
    failed = [gene for gene in _skipped_de_targets(adata) if gene and gene != control]
    used_mixscape = "class" in table.columns and set(table["class"].astype(str)) - {"all"}
    store.save(
        figure,
        name,
        "mixscale",
        table,
        "Mixscale_ScatterPlot; quantile bins of mixscale_score; "
        f"grouping={'mixscape_np_ko' if used_mixscape else 'by_target'}; "
        "NT mean dashed when available; y is target gene log-normalized expression; "
        f"n_targets={len(genes)}; nbin={nbin}; control={control}; "
        f"case_insensitive_y={','.join(case_insensitive) or 'none'}; "
        f"omitted_failed_de={','.join(sorted(set(failed))) or 'none'}; "
        f"omitted_fallback_score={','.join(fallback) or 'none'}; "
        f"omitted_no_finite_score={','.join(unscored) or 'none'}; "
        f"omitted_missing_gene={','.join(missing) or 'none'}",
    )


def _safe(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._+-" else "_" for character in value
    )


# Cell-cycle gene sets come from scib (Tirosh / Regev lab lists). Human symbols are
# uppercase; mouse symbols are title case. See scib.preprocessing.get_cell_cycle_genes.
_SUPPORTED_CELL_CYCLE_ORGANISMS = frozenset({"human", "mouse"})

_NT_COLOR = "#5C5C5C"
_PERTURBED_COLOR = "#C9A035"
_PHASE_COLORS = {"G1": "#E85D5D", "G2M": "#2EAE46", "S": "#3A6FD8"}


def _uniquify_symbols(symbols: pd.Index) -> tuple[pd.Index, list[str]]:
    seen: dict[str, int] = {}
    names: list[str] = []
    duplicated: list[str] = []
    for raw in symbols.astype(str):
        count = seen.get(raw, 0)
        seen[raw] = count + 1
        if count == 0:
            names.append(raw)
            continue
        names.append(f"{raw}__dup{count}")
        if count == 1:
            duplicated.append(raw)
    return pd.Index(names), duplicated


def _scib_phase_genes(gene_map: pd.DataFrame, phase: str, present: set[str]) -> list[str]:
    """Pick gene_name or gene_id with the largest overlap, matching scib's filter."""
    frame = gene_map.query("phase == @phase")
    best: list[str] = []
    for column in ("gene_name", "gene_id"):
        genes = [str(gene) for gene in frame[column] if str(gene) in present]
        if len(genes) > len(best):
            best = genes
    return best


def assign_cell_cycle(
    adata: AnnData,
    *,
    organism: str,
    symbol_column: str = "gene_symbol",
    random_state: int = 0,
) -> dict[str, Any]:
    """Score G1/S/G2M with scib cell-cycle genes after log-normalization.

    Gene lists come from ``scib.preprocessing.get_cell_cycle_genes`` (Tirosh /
    Regev lab markers packaged by scib). Scoring mirrors
    ``scib.preprocessing.score_cell_cycle`` via ``scanpy.tl.score_genes_cell_cycle``,
    but always uses ``X`` with ``use_raw=False``, does not regress out the cycle,
    and honors ``random_state``. Temporary unique symbols are used for matching,
    then original ``var_names`` are restored. Writes ``S_score``, ``G2M_score``,
    and ``phase``.
    """
    import scanpy as sc
    from scib.preprocessing import get_cell_cycle_genes

    species = str(organism).strip().lower()
    if species not in _SUPPORTED_CELL_CYCLE_ORGANISMS:
        raise ValueError(f"细胞周期只支持 human 或 mouse，收到 {organism}")
    if "log1p" not in adata.uns:
        raise ValueError("细胞周期打分必须在 log-normalization 之后，不能在原始 counts 上计算。")

    gene_map = get_cell_cycle_genes(species)
    original_names = pd.Index(adata.var_names.astype(str))
    if symbol_column in adata.var.columns:
        symbols = pd.Index(adata.var[symbol_column].astype(str))
    else:
        symbols = original_names
    unique_names, duplicated = _uniquify_symbols(symbols)
    present = set(unique_names)
    if "gene_id" in adata.var.columns:
        # Allow Ensembl ids when symbols miss; scoring still needs them on var_names.
        present |= set(adata.var["gene_id"].astype(str))

    s_genes = _scib_phase_genes(gene_map, "S", present)
    g2m_genes = _scib_phase_genes(gene_map, "G2/M", present)
    if not s_genes or not g2m_genes:
        raise ValueError(
            f"cell cycle genes not in adata\n organism: {species}\n"
            f"s_found={len(s_genes)} g2m_found={len(g2m_genes)}"
        )

    score_names = unique_names
    if any(gene not in set(unique_names) for gene in s_genes + g2m_genes):
        if "gene_id" not in adata.var.columns:
            raise ValueError(
                f"cell cycle genes not in adata\n organism: {species}\n"
                f"s_found={len(s_genes)} g2m_found={len(g2m_genes)}"
            )
        score_names, _ = _uniquify_symbols(pd.Index(adata.var["gene_id"].astype(str)))

    gene_pool = score_names
    if "feature_type" in adata.var.columns:
        gex = adata.var["feature_type"].astype(str).eq("Gene Expression").to_numpy()
        gene_pool = score_names[gex]

    try:
        adata.var_names = score_names
        present_names = set(score_names)
        s_query = [gene for gene in s_genes if gene in present_names]
        g2m_query = [gene for gene in g2m_genes if gene in present_names]
        if not s_query or not g2m_query:
            raise ValueError(
                f"cell cycle genes not in adata\n organism: {species}\n"
                f"s_found={len(s_query)} g2m_found={len(g2m_query)}"
            )
        sc.tl.score_genes_cell_cycle(
            adata,
            s_genes=s_query,
            g2m_genes=g2m_query,
            use_raw=False,
            random_state=random_state,
            gene_pool=gene_pool,
        )
    finally:
        adata.var_names = original_names

    phase = adata.obs["phase"].astype(str)
    unexpected = sorted(set(phase.unique()) - {"G1", "G2M", "S"})
    if unexpected:
        raise ValueError(f"细胞周期标签超出 G1/G2M/S: {unexpected}")
    s_set = gene_map.query("phase == 'S'")
    g2m_set = gene_map.query("phase == 'G2/M'")
    return {
        "s_set_size": int(len(s_set)),
        "s_found": len(s_query),
        "s_found_symbols": s_query,
        "s_missing_symbols": [
            str(name) for name in s_set["gene_name"] if str(name) not in present_names
        ],
        "g2m_set_size": int(len(g2m_set)),
        "g2m_found": len(g2m_query),
        "g2m_found_symbols": g2m_query,
        "g2m_missing_symbols": [
            str(name) for name in g2m_set["gene_name"] if str(name) not in present_names
        ],
        "duplicated_symbols_kept_once": [
            gene for gene in duplicated if gene in set(s_query) | set(g2m_query)
        ],
        "phase_counts": {label: int((phase == label).sum()) for label in ("G1", "G2M", "S")},
        "source": "scib.preprocessing.get_cell_cycle_genes",
    }


def score_cell_cycle_after_log1p(
    adata: AnnData,
    *,
    organism: str,
    random_state: int = 0,
) -> str:
    """Score human or mouse G1/S/G2M after log-normalization. Does not change ``X``."""
    species = str(organism).strip().lower()
    if species not in _SUPPORTED_CELL_CYCLE_ORGANISMS:
        return "skipped_unsupported_organism"
    if "phase" in adata.obs.columns:
        return "already_scored"
    if "log1p" not in adata.uns:
        return "skipped_missing_log1p"
    try:
        coverage = assign_cell_cycle(adata, organism=species, random_state=random_state)
    except ImportError:
        return "skipped_scib_unavailable"
    except ValueError as exc:
        message = str(exc)
        if "cell cycle genes not in adata" in message:
            return "skipped_missing_cell_cycle_genes"
        if "log-normalization" in message:
            return "skipped_missing_log1p"
        raise
    adata.uns["cell_cycle_s_found"] = int(coverage["s_found"])
    adata.uns["cell_cycle_g2m_found"] = int(coverage["g2m_found"])
    adata.uns["cell_cycle_source"] = str(coverage.get("source", "scib"))
    return "scored_after_log1p"


def _prepare_umap_axis(axis: Axes, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    axis.set_facecolor("white")
    axis.grid(False)
    axis.set_xlim(xlim)
    axis.set_ylim(ylim)
    axis.set_xlabel("UMAP 1")
    axis.set_ylabel("UMAP 2")
    axis.tick_params(labelsize=9)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)


def _proxy(color: str, label: str) -> Line2D:
    return Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor=color,
        markeredgecolor="none",
        markersize=7,
        label=label,
    )


def draw_nt_perturbed_cell_cycle_umap(
    adata: AnnData,
    output_dir: Path,
    *,
    status_column: str = "perturbation_status",
    phase_column: str = "phase",
    embedding: str = "X_umap",
    control_values: tuple[str, ...] | None = None,
    control_title: str | None = None,
    nt_values: tuple[str, ...] | None = None,
    perturbed_values: tuple[str, ...] = ("targeting",),
    point_size: float = 3.0,
    dpi: int = 300,
) -> dict[str, Any]:
    """Draw Perturbed, control, and cell-cycle panels on one shared UMAP.

    Order is left to right: Perturbed, Perturbation Status (designed NT or
    zero-guide cells), Cell Cycle Phase. ``control_values`` selects the middle
    panel; when omitted, ``nt_values`` is used, then ``("nt",)``. Designed NT
    and temporary zero-guide controls both use this middle panel.
    """
    output_dir = Path(output_dir)
    output_path = output_dir / "UMAP" / "umap_nt_perturbed_cell_cycle.png"
    if output_path.suffix.lower() != ".png":
        raise ValueError(f"只写 PNG: {output_path}")
    if embedding not in adata.obsm:
        raise ValueError(f"缺少嵌入 {embedding}")
    if status_column not in adata.obs:
        raise ValueError(f"缺少列 {status_column}")
    if phase_column not in adata.obs:
        raise ValueError(f"缺少列 {phase_column}")
    if control_values is None:
        control_values = nt_values if nt_values is not None else ("nt",)
    if not control_values:
        raise ValueError("对照面板缺少 control_values。")
    if control_title is None:
        if set(control_values) == {ZERO_GUIDE_LABEL}:
            control_title = "zero-guide cells"
        elif set(control_values) <= {"nt", ZERO_GUIDE_LABEL}:
            control_title = "NT / zero-guide"
        else:
            control_title = "NT" if set(control_values) == {"nt"} else " / ".join(control_values)
    coordinates = np.asarray(adata.obsm[embedding], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or coordinates.shape[0] != adata.n_obs:
        raise ValueError(f"{embedding} 必须是 n_obs x 2")
    if not np.isfinite(coordinates).all():
        raise ValueError(f"{embedding} 含有非有限坐标")

    status = adata.obs[status_column].astype(str)
    phase = adata.obs[phase_column].astype(str)
    control_mask = status.isin(control_values).to_numpy()
    perturbed_mask = status.isin(perturbed_values).to_numpy()
    omitted = ~(control_mask | perturbed_mask)
    phase_known = phase.isin(["G1", "G2M", "S"]).to_numpy()

    x_span = float(np.max(coordinates[:, 0]) - np.min(coordinates[:, 0])) or 1.0
    y_span = float(np.max(coordinates[:, 1]) - np.min(coordinates[:, 1])) or 1.0
    xlim = (
        float(np.min(coordinates[:, 0]) - 0.05 * x_span),
        float(np.max(coordinates[:, 0]) + 0.05 * x_span),
    )
    ylim = (
        float(np.min(coordinates[:, 1]) - 0.05 * y_span),
        float(np.max(coordinates[:, 1]) + 0.05 * y_span),
    )

    figure, axes = plt.subplots(1, 3, figsize=(16.8, 5.2))
    figure.patch.set_facecolor("white")
    perturbed_axis, control_axis, cycle_axis = axes
    for axis in (perturbed_axis, control_axis, cycle_axis):
        _prepare_umap_axis(axis, xlim, ylim)

    control_axis.scatter(
        coordinates[control_mask, 0],
        coordinates[control_mask, 1],
        c=_NT_COLOR,
        s=point_size,
        linewidths=0,
        rasterized=True,
    )
    control_axis.set_title(f"Perturbation Status\n{control_title}", fontsize=14, pad=8)
    perturbed_axis.scatter(
        coordinates[perturbed_mask, 0],
        coordinates[perturbed_mask, 1],
        c=_PERTURBED_COLOR,
        s=point_size,
        linewidths=0,
        rasterized=True,
    )
    perturbed_axis.set_title("Perturbed", fontsize=14, pad=8)
    for label in ("G1", "S", "G2M"):
        mask = phase_known & phase.eq(label).to_numpy()
        cycle_axis.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            c=_PHASE_COLORS[label],
            s=point_size,
            linewidths=0,
            rasterized=True,
        )
    cycle_axis.set_title("Cell Cycle Phase", fontsize=14, pad=8)
    cycle_axis.legend(
        handles=[_proxy(_PHASE_COLORS[label], label) for label in ("G1", "G2M", "S")],
        frameon=False,
        loc="lower left",
        fontsize=8,
    )
    figure.subplots_adjust(left=0.04, right=0.96, bottom=0.08, top=0.90, wspace=0.55)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    if output_path.stat().st_size <= 0:
        raise RuntimeError(f"图文件为空: {output_path}")
    stale_pdf = output_path.with_suffix(".pdf")
    if stale_pdf.exists():
        stale_pdf.unlink()
    omitted_counts = status[omitted].value_counts().astype(int).to_dict()
    summary = {
        "n_cells": int(adata.n_obs),
        "n_control": int(control_mask.sum()),
        "n_nt": int(control_mask.sum()),
        "n_perturbed": int(perturbed_mask.sum()),
        "n_status_omitted": int(omitted.sum()),
        "omitted_status_counts": {str(key): int(value) for key, value in omitted_counts.items()},
        "n_phase_omitted": int((~phase_known).sum()),
        "phase_counts": {label: int((phase == label).sum()) for label in ("G1", "G2M", "S")},
        "status_column": status_column,
        "control_values": list(control_values),
        "control_title": control_title,
        "nt_values": list(control_values),
        "perturbed_values": list(perturbed_values),
        "png_path": str(output_path),
        "panel_order": ["Perturbed", "Perturbation Status", "Cell Cycle Phase"],
    }
    stale_txt = output_path.with_suffix(".txt")
    if stale_txt.exists():
        stale_txt.unlink()
    return summary
