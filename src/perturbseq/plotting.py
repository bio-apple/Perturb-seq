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
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# mRNA and gRNA statistics. Section "plot" is shared with UMAPs, so routing is by name.
_QC_FIGURES = frozenset(
    {
        "qc_violin",
        "qc_scatter",
        "guides_per_cell",
        "assignment_status_bar",
        "cells_per_guide_bar",
        "cells_per_target_bar",
        "multi_guide_composition_bar",
        "target_cell_counts",
        "cell_filter_waterfall",
        "variable_features",
        "pca_elbow",
        "pca_scatter",
    }
)
# Gene-level perturbation figures. Descriptive log2FC lives here too; titles say it is not inference.
_DGE_FIGURES = frozenset(
    {
        "target_de_counts",
        "target_in_de_gene_counts",
        "target_in_de",
        "volcano",
        "mixscape_fraction_bar",
        "mixscape_target_panels",
        "mixscale_score_ridge",
        "mixscale_scatter",
        "mixscape_class_global_bar",
        "de_gene_heatmap",
        "target_expression_effect",
        "guide_consistency",
        "pseudobulk_library_size",
        "pseudobulk_pca",
        "augur_score_bar",
        "enrichment_dotplot",
        "perturbation_summary",
    }
)
_DGE_SECTIONS = frozenset({"mixscape", "mixscale", "target_in_de", "dge"})
# Pairwise guide-pseudobulk distances. The E-test scatter is written separately.
_DISTANCE_FIGURES = frozenset({"perturbation_distance_heatmap", "perturbation_embedding"})


def _output_bucket(name: str, section: str) -> str:
    # Signature UMAPs are saved with section "mixscape"; the image still goes to UMAP/.
    if name.startswith("umap_"):
        return "umap"
    if name in _DISTANCE_FIGURES:
        return "edistance"
    if section in _DGE_SECTIONS or name in _DGE_FIGURES or name.endswith("_mixscape_2x2"):
        return "dge"
    if name in _QC_FIGURES:
        return "qc"
    return "main"


class FigureStore:
    def __init__(self, output_dir: Path, config: dict[str, Any], context: dict[str, str]) -> None:
        self.output_dir = output_dir
        self.plotting = config["plotting"]
        self.context = context
        self.plot_data = output_dir / "plot_data"
        self.plot_data.mkdir(parents=True, exist_ok=True)
        self.manifest: list[dict[str, str]] = []
        sns.set_theme(style=str(self.plotting.get("theme", "whitegrid")))
        matplotlib.rcParams["font.sans-serif"] = [
            "Heiti SC",
            "Arial Unicode MS",
            "DejaVu Sans",
        ]
        matplotlib.rcParams["axes.unicode_minus"] = False

    def figure_dir(self, name: str, section: str) -> Path:
        bucket = _output_bucket(name, section)
        if bucket == "umap":
            return self.output_dir / "UMAP" / "figures"
        if bucket == "qc":
            return self.output_dir / "qc" / "figures"
        if bucket == "dge":
            return self.output_dir / "DGE" / "figures" / (section or "plot")
        if bucket == "edistance":
            return self.output_dir / "E-distance" / "figures"
        return self.output_dir / "figures" / section

    def data_dir(self, name: str, section: str = "") -> Path:
        bucket = _output_bucket(name, section)
        if bucket == "umap":
            return self.output_dir / "UMAP" / "tables"
        if bucket == "qc":
            return self.output_dir / "qc" / "tables"
        if bucket == "dge":
            return self.output_dir / "DGE" / "plot_data"
        if bucket == "edistance":
            return self.output_dir / "E-distance" / "tables"
        return self.plot_data

    def save(
        self, figure: Figure, name: str, section: str, data: pd.DataFrame, method: str
    ) -> None:
        folder = self.figure_dir(name, section)
        folder.mkdir(parents=True, exist_ok=True)
        figure.text(0.01, 0.005, self._caption(method), fontsize=7, ha="left", va="bottom")
        if getattr(figure, "perturbseq_manual_layout", False):
            figure.subplots_adjust(left=0.22, right=0.9, bottom=0.14, top=0.96)
        else:
            figure.tight_layout(rect=(0, 0.05, 1, 1))
        path = folder / f"{name}.png"
        figure.savefig(path, dpi=int(self.plotting["dpi"]), bbox_inches="tight")
        if path.stat().st_size <= 0:
            raise RuntimeError(f"图文件为空: {path}")
        stale_pdf = folder / f"{name}.pdf"
        if stale_pdf.exists():
            stale_pdf.unlink()
        plt.close(figure)
        data_folder = self.data_dir(name, section)
        data_folder.mkdir(parents=True, exist_ok=True)
        stale = data_folder / f"{name}.skipped.txt"
        if stale.exists():
            stale.unlink()
        legacy = self.plot_data / f"{name}.skipped.txt"
        if legacy != stale and legacy.exists():
            legacy.unlink()
        data_path = ""
        if self.plotting.get("save_plot_data", True):
            data_path = str(data_folder / f"{name}.tsv")
            data.to_csv(data_path, sep="\t", index=False)
        bucket = _output_bucket(name, section)
        if bucket == "umap":
            self._retire_old_umap(name)
        self.manifest.append(
            {
                "plot_name": name,
                "branch": "UMAP" if bucket == "umap" else section,
                "status": "written",
                "reason": "",
                "data_path": data_path,
                "png_path": str(path),
                "pdf_path": "",
            }
        )

    def skip(self, name: str, reason: str) -> None:
        data_folder = self.data_dir(name)
        data_folder.mkdir(parents=True, exist_ok=True)
        skipped = data_folder / f"{name}.skipped.txt"
        skipped.write_text(reason.strip() + "\n")
        if name.startswith("umap_"):
            for folder in (self.plot_data, self.output_dir / "DGE" / "plot_data"):
                stale = folder / f"{name}.skipped.txt"
                if stale != skipped and stale.is_file():
                    stale.unlink()
        self.manifest.append(
            {
                "plot_name": name,
                "branch": "UMAP" if name.startswith("umap_") else "",
                "status": "skipped",
                "reason": reason.strip(),
                "data_path": str(skipped),
                "png_path": "",
                "pdf_path": "",
            }
        )

    def _retire_old_umap(self, name: str) -> None:
        """Remove a UMAP product left under the pre-UMAP figure or table paths."""
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
            Path("DGE") / "figures" / "mixscape" / f"{name}.png",
            Path("DGE") / "figures" / "mixscape" / f"{name}.txt",
            Path("DGE") / "figures" / "mixscape" / f"{name}.pdf",
            Path("plot_data") / f"{name}.tsv",
            Path("plot_data") / f"{name}.skipped.txt",
            Path("DGE") / "plot_data" / f"{name}.tsv",
            Path("DGE") / "plot_data" / f"{name}.skipped.txt",
        ]
        for relative in relatives:
            stale = self.output_dir / relative
            if stale.is_file():
                stale.unlink()

    def write_manifest(self, path: Path) -> None:
        pd.DataFrame(self.manifest).to_csv(path, sep="\t", index=False)

    def _caption(self, method: str) -> str:
        parts = [f"{key}={value}" for key, value in self.context.items()]
        parts.append(f"method={method}")
        return "; ".join(parts)

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


def _main_branch(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep the main contrast when a table also carries the sensitivity branch."""
    if "branch" not in frame.columns:
        return frame
    return frame[frame["branch"].astype(str).eq("main")].copy()


def draw_descriptive_heatmap(store: FigureStore, effects: pd.DataFrame) -> None:
    if effects.empty or "log2fc_descriptive" not in effects.columns:
        store.skip("de_gene_heatmap", "missing_fields")
        return
    usable = effects[effects["effect_status"].astype(str) == "descriptive"].copy()
    if usable.empty:
        store.skip("de_gene_heatmap", "no_descriptive_effects")
        return
    values = pd.to_numeric(usable["log2fc_descriptive"], errors="coerce")
    usable = usable.assign(abs_effect=values.abs())
    if "gene_id" in usable.columns:
        usable["gene_key"] = usable["gene_id"].astype(str)
    else:
        usable["gene_key"] = usable["gene_symbol"].astype(str)
    label_by_key = (
        usable.drop_duplicates("gene_key").set_index("gene_key")["gene_symbol"].astype(str)
    )
    scores = usable.groupby("gene_key")["abs_effect"].max().sort_values(ascending=False)
    genes = list(scores.head(int(store.plotting["heatmap_max_genes"])).index)
    pivot = usable.pivot_table(
        index="gene_key", columns="target_gene", values="log2fc_descriptive", aggfunc="mean"
    )
    pivot = pivot.loc[[gene for gene in genes if gene in pivot.index]]
    labels = label_by_key.reindex(pivot.index)
    repeated = labels.value_counts()
    pivot.index = pd.Index(
        [
            f"{name} ({key})" if int(repeated.get(name, 0)) > 1 else str(name)
            for key, name in zip(pivot.index, labels, strict=True)
        ],
        name="gene",
    )
    if pivot.shape[0] < 2:
        store.skip("de_gene_heatmap", "no_descriptive_effects")
        return
    matrix = pivot.to_numpy(dtype=float)
    center = matrix - np.nanmean(matrix, axis=1, keepdims=True)
    scale = np.nanstd(matrix, axis=1, keepdims=True)
    constant = scale[:, 0] == 0
    scale[scale == 0] = 1
    zscore = center / scale
    shown = pd.DataFrame(zscore, index=pivot.index, columns=pivot.columns)
    figure, axis = plt.subplots(
        figsize=(max(8, 0.2 * shown.shape[1]), max(6, 0.2 * shown.shape[0]))
    )
    sns.heatmap(
        shown, cmap="vlag", center=0, xticklabels=shown.shape[1] <= 40, yticklabels=True, ax=axis
    )
    note = "constant rows kept" if constant.any() else "row z-score of descriptive log2FC"
    phrase = _contrast_name(str(store.context.get("control_label", "")))
    axis.set_title(f"descriptive log2FC vs {phrase}, not a replicate test")
    axis.set_xlabel("target")
    axis.set_ylabel("gene")
    long = shown.reset_index(names="gene").melt(
        id_vars="gene", var_name="target_gene", value_name="zscore"
    )
    store.save(figure, "de_gene_heatmap", "de", long, note)


def draw_assignment_status(store: FigureStore, status: pd.Series) -> None:
    counts = status.astype(str).value_counts()
    table = pd.DataFrame({"status": counts.index.astype(str), "n_cells": counts.to_numpy()})
    table["fraction"] = table["n_cells"] / table["n_cells"].sum()
    figure, axis = store._figure()
    colors = [store.palette(name, "#4C78A8") for name in table["status"]]
    axis.bar(table["status"], table["n_cells"], color=colors)
    _label_count_bars(axis, table["n_cells"])
    axis.set_ylabel("cells")
    axis.set_title("Guide assignment status")
    axis.tick_params(axis="x", rotation=45)
    store.save(figure, "assignment_status_bar", "qc", table, "cell counts by assignment status")


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


def draw_cells_per_guide(store: FigureStore, guides: pd.DataFrame, min_cells: int) -> None:
    if guides.empty:
        store.skip("cells_per_guide_bar", "no_guides")
        return
    top_n = int(store.plotting["top_n_guides"])
    shown = guides.sort_values("n_cells", ascending=False).head(top_n)
    figure, axis = store._figure(height=max(6, 0.28 * len(shown)))
    axis.barh(shown["guide_name"].astype(str), shown["n_cells"], color=store.palette("targeting"))
    axis.axvline(min_cells, color=store.palette("low_power"), linestyle="--")
    axis.set_xlabel("cells")
    axis.invert_yaxis()
    _label_count_bars(axis, shown["n_cells"], horizontal=True)
    store.save(
        figure,
        "cells_per_guide_bar",
        "qc",
        guides,
        f"top {top_n} guides by cell count; line is min_cells_per_guide",
    )


def draw_cells_per_target(store: FigureStore, targets: pd.DataFrame, min_cells: int) -> None:
    if targets.empty:
        store.skip("cells_per_target_bar", "no_targets")
        return
    shown = targets.sort_values("n_cells", ascending=False)
    figure, axis = store._figure(width=max(8, min(24, 0.12 * len(shown))))
    axis.bar(np.arange(len(shown)), shown["n_cells"], color=store.palette("targeting"))
    axis.axhline(min_cells, color=store.palette("low_power"), linestyle="--")
    axis.set_ylabel("cells")
    if len(shown) <= 40:
        axis.set_xticks(np.arange(len(shown)), shown["target_gene"].astype(str), rotation=90)
        _label_count_bars(axis, shown["n_cells"])
    else:
        axis.set_xticks([])
        axis.set_xlabel(f"{len(shown)} targets; names and counts are in the table")
    store.save(
        figure,
        "cells_per_target_bar",
        "qc",
        targets,
        "cells per target; line is min_cells_per_target",
    )


def draw_multi_composition(store: FigureStore, composition: pd.DataFrame) -> None:
    if composition.empty or int(composition["n_cells"].sum()) == 0:
        store.skip("multi_guide_composition_bar", "no_multi_guide_cells")
        return
    figure, axis = store._figure()
    axis.bar(
        composition["composition"], composition["n_cells"], color=store.palette("multi_excluded")
    )
    axis.tick_params(axis="x", rotation=30)
    axis.set_ylabel("cells")
    _label_count_bars(axis, composition["n_cells"])
    store.save(
        figure, "multi_guide_composition_bar", "qc", composition, "multi-guide cell composition"
    )


def _mito_column(metrics: pd.DataFrame) -> str | None:
    if "mito_frac" in metrics.columns:
        return "mito_frac"
    if "pct_counts_mito" in metrics.columns:
        return "pct_counts_mito"
    return None


def percent_mt_values(metrics: pd.DataFrame) -> tuple[pd.Series, str]:
    """Return mitochondrial percentage on a 0–100 scale."""
    column = _mito_column(metrics)
    if column is None:
        raise KeyError("mito_frac")
    values = metrics[column].astype(float)
    finite = values[np.isfinite(values)]
    if finite.empty or float(finite.max()) <= 1:
        return values * 100, "fraction_times_100"
    return values, "already_percent"


def qc_identity_violin(metrics: pd.DataFrame, identity: str) -> Figure:
    """Three Seurat-style violins: nFeature_RNA, nCount_RNA, percent.mt."""
    scaled, _note = percent_mt_values(metrics)
    shown = metrics.copy()
    shown["percent.mt"] = scaled
    figure, axes = plt.subplots(1, 3, figsize=(12, 4))
    panels = (
        ("n_genes_by_counts", "nFeature_RNA"),
        ("total_counts", "nCount_RNA"),
        ("percent.mt", "percent.mt"),
    )
    identity_labels = [identity] * len(shown)
    for axis, (column, title) in zip(axes, panels, strict=True):
        sns.violinplot(
            x=identity_labels, y=shown[column], ax=axis, inner=None, color="#F4A6A6", cut=0
        )
        sns.stripplot(
            x=identity_labels, y=shown[column], ax=axis, color="black", size=2, jitter=0.25
        )
        axis.set_title(title)
        axis.set_xlabel("Identity")
        axis.set_ylabel("")
    return figure


def draw_qc_violin(
    store: FigureStore, metrics: pd.DataFrame, *, population: str = "qc metrics"
) -> None:
    """Seurat VlnPlot: nFeature_RNA, nCount_RNA, percent.mt on a 0–100 scale."""
    mito = _mito_column(metrics)
    needed = {"total_counts", "n_genes_by_counts"}
    if metrics.empty or mito is None or not needed.issubset(metrics.columns):
        store.skip("qc_violin", "missing_fields")
        return
    identity = str(store.context.get("sample", "sample"))
    scaled, note = percent_mt_values(metrics)
    shown = metrics.copy()
    shown["percent.mt"] = scaled
    long = shown.melt(
        id_vars=[column for column in ("cell_barcode",) if column in shown.columns],
        value_vars=["total_counts", "n_genes_by_counts", "percent.mt"],
        var_name="metric",
        value_name="value",
    )
    figure = qc_identity_violin(metrics, identity)
    store.save(
        figure,
        "qc_violin",
        "qc",
        long,
        f"VlnPlot nFeature, nCount, percent.mt ({note}); {population}",
    )


_SINGLER_REFERENCES = {
    "mouse": "MouseRNAseqData",
    "human": "HumanPrimaryCellAtlasData",
}
SINGLER_IMPLEMENTATION = "R SingleR"

_SINGLER_RSCRIPT = r"""
options(timeout = 600)
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


def _run_singler(profiles: pd.DataFrame, reference: SingleRReference, organism: str) -> pd.DataFrame:
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
            if len(reference.features) != matrix.shape[0] or len(reference.labels) != matrix.shape[1]:
                raise ValueError("参考矩阵和基因或标签长度不一致。")
            ref_path = root / "reference.tsv"
            pd.DataFrame(matrix, index=list(reference.features)).reset_index(names="gene").to_csv(
                ref_path, sep="\t", index=False
            )
            label_path = root / "reference_labels.tsv"
            pd.DataFrame({"label": list(reference.labels)}).to_csv(label_path, sep="\t", index=False)
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


def draw_cluster_annotation(store: FigureStore, adata: AnnData) -> None:
    """One figure: Leiden clusters and their SingleR labels on the same UMAP."""
    name = "umap_cluster_annotation"
    if "leiden" not in adata.obs or "X_umap" not in adata.obsm:
        store.skip(name, "missing_fields")
        return
    organism = str(store.context.get("organism", ""))
    try:
        labels, summary = cluster_annotations(adata, organism)
    except ValueError as exc:
        store.skip(name, str(exc))
        return
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
    coordinates = np.asarray(adata.obsm["X_umap"], dtype=float)
    table = pd.DataFrame(
        {
            "UMAP1": coordinates[:, 0],
            "UMAP2": coordinates[:, 1],
            "leiden": adata.obs["leiden"].astype(str).to_numpy(),
            "singler_label": labels.to_numpy(),
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.4), sharex=True, sharey=True)
    sns.scatterplot(
        data=table, x="UMAP1", y="UMAP2", hue="leiden", s=4, linewidth=0, ax=axes[0], legend=False
    )
    sns.scatterplot(
        data=table,
        x="UMAP1",
        y="UMAP2",
        hue="singler_label",
        s=4,
        linewidth=0,
        ax=axes[1],
    )
    axes[0].set_title(f"Leiden (igraph)\nresolution = {resolution_text}")
    axes[1].set_title(f"SingleR ({reference_name})\n{label_field}")
    for axis in axes:
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
    for cluster, group in table.groupby("leiden"):
        axes[1].text(
            group["UMAP1"].median(),
            group["UMAP2"].median(),
            str(cluster),
            fontsize=8,
            ha="center",
            va="center",
        )
    axes[1].legend(title="singler_label", frameon=False, fontsize=7, loc="best")
    method = (
        f"clustering=Leiden igraph, resolution={resolution_text}; "
        f"annotation={summary['annotation_method'].iloc[0]}; "
        "one SingleR label per Leiden cluster; obs['cell_type'] unchanged"
    )
    store.save(figure, name, "plot", summary, method)
    note = store.figure_dir(name, "plot") / f"{name}.txt"
    note.parent.mkdir(parents=True, exist_ok=True)
    note.write_text(
        "\n".join(
            [
                f"聚类：Leiden，igraph，resolution = {resolution_text}",
                (
                    f"细胞注释：{SINGLER_IMPLEMENTATION} {version}，参考 {reference_name}，"
                    f"标签字段 {label_field}。"
                ),
                "先按 Leiden 簇对 log 归一化表达取平均，再对每个簇做 SingleR。",
                "标签写入 obs['singler_label']。不修改 obs['cell_type']。",
                "",
            ]
        )
    )


def draw_qc_scatter(
    store: FigureStore,
    metrics: pd.DataFrame,
    seed: int,
    *,
    population: str = "qc metrics",
) -> None:
    """Seurat FeatureScatter analogue: nCount vs percent.mt, and nCount vs nFeature."""
    mito = _mito_column(metrics)
    needed = {"total_counts", "n_genes_by_counts"}
    if metrics.empty or mito is None or not needed.issubset(metrics.columns):
        store.skip("qc_scatter", "missing_fields")
        return
    limit = int(store.plotting.get("scatter_max_points", 50000))
    shown = metrics
    if len(metrics) > limit:
        shown = metrics.sample(limit, random_state=seed)
    scaled_all, note = percent_mt_values(metrics)
    scaled = scaled_all.loc[shown.index]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(
        shown["total_counts"], scaled, s=8, linewidths=0, c="#4C78A8", rasterized=True
    )
    axes[0].set_xlabel("nCount_RNA")
    axes[0].set_ylabel("percent.mt")
    axes[1].scatter(
        shown["total_counts"],
        shown["n_genes_by_counts"],
        s=8,
        linewidths=0,
        c="#4C78A8",
        rasterized=True,
    )
    axes[1].set_xlabel("nCount_RNA")
    axes[1].set_ylabel("nFeature_RNA")
    store.save(
        figure,
        "qc_scatter",
        "qc",
        shown,
        f"FeatureScatter nCount vs percent.mt ({note}) and nFeature; {population}; seed={seed}",
    )


def draw_waterfall(store: FigureStore, waterfall: pd.DataFrame) -> None:
    if waterfall.empty:
        store.skip("cell_filter_waterfall", "missing_fields")
        return
    figure, axis = store._figure()
    axis.bar(waterfall["step"], waterfall["n_cells"], color=store.palette("targeting"))
    _label_count_bars(axis, waterfall["n_cells"])
    axis.tick_params(axis="x", rotation=45)
    axis.set_ylabel("cells remaining")
    store.save(
        figure, "cell_filter_waterfall", "qc", waterfall, "exclusive filters from cell_status"
    )


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
    store.save(
        figure,
        name,
        "embedding",
        table,
        f"descriptive UMAP from {embedding}; not used as a statistical distance",
    )


def draw_mixscape_fraction(store: FigureStore, fractions: pd.DataFrame) -> None:
    if fractions.empty:
        store.skip("mixscape_fraction_bar", "no_mixscape_cells")
        return
    melted = fractions.melt(
        id_vars=["target_gene", "n_cells"],
        value_vars=["fraction_perturbed", "fraction_non_perturbed", "fraction_ambiguous"],
        var_name="mixscape_status",
        value_name="fraction",
    )
    figure, axis = store._figure(width=max(8, min(24, 0.15 * fractions["target_gene"].nunique())))
    sns.barplot(data=melted, x="target_gene", y="fraction", hue="mixscape_status", ax=axis)
    if fractions["target_gene"].nunique() > 40:
        axis.set_xticks([])
    else:
        axis.tick_params(axis="x", rotation=90)
    store.save(
        figure,
        "mixscape_fraction_bar",
        "mixscape",
        melted,
        "Mixscape class fractions; n_cells in plot data",
    )


def draw_target_expression(store: FigureStore, effects: pd.DataFrame) -> None:
    """One point per target: descriptive CPM log2FC on the main contrast.

    The value is the descriptive CPM log2FC already stored on the effect table.
    The axis names the configured control and the table's pseudocount. It is not
    a mean of log1p expression, and the sensitivity branch is not drawn here.
    """
    effect_column = (
        "log2fc_descriptive" if "log2fc_descriptive" in effects.columns else "effect_log1p"
    )
    if effects.empty or effect_column not in effects.columns:
        store.skip("target_expression_effect", "missing_fields")
        return
    shown = _main_branch(effects).dropna(subset=[effect_column]).sort_values(effect_column)
    if shown.empty:
        store.skip("target_expression_effect", "no_target_expression")
        return
    direction = ""
    for column in ("direction_text", "direction_label"):
        if column in shown.columns and shown[column].notna().any():
            direction = str(shown[column].iloc[0])
            break
    figure, axis = store._figure(width=max(8, min(24, 0.12 * len(shown))))
    axis.axhline(0, color="black", linewidth=0.8)
    axis.scatter(
        np.arange(len(shown)), shown[effect_column], s=12, color=store.palette("targeting")
    )
    if len(shown) <= 40:
        axis.set_xticks(np.arange(len(shown)), shown["target_gene"].astype(str), rotation=90)
    else:
        axis.set_xticks([])
    phrase = _contrast_name(str(store.context.get("control_label", "")))
    alpha = 1.0
    if "pseudocount_cpm" in shown.columns and shown["pseudocount_cpm"].notna().any():
        alpha = float(shown["pseudocount_cpm"].iloc[0])
    axis.set_ylabel(f"log2((CPM_target+{alpha:g})/(CPM of {phrase}+{alpha:g}))")
    axis.set_title(f"target vs {phrase}, not a replicate test")
    store.save(
        figure,
        "target_expression_effect",
        "embedding",
        shown,
        (
            f"descriptive CPM log2FC of the target gene versus {phrase} in this one sample; "
            "cells and guides are not replicates; power and cross-sample reproducibility "
            f"are not estimable; {direction}"
        ),
    )


def draw_pseudobulk_pca(store: FigureStore, pseudobulk: pd.DataFrame, seed: int) -> None:
    if pseudobulk.empty or len(pseudobulk) < 3:
        store.skip("pseudobulk_pca", "fewer_than_3_guides")
        return
    from sklearn.decomposition import PCA

    matrix = pseudobulk.drop(
        columns=[column for column in ("guide_id", "target_gene") if column in pseudobulk.columns]
    )
    values = matrix.to_numpy(dtype=float)
    library = values.sum(axis=1)
    keep = np.isfinite(library) & (library > 0)
    if int(keep.sum()) < 3 or values.shape[1] < 2:
        store.skip("pseudobulk_pca", "fewer_than_3_guides")
        return
    pseudobulk = pseudobulk.iloc[np.flatnonzero(keep)].reset_index(drop=True)
    values = values[keep]
    library = library[keep]
    components = PCA(n_components=2, random_state=seed).fit_transform(
        np.log1p(1e4 * values / library[:, None])
    )
    table = pd.DataFrame(
        {
            "guide_id": pseudobulk["guide_id"].astype(str).to_numpy()
            if "guide_id" in pseudobulk
            else np.arange(len(pseudobulk)),
            "target_gene": pseudobulk["target_gene"].astype(str).to_numpy(),
            "pc1": components[:, 0],
            "pc2": components[:, 1],
        }
    )
    labels = table["target_gene"]
    limit = int(store.plotting["max_legend_categories"])
    if labels.nunique() > limit:
        keep = set(labels.value_counts().head(limit - 1).index)
        table["label"] = labels.where(labels.isin(keep), "Other")
    else:
        table["label"] = labels
    figure, axis = store._figure()
    sns.scatterplot(data=table, x="pc1", y="pc2", hue="label", s=18, linewidth=0, ax=axis)
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    axis.set_title("guide pseudobulk, not replicates")
    store.save(
        figure,
        "pseudobulk_pca",
        "de",
        table,
        "PCA of log1p guide-level CPM; library size is normalized out; "
        "guides are not biological replicates",
    )


def draw_library_size(store: FigureStore, guides: pd.DataFrame) -> None:
    needed = {"guide_id", "n_cells", "total_counts"}
    if guides.empty or not needed.issubset(guides.columns):
        store.skip("pseudobulk_library_size", "missing_fields")
        return
    figure, axis = store._figure()
    sns.scatterplot(
        data=guides,
        x="n_cells",
        y="total_counts",
        hue="target_gene" if guides["target_gene"].nunique() <= 20 else None,
        s=18,
        linewidth=0,
        ax=axis,
    )
    axis.set_xlabel("cells per guide")
    axis.set_ylabel("pseudobulk total counts")
    axis.set_title("guide totals, not replicates")
    store.save(
        figure,
        "pseudobulk_library_size",
        "de",
        guides,
        "guide-level total counts; a guide is not a biological replicate",
    )


def draw_guide_consistency(store: FigureStore, consistency: pd.DataFrame) -> None:
    """Per-guide CPM log2FC of the target gene on the main contrast.

    ``target_gene_log2fc`` is log2((CPM_guide+1)/(CPM_NT+1)). A guide is not a
    biological replicate, so the spread is not a reproducibility estimate.
    """
    effect_column = (
        "target_gene_log2fc" if "target_gene_log2fc" in consistency.columns else "effect_log1p"
    )
    if consistency.empty or effect_column not in consistency.columns:
        store.skip("guide_consistency", "missing_fields")
        return
    shown = _main_branch(consistency).dropna(subset=[effect_column])
    if shown.empty:
        store.skip("guide_consistency", "missing_fields")
        return
    figure, axis = store._figure(width=max(8, min(24, 0.12 * shown["target_gene"].nunique())))
    order = (
        shown.groupby("target_gene")[effect_column]
        .median()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    sns.stripplot(
        data=shown,
        x="target_gene",
        y=effect_column,
        order=order,
        ax=axis,
        color=store.palette("targeting"),
        size=3,
    )
    axis.axhline(0, color="black", linewidth=0.8)
    if len(order) > 40:
        axis.set_xticks([])
        axis.set_xlabel(f"{len(order)} targets; names are in the table")
    else:
        axis.tick_params(axis="x", rotation=90)
    phrase = _contrast_name(str(store.context.get("control_label", "")))
    axis.set_ylabel(f"target-gene CPM log2FC vs {phrase}")
    axis.set_title(f"per guide vs {phrase}, not replicates")
    store.save(
        figure,
        "guide_consistency",
        "de",
        shown,
        f"per-guide target-gene CPM log2FC versus the same {phrase} pool; "
        "guides are not biological replicates; this spread is not a reproducibility estimate",
    )


def draw_distance_heatmap(store: FigureStore, distance: pd.DataFrame) -> None:
    if distance.empty or list(distance.index) != list(distance.columns):
        store.skip("perturbation_distance_heatmap", "missing_fields")
        return
    values = distance.to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or not np.allclose(values, values.T)
        or not np.allclose(np.diag(values), 0)
    ):
        store.skip("perturbation_distance_heatmap", "invalid_distance_matrix")
        return
    figure, axis = plt.subplots(
        figsize=(max(8, min(20, 0.25 * len(distance))), max(6, min(20, 0.25 * len(distance))))
    )
    sns.heatmap(
        distance,
        cmap="mako",
        xticklabels=len(distance) <= 40,
        yticklabels=len(distance) <= 40,
        ax=axis,
    )
    axis.set_title("pairwise energy distance, not the E-test")
    axis.set_xlabel("target")
    axis.set_ylabel("target")
    long = distance.reset_index(names="target_a").melt(
        id_vars="target_a", var_name="target_b", value_name="distance"
    )
    store.save(
        figure,
        "perturbation_distance_heatmap",
        "summary",
        long,
        "pairwise energy distance among guide pseudobulks; "
        "not the contrast against the configured control",
    )


def draw_distance_embedding(store: FigureStore, embedding: pd.DataFrame) -> None:
    needed = {"target_gene", "dimension_1", "dimension_2"}
    if embedding.empty or not needed.issubset(embedding.columns):
        store.skip("perturbation_embedding", "guide_distance_not_estimated")
        return
    if not np.isfinite(embedding[["dimension_1", "dimension_2"]].to_numpy(dtype=float)).all():
        store.skip("perturbation_embedding", "non_finite_coordinates")
        return
    figure, axis = store._figure()
    axis.scatter(
        embedding["dimension_1"], embedding["dimension_2"], s=18, color=store.palette("targeting")
    )
    if len(embedding) <= 40:
        for record in embedding.itertuples(index=False):
            axis.text(record.dimension_1, record.dimension_2, str(record.target_gene), fontsize=7)
    axis.set_xlabel("dimension 1")
    axis.set_ylabel("dimension 2")
    axis.set_title("MDS of pairwise distances, not the E-test")
    store.save(
        figure,
        "perturbation_embedding",
        "summary",
        embedding,
        "MDS of pairwise perturbation distances; not the E-test against the configured control; "
        "axes are not biological variables",
    )


def draw_augur(store: FigureStore, scores: pd.DataFrame) -> None:
    if scores.empty or "auc" not in scores.columns:
        store.skip("augur_score_bar", "missing_fields")
        return
    tested = scores.dropna(subset=["auc"]).sort_values("auc", ascending=False)
    if tested.empty:
        store.skip("augur_score_bar", "insufficient_cells")
        return
    figure, axis = store._figure(width=max(8, min(24, 0.12 * len(tested))))
    x = np.arange(len(tested))
    axis.bar(x, tested["auc"], color=store.palette("targeting"))
    if "auc_sd" in tested and tested["auc_sd"].notna().any():
        axis.errorbar(
            x, tested["auc"], yerr=tested["auc_sd"], fmt="none", ecolor="black", elinewidth=0.6
        )
    axis.axhline(0.5, color="black", linestyle="--")
    if len(tested) <= 40:
        if "cell_type" in tested.columns and tested["cell_type"].astype(str).nunique() > 1:
            tick = tested["target_gene"].astype(str) + "|" + tested["cell_type"].astype(str)
        else:
            tick = tested["target_gene"].astype(str)
        axis.set_xticks(x, tick, rotation=90)
    else:
        axis.set_xticks([])
    axis.set_ylabel("AUC")
    store.save(
        figure,
        "augur_score_bar",
        "summary",
        scores,
        "Augur AUC is label separability, not editing efficiency",
    )


def draw_enrichment(store: FigureStore, enrichment: pd.DataFrame) -> None:
    top_n = int(store.plotting["enrichment_top_n"])
    if (
        enrichment.empty
        or "direction" not in enrichment.columns
        or "term" not in enrichment.columns
    ):
        store.skip("enrichment_dotplot", "no_candidate_genes")
        return
    tested = (
        enrichment[enrichment["reason"] == "tested"]
        if "reason" in enrichment.columns
        else enrichment
    )
    fdr_column = "enrichment_fdr" if "enrichment_fdr" in tested.columns else "fdr"
    if fdr_column not in tested.columns:
        store.skip("enrichment_dotplot", "no_candidate_genes")
        return
    tested = tested.copy()
    tested["adjusted_p"] = pd.to_numeric(tested[fdr_column], errors="coerce")
    tested = tested.dropna(subset=["adjusted_p"])
    if tested.empty:
        store.skip("enrichment_dotplot", "no_candidate_genes")
        return
    best = tested.sort_values("adjusted_p").groupby(["direction", "term"], sort=False).head(1)
    frames = []
    for direction in ("up", "down"):
        part = best[best["direction"] == direction].nsmallest(top_n, "adjusted_p")
        frames.append(part)
    shown = pd.concat(frames, ignore_index=True)
    if shown.empty:
        store.skip("enrichment_dotplot", "no_candidate_genes")
        return
    figure, axes = plt.subplots(1, 2, figsize=(14, 6))
    for axis, direction in zip(axes, ("up", "down"), strict=True):
        part = shown[shown["direction"] == direction]
        if part.empty:
            axis.set_axis_off()
            continue
        axis.scatter(
            part["adjusted_p"],
            np.arange(len(part)),
            s=20 + part["overlap"].to_numpy(dtype=float) * 8,
            c=-np.log10(part["adjusted_p"].clip(lower=1e-300)),
        )
        axis.set_yticks(np.arange(len(part)), part["term"].astype(str))
        axis.set_xlabel("Enrichr adjusted P")
        axis.set_title(f"{direction}; not a replicate test")
    store.save(
        figure,
        "enrichment_dotplot",
        "summary",
        shown,
        "Enrichr overrepresentation of genes with |descriptive CPM log2FC| at or above "
        "volcano_lfc; adjusted P is a gene-set statistic on that fixed list, not a "
        "biological-replicate FDR; power is not estimable from one sample",
    )


def draw_summary(store: FigureStore, summary: pd.DataFrame) -> None:
    metrics = [
        column
        for column in (
            "n_cells",
            "n_guides_passing",
            "fraction_perturbed",
            "effect_log1p",
            "n_significant",
            "auc",
        )
        if column in summary.columns
        and pd.to_numeric(summary[column], errors="coerce").notna().any()
    ]
    if summary.empty or not metrics:
        store.skip("perturbation_summary", "missing_fields")
        return
    figure, axes = plt.subplots(len(metrics), 1, figsize=(8, 2.2 * len(metrics)), sharex=False)
    if len(metrics) == 1:
        axes = [axes]
    order = summary.sort_values(metrics[-1], na_position="last")
    x = np.arange(len(order))
    for axis, metric in zip(axes, metrics, strict=True):
        axis.bar(x, order[metric], color=store.palette("targeting"))
        axis.set_ylabel(metric)
        axis.set_xticks([])
    store.save(
        figure,
        "perturbation_summary",
        "summary",
        summary,
        "one panel per quantity; axes are not shared",
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


def draw_pca_scatter(store: FigureStore, adata: AnnData) -> None:
    """PC1 vs PC2. Leiden is used when it exists; cell type is not invented."""
    if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] < 2:
        store.skip("pca_scatter", "missing_fields")
        return
    coordinates = np.asarray(adata.obsm["X_pca"][:, :2], dtype=float)
    if not np.isfinite(coordinates).all():
        store.skip("pca_scatter", "non_finite_coordinates")
        return
    table = pd.DataFrame({"PC1": coordinates[:, 0], "PC2": coordinates[:, 1]})
    figure, axis = store._figure()
    if "leiden" in adata.obs:
        table["leiden"] = adata.obs["leiden"].astype(str).to_numpy()
        sns.scatterplot(data=table, x="PC1", y="PC2", hue="leiden", s=8, linewidth=0, ax=axis)
        note = "colored by Leiden"
    else:
        axis.scatter(table["PC1"], table["PC2"], s=8, c="#4C78A8", linewidths=0)
        note = "uncolored; no invented cell type"
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    store.save(figure, "pca_scatter", "embedding", table, f"DimPlot PCA; {note}")


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
    axis.set_title("mixscape_class_global")
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    counts = labels.value_counts()
    summary = ", ".join(f"{label}={int(counts.get(label, 0))}" for label in order)
    store.save(
        figure,
        name,
        "mixscape",
        table,
        f"X_umap_mixscape signature space; classes kept: {summary}",
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
    labels = (
        ["NP"] * int(np_index.size) + ["NT"] * int(nt_index.size) + ["KO"] * int(ko_index.size)
    )
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
    """Label NP, NT and KO inside the axes. Heatmap color is only the z-score."""
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
        label = key if key in {"NP", "NT", "KO"} else str(display_names.get(key, key))
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
    adj_p = np.asarray(result["pvals_adj"][target], dtype=float)
    table = pd.DataFrame({"gene": genes, "log_fc": log_fc, "adj_p_value": adj_p})
    up = (table["log_fc"] >= lfc_line) & (table["adj_p_value"] < pval_cutoff)
    down = (table["log_fc"] <= -lfc_line) & (table["adj_p_value"] < pval_cutoff)
    table["direction"] = np.where(up, "Up", np.where(down, "Down", "not DE"))
    positive = table.loc[table["adj_p_value"] > 0, "adj_p_value"]
    floor = float(positive.min() / 10) if len(positive) else 1e-300
    plotted = table["adj_p_value"].where(table["adj_p_value"] > 0, floor)
    table["neg_log10_adj_p"] = -np.log10(plotted)
    return table


def _draw_posterior(axis: Axes, adata: AnnData, target: str, perturbation_label: str) -> None:
    status = _status_labels(adata.obs["mixscape_status"])
    targets = adata.obs["target_gene"].astype(str)
    probability = adata.obs["mixscape_class_p"].astype(float)
    control_name = (
        "zero-guide cells" if (status == "zero_guide").any() else "NT"
    )
    groups = {
        control_name: probability[status.isin(["nt", "zero_guide"])],
        f"{target} NP": probability[(targets == target) & (status == "non_perturbed")],
        f"{target} {perturbation_label}": probability[
            (targets == target) & (status == "perturbed")
        ],
    }
    order = list(groups)
    long = pd.DataFrame(
        {
            "Identity": np.concatenate([[name] * len(values) for name, values in groups.items()]),
            "probability": np.concatenate([values.to_numpy() for values in groups.values()]),
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
    axis.set_title("mixscape posterior probabilities")
    axis.set_xlabel("Identity")
    axis.set_ylabel("")
    axis.tick_params(axis="x", labelrotation=20)


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
    axis.set_xlabel("perturbation score")
    axis.set_ylabel("Cell density")
    axis.set_title(target)


def _draw_volcano_axis(
    axis: Axes,
    table: pd.DataFrame,
    lfc_line: float,
    pval_cutoff: float,
    *,
    target: str = "",
    control: str = "",
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
    labeled = table[table["direction"] != "not DE"].nsmallest(8, "adj_p_value")
    for record in labeled.itertuples(index=False):
        axis.text(record.log_fc, record.neg_log10_adj_p, str(record.gene), fontsize=6)
    if target and control:
        axis.set_xlabel(f"logFC ({target} vs {_contrast_name(control)})")
        axis.set_title("cell-level Wilcoxon, not a replicate test")
    else:
        axis.set_xlabel("log_fc")
    axis.set_ylabel("-log10 adj_p_value")
    axis.legend(frameon=False, fontsize=7)


def _heatmap_panel(
    axis: Axes,
    adata: AnnData,
    target: str,
    perturbation_label: str,
    *,
    seed: int,
    max_genes: int,
    max_cells: int,
) -> int:
    from scipy import sparse

    status = _status_labels(adata.obs["mixscape_status"]).to_numpy()
    targets = adata.obs["target_gene"].astype(str).to_numpy()
    nt_index = _control_positions(status)
    control_block = "zero-guide cells" if np.any(status == "zero_guide") else "NT"
    ko = np.flatnonzero((targets == target) & (status == "perturbed"))
    np_cells = np.flatnonzero((targets == target) & (status == "non_perturbed"))
    layer = "X_pert" if "X_pert" in adata.layers else None
    source = adata.layers[layer] if layer else adata.X
    ko_block = source[ko]
    np_block = source[np_cells]
    if sparse.issparse(ko_block):
        difference = np.abs(
            np.asarray(ko_block.mean(axis=0)).ravel() - np.asarray(np_block.mean(axis=0)).ravel()
        )
    else:
        difference = np.abs(
            np.asarray(ko_block, dtype=float).mean(axis=0)
            - np.asarray(np_block, dtype=float).mean(axis=0)
        )
    genes = np.asarray(
        [index for index in np.argsort(-difference) if difference[index] > 0][:max_genes]
    )
    if genes.size == 0 or min(nt_index.size, np_cells.size, ko.size) == 0:
        axis.set_axis_off()
        axis.set_title(f"{target}: missing NP, {control_block}, or KO")
        return 0
    columns, group_labels = mixscape_heatmap_cell_order(
        _group_index(np_cells, max_cells, seed + 1),
        _group_index(nt_index, max_cells, seed + 2),
        _group_index(ko, max_cells, seed),
    )
    if control_block != "NT":
        group_labels = [control_block if label == "NT" else label for label in group_labels]
    values = _row_matrix(adata, columns, genes, layer).T
    center = values.mean(axis=1, keepdims=True)
    scale = values.std(axis=1, keepdims=True)
    scale[scale == 0] = 1
    scaled = (values - center) / scale
    largest = max(group_labels.count("NP"), group_labels.count("NT"), group_labels.count("KO"))
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
    symbols = (
        adata.var["gene_symbol"].astype(str).to_numpy()
        if "gene_symbol" in adata.var.columns
        else adata.var_names.astype(str).to_numpy()
    )
    axis.set_yticks(np.arange(genes.size))
    axis.set_yticklabels(symbols[genes], fontsize=6)
    axis.set_xticks([])
    _annotate_heatmap_groups(
        axis,
        spans,
        int(genes.size),
        {"NP": f"{target} NP", "NT": "NT", "KO": f"{target} {perturbation_label}"},
    )
    axis.set_title(f"{target} Top {genes.size}")
    return int(genes.size)


def target_in_de_table(adata: AnnData, *, control: str) -> pd.DataFrame:
    """Passed targets, and whether each target symbol is exactly in its own DE list."""
    columns = ["target_gene", "target_in_de", "n_de_genes", "logfc", "adj_p"]
    gate = adata.uns.get("mixscape_de_gate")
    genes = adata.uns.get("mixscape_de_genes")
    if not isinstance(gate, pd.DataFrame) or gate.empty or "passed" not in gate.columns:
        return pd.DataFrame(columns=columns)
    if not isinstance(genes, pd.DataFrame):
        genes = pd.DataFrame(columns=["target_gene", "gene_symbol", "logfc", "adj_p"])
    symbol_col = "gene_symbol" if "gene_symbol" in genes.columns else "gene"
    rows = []
    passed = gate[gate["passed"].astype(bool)]
    for record in passed.itertuples(index=False):
        target = str(record.target_gene)
        if target == control or target == "":
            continue
        if "target_gene" in genes.columns and symbol_col in genes.columns:
            subset = genes[genes["target_gene"].astype(str) == target]
            symbols = subset[symbol_col].astype(str)
            matched = subset[symbols.eq(target)]
        else:
            subset = genes.iloc[0:0]
            matched = subset
        hit = not matched.empty
        logfc = float(matched["logfc"].iloc[0]) if hit and "logfc" in matched.columns else pd.NA
        adj_p = float(matched["adj_p"].iloc[0]) if hit and "adj_p" in matched.columns else pd.NA
        rows.append(
            {
                "target_gene": target,
                "target_in_de": hit,
                "n_de_genes": int(record.n_de_genes),
                "logfc": logfc,
                "adj_p": adj_p,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def export_target_in_de(store: FigureStore, adata: AnnData, *, control: str) -> pd.DataFrame:
    """Copy 2x2 figures whose target symbol is in that target's DE list.

    Every passing target stays in ``DGE/figures/mixscape``. This folder is the subset.
    """
    table = target_in_de_table(adata, control=control)
    data = store.data_dir("target_in_de", "target_in_de")
    data.mkdir(parents=True, exist_ok=True)
    table.to_csv(data / "target_in_de.tsv", sep="\t", index=False)
    destination = store.figure_dir("target_in_de", "target_in_de")
    destination.mkdir(parents=True, exist_ok=True)
    source = store.figure_dir("panel_mixscape_2x2", "mixscape")
    hits = set(table.loc[table["target_in_de"].astype(bool), "target_gene"].astype(str))
    keep = {f"{_safe(target)}_mixscape_2x2" for target in hits}
    for stem in keep:
        for suffix in (".png", ".txt"):
            src = source / f"{stem}{suffix}"
            if src.exists():
                shutil.copy2(src, destination / src.name)
    for path in destination.glob("*_mixscape_2x2.*"):
        if path.stem not in keep:
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
    adata: AnnData, *, control: str, min_de_genes: int = 5
) -> pd.DataFrame:
    """Three counts: all targets, targets with enough DE genes, targets inside their DE list."""
    targets = analysis_target_genes(adata, control)
    target_set = set(targets)
    gate = adata.uns.get("mixscape_de_gate")
    genes = adata.uns.get("mixscape_de_genes")
    with_de: list[str] = []
    if isinstance(gate, pd.DataFrame) and not gate.empty and "passed" in gate.columns:
        enough = gate["passed"].astype(bool) & (gate["n_de_genes"].astype(int) >= int(min_de_genes))
        with_de = [
            gene
            for gene in gate.loc[enough, "target_gene"].astype(str)
            if gene in target_set
        ]
    symbol_col = "gene"
    if isinstance(genes, pd.DataFrame) and "gene_symbol" in genes.columns:
        symbol_col = "gene_symbol"
    in_own_list: list[str] = []
    if isinstance(genes, pd.DataFrame) and not genes.empty and "target_gene" in genes.columns:
        for gene in with_de:
            symbols = genes.loc[genes["target_gene"].astype(str).eq(gene), symbol_col].astype(str)
            if gene in set(symbols):
                in_own_list.append(gene)
    return pd.DataFrame(
        {
            "key": ["n_targets", "n_targets_with_de", "n_targets_in_own_de"],
            "label": ["靶基因总数", "有差异表达基因的靶基因", "差异基因包含对应靶基因"],
            "n": [len(targets), len(with_de), len(in_own_list)],
        }
    )


def target_de_count_bars(counts: pd.DataFrame) -> tuple[Figure, pd.DataFrame]:
    """Three blue bars with integer labels. Order is total, with DE, target-in-DE."""
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
    axis.set_ylabel("靶基因数")
    axis.set_title("靶基因差异表达统计")
    return figure, frame


def target_de_count_note(control: str, min_de_genes: int = 5) -> str:
    return (
        "口径\n"
        "靶基因总数 = 分析里出现的 targeting 靶基因数，"
        f"不含 NT / {control}。\n"
        "有差异表达基因 = 按 Mixscape 门槛筛完后，差异基因数 ≥ "
        f"min_de_genes（{min_de_genes}）的靶基因。"
        "门槛：min_de_genes=5，logfc_threshold=0.25，pval_cutoff=0.05，max_de_genes=100。"
        "不够 5 个的不算「存在差异表达基因」。\n"
        "差异基因包含对应靶基因 = 上一类的子集："
        "该靶基因自己的 symbol 精确出现在它的差异基因名单里。\n"
        "柱子顺序：靶基因总数；有差异表达基因的靶基因；差异基因包含对应靶基因。\n"
    )


def draw_target_de_counts(
    store: FigureStore, adata: AnnData, *, control: str, min_de_genes: int = 5
) -> None:
    name = "target_de_counts"
    if adata is None or "target_gene" not in adata.obs:
        store.skip(name, "missing_fields")
        return
    if not isinstance(adata.uns.get("mixscape_de_gate"), pd.DataFrame):
        store.skip(name, "de_gate_not_run")
        return
    frame = target_de_count_table(adata, control=control, min_de_genes=min_de_genes)
    figure, frame = target_de_count_bars(frame)
    store.save(
        figure,
        name,
        "plot",
        frame,
        "three target-DE counts; integer labels; blue bars; no flag line",
    )
    sidecar = store.figure_dir(name, "plot") / "target_de_counts.txt"
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(target_de_count_note(control, min_de_genes))


def target_in_de_gene_count_table(
    adata: AnnData, *, control: str, min_de_genes: int = 5
) -> pd.DataFrame:
    """Targets that passed the DE gate and contain their own symbol in the stored DE list.

    ``n_de_genes`` is the full gate count. Membership uses that stored list, which
    is capped at ``max_de_genes`` and ranked by raw p. The match is exact.
    """
    columns = ["target_gene", "n_de_genes"]
    table = target_in_de_table(adata, control=control)
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


def target_in_de_gene_count_bars(counts: pd.DataFrame) -> tuple[Figure, pd.DataFrame]:
    """One bar per qualifying target. Height and the integer label are full n_de_genes."""
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
    axis.set_ylabel("差异基因数")
    axis.set_title("差异基因名单包含自身的靶基因")
    if len(frame) > 40:
        axis.set_xticks([])
        axis.set_xlabel(f"{len(frame)} 个靶基因，名称在表里")
    else:
        axis.set_xlabel("靶基因")
        plt.setp(axis.get_xticklabels(), rotation=40, ha="right", fontsize=10)
    legend = axis.get_legend()
    if legend is not None:
        legend.remove()
    return figure, frame


def target_in_de_gene_count_note(control: str, min_de_genes: int = 5) -> str:
    excluded: list[str] = []
    for label in (str(control), "zero_guide"):
        if label and label not in excluded:
            excluded.append(label)
    return (
        "口径\n"
        "只画同时满足下面两条的靶基因，每根柱是一个靶基因。\n"
        "1. 有差异基因：n_de_genes ≥ "
        f"min_de_genes（{min_de_genes}）。"
        "一条基因算差异基因，当 abs(log fold change) ≥ logfc_threshold "
        "且 BH 校正后的 adj_p < 0.05。"
        "n_de_genes 是完整计数，不按 max_de_genes 截断。\n"
        "2. 差异基因名单包含该靶基因：该靶基因的 symbol 与名单里的基因名完全一致"
        "（区分大小写，不是子串匹配）。"
        "名单是打分时保存的差异基因，按 raw p 排序，长度不超过 max_de_genes。\n"
        f"不含 {' / '.join(excluded)}。\n"
        "柱高和柱顶整数都是未截断的 n_de_genes，不是截断后的名单长度。\n"
        "按 n_de_genes 从高到低排列。\n"
    )


def draw_target_in_de_gene_counts(
    store: FigureStore, adata: AnnData, *, control: str, min_de_genes: int = 5
) -> None:
    name = "target_in_de_gene_counts"
    if adata is None or "target_gene" not in adata.obs:
        store.skip(name, "missing_fields")
        return
    if not isinstance(adata.uns.get("mixscape_de_gate"), pd.DataFrame):
        store.skip(name, "de_gate_not_run")
        return
    frame = target_in_de_gene_count_table(adata, control=control, min_de_genes=min_de_genes)
    if frame.empty:
        store.skip(name, "no_targets_in_own_de")
        return
    figure, frame = target_in_de_gene_count_bars(frame)
    store.save(
        figure,
        name,
        "plot",
        frame,
        "targets whose stored DE list contains their own symbol; label is full n_de_genes",
    )
    sidecar = store.figure_dir(name, "plot") / f"{name}.txt"
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(target_in_de_gene_count_note(control, min_de_genes))


_EXCLUDED_TARGET_LABELS = {"", "nt", "zero_guide", "no_guide", "nan", "None", "<NA>", "NA"}
_NOT_SINGLE_TARGET_STATUS = {
    "no_guide",
    "nt",
    "single_nt",
    "multi_guide",
    "multi_excluded",
    "flagged",
    "flagged_high_guides",
}


def single_targeting_cell_counts(
    cells: pd.DataFrame,
    *,
    catalog: list[str] | None = None,
    control: str = "",
) -> pd.DataFrame:
    """Count cells with exactly one guide assigned to each targeting gene.

    NT, no-guide, and multi-guide cells are not added to a target. Targets listed
    in ``catalog`` but absent here stay at 0 so the zero count can be reported.
    """
    excluded = set(_EXCLUDED_TARGET_LABELS)
    excluded.add(control)
    frame = cells
    if "assignment_status" in frame.columns:
        frame = frame[frame["assignment_status"].astype(str).eq("single_targeting")]
    else:
        if "n_guides" in frame.columns:
            frame = frame[pd.to_numeric(frame["n_guides"], errors="coerce").fillna(0).astype(int).eq(1)]
        if "perturbation_status" in frame.columns:
            frame = frame[~frame["perturbation_status"].astype(str).isin(_NOT_SINGLE_TARGET_STATUS)]
    if frame.empty or "target_gene" not in frame.columns:
        observed = pd.Series(dtype=int)
    else:
        genes = frame["target_gene"].astype(str)
        observed = genes[~genes.isin(excluded)].value_counts()
    names = {str(gene) for gene in observed.index if str(gene) not in excluded}
    if catalog is not None:
        names.update(str(gene) for gene in catalog if str(gene) not in excluded)
    rows = [
        {"target_gene": gene, "n_cells": int(observed.get(gene, 0))}
        for gene in sorted(names, key=lambda gene: (-int(observed.get(gene, 0)), gene))
    ]
    return pd.DataFrame(rows, columns=["target_gene", "n_cells"])


def target_cell_count_summary(table: pd.DataFrame) -> dict[str, int]:
    counts = table["n_cells"].astype(int)
    return {
        "n_targets_zero": int((counts == 0).sum()),
        "n_targets_ge_20": int((counts >= 20).sum()),
        "n_targets_ge_40": int((counts >= 40).sum()),
    }


def target_cell_count_bars(table: pd.DataFrame) -> tuple[Figure, pd.DataFrame]:
    """One blue bar per target with cells, highest first. Zeros are counted, not drawn."""
    summary = target_cell_count_summary(table)
    drawn = table[table["n_cells"].astype(int) > 0].reset_index(drop=True)
    figure, axis = plt.subplots(figsize=(max(12.0, min(len(drawn), 800) * 0.02), 5.5))
    axis.bar(
        np.arange(len(drawn)),
        drawn["n_cells"].astype(int),
        color="#4C78A8",
        width=1.0,
        linewidth=0,
    )
    if len(drawn) <= 40:
        axis.set_xticks(np.arange(len(drawn)))
        axis.set_xticklabels(drawn["target_gene"].astype(str), rotation=90)
    else:
        axis.set_xticks([])
        axis.set_xlabel(f"{len(drawn)} 个靶基因，按细胞数从高到低")
    ceiling = float(drawn["n_cells"].max()) if len(drawn) else 1.0
    axis.set_ylim(0, max(ceiling, 40.0) * 1.18)
    axis.axhline(20, color="red", linewidth=1.0)
    axis.axhline(40, color="red", linewidth=1.0)
    axis.text(0, 20, ">=20", color="red", va="bottom", ha="left", fontsize=9)
    axis.text(0, 40, ">=40", color="red", va="bottom", ha="left", fontsize=9)
    axis.text(
        0.98,
        0.98,
        (
            f"细胞数 = 0 的靶基因：{summary['n_targets_zero']}\n"
            f"细胞数 ≥ 20 的靶基因：{summary['n_targets_ge_20']}\n"
            f"细胞数 ≥ 40 的靶基因：{summary['n_targets_ge_40']}"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=10,
    )
    axis.set_ylabel("细胞数")
    return figure, drawn


def target_cell_count_note(summary: dict[str, int], *, source: str) -> str:
    return (
        "口径\n"
        "细胞数 = 该靶基因的 single targeting 细胞数：一个细胞恰好 1 条 guide，"
        "并且分配到这个靶基因。\n"
        "不含配置里的对照标签，不含 no_guide，不含 multi_guide。\n"
        "细胞数为 0 的靶基因不画柱，数量写在图上。\n"
        "细胞数 ≥ 20 与细胞数 ≥ 40 分别计数；≥ 40 是 ≥ 20 的子集，不是互斥区间。\n"
        "柱按细胞数从高到低排列。靶基因过多时不标注每个基因名。\n"
        f"来源：{source}\n"
        f"细胞数 = 0 的靶基因：{summary['n_targets_zero']}\n"
        f"细胞数 ≥ 20 的靶基因：{summary['n_targets_ge_20']}\n"
        f"细胞数 ≥ 40 的靶基因：{summary['n_targets_ge_40']}\n"
    )


def draw_target_cell_counts(
    store: FigureStore,
    table: pd.DataFrame,
    *,
    source: str,
) -> dict[str, int]:
    summary = target_cell_count_summary(table)
    figure, _drawn = target_cell_count_bars(table)
    folder = store.figure_dir("target_cell_counts", "qc")
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "target_cell_counts.png"
    figure.savefig(path, dpi=int(store.plotting.get("dpi", 300)), bbox_inches="tight")
    plt.close(figure)
    stale = path.with_suffix(".pdf")
    if stale.exists():
        stale.unlink()
    (folder / "target_cell_counts.txt").write_text(target_cell_count_note(summary, source=source))
    counts = store.data_dir("target_cell_counts")
    counts.mkdir(parents=True, exist_ok=True)
    table.to_csv(counts / "target_cell_counts.tsv", sep="\t", index=False)
    return summary


def log2fc_bar_genes(frame: pd.DataFrame, target: str, *, n_top: int = 15) -> pd.DataFrame:
    """Highest and lowest log2FC genes. The target symbol is kept even outside that window.

    ``n_top`` matches pertpy ``plot_fold_change(..., n_top_vars=15)``: up to 15
    increased genes and 15 decreased genes.
    """
    columns = ["gene_symbol", "log2fc"]
    if frame.empty or "gene_symbol" not in frame.columns or "log2fc" not in frame.columns:
        return pd.DataFrame(columns=columns)
    work = frame.loc[:, columns].copy()
    work["gene_symbol"] = work["gene_symbol"].astype(str)
    work["log2fc"] = pd.to_numeric(work["log2fc"], errors="coerce")
    work = work.dropna(subset=["log2fc"]).drop_duplicates("gene_symbol", keep="first")
    ordered = work.sort_values(["log2fc", "gene_symbol"], ascending=[False, True])
    increased = ordered[ordered["log2fc"] > 0].head(n_top)
    decreased = ordered[ordered["log2fc"] < 0].nsmallest(n_top, "log2fc")
    chosen = pd.concat([increased, decreased], ignore_index=True)
    hit = ordered[ordered["gene_symbol"].eq(target)]
    if not hit.empty and target not in set(chosen["gene_symbol"]):
        chosen = pd.concat([chosen, hit.iloc[:1]], ignore_index=True)
    chosen = chosen.drop_duplicates("gene_symbol").sort_values(
        ["log2fc", "gene_symbol"], ascending=[False, True]
    )
    return chosen.reset_index(drop=True)


def log2fc_bar_figure(genes: pd.DataFrame) -> Figure:
    """Vertical log2FC bars. High values are dark red and low values are dark blue."""
    from matplotlib.colors import TwoSlopeNorm

    labels = genes["gene_symbol"].astype(str).tolist()
    values = genes["log2fc"].to_numpy(dtype=float)
    figure, axis = plt.subplots(figsize=(max(8.0, 0.32 * max(len(labels), 1)), 5))
    limit = float(np.max(np.abs(values))) if len(values) else 1.0
    if not np.isfinite(limit) or limit == 0:
        limit = 1.0
    colors = plt.get_cmap("RdBu_r")(TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)(values))
    axis.bar(labels, values, color=colors)
    axis.set_ylabel("Log2 fold change")
    axis.set_xlabel("")
    axis.tick_params(axis="x", labelrotation=90)
    axis.set_title("")
    return figure


def target_log2fc_note(target: str, genes: pd.DataFrame, *, definition: str) -> str:
    listed = ", ".join(
        f"{record.gene_symbol}:{record.log2fc:.4g}" for record in genes.itertuples(index=False)
    )
    return (
        f"target={target}\n"
        f"log2fc={definition}\n"
        "layout=pertpy plot_fold_change; n_top_vars=15 up and 15 down; "
        "target symbol is always included; genes ordered by log2FC descending\n"
        "y=Log2 fold change\n"
        f"genes={listed}\n"
    )


def save_target_log2fc_bar(
    folder: Path,
    target: str,
    genes: pd.DataFrame,
    *,
    definition: str,
) -> Path:
    selected = log2fc_bar_genes(genes, target, n_top=15)
    figure = log2fc_bar_figure(selected)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{_safe(target)}_log2fc.png"
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    stale = path.with_suffix(".pdf")
    if stale.exists():
        stale.unlink()
    path.with_suffix(".txt").write_text(target_log2fc_note(target, selected, definition=definition))
    return path


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
    lfc_line: float = 0.25,
) -> None:
    """One 2x2 PNG per target: posterior, perturbation score, cell-level volcano, Top 20 heatmap.

    ``lfc_line`` is MixscapeHeatmap / RunMixscape logfc.threshold (0.25). It is not
    ``plotting.volcano_lfc``. The heatmap stays in the bottom-right panel.
    """
    if "mixscape_status" not in adata.obs or "mixscape_class_p" not in adata.obs:
        store.skip("mixscape_target_panels", "mixscape_not_run")
        return
    if "mixscape" not in adata.uns:
        store.skip("mixscape_target_panels", "perturbation_scores_missing")
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
        record = {"target_gene": target, "n_ko": n_ko, "n_np": n_np, "n_nt": n_nt}
        if min(n_ko, n_np, n_nt) == 0 or _score_frame(adata, str(target)) is None:
            record.update(status="skipped", reason="missing_ko_np_nt_or_score")
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
        figure, axes = plt.subplots(2, 2, figsize=(14, 10))
        figure.perturbseq_manual_layout = True
        _draw_posterior(axes[0, 0], adata, str(target), perturbation_label)
        _draw_perturb_score(axes[0, 1], adata, str(target), perturbation_label, seed)
        _draw_volcano_axis(
            axes[1, 0],
            volcano,
            lfc_line,
            pval_cutoff,
            target=str(target),
            control=str(control),
        )
        n_genes = _heatmap_panel(
            axes[1, 1],
            adata,
            str(target),
            perturbation_label,
            seed=seed,
            max_genes=max_genes,
            max_cells=max_cells,
        )
        figure.suptitle(f"{target} vs {control}")
        store.save(
            figure,
            f"{_safe(str(target))}_mixscape_2x2",
            "mixscape",
            volcano,
            (
                f"cell-level Wilcoxon {target} vs {control}; not a biological-replicate test; "
                f"heatmap genes={n_genes}; lfc={lfc_line}; pval_cutoff={pval_cutoff}"
            ),
        )
        record.update(status="plotted", reason="", n_heatmap_genes=n_genes)
        rows.append(record)
    selection = pd.DataFrame(rows)
    panels = store.data_dir("mixscape_target_panels", "mixscape")
    panels.mkdir(parents=True, exist_ok=True)
    selection.to_csv(panels / "mixscape_target_panels.tsv", sep="\t", index=False)
    export_target_in_de(store, adata, control=control)
    if selection.empty or (selection["status"] != "plotted").all():
        if selection.empty or not (selection["status"] == "plotted").any():
            store.skip("mixscape_target_panels", "no_ko_np_nt_cells")


def _mixscale_gate_frame(adata: AnnData) -> pd.DataFrame | None:
    """Prefer the Mixscale gate. Fall back only when Mixscale has not stored one."""
    for key in ("mixscale_de_gate", "mixscape_de_gate"):
        gate = adata.uns.get(key)
        if (
            isinstance(gate, pd.DataFrame)
            and "passed" in gate.columns
            and "target_gene" in gate.columns
        ):
            return gate
    return None


def _passed_de_targets(adata: AnnData, control: str) -> set[str] | None:
    """Targets that cleared the Mixscale DE gate. None means the gate was not run."""
    gate = _mixscale_gate_frame(adata)
    if gate is None:
        return None
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    passed.discard(str(control))
    passed.discard("")
    return passed


def _skipped_de_targets(adata: AnnData) -> list[str]:
    gate = _mixscale_gate_frame(adata)
    if gate is not None:
        skipped = gate.loc[~gate["passed"].astype(bool), "target_gene"].astype(str)
        return skipped.tolist()
    return [str(gene) for gene in adata.uns.get("mixscale_unscored_targets", [])]


def draw_mixscale_ridge(store: FigureStore, adata: AnnData, control: str) -> None:
    """Ridge plot of mixscale_score for the control and targets that passed the DE gate."""
    name = "mixscale_score_ridge"
    if adata is None or "mixscale_score" not in adata.obs or "target_gene" not in adata.obs:
        store.skip(name, "mixscale_not_run")
        return
    from scipy.stats import gaussian_kde

    passed = _passed_de_targets(adata, control)
    if passed is None:
        store.skip(name, "de_gate_not_run")
        return
    labels = adata.obs["target_gene"].astype(str)
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce")
    identities = [control] if np.isfinite(scores[labels == control]).any() else []
    for identity in sorted(passed):
        if np.isfinite(scores[labels == identity]).any():
            identities.append(identity)
    finite = scores[labels.isin(identities)]
    finite = finite[np.isfinite(finite)]
    if finite.empty or len(identities) < 2:
        store.skip(name, "no_targets_with_de_genes")
        return
    grid = np.linspace(float(finite.min()) - 0.5, float(finite.max()) + 0.5, 200)
    figure, axis = plt.subplots(figsize=(8, max(4, 0.32 * len(identities))))
    figure.perturbseq_manual_layout = True
    cmap = plt.get_cmap("turbo")
    rows = []
    for index, identity in enumerate(identities):
        values = scores[labels == identity].to_numpy(dtype=float)
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
    axis.set_yticklabels([_contrast_name(str(identity)) for identity in identities], fontsize=7)
    axis.set_xlabel("mixscale_score")
    axis.set_ylabel("Identity")
    axis.set_title(f"mixscale_score; control = {_contrast_name(control)}")
    axis.set_ylim(-0.2, len(identities))
    table = pd.DataFrame(rows)
    skipped = _skipped_de_targets(adata)
    store.save(
        figure,
        name,
        "mixscale",
        table,
        "Mixscale ridge of targets that passed the DE gate; "
        "score is perturbation strength, not editing success; no legend",
    )
    sidecar = store.figure_dir(name, "mixscale") / "mixscale_score_ridge.txt"
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    skipped_text = ",".join(skipped) if skipped else "none"
    unscored_passed = sorted(set(passed) - set(identities))
    unscored_text = ",".join(unscored_passed) if unscored_passed else "none"
    sidecar.write_text(
        f"n_cells={int(finite.size)}\n"
        f"n_identities={len(identities)}\n"
        f"control={control}\n"
        f"skipped_targets={skipped_text}\n"
        f"unscored_passed_targets={unscored_text}\n"
        "omitted_reason=not enough DE genes; a score without a passed DE gate is not drawn\n"
        "unscored_passed_reason=expression DE gate passed, but fewer than min_de_genes "
        "could be projected, so no score is drawn\n"
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


def draw_mixscale_scatter(store: FigureStore, adata: AnnData, control: str, nbin: int = 10) -> None:
    """Connected scatter of target-gene expression across mixscale score bins.

    One facet per target that passed the DE gate and has a real finite score.
    NT is not a facet. A target whose scores are all the software fallback 1 is
    omitted. ``nbin`` is 10, matching the Mixscale tutorial, and every qualifying
    target is drawn.
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
        store.skip(name, "de_gate_not_run")
        return
    labels = adata.obs["target_gene"].astype(str).to_numpy()
    scores = pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").to_numpy(dtype=float)
    control = str(control)
    fallback: list[str] = []
    unscored: list[str] = []
    missing: list[str] = []
    case_insensitive: list[str] = []
    frames: list[pd.DataFrame] = []
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
        use = chosen & np.isfinite(scores)
        binned = _mixscale_expression_bins(scores[use], expression[use], nbin)
        binned.insert(0, "target_gene", gene)
        frames.append(binned)
    if not frames:
        store.skip(name, "no_scored_de_targets")
        return
    table = pd.concat(frames, ignore_index=True)
    genes = table["target_gene"].drop_duplicates().tolist()
    ncol = int(np.ceil(np.sqrt(len(genes))))
    nrow = int(np.ceil(len(genes) / ncol))
    figure, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.45 * nrow), squeeze=False)
    colors = sns.color_palette("husl", n_colors=len(genes))
    for index, gene in enumerate(genes):
        axis = axes[index // ncol, index % ncol]
        piece = table[table["target_gene"] == gene]
        axis.plot(
            piece["bin"],
            piece["mean_expression"],
            color=colors[index],
            marker="o",
            linewidth=1.4,
            markersize=4.5,
        )
        axis.set_title(gene, fontsize=10)
        axis.set_xlim(0.5, nbin + 0.5)
        axis.tick_params(labelsize=8)
    for index in range(len(genes), nrow * ncol):
        axes[index // ncol, index % ncol].axis("off")
    figure.supxlabel("bins for perturbation scores")
    figure.supylabel("expression levels of target")
    failed = [gene for gene in _skipped_de_targets(adata) if gene and gene != control]
    store.save(
        figure,
        name,
        "mixscale",
        table,
        "Mixscale_ScatterPlot; quantile bins of mixscale_score within each target; "
        "NT excluded; y is that target gene's own log-normalized expression; no legend",
    )
    sidecar = store.figure_dir(name, "mixscale") / f"{name}.txt"
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(
        f"n_targets={len(genes)}\n"
        f"nbin={nbin}\n"
        f"control={control}\n"
        "nt_in_facets=false\n"
        "binning=within each target, quantile bins of finite mixscale_score "
        "(R quantile type 7, probs=seq(0, 1, 1/nbin)); "
        "each bin is a closed interval, so a score on a cutoff is counted in both "
        "adjacent bins; NT cells are excluded before binning\n"
        "y=mean log-normalized expression of the target gene itself "
        "(adata.X, Seurat RNA data slot)\n"
        "y_axis_symbol_rule=exact gene_symbol or var name first; if both miss, "
        "exactly one case-insensitive gene_symbol is used; two symbols that differ "
        "only by case stay missing; DE-list membership is not changed and stays "
        "case-sensitive\n"
        f"case_insensitive_y={','.join(case_insensitive) or 'none'}\n"
        f"omitted_failed_de={','.join(sorted(set(failed))) or 'none'}\n"
        f"omitted_fallback_score={','.join(fallback) or 'none'}\n"
        f"omitted_no_finite_score={','.join(unscored) or 'none'}\n"
        f"omitted_missing_gene={','.join(missing) or 'none'}\n"
        "fallback_rule=a target whose finite mixscale scores are all 1 is omitted\n"
    )


def _safe(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._+-" else "_" for character in value
    )


# Seurat cc.genes.updated.2019. The original Tirosh 2016 symbols are aliases.
SEURAT_S_GENES = (
    "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU "
    "HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 "
    "RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8"
).split()
SEURAT_G2M_GENES = (
    "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF "
    "TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B "
    "HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 "
    "KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA"
).split()
SEURAT_CC_ALIASES = {"CENPU": "MLF1IP", "PIMREG": "FAM64A", "JPT1": "HN1"}

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


def _resolve_cycle_genes(
    genes: list[str], present: set[str]
) -> tuple[list[str], list[str], list[str], dict[str, str]]:
    """Map one Seurat set onto symbols stored in the object.

    A renamed gene is scored only under the symbol that is present. A symbol
    that is absent is kept in the query so Scanpy drops it, and is listed as
    missing. No stand-in gene is added.
    """
    query: list[str] = []
    found: list[str] = []
    missing: list[str] = []
    aliases: dict[str, str] = {}
    for gene in genes:
        alias = SEURAT_CC_ALIASES.get(gene)
        if gene in present:
            query.append(gene)
            found.append(gene)
        elif alias and alias in present:
            query.append(alias)
            found.append(alias)
            aliases[gene] = alias
        else:
            query.append(gene)
            missing.append(gene)
    return query, found, missing, aliases


def assign_human_cell_cycle(
    adata: AnnData,
    *,
    symbol_column: str = "gene_symbol",
    random_state: int = 0,
) -> dict[str, Any]:
    """Score G1/S/G2M with Scanpy's Seurat/Tirosh cell-cycle scorer.

    Uses the log-normalized matrix in ``X`` (``use_raw=False``). Refuses to score
    when Scanpy has not recorded ``log1p``, so raw counts are not treated as
    log expression. Genes from the Seurat updated-2019 lists that are not in the
    object are omitted from the score and counted in the returned coverage.
    """
    import scanpy as sc

    if "log1p" not in adata.uns:
        raise ValueError("细胞周期打分必须在 log-normalization 之后，不能在原始 counts 上计算。")
    original_names = pd.Index(adata.var_names.astype(str))
    if symbol_column in adata.var.columns:
        symbols = pd.Index(adata.var[symbol_column].astype(str))
    else:
        symbols = original_names
    unique_names, duplicated = _uniquify_symbols(symbols)
    present = set(unique_names)
    s_query, s_found, s_missing, s_alias = _resolve_cycle_genes(SEURAT_S_GENES, present)
    g2m_query, g2m_found, g2m_missing, g2m_alias = _resolve_cycle_genes(SEURAT_G2M_GENES, present)
    gene_pool = unique_names
    if "feature_type" in adata.var.columns:
        gex = adata.var["feature_type"].astype(str).eq("Gene Expression").to_numpy()
        gene_pool = unique_names[gex]
    try:
        adata.var_names = unique_names
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
    duplicated_cycle = [
        gene
        for gene in duplicated
        if gene in set(s_found) | set(g2m_found) | set(SEURAT_CC_ALIASES)
    ]
    return {
        "s_set_size": len(SEURAT_S_GENES),
        "s_found": len(s_found),
        "s_found_symbols": s_found,
        "s_missing_symbols": s_missing,
        "s_alias_used": s_alias,
        "g2m_set_size": len(SEURAT_G2M_GENES),
        "g2m_found": len(g2m_found),
        "g2m_found_symbols": g2m_found,
        "g2m_missing_symbols": g2m_missing,
        "g2m_alias_used": g2m_alias,
        "duplicated_symbols_kept_once": duplicated_cycle,
        "phase_counts": {label: int((phase == label).sum()) for label in ("G1", "G2M", "S")},
    }


def score_human_cell_cycle_after_log1p(
    adata: AnnData,
    *,
    organism: str,
    random_state: int = 0,
) -> str:
    """Score human G1/S/G2M only after log-normalization. Does not change ``X``."""
    if organism != "human":
        return "skipped_not_human"
    if "phase" in adata.obs.columns:
        return "already_scored"
    coverage = assign_human_cell_cycle(adata, random_state=random_state)
    adata.uns["cell_cycle_s_found"] = int(coverage["s_found"])
    adata.uns["cell_cycle_g2m_found"] = int(coverage["g2m_found"])
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
    nt_values: tuple[str, ...] = ("nt",),
    perturbed_values: tuple[str, ...] = ("targeting",),
    point_size: float = 3.0,
    dpi: int = 300,
) -> dict[str, Any]:
    """Draw Perturbed, NT, and cell-cycle panels on one shared UMAP.

    Order is left to right: Perturbed, Perturbation Status (NT), Cell Cycle
    Phase. ``perturbed_values`` selects the Perturbed panel. Mixscape
    ``non_perturbed`` is not included unless the caller puts that label in
    ``perturbed_values``. The PNG and its sidecar txt are written under
    ``output_dir/UMAP/figures/``.
    """
    output_dir = Path(output_dir)
    output_path = output_dir / "UMAP" / "figures" / "umap_nt_perturbed_cell_cycle.png"
    if output_path.suffix.lower() != ".png":
        raise ValueError(f"只写 PNG: {output_path}")
    if embedding not in adata.obsm:
        raise ValueError(f"缺少嵌入 {embedding}")
    if status_column not in adata.obs:
        raise ValueError(f"缺少列 {status_column}")
    if phase_column not in adata.obs:
        raise ValueError(f"缺少列 {phase_column}")
    coordinates = np.asarray(adata.obsm[embedding], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or coordinates.shape[0] != adata.n_obs:
        raise ValueError(f"{embedding} 必须是 n_obs x 2")
    if not np.isfinite(coordinates).all():
        raise ValueError(f"{embedding} 含有非有限坐标")

    status = adata.obs[status_column].astype(str)
    phase = adata.obs[phase_column].astype(str)
    nt_mask = status.isin(nt_values).to_numpy()
    perturbed_mask = status.isin(perturbed_values).to_numpy()
    omitted = ~(nt_mask | perturbed_mask)
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
    perturbed_axis, nt_axis, cycle_axis = axes
    for axis in (perturbed_axis, nt_axis, cycle_axis):
        _prepare_umap_axis(axis, xlim, ylim)

    nt_axis.scatter(
        coordinates[nt_mask, 0],
        coordinates[nt_mask, 1],
        c=_NT_COLOR,
        s=point_size,
        linewidths=0,
        rasterized=True,
    )
    nt_axis.set_title("Perturbation Status\nNT", fontsize=14, pad=8)
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
        "n_nt": int(nt_mask.sum()),
        "n_perturbed": int(perturbed_mask.sum()),
        "n_status_omitted": int(omitted.sum()),
        "omitted_status_counts": {str(key): int(value) for key, value in omitted_counts.items()},
        "n_phase_omitted": int((~phase_known).sum()),
        "phase_counts": {label: int((phase == label).sum()) for label in ("G1", "G2M", "S")},
        "status_column": status_column,
        "nt_values": list(nt_values),
        "perturbed_values": list(perturbed_values),
        "png_path": str(output_path),
        "panel_order": ["Perturbed", "Perturbation Status", "Cell Cycle Phase"],
    }
    note = output_path.with_suffix(".txt")
    omitted_text = (
        ", ".join(f"{key}={value}" for key, value in summary["omitted_status_counts"].items())
        or "无"
    )
    note.write_text(
        "\n".join(
            [
                "panel_order: Perturbed, Perturbation Status (NT), Cell Cycle Phase",
                f"嵌入: obsm['{embedding}']",
                f"NT 面板: {status_column} in {list(nt_values)}",
                (f"Perturbed 面板: {status_column} in {list(perturbed_values)}，不是 Mixscape NP"),
                f"细胞数: {summary['n_cells']}",
                f"NT: {summary['n_nt']}",
                f"Perturbed: {summary['n_perturbed']}",
                f"省略: {summary['n_status_omitted']} ({omitted_text})",
                f"细胞周期: {summary['phase_counts']}",
                f"png: {output_path}",
                "",
            ]
        )
    )
    summary["txt_path"] = str(note)
    return summary
