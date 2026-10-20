"""Assemble tables and figures."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from perturbseq.config import control_spec, expected_expression, mito_prefix
from perturbseq.distance_tests import refresh_edistance_figure
from perturbseq.mixscape import perturbation_type
from perturbseq.plotting import (
    FigureStore,
    draw_assignment_status,
    draw_cells_per_target,
    draw_cluster_annotation,
    draw_descriptive_heatmap,
    draw_distance_embedding,
    draw_distance_heatmap,
    draw_guide_consistency,
    draw_guides_per_cell,
    draw_library_size,
    draw_mixscale_ridge,
    draw_mixscale_scatter,
    draw_mixscape_class_counts,
    draw_mixscape_class_global_umap,
    draw_mixscape_fraction,
    draw_mixscape_target_panels,
    draw_pca_elbow,
    draw_pca_scatter,
    draw_pseudobulk_pca,
    draw_qc_scatter,
    draw_qc_violin,
    draw_summary,
    draw_target_de_counts,
    draw_target_expression,
    draw_target_in_de_gene_counts,
    draw_temporary_control_qc,
    draw_umap,
    draw_variable_features,
    remove_temporary_control_qc,
)


def write_report(output_dir: Path, config: dict[str, Any], adata: AnnData | None) -> None:
    tables = output_dir / "tables"
    _build_summary(output_dir, config).to_csv(output_dir / "qc_summary.tsv", sep="\t", index=False)
    context = _context(output_dir, config, adata)
    store = FigureStore(output_dir, config, context)
    ledger_path = output_dir / "cell_status.tsv"
    cell_status = _read(ledger_path if ledger_path.exists() else tables / "cell_status.tsv")
    if not cell_status.empty and "assignment_status" in cell_status.columns:
        draw_assignment_status(store, cell_status["assignment_status"])
        draw_guides_per_cell(store, cell_status["n_guides"])
    else:
        for name in ("assignment_status_bar", "guides_per_cell"):
            store.skip(name, "missing_fields")
    qc_tables = output_dir / "qc" / "tables"
    targets = _read(qc_tables / "cells_per_target.tsv")
    draw_cells_per_target(store, targets, int(config["qc"]["min_cells_per_target"]))
    qc_table, qc_note = _qc_population(cell_status, adata)
    if qc_table.empty:
        store.skip("qc_violin", "missing_fields")
        store.skip("qc_scatter", "missing_fields")
    else:
        draw_qc_violin(store, qc_table, population=qc_note)
        draw_qc_scatter(store, qc_table, int(config["runtime"]["seed"]), population=qc_note)
    if control_spec(config)["temporary"]:
        draw_temporary_control_qc(store, cell_status)
    else:
        remove_temporary_control_qc(output_dir)

    if adata is not None:
        draw_variable_features(store, adata)
        draw_pca_elbow(store, adata)
        draw_pca_scatter(store, adata)
    else:
        for name in ("variable_features", "pca_elbow", "pca_scatter"):
            store.skip(name, "missing_fields")
    _draw_cell_embeddings(store, adata)
    _draw_cell_cycle(output_dir, adata, config)
    if adata is not None:
        draw_cluster_annotation(store, adata)
    else:
        store.skip("umap_cluster_annotation", "missing_fields")
    label = perturbation_type(config)
    contrast = control_spec(config)["control_label"]
    if adata is not None:
        draw_mixscape_class_global_umap(store, adata, label)
        draw_mixscape_class_counts(store, adata, label)
        plotting = config.get("plotting", {})
        draw_mixscape_target_panels(
            store,
            adata,
            label,
            control=contrast,
            seed=int(config["runtime"]["seed"]),
            max_genes=int(plotting.get("mixscape_heatmap_max_genes", 20)),
            max_cells=int(plotting.get("mixscape_heatmap_max_cells", 300)),
            pval_cutoff=float(config.get("mixscape", {}).get("pval_cutoff", 0.05)),
            lfc_line=float(config.get("mixscape", {}).get("logfc_threshold", 0.25)),
        )
        draw_mixscale_ridge(store, adata, contrast)
        draw_mixscale_scatter(store, adata, contrast)
        draw_target_de_counts(
            store,
            adata,
            control=contrast,
            min_de_genes=int(config.get("mixscape", {}).get("min_de_genes", 5)),
        )
        draw_target_in_de_gene_counts(
            store,
            adata,
            control=contrast,
            min_de_genes=int(config.get("mixscape", {}).get("min_de_genes", 5)),
        )
    else:
        store.skip("umap_mixscape_class_global", "mixscape_not_run")
        store.skip("mixscape_class_global_bar", "mixscape_not_run")
        store.skip("mixscape_target_panels", "mixscape_not_run")
        store.skip("mixscale_score_ridge", "mixscale_not_run")
        store.skip("mixscale_scatter", "mixscale_not_run")
        store.skip("target_de_counts", "de_gate_not_run")
        store.skip("target_in_de_gene_counts", "de_gate_not_run")
    draw_mixscape_fraction(
        store,
        _read_first(
            output_dir / "DGE" / "Mixscape" / "tables" / "mixscape_summary.tsv",
            output_dir / "DGE" / "mixscape" / "tables" / "mixscape_summary.tsv",
            output_dir / "DGE" / "tables" / "mixscape_summary.tsv",
            tables / "mixscape_fraction.tsv",
        ),
    )
    direction = _read_first(
        output_dir / "DGE" / "descriptive" / "tables" / "target_expression_effect.tsv",
        output_dir / "target_expression_effect.tsv",
        tables / "target_expression.tsv",
    )
    draw_target_expression(store, direction)

    store.skip("volcano", "inference_not_available")
    store.skip("target_effect_forest", "inference_not_available")
    library_path = output_dir / "pseudobulk" / "guide_metadata.tsv"
    library = _read(library_path if library_path.exists() else tables / "pseudobulk_guides.tsv")
    draw_library_size(store, library)
    pseudo_path = output_dir / "pseudobulk" / "guide_counts.h5ad"
    if not pseudo_path.exists():
        pseudo_path = output_dir / "checkpoints" / "pseudobulk.h5ad"
    pca_table = _pseudobulk_matrix(pseudo_path)
    draw_pseudobulk_pca(store, pca_table, int(config["runtime"]["seed"]))
    draw_descriptive_heatmap(
        store,
        _read_first(
            output_dir / "DGE" / "descriptive" / "tables" / "effect_main.tsv",
            output_dir / "effect_main.tsv",
        ),
    )
    guides = _read_first(
        output_dir / "DGE" / "descriptive" / "tables" / "guide_effects.tsv",
        output_dir / "guide_effects.tsv",
    )
    draw_guide_consistency(store, guides)
    distance = _read_matrix_first(
        output_dir / "E-distance" / "tables" / "perturbation_distance.tsv",
        output_dir / "perturbation_distance.tsv",
    )
    draw_distance_heatmap(store, distance)
    draw_distance_embedding(
        store,
        _read_first(
            output_dir / "E-distance" / "tables" / "perturbation_embedding.tsv",
            output_dir / "perturbation_embedding.tsv",
        ),
    )
    summary = _read(output_dir / "qc_summary.tsv")
    if not summary.empty:
        summary["n_significant_de_genes"] = pd.NA
        if "n_significant" in summary.columns:
            summary["n_significant"] = pd.NA
        summary.to_csv(output_dir / "qc_summary.tsv", sep="\t", index=False)
        summary.to_csv(output_dir / "perturbation_summary.tsv", sep="\t", index=False)
    draw_summary(store, summary)
    refresh_edistance_figure(output_dir, adata, config)
    store.write_manifest(output_dir / "plot_manifest.tsv")
    if adata is not None:
        adata.write_h5ad(output_dir / "adata_processed.h5ad", compression="lzf")


def _qc_population(cell_status: pd.DataFrame, adata: AnnData | None) -> tuple[pd.DataFrame, str]:
    """Prefer the full annotated ledger, which is recorded before analysis filtering."""
    mito = {"mito_frac", "pct_counts_mito"}
    needed = {"total_counts", "n_genes_by_counts"}
    if (
        not cell_status.empty
        and needed.issubset(cell_status.columns)
        and mito & set(cell_status.columns)
    ):
        return cell_status, "before filtering; full annotated cells"
    if adata is not None:
        frame = adata.obs.reset_index(drop=False)
        if needed.issubset(frame.columns) and mito & set(frame.columns):
            return frame, "analysis object; full annotated QC columns were absent"
    return pd.DataFrame(), ""


def _draw_cell_cycle(output_dir: Path, adata: AnnData | None, config: dict[str, Any]) -> None:
    """Draw the cell-cycle UMAP after log-normalized scoring.

    Human and mouse use scib gene lists. A matrix that has not been log-normalized
    is not scored.
    """
    if adata is None or "X_umap" not in adata.obsm:
        return
    from perturbseq.plotting import (
        draw_nt_perturbed_cell_cycle_umap,
        score_cell_cycle_after_log1p,
    )

    organism = str(config["experiment"].get("organism", ""))
    if "phase" not in adata.obs.columns and "log1p" in adata.uns:
        try:
            status = score_cell_cycle_after_log1p(
                adata,
                organism=organism,
                random_state=int(config["runtime"]["seed"]),
            )
        except ValueError as exc:
            print(f"  cell cycle skipped: {exc}", flush=True)
            return
        adata.uns["cell_cycle_status"] = status
        if status != "scored_after_log1p" and status != "already_scored":
            return
    if "phase" not in adata.obs.columns:
        return
    spec = control_spec(config)
    control_label = str(spec["control_label"])
    if spec["temporary"]:
        control_values = (control_label,)
        control_title = "zero-guide cells"
    else:
        control_values = tuple(
            dict.fromkeys(
                value for value in ("nt", control_label) if value and value not in {"nan", "None"}
            )
        )
        control_title = "NT"
    draw_nt_perturbed_cell_cycle_umap(
        adata,
        output_dir,
        control_values=control_values,
        control_title=control_title,
        dpi=int(config.get("plotting", {}).get("dpi", 300)),
    )


def _draw_cell_embeddings(store: FigureStore, adata: AnnData | None) -> None:
    # Leiden is drawn in umap_cluster_annotation (left panel); do not duplicate.
    _retire_removed_umaps(store)
    if adata is None or "X_umap" not in adata.obsm:
        store.skip("umap_mixscape_class", "missing_fields")
        store.skip("umap_cell_type", "unstratified_cell_type")
        store.skip("umap_mixscape", "signature_umap_not_run")
        return
    if "mixscape_class" in adata.obs:
        draw_umap(store, adata, "mixscape_class", "umap_mixscape_class")
    else:
        store.skip("umap_mixscape_class", "missing_fields")
    scope = str(adata.uns.get("cell_type_scope", "unstratified"))
    labels = set(adata.obs["cell_type"].astype(str).unique()) if "cell_type" in adata.obs else set()
    if scope == "annotated" and labels - {"ALL", "nan", "None", "<NA>"}:
        draw_umap(store, adata, "cell_type", "umap_cell_type")
    else:
        store.skip("umap_cell_type", "unstratified_cell_type")
    if "mixscape_class" in adata.obs and "X_umap_mixscape" in adata.obsm:
        draw_umap(store, adata, "mixscape_class", "umap_mixscape", embedding="X_umap_mixscape")
    else:
        reason = (
            "mixscape_not_run" if "mixscape_class" not in adata.obs else "signature_umap_not_run"
        )
        store.skip("umap_mixscape", reason)


def _retire_removed_umaps(store: FigureStore) -> None:
    """Drop UMAPs retired from the report (Leiden / assignment / target / status)."""
    names = ("umap_leiden", "umap_assignment", "umap_target", "umap_status")
    relatives: list[Path] = []
    for name in names:
        relatives.extend(
            (
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
        stale = store.output_dir / relative
        if stale.is_file():
            stale.unlink()


def _build_summary(output_dir: Path, config: dict[str, Any]) -> pd.DataFrame:
    tables = output_dir / "tables"
    targets = _read(output_dir / "qc" / "tables" / "cells_per_target.tsv")
    if targets.empty:
        return pd.DataFrame()
    summary = targets.rename(columns={"n_cells": "n_cells"}).copy()
    mix = _read_first(
        output_dir / "DGE" / "Mixscape" / "tables" / "mixscape_summary.tsv",
        output_dir / "DGE" / "mixscape" / "tables" / "mixscape_summary.tsv",
        output_dir / "DGE" / "tables" / "mixscape_summary.tsv",
        tables / "mixscape_fraction.tsv",
    )
    if not mix.empty:
        summary = summary.merge(
            mix[["target_gene", "fraction_perturbed", "n_perturbed"]],
            on="target_gene",
            how="left",
        )
    effects = _read_first(
        output_dir / "DGE" / "descriptive" / "tables" / "target_expression_effect.tsv",
        output_dir / "target_expression_effect.tsv",
        tables / "target_expression.tsv",
    )
    if not effects.empty and "effect_log1p" in effects.columns:
        summary = summary.merge(
            effects[["target_gene", "effect_log1p"]], on="target_gene", how="left"
        )
    de_notes = _read(tables / "de_notes.tsv")
    if not de_notes.empty and "n_fdr_significant" in de_notes.columns:
        summary = summary.merge(
            de_notes[["target_gene", "n_fdr_significant"]].rename(
                columns={"n_fdr_significant": "n_significant"}
            ),
            on="target_gene",
            how="left",
        )
    statistics = config.get("statistics", {})
    if isinstance(statistics, dict) and "fdr_threshold" in statistics:
        summary["fdr_threshold"] = float(statistics["fdr_threshold"])
    return summary


def _context(output_dir: Path, config: dict[str, Any], adata: AnnData | None) -> dict[str, str]:
    _direction, sentence = expected_expression(config)
    spec = control_spec(config)
    thresholds = {}
    threshold_path = output_dir / "tables" / "qc_thresholds.json"
    if threshold_path.exists():
        thresholds = json.loads(threshold_path.read_text())
    n_cells = int(adata.n_obs) if adata is not None else 0
    return {
        "sample": str(config["experiment"]["sample_id"]),
        "organism": str(config["experiment"]["organism"]),
        "crispr_mode": str(config["experiment"]["crispr_mode"]),
        "n_cells": str(n_cells),
        "mito_prefix": mito_prefix(config),
        "qc": json.dumps(thresholds, ensure_ascii=False),
        "direction": sentence,
        "seed": str(config["runtime"]["seed"]),
        "control_source": spec["control_source"],
        "control_label": spec["control_label"],
        "control_phrase": spec["control_phrase"],
    }


def _read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _read_first(*paths: Path) -> pd.DataFrame:
    for path in paths:
        frame = _read(path)
        if not frame.empty or path.exists():
            return frame
    return pd.DataFrame()


def _read_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", index_col=0)


def _read_matrix_first(*paths: Path) -> pd.DataFrame:
    for path in paths:
        frame = _read_matrix(path)
        if not frame.empty or path.exists():
            return frame
    return pd.DataFrame()


def _pseudobulk_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    adata = sc.read_h5ad(path)
    if adata.n_obs < 3:
        return pd.DataFrame()
    if "highly_variable" in adata.var.columns:
        adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    values = adata.layers["counts_exact"] if "counts_exact" in adata.layers else adata.X
    frame = pd.DataFrame(np.asarray(values), index=adata.obs_names.astype(str))
    frame.insert(0, "guide_id", adata.obs_names.astype(str))
    frame.insert(1, "target_gene", adata.obs["target_gene"].astype(str).to_numpy())
    return frame.reset_index(drop=True)
