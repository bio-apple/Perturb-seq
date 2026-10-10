"""Tests for biologist-module figure / table helpers (QC, gRNA, composition, DE heatmap)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq.composition import (
    composition_by_perturbation,
    guide_assignment_summary,
    write_composition_by_perturbation,
    write_guide_assignment_summary,
)
from perturbseq.plots import (
    plot_composition_by_cluster,
    plot_de_heatmap,
    plot_guide_cell_counts,
    plot_qc,
    plot_qc_cell_counts,
)


def _tiny_adata(n_cells: int = 40) -> AnnData:
    rng = np.random.default_rng(0)
    X = rng.poisson(5, size=(n_cells, 8)).astype(float)
    guides = [f"g{i % 5}" for i in range(n_cells)]
    targets = ["NT" if i % 5 == 0 else f"GENE{i % 5}" for i in range(n_cells)]
    obs = pd.DataFrame(
        {
            "n_counts": X.sum(axis=1),
            "n_genes": (X > 0).sum(axis=1),
            "pct_counts_mt": rng.uniform(0, 10, size=n_cells),
            "guide_id": guides,
            "gene_target": targets,
            "num_features": np.ones(n_cells, dtype=int),
            "leiden": [str(i % 3) for i in range(n_cells)],
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    return AnnData(X=X, obs=obs)


def test_plot_qc_writes_violin_and_scatter(tmp_path: Path) -> None:
    adata = _tiny_adata()
    plot_qc(adata, tmp_path)
    assert (tmp_path / "qc_histograms.png").exists()
    assert (tmp_path / "qc_violin.png").exists()
    assert (tmp_path / "qc_scatter_counts_genes.png").exists()
    assert (tmp_path / "qc_scatter_counts_mt.png").exists()


def test_plot_qc_cell_counts_and_guide_bars(tmp_path: Path) -> None:
    adata = _tiny_adata()
    counts = pd.DataFrame(
        [{"stage": "loaded", "n_cells": 50}, {"stage": "after_qc", "n_cells": 40}, {"stage": "after_singlet", "n_cells": 38}]
    )
    plot_qc_cell_counts(counts, tmp_path)
    plot_guide_cell_counts(adata, tmp_path)
    assert (tmp_path / "qc_cell_counts.png").exists()
    assert (tmp_path / "guide_cell_counts.png").exists()


def test_composition_by_perturbation_table_and_plot(tmp_path: Path) -> None:
    adata = _tiny_adata()
    table = composition_by_perturbation(adata)
    assert not table.empty
    assert {"cluster", "perturbation", "n_cells", "fraction_in_cluster"}.issubset(table.columns)
    path = tmp_path / "composition_by_perturbation.csv"
    write_composition_by_perturbation(table, path)
    assert path.exists()
    plot_composition_by_cluster(table, tmp_path)
    assert (tmp_path / "composition_by_perturbation.png").exists()


def test_guide_assignment_summary(tmp_path: Path) -> None:
    adata = _tiny_adata()
    summary = guide_assignment_summary(adata, control="NT")
    assert "n_cells" in set(summary["metric"])
    assert "ntc_cells" in set(summary["metric"])
    out = tmp_path / "guide_assignment_summary.csv"
    write_guide_assignment_summary(summary, out)
    assert out.exists()


def test_plot_de_heatmap(tmp_path: Path) -> None:
    tables = [
        (
            "GENE1",
            pd.DataFrame(
                {
                    "gene": [f"g{i}" for i in range(30)],
                    "log2FoldChange": np.linspace(-2, 2, 30),
                    "padj": np.linspace(1e-5, 0.2, 30),
                }
            ),
        ),
        (
            "GENE2",
            pd.DataFrame(
                {
                    "gene": [f"g{i}" for i in range(30)],
                    "log2FoldChange": np.linspace(1.5, -1.5, 30),
                    "padj": np.linspace(1e-4, 0.3, 30),
                }
            ),
        ),
    ]
    assert plot_de_heatmap(tables, tmp_path) is True
    assert (tmp_path / "de_heatmap.png").exists()
    assert plot_de_heatmap([], tmp_path) is False


@pytest.mark.parametrize("empty", [pd.DataFrame(), None])
def test_plot_composition_empty_noop(tmp_path: Path, empty) -> None:
    plot_composition_by_cluster(empty, tmp_path)
    assert not (tmp_path / "composition_by_perturbation.png").exists()
