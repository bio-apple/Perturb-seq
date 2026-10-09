"""Tests for biologist-facing report PNG helpers."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd

from perturbseq.report import write_analysis_report
from perturbseq.report_plots import (
    EFFECT_SUMMARY_NAME,
    generate_report_plots,
    plot_guide_consistency_matrix,
    plot_perturbation_effect_summary,
    plot_target_gene_validation,
)


def _tiny_adata(rng: np.random.Generator | None = None) -> ad.AnnData:
    rng = rng or np.random.default_rng(0)
    genes = [f"G{i}" for i in range(30)] + ["GENEA", "GENEB"]
    n_nt, n_a1, n_a2, n_b = 20, 12, 12, 15
    n = n_nt + n_a1 + n_a2 + n_b
    X = rng.poisson(5, size=(n, len(genes))).astype(np.float32)
    # Knock down GENEA in targeting cells
    ga = genes.index("GENEA")
    X[n_nt : n_nt + n_a1 + n_a2, ga] = rng.poisson(0.5, size=n_a1 + n_a2).astype(np.float32)
    X[:n_nt, ga] = rng.poisson(8, size=n_nt).astype(np.float32)
    obs = pd.DataFrame(
        {
            "gene_target": (["NT"] * n_nt)
            + (["GENEA"] * n_a1)
            + (["GENEA"] * n_a2)
            + (["GENEB"] * n_b),
            "guide_id": (["NT|1"] * n_nt)
            + (["GENEA|g1"] * n_a1)
            + (["GENEA|g2"] * n_a2)
            + (["GENEB|g1"] * n_b),
        },
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))


def _fake_tables(tmp_path, *, with_guide_qc: bool = True):
    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    pd.DataFrame(
        {
            "gene_target": ["GENEA", "GENEB"],
            "edistance": [3.0, 1.2],
            "n_cells": [24, 15],
        }
    ).to_csv(tables / "edistance.csv", index=False)
    pd.DataFrame(
        {
            "distance": [3.0, 1.2],
            "pvalue": [0.01, 0.2],
            "pvalue_adj": [0.02, 0.4],
            "significant_adj": [True, False],
            "significant_adj_reported": [True, False],
        },
        index=["GENEA", "GENEB"],
    ).to_csv(tables / "etest.csv")
    if with_guide_qc:
        pd.DataFrame(
            {
                "guide_id": ["GENEA|g1", "GENEA|g2", "GENEB|g1"],
                "gene_target": ["GENEA", "GENEA", "GENEB"],
                "is_control": [False, False, False],
                "n_cells": [12, 12, 15],
                "target_expr_log2fc_vs_control": [-1.0, -0.8, 0.1],
            }
        ).to_csv(tables / "guide_qc.csv", index=False)
    (tmp_path / "report.json").write_text(
        json.dumps(
            {
                "sample_id": "tiny",
                "config": {"control": "NT", "min_cells_per_pert": 10, "report_plot_top_n": 5},
            }
        )
    )
    return tmp_path


def test_effect_summary_writes_png(tmp_path):
    out = _fake_tables(tmp_path)
    ed = pd.read_csv(out / "tables" / "edistance.csv")
    et = pd.read_csv(out / "tables" / "etest.csv", index_col=0)
    path = plot_perturbation_effect_summary(ed, et, out / "figures")
    assert path is not None
    assert path.name == EFFECT_SUMMARY_NAME
    assert path.exists() and path.stat().st_size > 0


def test_guide_consistency_and_validation_on_synthetic(tmp_path):
    figures = tmp_path / "figures"
    figures.mkdir()
    adata = _tiny_adata()
    cpath = plot_guide_consistency_matrix(adata, "GENEA", figures, control="NT", min_cells=5, n_genes=20)
    assert cpath is not None and cpath.exists()
    assert "guide_consistency_GENEA" in cpath.name

    vpath = plot_target_gene_validation(adata, "GENEA", figures, control="NT")
    assert vpath is not None and vpath.exists()
    assert "target_validation_GENEA" in vpath.name


def test_generate_skips_h5ad_dependent_when_absent(tmp_path):
    out = _fake_tables(tmp_path)
    status = generate_report_plots(out, top_n=5, min_guide_cells=5)
    assert status["effect_summary"] == f"figures/{EFFECT_SUMMARY_NAME}"
    assert (out / "figures" / EFFECT_SUMMARY_NAME).exists()
    reasons = {s["reason"] for s in status["skipped_plots"]}
    assert "h5ad_missing" in reasons
    assert status["guide_consistency"] == {}
    assert status["target_validation"] == {}


def test_write_analysis_report_embeds_figure_paths(tmp_path):
    out = _fake_tables(tmp_path)
    adata = _tiny_adata()
    h5ad = out / "tiny.tertiary.h5ad"
    adata.write_h5ad(h5ad)

    report = write_analysis_report(out, control="NT", min_cells=5, adata=adata, report_plot_top_n=5)
    html = (out / "report.html").read_text()

    assert "figures/perturbation_effect_summary.png" in html
    assert 'id="key-visualizations"' in html
    assert "report_plots" in report
    assert report["report_plots"]["effect_summary"] == "figures/perturbation_effect_summary.png"

    # GENEA has two guides → consistency + validation expected
    assert "GENEA" in report["report_plots"]["guide_consistency"]
    assert "GENEA" in report["report_plots"]["target_validation"]
    assert "figures/guide_consistency_GENEA.png" in html
    assert "figures/target_validation_GENEA.png" in html

    genea = report["perturbations"]["GENEA"]
    assert genea["figures"]["guide_consistency"] == "figures/guide_consistency_GENEA.png"
    assert genea["figures"]["target_validation"] == "figures/target_validation_GENEA.png"
