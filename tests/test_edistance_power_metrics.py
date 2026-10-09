"""Unit tests for E-test low_power flags, provenance keys, and secondary-distance helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbseq.perturbation import (
    annotate_edistance_bootstrap_ci,
    annotate_etest_power,
    combine_distance_tables,
    energy_distance,
    run_secondary_distances,
)
from perturbseq.pipeline import PipelineConfig, _embedding_provenance, pipeline_config_from_mapping


def test_annotate_etest_power_low_power_suppresses_reported_sig():
    etest = pd.DataFrame(
        {
            "pvalue_adj": [0.01, 0.01, 0.4],
            "significant_adj": [True, True, False],
        },
        index=["SMALL", "LARGE", "NULL"],
    )
    n_cells = pd.Series({"SMALL": 20, "LARGE": 80, "NULL": 15})
    out = annotate_etest_power(etest, n_cells, power_min_cells=50)

    assert bool(out.loc["SMALL", "low_power"]) is True
    assert bool(out.loc["LARGE", "low_power"]) is False
    assert bool(out.loc["NULL", "low_power"]) is True
    assert bool(out.loc["SMALL", "significant_adj"]) is True
    assert bool(out.loc["SMALL", "significant_adj_reported"]) is False
    assert bool(out.loc["LARGE", "significant_adj_reported"]) is True
    assert bool(out.loc["NULL", "significant_adj_reported"]) is False
    assert int(out.loc["SMALL", "n_cells"]) == 20


def test_annotate_etest_power_derives_significant_adj_from_padj():
    etest = pd.DataFrame({"pvalue_adj": [0.01, 0.2]}, index=["A", "B"])
    out = annotate_etest_power(etest, {"A": 100, "B": 100}, power_min_cells=50)
    assert bool(out.loc["A", "significant_adj"]) is True
    assert bool(out.loc["A", "significant_adj_reported"]) is True
    assert bool(out.loc["B", "significant_adj_reported"]) is False


def test_bootstrap_edistance_ci_smoke():
    """Tiny cell-bootstrap CI: finite, ordered, and contains the point estimate."""
    rng = np.random.default_rng(0)
    n_nt, n_ko = 30, 25
    pca = np.zeros((n_nt + n_ko, 6), dtype=float)
    pca[:n_nt] = rng.normal(0, 0.3, size=(n_nt, 6))
    pca[n_nt:] = rng.normal(0, 0.3, size=(n_ko, 6))
    pca[n_nt:, 0] += 3.0
    obs = pd.DataFrame({"gene_target": ["NT"] * n_nt + ["KO"] * n_ko})
    adata = AnnData(X=rng.normal(size=(n_nt + n_ko, 8)), obs=obs)
    adata.obsm["X_pca"] = pca

    point = energy_distance(pca[n_nt:], pca[:n_nt])
    edist = pd.DataFrame({"edistance": [point], "n_cells": [n_ko]}, index=["KO"])
    out = annotate_edistance_bootstrap_ci(
        edist,
        adata,
        groupby="gene_target",
        contrast="NT",
        n_bootstrap=20,
        ci_level=0.95,
        random_state=0,
        n_jobs=1,
    )
    lo = float(out.loc["KO", "edistance_ci_low"])
    hi = float(out.loc["KO", "edistance_ci_high"])
    assert np.isfinite(lo) and np.isfinite(hi)
    assert lo <= hi
    assert lo <= point <= hi
    assert int(out.loc["KO", "n_bootstrap"]) == 20
    assert float(out.loc["KO", "ci_level"]) == 0.95


def test_bootstrap_ci_skipped_when_n_bootstrap_zero():
    rng = np.random.default_rng(1)
    adata = AnnData(X=rng.normal(size=(20, 4)))
    adata.obs["gene_target"] = ["NT"] * 10 + ["KO"] * 10
    adata.obsm["X_pca"] = rng.normal(size=(20, 3))
    edist = pd.DataFrame({"edistance": [1.0], "n_cells": [10]}, index=["KO"])
    out = annotate_edistance_bootstrap_ci(edist, adata, n_bootstrap=0)
    assert "edistance_ci_low" not in out.columns


def test_embedding_provenance_always_has_required_keys(tmp_path: Path):
    rng = np.random.default_rng(0)
    adata = AnnData(X=rng.normal(size=(30, 40)))
    adata.var["highly_variable"] = [True] * 12 + [False] * 28
    adata.obsm["X_pca"] = rng.normal(size=(30, 8))
    cfg = PipelineConfig(input_dir=tmp_path / "in", output_dir=tmp_path / "out", n_top_genes=2000, n_pcs=30)
    prov = _embedding_provenance(adata, "log1p_hvg", cfg)
    assert set(prov) >= {"pca_source", "n_hvg", "n_pcs"}
    assert prov["pca_source"] == "log1p_hvg"
    assert prov["n_hvg"] == 12
    assert prov["n_pcs"] == 8


def test_embedding_provenance_fallbacks_without_hvg_pca(tmp_path: Path):
    adata = AnnData(X=np.ones((5, 3)))
    cfg = PipelineConfig(input_dir=tmp_path / "in", output_dir=tmp_path / "out", n_top_genes=99, n_pcs=7)
    prov = _embedding_provenance(adata, "X_pert", cfg)
    assert prov == {"pca_source": "X_pert", "n_hvg": 99, "n_pcs": 7}


def test_combine_distance_tables_merges_secondary():
    edist = pd.DataFrame(
        {"edistance": [3.0, 1.0], "n_cells": [40, 20]},
        index=["G1", "G2"],
    )
    mmd = pd.DataFrame({"mmd": [0.5, 0.1], "n_cells": [40, 20]}, index=["G1", "G2"])
    out = combine_distance_tables(edist, {"mmd": mmd})
    assert float(out.loc["G1", "edistance"]) == 3.0
    assert float(out.loc["G1", "mmd"]) == 0.5
    assert "n_cells" in out.columns


def test_run_secondary_distances_skips_with_reason(tmp_path: Path):
    rng = np.random.default_rng(1)
    n = 40
    adata = AnnData(X=rng.normal(size=(n, 10)))
    adata.obs["gene_target"] = ["NT"] * 20 + ["P1"] * 20
    adata.obsm["X_pca"] = rng.normal(size=(n, 4))

    def _fake_frame(*_a, **_k):
        raise ImportError("No module named 'jax'")

    with patch("perturbseq.perturbation._onesided_distance_frame", side_effect=_fake_frame):
        tables, statuses = run_secondary_distances(
            adata,
            metrics=("mmd", "wasserstein"),
            min_cells=5,
        )
    assert tables == {}
    assert len(statuses) == 2
    assert all(s["skipped"] is True for s in statuses)
    assert statuses[0]["metric"] == "mmd"
    assert statuses[1]["metric"] == "wasserstein"
    assert statuses[1]["reason"] == "missing_jax"


def test_pipeline_config_secondary_metrics_and_power(tmp_path: Path):
    cfg = pipeline_config_from_mapping(
        {
            "input_dir": tmp_path / "in",
            "output_dir": tmp_path / "out",
            "etest_power_min_cells": 80,
            "secondary_distance_metrics": ["mmd"],
        }
    )
    assert cfg.etest_power_min_cells == 80
    assert cfg.secondary_distance_metrics == ("mmd",)


def test_smoke_report_matrix_provenance_keys(tmp_path: Path):
    """Pipeline smoke with skip-distance still records top-level provenance keys."""
    from perturbseq.cli import main

    demo = tmp_path / "demo"
    out = tmp_path / "results"
    assert main(["write-demo", "--output-dir", str(demo)]) == 0
    rc = main(
        [
            "run",
            "--input-dir",
            str(demo),
            "--output-dir",
            str(out),
            "--sample-id",
            "sample1",
            "--skip-mixscape",
            "--skip-distance",
            "--n-perms",
            "10",
        ]
    )
    assert rc == 0
    import json

    report = json.loads((out / "report.json").read_text())
    prov = report["matrix_provenance"]
    assert prov["pca_source"] == "log1p_hvg"
    assert isinstance(prov["n_hvg"], int) and prov["n_hvg"] > 0
    assert isinstance(prov["n_pcs"], int) and prov["n_pcs"] > 0
