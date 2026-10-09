"""Perturbation-aware MAD QC: thresholds from NT; cytotoxicity when guides are depleted."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.guide_qc import merge_qc_filter_stats, run_guide_qc
from perturbseq.qc import add_qc_metrics, filter_cells


def _synthetic_screen() -> AnnData:
    """Controlled QC metrics: NT healthy, MILD within NT MAD, TOXIC far below."""
    rng = np.random.default_rng(0)
    n_nt, n_mild, n_toxic = 40, 20, 20
    n = n_nt + n_mild + n_toxic
    n_genes = 30
    X = rng.poisson(50, size=(n, n_genes)).astype(np.float32)
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_genes - 1)] + ["MT-CO1"])
    obs = pd.DataFrame(
        {
            "guide_id": (
                ["NT_ctrl"] * n_nt + ["MILD|g1"] * n_mild + ["TOXIC|g1"] * n_toxic
            ),
            "gene_target": ["NT"] * n_nt + ["MILD"] * n_mild + ["TOXIC"] * n_toxic,
            "perturbation": ["NT"] * n_nt + ["perturbed"] * (n_mild + n_toxic),
            "num_features": 1,
            "sample_id": "s1",
            "guide_umi": 30.0,
        }
    )
    adata = AnnData(X=sparse.csr_matrix(X), obs=obs, var=var)
    adata = add_qc_metrics(adata)
    # Overwrite metrics so thresholds are deterministic (expression matrix unused for MAD).
    # Mild is within ~5 MAD of NT (spread ~200); toxic is far below → QC drop + cyto flag.
    adata.obs["n_counts"] = np.array(
        [10_000.0] * n_nt + [9_200.0] * n_mild + [200.0] * n_toxic, dtype=float
    )
    adata.obs["n_genes"] = np.array(
        [2_000.0] * n_nt + [1_900.0] * n_mild + [50.0] * n_toxic, dtype=float
    )
    adata.obs["pct_counts_mt"] = np.array(
        [2.0] * n_nt + [3.0] * n_mild + [40.0] * n_toxic, dtype=float
    )
    # Add mild NT jitter so MAD > 0.
    adata.obs.loc[adata.obs["gene_target"] == "NT", "n_counts"] = 10_000 + rng.normal(
        0, 200, size=n_nt
    )
    adata.obs.loc[adata.obs["gene_target"] == "NT", "n_genes"] = 2_000 + rng.normal(
        0, 40, size=n_nt
    )
    return adata


def test_perturbation_aware_retains_mild_phenotype():
    """Mild low-count targeting cells survive when thresholds are fit on NT."""
    adata = _synthetic_screen()
    filtered_aware, log_aware = filter_cells(
        adata,
        n_mads=5.0,
        min_cells=1,
        perturbation_aware=True,
        control="NT",
    )
    assert log_aware["qc_threshold_source"] == "control:NT"
    assert log_aware["perturbation_aware_qc"] is True
    mild_kept = int((filtered_aware.obs["guide_id"] == "MILD|g1").sum())
    assert mild_kept == 20, f"expected all mild cells retained, got {mild_kept}"

    by_guide = log_aware["qc_filter_by_guide"]
    mild_row = by_guide.loc[
        (by_guide["category"] == "guide_id") & (by_guide["value"] == "MILD|g1")
    ].iloc[0]
    assert mild_row["frac_removed_qc"] == 0.0
    assert not bool(mild_row["cytotoxicity_qc_flag"])


def test_cytotoxicity_flag_when_guide_mostly_removed():
    """≥90% QC removal → cytotoxicity_qc_flag + guide_qc warning / stub row."""
    adata = _synthetic_screen()
    filtered, log = filter_cells(
        adata,
        n_mads=5.0,
        min_cells=1,
        perturbation_aware=True,
        control="NT",
    )
    by_guide = log["qc_filter_by_guide"]
    toxic = by_guide.loc[
        (by_guide["category"] == "guide_id") & (by_guide["value"] == "TOXIC|g1")
    ].iloc[0]
    assert float(toxic["frac_removed_qc"]) >= 0.9
    assert bool(toxic["cytotoxicity_qc_flag"])
    assert toxic["gene_target"] == "TOXIC"
    assert int((filtered.obs["guide_id"] == "TOXIC|g1").sum()) == 0

    guide_df, _cons, warnings_df, _summary = run_guide_qc(
        filtered, control="NT", min_cells=5, qc_filter_df=by_guide
    )
    toxic_guides = guide_df.loc[guide_df["guide_id"] == "TOXIC|g1"]
    assert not toxic_guides.empty, "fully depleted guide must appear as stub in guide_qc"
    assert bool(toxic_guides.iloc[0]["cytotoxicity_qc_flag"])
    assert toxic_guides.iloc[0]["interpretation"] == "cytotoxicity_qc_depletion"
    assert (warnings_df["warning"] == "cytotoxicity_qc_depletion").any()


def test_default_global_mad_still_available():
    """Flag off → thresholds from all cells (backward compatible)."""
    adata = _synthetic_screen()
    _filtered, log = filter_cells(adata, n_mads=5.0, min_cells=1)
    assert log["perturbation_aware_qc"] is False
    assert log["qc_threshold_source"] == "all_cells"
    assert isinstance(log["qc_filter_by_guide"], pd.DataFrame)
    assert not log["qc_filter_by_guide"].empty


def test_merge_qc_filter_stats_standalone():
    guide_df = pd.DataFrame(
        {
            "guide_id": ["MILD|g1"],
            "gene_target": ["MILD"],
            "is_control": [False],
            "n_cells": [10],
            "interpretation": ["pending"],
        }
    )
    filter_df = pd.DataFrame(
        [
            {
                "category": "guide_id",
                "value": "MILD|g1",
                "gene_target": "MILD",
                "n_cells_before_qc": 20,
                "n_cells_after_qc": 10,
                "n_cells_removed_qc": 10,
                "frac_removed_qc": 0.5,
                "cytotoxicity_qc_flag": False,
            },
            {
                "category": "guide_id",
                "value": "TOXIC|g1",
                "gene_target": "TOXIC",
                "n_cells_before_qc": 20,
                "n_cells_after_qc": 0,
                "n_cells_removed_qc": 20,
                "frac_removed_qc": 1.0,
                "cytotoxicity_qc_flag": True,
            },
        ]
    )
    merged = merge_qc_filter_stats(guide_df, filter_df, control="NT")
    assert set(merged["guide_id"]) == {"MILD|g1", "TOXIC|g1"}
    toxic = merged.loc[merged["guide_id"] == "TOXIC|g1"].iloc[0]
    assert bool(toxic["cytotoxicity_qc_flag"])
    assert toxic["interpretation"] == "cytotoxicity_qc_depletion"
