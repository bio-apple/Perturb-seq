"""Toy-data tests for optional guide reassignment helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbseq.guide_reassignment import (
    GUIDE_REASSIGN_MODES,
    apply_guide_ids,
    assign_gmm_per_guide,
    assign_max_count,
    compare_guide_assignments,
    run_guide_reassignment,
    summarize_guide_reassignment,
)


def _toy_crispr(seed: int = 0) -> tuple[AnnData, AnnData]:
    """RNA + CRISPR AnnData with clear max-count winners."""
    rng = np.random.default_rng(seed)
    n_cells = 40
    guides = ["GENE1_g1", "GENE2_g1", "NTC_01"]
    # Each third of cells dominated by one guide.
    mat = rng.poisson(0.3, size=(n_cells, 3)).astype(float)
    mat[0:12, 0] += 20
    mat[12:24, 1] += 20
    mat[24:40, 2] += 20
    crispr = AnnData(
        X=mat,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=guides),
    )
    rna = AnnData(
        X=rng.normal(size=(n_cells, 5)).astype(np.float32),
        obs=pd.DataFrame(
            {
                "guide_id": ["GENE1_g1"] * 12 + ["GENE2_g1"] * 12 + ["NTC_01"] * 16,
                "gene_target": ["GENE1"] * 12 + ["GENE2"] * 12 + ["NT"] * 16,
                "perturbation": ["perturbed"] * 24 + ["NT"] * 16,
                "num_features": 1,
            },
            index=crispr.obs_names,
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(5)]),
    )
    return rna, crispr


def test_assign_max_count_picks_dominant_guide():
    _, crispr = _toy_crispr()
    assigned = assign_max_count(crispr, min_umi=1.0)
    assert (assigned.iloc[:12] == "GENE1_g1").all()
    assert (assigned.iloc[12:24] == "GENE2_g1").all()
    assert (assigned.iloc[24:] == "NTC_01").all()


def test_assign_gmm_runs_and_returns_series():
    _, crispr = _toy_crispr()
    assigned = assign_gmm_per_guide(crispr, min_umi=1.0, random_state=0)
    assert len(assigned) == crispr.n_obs
    assert assigned.astype(str).ne("").mean() > 0.8


def test_compare_and_summarize_agreement():
    rna, crispr = _toy_crispr()
    max_ids = assign_max_count(crispr)
    comparison = compare_guide_assignments(rna.obs["guide_id"], {"max": max_ids})
    summary = summarize_guide_reassignment(comparison)
    assert summary["n_cells"] == 40
    assert summary["methods"]["max"]["frac_agree_with_dragen"] == 1.0


def test_apply_guide_ids_keeps_dragen_backup():
    rna, crispr = _toy_crispr()
    # Flip a few DRAGEN labels so apply changes something.
    rna.obs.loc[rna.obs_names[:3], "guide_id"] = "WRONG"
    max_ids = assign_max_count(crispr)
    apply_guide_ids(rna, max_ids)
    assert "guide_id_dragen" in rna.obs.columns
    assert (rna.obs["guide_id"].iloc[:12] == "GENE1_g1").all()


def test_run_guide_reassignment_compare_and_off(tmp_path):
    rna, crispr = _toy_crispr()
    off = run_guide_reassignment(rna, crispr, mode="off")
    assert off["skipped"] is True
    summary = run_guide_reassignment(
        rna,
        crispr,
        mode="compare",
        tables_dir=tmp_path / "tables",
        figures_dir=tmp_path / "figures",
    )
    assert summary["skipped"] is False
    assert summary["applied"] is False
    assert (tmp_path / "tables" / "guide_reassignment_comparison.csv").is_file()
    assert set(GUIDE_REASSIGN_MODES) >= {"off", "compare", "apply_max", "apply_gmm"}
