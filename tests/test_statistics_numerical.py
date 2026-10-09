"""Numerical assertions on synthetic planted-signal DE (no pertpy required)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq.statistics import (
    build_deseq2_design,
    check_experimental_design,
    exploratory_wilcoxon,
    resolve_de_covariates,
    run_de_contrasts,
    run_deseq2_or_wilcoxon,
)


def _planted_wilcoxon_adata(seed: int = 42) -> AnnData:
    """NT vs KO with a strong gene-0 up, mild gene-1 up, null gene-2."""
    rng = np.random.default_rng(seed)
    n_nt, n_ko, n_genes = 50, 50, 40
    X = rng.normal(0.0, 1.0, size=(n_nt + n_ko, n_genes)).astype(np.float32)
    X[n_nt:, 0] += 4.0
    X[n_nt:, 1] += 1.2
    obs = pd.DataFrame(
        {
            "gene_target": ["NT"] * n_nt + ["KOGENE"] * n_ko,
            "de_group": ["NT"] * n_nt + ["KOGENE"] * n_ko,
        }
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return AnnData(X=X, obs=obs, var=var)


def test_experimental_design_flags_no_replicates():
    adata = _planted_wilcoxon_adata()
    design = check_experimental_design(adata, replicate_col=None, groupby="gene_target")
    assert design["replicate_aware"] is False
    assert design["recommended_method"] == "wilcoxon_cell_level_exploratory"
    assert design["evidence_level"] == "exploratory"
    assert design["n_groups"] == 2
    assert design["note"] is not None


def test_experimental_design_replicate_aware():
    adata = _planted_wilcoxon_adata()
    adata.obs["replicate"] = (["r1"] * 25 + ["r2"] * 25) * 2
    design = check_experimental_design(adata, replicate_col="replicate", groupby="gene_target")
    assert design["replicate_aware"] is True
    assert design["n_replicates"] == 2
    assert design["recommended_method"] == "pydeseq2_pseudobulk"
    assert design["evidence_level"] == "inferential"


def test_wilcoxon_planted_signal_ranks_and_pvalues():
    adata = _planted_wilcoxon_adata(seed=42)
    table = exploratory_wilcoxon(adata, groupby="gene_target", group="KOGENE", reference="NT")
    assert not table.empty
    assert table["method"].iloc[0] == "wilcoxon_cell_level_exploratory"
    assert table["evidence_level"].iloc[0] == "exploratory"
    ranked = table.set_index("names")
    assert ranked.index[0] == "g0"
    p0 = float(ranked.loc["g0", "pvals"])
    p1 = float(ranked.loc["g1", "pvals"])
    p2 = float(ranked.loc["g2", "pvals"])
    assert p0 < 1e-6
    assert p0 < p1 < p2
    assert float(ranked.loc["g0", "scores"]) > float(ranked.loc["g1", "scores"])


def test_wilcoxon_empty_when_single_group():
    adata = _planted_wilcoxon_adata()
    adata = adata[adata.obs["gene_target"] == "NT"].copy()
    table = exploratory_wilcoxon(adata, groupby="gene_target", group="KOGENE", reference="NT")
    assert table.empty


def test_run_deseq2_or_wilcoxon_falls_back_without_replicates():
    adata = _planted_wilcoxon_adata()
    table = run_deseq2_or_wilcoxon(
        adata, group="KOGENE", reference="NT", replicate_col=None, groupby="gene_target"
    )
    assert table["method"].iloc[0] == "wilcoxon_cell_level_exploratory"
    assert table["evidence_level"].iloc[0] == "exploratory"
    assert table.set_index("names").index[0] == "g0"


def test_evidence_level_assignment_exploratory_vs_inferential_helper():
    from perturbseq.statistics import de_evidence_level

    assert de_evidence_level(replicate_aware=False) == "exploratory"
    assert de_evidence_level(replicate_aware=True) == "inferential"


def test_de_contrasts_skips_below_min_cells():
    adata = _planted_wilcoxon_adata()
    # One-cell "RARE" group must be skipped at min_cells=10.
    rare = adata[:1].copy()
    rare.obs["de_group"] = "RARE"
    rare.obs["gene_target"] = "RARE"
    combined = AnnData(
        X=np.vstack([adata.X, rare.X]),
        obs=pd.concat([adata.obs, rare.obs], ignore_index=True),
        var=adata.var.copy(),
    )
    tables, errors = run_de_contrasts(
        combined,
        groups=["KOGENE", "RARE"],
        reference="NT",
        groupby="de_group",
        min_cells=10,
        n_jobs=1,
    )
    assert not errors
    assert [g for g, _ in tables] == ["KOGENE"]
    assert tables[0][1].set_index("names").index[0] == "g0"


@pytest.mark.pertpy
def test_edistance_ranking_with_pertpy_if_available():
    """Optional: strong PCA shift ranks above null when pertpy is installed.

    Mixscape post≥pre is *not* asserted here: Mixscape needs a large NT pool,
    is stochastic, and pertpy is an optional heavy extra — unsuitable for the
    default CI job. Ranking of planted PCA clusters is the guaranteed invariant.
    """
    pytest.importorskip("pertpy")
    from perturbseq.perturbation import run_edistance

    rng = np.random.default_rng(0)
    n_nt, n_strong, n_null = 40, 30, 30
    pca = np.zeros((n_nt + n_strong + n_null, 8), dtype=float)
    pca[:n_nt] = rng.normal(0, 0.3, size=(n_nt, 8))
    pca[n_nt : n_nt + n_strong] = rng.normal(0, 0.3, size=(n_strong, 8))
    pca[n_nt : n_nt + n_strong, 0] += 5.0
    pca[n_nt + n_strong :] = rng.normal(0, 0.3, size=(n_null, 8))
    obs = pd.DataFrame(
        {
            "gene_target": ["NT"] * n_nt + ["STRONG"] * n_strong + ["NULL"] * n_null,
        }
    )
    adata = AnnData(X=rng.normal(0, 1, size=(pca.shape[0], 20)), obs=obs)
    adata.obsm["X_pca"] = pca
    edist, _ = run_edistance(
        adata,
        groupby="gene_target",
        contrast="NT",
        n_perms=20,
        min_cells=10,
        etest_top_n=2,
        etest_random_n=0,
        random_state=0,
        n_jobs=1,
        n_bootstrap=0,
    )
    assert float(edist.loc["STRONG", "edistance"]) > float(edist.loc["NULL", "edistance"])


def test_build_deseq2_design_order_and_dedupe():
    assert build_deseq2_design("de_group") == "~ de_group"
    assert (
        build_deseq2_design("de_group", replicate_col="replicate", covariates=["phase", "log_n_counts"])
        == "~ replicate + phase + log_n_counts + de_group"
    )
    # Duplicate groupby / covariate names dropped.
    assert build_deseq2_design("de_group", covariates=["de_group", "phase"]) == "~ phase + de_group"


def test_resolve_de_covariates_true_false_and_list():
    adata = _planted_wilcoxon_adata()
    adata.obs["n_counts"] = 1000.0
    adata.obs["pct_counts_mt"] = 2.0
    adata.obs["phase"] = "G1"
    assert resolve_de_covariates(adata, False) == []
    resolved = resolve_de_covariates(adata, True)
    assert "log_n_counts" in resolved
    assert "pct_counts_mt" in resolved
    assert "phase" in resolved
    assert "log_n_counts" in adata.obs.columns
    assert resolve_de_covariates(adata, ["phase", "missing_col"]) == ["phase"]


def test_check_experimental_design_records_covariates_for_wilcoxon():
    adata = _planted_wilcoxon_adata()
    adata.obs["phase"] = "G1"
    design = check_experimental_design(
        adata, replicate_col=None, groupby="gene_target", covariates=["phase"]
    )
    assert design["replicate_aware"] is False
    assert design["covariates"] == ["phase"]
    assert design["design_formula"] is None
    assert "covariates" in (design["note"] or "").lower() or "phase" in (design["note"] or "")


def test_wilcoxon_records_ignored_covariates():
    adata = _planted_wilcoxon_adata()
    adata.obs["phase"] = "G1"
    table = run_deseq2_or_wilcoxon(
        adata,
        group="KOGENE",
        reference="NT",
        replicate_col=None,
        groupby="gene_target",
        covariates=["phase"],
    )
    assert "covariates_ignored" in table.columns
    assert "phase" in str(table["covariates_ignored"].iloc[0])


def test_pipeline_config_parses_guide_reassign_and_de_covariates():
    from perturbseq.pipeline import pipeline_config_from_mapping

    cfg = pipeline_config_from_mapping(
        {
            "input_dir": "/tmp/in",
            "output_dir": "/tmp/out",
            "guide_reassign": "compare",
            "de_covariates": ["phase", "log_n_counts"],
        }
    )
    assert cfg.guide_reassign == "compare"
    assert cfg.de_covariates == ("phase", "log_n_counts")
