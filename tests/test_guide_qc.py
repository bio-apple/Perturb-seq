import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.guide_qc import (
    build_qc_warnings,
    compute_gene_guide_consistency,
    compute_guide_metrics,
    run_guide_qc,
    _effect_direction,
    _interpretation,
)


def _make_adata() -> AnnData:
    """Synthetic screen: GENEA has consistent KD; GENEB has opposing guides; GENEC too few cells."""
    genes = ["GENEA", "GENEB", "GENEC", "OTHER"]
    # cells: NT×8, GENEA_g1×12, GENEA_g2×12, GENEB_g1×12, GENEB_g2×12, GENEC_g1×3
    n_nt, n_a, n_b, n_c = 8, 12, 12, 3
    n = n_nt + 2 * n_a + 2 * n_b + n_c
    rng = np.random.default_rng(0)
    # columns = genes
    X = rng.poisson(5, size=(n, len(genes))).astype(np.float32)
    # NT baseline for GENEA/GENEB ~ 10
    X[:n_nt, 0] = 10
    X[:n_nt, 1] = 10
    # GENEA guides both knockdown
    X[n_nt : n_nt + n_a, 0] = 2
    X[n_nt + n_a : n_nt + 2 * n_a, 0] = 3
    # GENEB g1 down, g2 up
    X[n_nt + 2 * n_a : n_nt + 2 * n_a + n_b, 1] = 2
    X[n_nt + 2 * n_a + n_b : n_nt + 2 * n_a + 2 * n_b, 1] = 20
    # GENEC few cells mild KD
    start_c = n_nt + 2 * n_a + 2 * n_b
    X[start_c:, 2] = 4

    guide_ids = (
        ["NT_ctrl"] * n_nt
        + ["GENEA|design_1"] * n_a
        + ["GENEA|design_2"] * n_a
        + ["GENEB|design_1"] * n_b
        + ["GENEB|design_2"] * n_b
        + ["GENEC|design_1"] * n_c
    )
    gene_targets = (
        ["NT"] * n_nt
        + ["GENEA"] * n_a
        + ["GENEA"] * n_a
        + ["GENEB"] * n_b
        + ["GENEB"] * n_b
        + ["GENEC"] * n_c
    )
    samples = ["s1"] * (n // 2) + ["s2"] * (n - n // 2)
    umi = (
        [50] * n_nt
        + [40] * n_a
        + [35] * n_a
        + [30] * n_b
        + [2] * n_b  # low UMI on GENEB design_2
        + [20] * n_c
    )
    obs = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "gene_target": gene_targets,
            "sample_id": samples,
            "guide_umi": umi,
            "num_features": 1,
            "n_counts": rng.integers(1000, 2000, size=n),
            "n_genes": rng.integers(200, 400, size=n),
            "pct_counts_mt": rng.uniform(1, 5, size=n),
            "perturbation": ["NT"] * n_nt + ["perturbed"] * (n - n_nt),
        }
    )
    return AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))


def test_effect_direction_thresholds():
    assert _effect_direction(-0.5) == "down"
    assert _effect_direction(0.5) == "up"
    assert _effect_direction(0.1) == "none"
    assert _effect_direction(float("nan")) == "unknown"


def test_interpretation_layers():
    assert (
        _interpretation(
            n_cells=3,
            min_cells=10,
            median_umi=20,
            min_umi=5,
            detection_rate=1.0,
            min_detection=0.5,
            target_in_matrix=True,
            direction="none",
            gene_consistent=True,
            n_guides_for_gene=2,
        )
        == "too_few_cells"
    )
    assert (
        _interpretation(
            n_cells=20,
            min_cells=10,
            median_umi=20,
            min_umi=5,
            detection_rate=1.0,
            min_detection=0.5,
            target_in_matrix=True,
            direction="none",
            gene_consistent=True,
            n_guides_for_gene=2,
        )
        == "no_target_effect_with_adequate_guides"
    )


def test_guide_metrics_and_consistency():
    adata = _make_adata()
    guide_df = compute_guide_metrics(adata, control="NT", min_cells=10, min_median_umi=5)
    assert "GENEA|design_1" in set(guide_df["guide_id"])
    genea = guide_df.loc[guide_df["gene_target"] == "GENEA"]
    assert (genea["target_effect_direction"] == "down").all()
    assert genea["n_cells"].min() >= 10

    consistency = compute_gene_guide_consistency(guide_df)
    genea_c = consistency.loc[consistency["gene_target"] == "GENEA"].iloc[0]
    geneb_c = consistency.loc[consistency["gene_target"] == "GENEB"].iloc[0]
    assert bool(genea_c["guides_consistent"]) is True
    assert bool(geneb_c["guides_consistent"]) is False


def test_warnings_flag_insufficient_low_umi_inconsistent():
    adata = _make_adata()
    guide_df, consistency_df, warnings_df, summary = run_guide_qc(
        adata, control="NT", min_cells=10, min_median_umi=5
    )
    kinds = set(warnings_df["warning"])
    assert "insufficient_cells" in kinds  # GENEC
    assert "low_guide_umi" in kinds  # GENEB|design_2
    assert "inconsistent_guides" in kinds  # GENEB
    assert summary["n_genes_inconsistent_guides"] >= 1
    assert summary["n_genes_consistent_guides"] >= 1
    # interpretation distinguishes biology vs technical for GENEA (KD detected)
    genea = guide_df.loc[guide_df["guide_id"] == "GENEA|design_1"].iloc[0]
    assert genea["interpretation"] == "target_knockdown_detected"


def test_detection_rate_across_samples():
    adata = _make_adata()
    # Put GENEC only in s1
    mask = adata.obs["guide_id"].astype(str) == "GENEC|design_1"
    adata.obs.loc[mask, "sample_id"] = "s1"
    guide_df = compute_guide_metrics(adata, control="NT", min_cells=1)
    genec = guide_df.loc[guide_df["guide_id"] == "GENEC|design_1"].iloc[0]
    assert genec["n_samples_total"] == 2
    assert genec["n_samples_detected"] == 1
    assert genec["detection_rate"] == 0.5


def test_build_qc_warnings_empty_inputs():
    warn = build_qc_warnings(pd.DataFrame(), pd.DataFrame())
    assert list(warn.columns) == ["level", "entity", "warning", "detail"]
