import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.guide_qc import (
    _effect_direction,
    _interpretation,
    build_qc_warnings,
    compute_gene_guide_consistency,
    compute_guide_metrics,
    compute_weighted_guide_summary,
    expected_on_target_direction,
    on_target_effect_score,
    on_target_pass,
    run_guide_qc,
    write_guide_qc_tables,
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


def test_weighted_merge_umi_math():
    """UMI-weighted mean LFC; inconsistency flags retained under merge modes."""
    adata = _make_adata()
    guide_df, consistency_df, _warnings, summary = run_guide_qc(
        adata, control="NT", min_cells=10, min_median_umi=5, guide_merge="umi"
    )
    assert summary["guide_merge"] == "umi"
    assert summary["weighted"]["enabled"] is True
    geneb = consistency_df.loc[consistency_df["gene_target"] == "GENEB"].iloc[0]
    # GENEB remains inconsistent even with weighted summary
    assert bool(geneb["guides_consistent"]) is False
    assert geneb["guide_merge"] == "umi"
    assert np.isfinite(geneb["weighted_target_log2fc"])

    # Manual UMI-weighted mean for GENEB guides
    sub = guide_df.loc[guide_df["gene_target"] == "GENEB"]
    umi = sub["median_guide_umi"].astype(float).to_numpy()
    lfc = sub["target_expr_log2fc_vs_control"].astype(float).to_numpy()
    expected = float(np.dot(lfc, umi) / umi.sum())
    assert abs(float(geneb["weighted_target_log2fc"]) - expected) < 1e-6
    # High-UMI guide (design_1, down) should dominate direction vote
    assert geneb["weighted_direction_vote"] == "down"

    none_df = compute_weighted_guide_summary(guide_df, guide_merge="none")
    assert none_df.empty


def test_weighted_merge_equal_and_confidence(tmp_path):
    adata = _make_adata()
    guide_df = compute_guide_metrics(adata, control="NT", min_cells=10, min_median_umi=5)
    eq = compute_weighted_guide_summary(guide_df, guide_merge="equal")
    genea = eq.loc[eq["gene_target"] == "GENEA"].iloc[0]
    lfcs = guide_df.loc[guide_df["gene_target"] == "GENEA", "target_expr_log2fc_vs_control"].astype(float)
    assert abs(float(genea["weighted_target_log2fc"]) - float(lfcs.mean())) < 1e-6

    conf = compute_weighted_guide_summary(guide_df, guide_merge="confidence")
    assert not conf.empty
    umi_conf = compute_weighted_guide_summary(guide_df, guide_merge="umi_confidence")
    assert "guide_weights" in umi_conf.columns

    guide_df2, cons, warns, _ = run_guide_qc(
        adata, control="NT", min_cells=10, guide_merge="umi_confidence"
    )
    write_guide_qc_tables(guide_df2, cons, warns, tmp_path)
    assert (tmp_path / "gene_guide_weighted.csv").exists()
    assert (tmp_path / "gene_guide_consistency.csv").exists()
    # Default none does not write companion weighted table
    guide_df3, cons3, warns3, _ = run_guide_qc(adata, control="NT", guide_merge="none")
    out2 = tmp_path / "none"
    write_guide_qc_tables(guide_df3, cons3, warns3, out2)
    assert not (out2 / "gene_guide_weighted.csv").exists()
    assert (out2 / "gene_guide_consistency.csv").exists()
    assert "guide_merge" in pd.read_csv(out2 / "gene_guide_consistency.csv").columns


def test_on_target_score_direction_by_perturbation_type():
    assert expected_on_target_direction("KO") == "down"
    assert expected_on_target_direction("KD") == "down"
    assert expected_on_target_direction("CRISPRi") == "down"
    assert expected_on_target_direction("CRISPRa") == "up"
    assert on_target_pass(-0.5, "KO", 0.25) is True
    assert on_target_pass(-0.1, "KO", 0.25) is False
    assert on_target_pass(0.5, "CRISPRa", 0.25) is True
    assert on_target_pass(-0.5, "CRISPRa", 0.25) is False
    assert on_target_effect_score(-1.0, "KO") == 1.0
    assert on_target_effect_score(1.0, "CRISPRa") == 1.0


def _make_low_efficiency_adata() -> AnnData:
    """GENEA: 3 adequate guides, 2 fail KO downregulation → potential_low_efficiency."""
    genes = ["GENEA", "OTHER"]
    n_nt, n_g = 10, 12
    n = n_nt + 3 * n_g
    rng = np.random.default_rng(1)
    X = rng.poisson(5, size=(n, len(genes))).astype(np.float32)
    X[:n_nt, 0] = 10
    # g1: strong KD; g2/g3: no change (low efficiency)
    X[n_nt : n_nt + n_g, 0] = 2
    X[n_nt + n_g : n_nt + 2 * n_g, 0] = 10
    X[n_nt + 2 * n_g :, 0] = 11
    guide_ids = (
        ["NT_ctrl"] * n_nt
        + ["GENEA|g1"] * n_g
        + ["GENEA|g2"] * n_g
        + ["GENEA|g3"] * n_g
    )
    gene_targets = ["NT"] * n_nt + ["GENEA"] * (3 * n_g)
    obs = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "gene_target": gene_targets,
            "sample_id": ["s1"] * n,
            "guide_umi": [40] * n,
            "num_features": 1,
            "n_counts": rng.integers(1000, 2000, size=n),
            "n_genes": rng.integers(200, 400, size=n),
            "pct_counts_mt": rng.uniform(1, 5, size=n),
        }
    )
    return AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))


def test_low_efficiency_guides_not_no_phenotype():
    adata = _make_low_efficiency_adata()
    guide_df, consistency_df, warnings_df, summary = run_guide_qc(
        adata,
        control="NT",
        min_cells=10,
        min_median_umi=5,
        perturbation_type="KO",
        on_target_lfc_cutoff=0.25,
        on_target_min_fail_guides=2,
    )
    genea = guide_df.loc[guide_df["gene_target"] == "GENEA"]
    assert set(genea.columns) >= {
        "on_target_score",
        "on_target_pass",
        "low_efficiency_guide",
    }
    assert int((~genea["on_target_pass"]).sum()) == 2
    assert int(genea["low_efficiency_guide"].sum()) == 2
    assert (genea.loc[genea["low_efficiency_guide"], "interpretation"] == "potential_low_efficiency_guide").all()
    # Passing guide keeps knockdown label; gene-level table carries potential_low_efficiency.
    assert (
        genea.loc[genea["guide_id"] == "GENEA|g1", "interpretation"].iloc[0]
        == "target_knockdown_detected"
    )

    row = consistency_df.loc[consistency_df["gene_target"] == "GENEA"].iloc[0]
    assert int(row["n_low_efficiency_guides"]) == 2
    assert bool(row["potential_low_efficiency"]) is True
    assert row["interpretation"] == "potential_low_efficiency_guides"
    assert "no phenotype" not in str(row["interpretation"]).lower()

    kinds = set(warnings_df["warning"])
    assert "potential_low_efficiency_guides" in kinds
    assert "low_efficiency_guide" in kinds
    assert "no_clear_effect_adequate_guides" not in kinds
    assert summary["n_genes_potential_low_efficiency"] == 1


def test_crispr_a_on_target_expected_up():
    """CRISPRa: upregulation passes; downregulation fails."""
    genes = ["GENEA"]
    n_nt, n_g = 10, 12
    n = n_nt + 2 * n_g
    X = np.ones((n, 1), dtype=np.float32) * 5
    X[:n_nt, 0] = 5
    X[n_nt : n_nt + n_g, 0] = 20  # up → pass for a
    X[n_nt + n_g :, 0] = 2  # down → fail for a
    obs = pd.DataFrame(
        {
            "guide_id": ["NT_ctrl"] * n_nt + ["GENEA|up"] * n_g + ["GENEA|down"] * n_g,
            "gene_target": ["NT"] * n_nt + ["GENEA"] * (2 * n_g),
            "sample_id": ["s1"] * n,
            "guide_umi": [40] * n,
            "num_features": 1,
            "n_counts": [1500] * n,
            "n_genes": [300] * n,
            "pct_counts_mt": [2.0] * n,
        }
    )
    adata = AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    guide_df, _, _, _ = run_guide_qc(
        adata, control="NT", min_cells=10, perturbation_type="CRISPRa", on_target_lfc_cutoff=0.25
    )
    up = guide_df.loc[guide_df["guide_id"] == "GENEA|up"].iloc[0]
    down = guide_df.loc[guide_df["guide_id"] == "GENEA|down"].iloc[0]
    assert bool(up["on_target_pass"]) is True
    assert bool(down["on_target_pass"]) is False
    assert bool(down["low_efficiency_guide"]) is True
