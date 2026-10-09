"""Edge cases: all-NT pools, singleton guides, DRAGEN feature_call pipe formats."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.annotate import annotate_guides, parse_gene_target, split_feature_calls
from perturbseq.guide_qc import compute_guide_metrics, run_guide_qc
from perturbseq.statistics import exploratory_wilcoxon


def _feature_ref(rows: list[tuple[str, str]]) -> pd.DataFrame:
    """rows: (name, target_gene)."""
    return pd.DataFrame(
        {
            "id": [r[0] for r in rows],
            "name": [r[0] for r in rows],
            "target_gene": [r[1] for r in rows],
        }
    )


def test_split_feature_calls_without_known_names():
    """Heuristic path when feature_ref / known_names is empty."""
    assert split_feature_calls("NEK4|design_3") == ["NEK4|design_3"]
    assert split_feature_calls("NEK4|design_3|STRADB|design_1") == [
        "NEK4|design_3",
        "STRADB|design_1",
    ]
    assert split_feature_calls("Non-Targeting_Control|GENE|design_1") == [
        "Non-Targeting_Control",
        "GENE|design_1",
    ]
    assert split_feature_calls(
        "Non-Targeting_Control|Non-Targeting_Control|PIM1|design_3"
    ) == ["Non-Targeting_Control", "Non-Targeting_Control", "PIM1|design_3"]


def test_parse_gene_target_pipe_and_control_tokens():
    assert parse_gene_target("NEK4|design_3") == "NEK4"
    assert parse_gene_target("GENE|design_1") == "GENE"
    assert parse_gene_target("Non-Targeting_Control") == "NT"
    assert parse_gene_target("") == "unassigned"


def test_annotate_all_cells_nt():
    n, g = 20, 5
    X = sparse.csr_matrix(np.ones((n, g), dtype=np.float32))
    adata = AnnData(X=X, obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    assignments = pd.DataFrame(
        {
            "num_features": [1] * n,
            "feature_call": ["Non-Targeting_Control"] * n,
            "num_transcripts": ["20"] * n,
        },
        index=adata.obs_names,
    )
    ref = _feature_ref([("Non-Targeting_Control", "NT")])
    out = annotate_guides(adata, assignments, ref)
    assert set(out.obs["gene_target"]) == {"NT"}
    assert set(out.obs["perturbation"]) == {"NT"}
    assert (out.obs["guide_id"] == "Non-Targeting_Control").all()


def test_annotate_multi_guide_pipe_formats():
    barcodes = ["a", "b", "c"]
    X = sparse.csr_matrix(np.ones((3, 4), dtype=np.float32))
    adata = AnnData(X=X, obs=pd.DataFrame(index=barcodes))
    known = [
        ("NEK4|design_3", "NEK4"),
        ("STRADB|design_1", "STRADB"),
        ("Non-Targeting_Control", "NT"),
        ("PIM1|design_3", "PIM1"),
        ("GENE|design_1", "GENE"),
    ]
    assignments = pd.DataFrame(
        {
            "num_features": [1, 2, 3],
            "feature_call": [
                "NEK4|design_3",
                "STRADB|design_1|KSR2|design_1",  # KSR2 unknown → still split if in known? not in known
                "Non-Targeting_Control|GENE|design_1",
            ],
            "num_transcripts": ["10", "5|6", "8|9"],
        },
        index=barcodes,
    )
    # Include KSR2 so multi-guide parse is exact.
    known.append(("KSR2|design_1", "KSR2"))
    ref = _feature_ref(known)
    out = annotate_guides(adata, assignments, ref)
    assert out.obs.loc["a", "gene_target"] == "NEK4"
    assert out.obs.loc["a", "guide_id"] == "NEK4|design_3"
    # Multi-guide: unique gene targets joined
    assert set(out.obs.loc["b", "gene_target"].split("|")) == {"STRADB", "KSR2"}
    assert out.obs.loc["b", "num_features"] == 2
    # NT + targeting → gene_target keeps both until NT remap only when *all* NT
    # Here Non-Targeting_Control|GENE → targets NT|GENE then nt_mask is False for multi
    assert "GENE" in out.obs.loc["c", "gene_target"]


def test_guide_with_one_cell_flagged_too_few():
    genes = ["HIT", "OTHER"]
    # NT×15, HIT guide with exactly 1 cell
    n_nt = 15
    rng = np.random.default_rng(0)
    X = rng.poisson(5, size=(n_nt + 1, 2)).astype(np.float32)
    X[:n_nt, 0] = 10
    X[n_nt, 0] = 2
    obs = pd.DataFrame(
        {
            "guide_id": ["NT_ctrl"] * n_nt + ["HIT|design_1"],
            "gene_target": ["NT"] * n_nt + ["HIT"],
            "sample_id": ["s1"] * (n_nt + 1),
            "guide_umi": [40] * n_nt + [30],
            "num_features": 1,
            "n_counts": 1000,
            "n_genes": 200,
            "pct_counts_mt": 2.0,
            "perturbation": ["NT"] * n_nt + ["perturbed"],
        }
    )
    adata = AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    guide_df, _consistency, warnings_df, _summary = run_guide_qc(
        adata, control="NT", min_cells=10, min_median_umi=5
    )
    hit = guide_df.loc[guide_df["guide_id"] == "HIT|design_1"].iloc[0]
    assert int(hit["n_cells"]) == 1
    assert hit["interpretation"] == "too_few_cells"
    assert "insufficient_cells" in set(warnings_df["warning"])


def test_guide_metrics_all_nt_no_targeting_rows():
    n = 12
    X = sparse.csr_matrix(np.ones((n, 3), dtype=np.float32))
    obs = pd.DataFrame(
        {
            "guide_id": ["NT_ctrl"] * n,
            "gene_target": ["NT"] * n,
            "sample_id": ["s1"] * n,
            "guide_umi": [20] * n,
            "num_features": 1,
            "n_counts": 1000,
            "n_genes": 100,
            "pct_counts_mt": 1.0,
            "perturbation": ["NT"] * n,
        }
    )
    adata = AnnData(X=X, obs=obs, var=pd.DataFrame(index=["g0", "g1", "g2"]))
    metrics = compute_guide_metrics(adata, control="NT", min_cells=10)
    # Only NT guide present; no targeting gene rows expected (or only NT).
    targeting = metrics.loc[metrics["gene_target"].astype(str) != "NT"]
    assert targeting.empty


def test_wilcoxon_all_nt_returns_empty():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, 10)).astype(np.float32)
    obs = pd.DataFrame({"gene_target": ["NT"] * 30})
    adata = AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(10)]))
    table = exploratory_wilcoxon(adata, groupby="gene_target", group="MISSING", reference="NT")
    assert table.empty
