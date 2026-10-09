import anndata as ad
import numpy as np
import pandas as pd

from perturbseq.composition import (
    append_composition_audit,
    composition_snapshot,
    write_composition_audit,
)


def _toy_adata(n: int = 20) -> ad.AnnData:
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=rng.poisson(2, size=(n, 5)).astype(float))
    adata.obs["gene_target"] = (["NT"] * 5 + ["GENE1"] * 10 + ["GENE2"] * 5)[:n]
    adata.obs["guide_id"] = [f"g{i % 4}" for i in range(n)]
    adata.obs["perturbation"] = np.where(adata.obs["gene_target"].eq("NT"), "NT", "perturbed")
    adata.obs["num_features"] = 1
    return adata


def test_composition_snapshot_counts_and_fractions():
    adata = _toy_adata(20)
    snap = composition_snapshot(adata, "loaded", columns=("gene_target", "perturbation"))
    assert set(snap.columns) == {"stage", "category", "value", "n_cells", "fraction"}
    total = snap.loc[snap["category"].eq("n_cells"), "n_cells"].iloc[0]
    assert total == 20
    gene_rows = snap[snap["category"] == "gene_target"]
    assert int(gene_rows.loc[gene_rows["value"] == "GENE1", "n_cells"].iloc[0]) == 10
    assert abs(float(gene_rows.loc[gene_rows["value"] == "GENE1", "fraction"].iloc[0]) - 0.5) < 1e-9


def test_composition_audit_tracks_filter_change(tmp_path):
    adata = _toy_adata(20)
    rows: list[pd.DataFrame] = []
    append_composition_audit(rows, adata, "before")
    filtered = adata[adata.obs["gene_target"] != "GENE2"].copy()
    append_composition_audit(rows, filtered, "after_drop_gene2")
    table = write_composition_audit(rows, tmp_path / "composition_audit.csv")
    assert (tmp_path / "composition_audit.csv").exists()
    before_n = int(table.loc[(table["stage"] == "before") & (table["category"] == "n_cells"), "n_cells"].iloc[0])
    after_n = int(
        table.loc[(table["stage"] == "after_drop_gene2") & (table["category"] == "n_cells"), "n_cells"].iloc[0]
    )
    assert before_n == 20
    assert after_n == 15
    assert "GENE2" not in set(
        table.loc[(table["stage"] == "after_drop_gene2") & (table["category"] == "gene_target"), "value"]
    )
