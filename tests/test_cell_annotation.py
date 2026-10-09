import numpy as np
from anndata import AnnData
from scipy import sparse

from perturbseq.cell_annotation import annotate_cells, assign_cell_state, score_state_programs


def _toy_adata(n_cells: int = 60, n_genes: int = 40) -> AnnData:
    genes = (
        ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3"]
        + ["VIM", "CDH2", "FN1", "ZEB1", "SNAI1", "SNAI2"]
        + ["MKI67", "TOP2A", "PCNA", "CDK1", "UBE2C", "BIRC5"]
        + [f"GENE{i}" for i in range(n_genes - 18)]
    )
    x = np.random.default_rng(0).poisson(1.0, size=(n_cells, len(genes))).astype(float)
    x[:20, :6] += 8
    x[20:40, 6:12] += 8
    x[40:, 12:18] += 8
    adata = AnnData(X=sparse.csr_matrix(x))
    adata.var_names = genes
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.obs["leiden"] = ["0"] * 20 + ["1"] * 20 + ["2"] * 20
    adata.layers["counts"] = adata.X.copy()
    return adata


def test_state_assignment_separates_programs():
    adata = _toy_adata()
    score_state_programs(adata, ctrl_size=5)
    assign_cell_state(adata, min_score=0.0, min_margin=0.0)
    assert set(adata.obs["cell_state"]).issubset(
        {"epithelial", "EMT", "proliferative", "stress", "hypoxia", "interferon", "unassigned"}
    )
    assert adata.obs.loc[adata.obs["leiden"] == "0", "cell_state"].mode().iloc[0] == "epithelial"


def test_annotate_cells_adds_expected_columns():
    adata = _toy_adata()
    adata, markers = annotate_cells(adata, n_marker_genes=5)
    for col in ("phase", "cell_state", "cluster_annotation"):
        assert col in adata.obs
    assert not markers.empty
    assert "cluster" in markers.columns
