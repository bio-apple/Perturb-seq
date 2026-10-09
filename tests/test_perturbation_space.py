"""Perturbation-space construction: pca_silhouette / kmeans / lr_classifier."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq.perturbation import (
    PERTURBATION_SPACE_METHODS,
    build_perturbation_space,
    cluster_perturbations,
)
from perturbseq.pipeline import pipeline_config_from_mapping


def _tiny_adata(n_per_group: int = 12, n_genes: int = 16, n_pcs: int = 4) -> AnnData:
    """Planted PCA shifts so groups are separable; keep n_cells small for LR."""
    rng = np.random.default_rng(0)
    groups = ["NT", "A", "B", "C"]
    shifts = {
        "NT": np.zeros(n_pcs),
        "A": np.array([3.0, 0.0, 0.0, 0.0]),
        "B": np.array([0.0, 3.0, 0.0, 0.0]),
        "C": np.array([0.0, 0.0, 3.0, 0.0]),
    }
    rows = []
    labels = []
    for g in groups:
        for _ in range(n_per_group):
            rows.append(shifts[g] + rng.normal(0.0, 0.25, size=n_pcs))
            labels.append(g)
    emb = np.asarray(rows, dtype=float)
    # Pad expression so AnnData is valid; PCA embedding is what spaces use.
    x = rng.poisson(2.0, size=(len(labels), n_genes)).astype(float)
    adata = AnnData(X=x)
    adata.obs["gene_target"] = pd.Categorical(labels)
    adata.obsm["X_pca"] = emb
    adata.obs_names = [f"c{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


@pytest.mark.parametrize("method", list(PERTURBATION_SPACE_METHODS))
def test_build_perturbation_space_each_method(method: str):
    adata = _tiny_adata()
    result = build_perturbation_space(adata, method=method, random_state=0)
    assert result.method == method
    assert result.backend in {"local", "pertpy"}
    assert "pert_cluster" in result.clusters.columns
    assert set(result.clusters.index) == {"NT", "A", "B", "C"}
    assert result.embeddings is not None
    assert result.embeddings.shape[0] == 4
    assert result.clusters["pert_cluster"].nunique() >= 1


def test_cluster_perturbations_default_matches_pca_silhouette():
    adata = _tiny_adata()
    a = cluster_perturbations(adata, method="pca_silhouette", random_state=0)
    b = build_perturbation_space(adata, method="pca_silhouette", random_state=0).clusters
    pd.testing.assert_frame_equal(a, b)


def test_invalid_method_raises():
    adata = _tiny_adata()
    with pytest.raises(ValueError, match="perturbation_space"):
        build_perturbation_space(adata, method="not_a_method")


def test_pipeline_config_accepts_perturbation_space():
    cfg = pipeline_config_from_mapping(
        {
            "input_dir": "data/demo",
            "output_dir": "results/demo",
            "perturbation_space": "kmeans",
        }
    )
    assert cfg.perturbation_space == "kmeans"


def test_pipeline_config_rejects_bad_perturbation_space():
    with pytest.raises(ValueError, match="perturbation_space"):
        pipeline_config_from_mapping(
            {
                "input_dir": "data/demo",
                "output_dir": "results/demo",
                "perturbation_space": "embedding_space",
            }
        )


def test_kmeans_fixed_n_clusters():
    adata = _tiny_adata()
    result = build_perturbation_space(adata, method="kmeans", n_clusters=2, random_state=0)
    assert result.metadata.get("k_selection") == "user"
    assert int(result.clusters["pert_cluster"].nunique()) <= 2
