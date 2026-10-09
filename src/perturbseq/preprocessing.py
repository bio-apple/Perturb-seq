"""Normalization, PCA, neighbors, UMAP, and Leiden for visualization / downstream analysis."""

from __future__ import annotations

from scipy import sparse
from anndata import AnnData
import scanpy as sc


def preprocess_rna(
    adata: AnnData,
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    leiden_resolution: float = 0.5,
    random_state: int = 0,
) -> AnnData:
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    # Preserve raw UMIs; never replace an existing counts layer with normalized X.
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_top = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, subset=False)
    n_comps = min(n_pcs, max(2, adata.n_obs - 1), max(2, adata.n_vars - 1))
    sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)
    n_neigh = min(n_neighbors, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neigh, n_pcs=n_comps, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    try:
        sc.tl.leiden(
            adata,
            resolution=leiden_resolution,
            flavor="igraph",
            n_iterations=2,
            random_state=random_state,
            key_added="leiden",
        )
    except Exception:
        try:
            sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state, key_added="leiden")
        except Exception:
            adata.obs["leiden"] = "0"
    adata.uns["matrix_provenance"] = {
        "layers_counts": "raw_umi",
        "X": "log1p_normalized",
        "obsm_X_pca": "PCA_on_HVG_of_log1p_X",
        "obsm_X_umap": "visualization_only",
        "obs_leiden": "descriptive_cluster_not_perturbation_evidence",
    }
    return adata
