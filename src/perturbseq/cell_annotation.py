from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

# Regev cell-cycle gene lists (Tirosh et al.), as used by scanpy tutorials.
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
    "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
    "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
    "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# State panels for epithelial / cancer cell-line Perturb-seq (not tissue cell types).
DEFAULT_STATE_MARKERS: dict[str, list[str]] = {
    "epithelial": ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4"],
    "EMT": ["VIM", "CDH2", "FN1", "ZEB1", "SNAI1", "SNAI2", "TWIST1"],
    "proliferative": ["MKI67", "TOP2A", "PCNA", "CDK1", "UBE2C", "BIRC5"],
    "stress": ["HSPA5", "HSP90B1", "DDIT3", "ATF4", "XBP1", "HERPUD1"],
    # Avoid constitutive glycolysis genes (ENO1/LDHA) that dominate cancer cell lines.
    "hypoxia": ["HIF1A", "VEGFA", "BNIP3", "SLC2A1", "CA9", "NDRG1", "P4HA1", "ADM"],
    "interferon": ["ISG15", "IFIT1", "IFIT3", "MX1", "STAT1", "OAS1"],
}


def _present_genes(adata: AnnData, genes: list[str]) -> list[str]:
    available = set(map(str, adata.var_names))
    return [g for g in genes if g in available]


def score_cell_cycle(adata: AnnData) -> AnnData:
    s_genes = _present_genes(adata, S_GENES)
    g2m_genes = _present_genes(adata, G2M_GENES)
    if len(s_genes) < 3 or len(g2m_genes) < 3:
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"
        return adata
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata


def score_state_programs(
    adata: AnnData,
    marker_sets: dict[str, list[str]] | None = None,
    ctrl_size: int = 50,
) -> pd.DataFrame:
    marker_sets = marker_sets or DEFAULT_STATE_MARKERS
    used: dict[str, list[str]] = {}
    for name, genes in marker_sets.items():
        present = _present_genes(adata, genes)
        if len(present) < 2:
            continue
        score_name = f"score_{name}"
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name, ctrl_size=ctrl_size)
        used[name] = present
    adata.uns["cell_state_markers_used"] = used
    score_cols = [f"score_{name}" for name in used]
    return adata.obs[score_cols] if score_cols else pd.DataFrame(index=adata.obs_names)


def assign_cell_state(adata: AnnData, min_score: float = 0.1, min_margin: float = 0.05) -> AnnData:
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if not score_cols:
        adata.obs["cell_state"] = "unassigned"
        return adata
    scores = adata.obs[score_cols]
    # Median-center across cells so labels reflect relative enrichment, not baseline expression.
    centered = scores - scores.median(axis=0)
    values = centered.to_numpy(dtype=float)
    order = np.argsort(values, axis=1)
    best_idx = order[:, -1]
    second_idx = order[:, -2] if values.shape[1] > 1 else best_idx
    rows = np.arange(values.shape[0])
    best_vals = values[rows, best_idx]
    second_vals = values[rows, second_idx]
    labels = np.array([score_cols[i].removeprefix("score_") for i in best_idx], dtype=object)
    unsure = (best_vals < min_score) | ((best_vals - second_vals) < min_margin)
    adata.obs["cell_state"] = np.where(unsure, "unassigned", labels)
    adata.obs["cell_state_score"] = best_vals
    return adata


def ensure_leiden_clusters(adata: AnnData, cluster_key: str = "leiden", resolutions: tuple[float, ...] = (0.8, 1.2, 1.6)) -> AnnData:
    """Recluster if previous Leiden collapsed to a single community."""
    if cluster_key in adata.obs and adata.obs[cluster_key].nunique() > 1:
        return adata
    if "neighbors" not in adata.uns:
        n_pcs = min(30, adata.obsm["X_pca"].shape[1]) if "X_pca" in adata.obsm else None
        sc.pp.neighbors(adata, n_pcs=n_pcs)
    for resolution in resolutions:
        try:
            sc.tl.leiden(
                adata,
                resolution=resolution,
                flavor="igraph",
                n_iterations=2,
                key_added=cluster_key,
            )
        except Exception:
            try:
                sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)
            except Exception:
                continue
        if adata.obs[cluster_key].nunique() > 1:
            adata.uns["leiden_resolution_used"] = resolution
            return adata
    if cluster_key not in adata.obs:
        adata.obs[cluster_key] = "0"
    return adata


def annotate_clusters_by_majority(adata: AnnData, cluster_key: str = "leiden") -> AnnData:
    if cluster_key not in adata.obs or "cell_state" not in adata.obs:
        return adata

    def _majority(series: pd.Series) -> str:
        counts = series.astype(str).value_counts()
        assigned = counts.drop(labels=["unassigned"], errors="ignore")
        if not assigned.empty and assigned.iloc[0] >= max(3, 0.15 * len(series)):
            return str(assigned.index[0])
        return str(counts.index[0])

    mapping = adata.obs.groupby(cluster_key, observed=True)["cell_state"].agg(_majority).to_dict()
    adata.obs["cluster_annotation"] = adata.obs[cluster_key].astype(str).map(mapping)
    adata.uns["cluster_annotation_map"] = {str(k): str(v) for k, v in mapping.items()}
    return adata


def find_cluster_markers(
    adata: AnnData,
    cluster_key: str = "leiden",
    n_genes: int = 25,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    if cluster_key not in adata.obs or adata.obs[cluster_key].nunique() < 2:
        return pd.DataFrame()
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, n_genes=n_genes)
    frames = []
    for group in adata.obs[cluster_key].astype(str).unique():
        frame = sc.get.rank_genes_groups_df(adata, group=group)
        frame.insert(0, "cluster", group)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def annotate_cells(
    adata: AnnData,
    cluster_key: str = "leiden",
    marker_sets: dict[str, list[str]] | None = None,
    n_marker_genes: int = 25,
) -> tuple[AnnData, pd.DataFrame]:
    """Cell-cycle + state-program scoring + cluster labels.

    For cell-line Perturb-seq this is state annotation, not tissue cell-type calling.
    """
    ensure_leiden_clusters(adata, cluster_key=cluster_key)
    score_cell_cycle(adata)
    score_state_programs(adata, marker_sets=marker_sets)
    assign_cell_state(adata)
    annotate_clusters_by_majority(adata, cluster_key=cluster_key)
    markers = find_cluster_markers(adata, cluster_key=cluster_key, n_genes=n_marker_genes)
    return adata, markers


def annotation_summary(adata: AnnData) -> dict:
    summary: dict = {}
    if "phase" in adata.obs:
        summary["phase"] = adata.obs["phase"].astype(str).value_counts().to_dict()
    if "cell_state" in adata.obs:
        summary["cell_state"] = adata.obs["cell_state"].astype(str).value_counts().to_dict()
    if "cluster_annotation" in adata.obs and "leiden" in adata.obs:
        summary["cluster_annotation"] = (
            adata.obs.groupby("leiden", observed=True)["cluster_annotation"]
            .agg(lambda s: s.value_counts().index[0])
            .to_dict()
        )
    if "cell_state_markers_used" in adata.uns:
        summary["markers_used"] = {
            k: len(v) for k, v in adata.uns["cell_state_markers_used"].items()
        }
    return summary


def write_annotation_tables(adata: AnnData, markers: pd.DataFrame, tables: Path) -> None:
    tables.mkdir(parents=True, exist_ok=True)
    cols = [
        c
        for c in (
            "leiden",
            "phase",
            "cell_state",
            "cell_state_score",
            "cluster_annotation",
            "gene_target",
            "perturbation",
        )
        if c in adata.obs
    ]
    adata.obs[cols].to_csv(tables / "cell_annotations.csv")
    if not markers.empty:
        markers.to_csv(tables / "cluster_markers.csv", index=False)
    if "cluster_annotation" in adata.obs and "leiden" in adata.obs:
        (
            adata.obs.groupby(["leiden", "cluster_annotation"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .to_csv(tables / "cluster_annotation_summary.csv", index=False)
        )
