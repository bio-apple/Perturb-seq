"""ScType marker-gene cell-type annotation (Ianevski et al.).

Pure-Python port of the official scoring rule: gene-wise z-scores, positive /
negative marker sets, marker-sensitivity weights, then Leiden-cluster majority
labels. Marker DB defaults to the bundled ScTypeDB TSV under ``data/sctype/``.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

SCTYPE_OBS_COLUMN = "sctype_classification"
SCTYPE_IMPLEMENTATION = "python ScType (Ianevski scoring)"
_DB_COLUMNS = ("tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2")
_UNKNOWN = "Unknown"


def default_sctype_db_path() -> Path:
    """Bundled ScTypeDB (repo ``data/sctype/ScTypeDB_full.tsv``)."""
    return Path(__file__).resolve().parents[2] / "data" / "sctype" / "ScTypeDB_full.tsv"


def resolve_sctype_db_path(config: dict[str, Any]) -> Path:
    configured = config.get("sctype", {}).get("db_path")
    if configured is not None and str(configured).strip():
        return Path(str(configured)).expanduser()
    return default_sctype_db_path()


def list_sctype_tissues(db_path: Path | None = None) -> list[str]:
    """Tissue names available in the marker DB."""
    path = Path(db_path) if db_path is not None else default_sctype_db_path()
    table = _read_marker_db(path)
    return sorted(table["tissueType"].astype(str).unique().tolist())


def annotate_sctype(adata: AnnData, config: dict[str, Any]) -> pd.DataFrame:
    """Assign ScType labels to Leiden clusters; store ``obs['sctype_classification']``.

    Does not modify ``obs['cell_type']``. Low-confidence clusters (top score
    ``< n_cells * low_score_fraction``) become ``Unknown``.
    """
    settings = config.get("sctype", {})
    tissue = settings.get("tissue")
    if tissue is None or not str(tissue).strip():
        raise ValueError("sctype.tissue 未设置（例如 Immune system、Lung、Brain）。")
    tissue = str(tissue).strip()
    fraction = settings.get("low_score_fraction", 0.25)
    if not isinstance(fraction, (int, float)) or isinstance(fraction, bool) or not (
        0.0 <= float(fraction) <= 1.0
    ):
        raise ValueError("sctype.low_score_fraction 必须是 0–1 之间的数（默认 0.25）。")
    fraction = float(fraction)
    if "leiden" not in adata.obs:
        raise ValueError("缺少 leiden。")

    db_path = resolve_sctype_db_path(config)
    if not db_path.is_file():
        raise ValueError(f"ScType marker DB 不存在: {db_path}")
    positive, negative = gene_sets_prepare(db_path, tissue)
    if not positive:
        available = ", ".join(list_sctype_tissues(db_path))
        raise ValueError(f"ScType tissue={tissue!r} 在 DB 中无标记基因。可选: {available}")

    matrix = _expression_matrix(adata)
    scores = sctype_score(matrix, positive, negative, scaled=False)
    clusters = adata.obs["leiden"].astype(str).to_numpy()
    order = _cluster_order(clusters)
    rows: list[dict[str, Any]] = []
    mapping: dict[str, str] = {}
    for cluster in order:
        mask = clusters == cluster
        n_cells = int(mask.sum())
        if n_cells == 0 or scores.shape[1] == 0:
            mapping[cluster] = _UNKNOWN
            rows.append(
                {
                    "leiden": cluster,
                    SCTYPE_OBS_COLUMN: _UNKNOWN,
                    "sctype_score": 0.0,
                    "n_cells": n_cells,
                    "tissue": tissue,
                    "db_path": str(db_path),
                    "implementation": SCTYPE_IMPLEMENTATION,
                }
            )
            continue
        cluster_sum = scores.loc[:, mask].sum(axis=1).sort_values(ascending=False)
        top_type = str(cluster_sum.index[0])
        top_score = float(cluster_sum.iloc[0])
        label = top_type if top_score >= n_cells * fraction else _UNKNOWN
        mapping[cluster] = label
        rows.append(
            {
                "leiden": cluster,
                SCTYPE_OBS_COLUMN: label,
                "sctype_score": top_score,
                "n_cells": n_cells,
                "tissue": tissue,
                "db_path": str(db_path),
                "implementation": SCTYPE_IMPLEMENTATION,
            }
        )

    adata.obs[SCTYPE_OBS_COLUMN] = [mapping[cluster] for cluster in clusters]
    summary = pd.DataFrame(rows)
    adata.uns["sctype"] = {
        "tissue": tissue,
        "db_path": str(db_path),
        "low_score_fraction": fraction,
        "implementation": SCTYPE_IMPLEMENTATION,
        "n_cell_types_scored": int(len(positive)),
    }
    return summary


def gene_sets_prepare(
    db_path: Path | str, tissue: str
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Load positive / negative marker lists for one tissue from ScTypeDB."""
    table = _read_marker_db(Path(db_path))
    tissue_key = str(tissue).strip()
    block = table[table["tissueType"].astype(str) == tissue_key].copy()
    if block.empty:
        return {}, {}
    positive: dict[str, list[str]] = {}
    negative: dict[str, list[str]] = {}
    for _, row in block.iterrows():
        name = str(row["cellName"]).strip()
        if not name or name in {"nan", "None"}:
            continue
        pos = _parse_genes(row.get("geneSymbolmore1", ""))
        neg = _parse_genes(row.get("geneSymbolmore2", ""))
        # Multiple rows with the same cellName: union markers.
        positive[name] = sorted(set(positive.get(name, [])) | set(pos))
        negative[name] = sorted(set(negative.get(name, [])) | set(neg))
    return positive, negative


def sctype_score(
    expression: pd.DataFrame,
    gs_positive: dict[str, list[str]],
    gs_negative: dict[str, list[str]] | None = None,
    *,
    scaled: bool = True,
) -> pd.DataFrame:
    """ScType enrichment scores (cell types × cells).

    ``expression`` is genes × cells. When ``scaled`` is False, each gene is
    z-scored across cells (official ``scaled=FALSE`` path). Gene names are
    matched case-insensitively via uppercase.
    """
    if gs_negative is None:
        gs_negative = {name: [] for name in gs_positive}
    if expression.empty or not gs_positive:
        return pd.DataFrame()

    matrix = expression.copy()
    matrix.index = matrix.index.astype(str).str.upper()
    # Collapse duplicate symbols by mean.
    if matrix.index.duplicated().any():
        matrix = matrix.groupby(level=0, sort=False).mean()

    n_types = len(gs_positive)
    marker_counts = Counter(gene for genes in gs_positive.values() for gene in genes)
    if n_types <= 1:
        sensitivity = {gene: 1.0 for gene in marker_counts}
    else:
        # scales::rescale(count, to=c(0,1), from=c(n_types, 1))
        sensitivity = {
            gene: (n_types - count) / (n_types - 1) for gene, count in marker_counts.items()
        }

    present = set(matrix.index)
    pos_sets = {
        name: [gene for gene in genes if gene in present]
        for name, genes in gs_positive.items()
    }
    neg_sets = {
        name: [gene for gene in gs_negative.get(name, []) if gene in present]
        for name in gs_positive
    }
    used = sorted(
        {
            gene
            for genes in list(pos_sets.values()) + list(neg_sets.values())
            for gene in genes
        }
    )
    if not used:
        raise ValueError("ScType 标记基因与表达矩阵无交集。")

    values = matrix.loc[used].to_numpy(dtype=float, copy=True)
    if not scaled:
        # Gene-wise z-score across cells (rows).
        means = np.nanmean(values, axis=1, keepdims=True)
        stds = np.nanstd(values, axis=1, ddof=0, keepdims=True)
        stds = np.where(stds == 0, 1.0, stds)
        values = (values - means) / stds
        values = np.nan_to_num(values, nan=0.0)

    weights = np.array([sensitivity.get(gene, 1.0) for gene in used], dtype=float)
    # Only positive markers receive sensitivity weights in the R reference.
    pos_genes = {gene for genes in pos_sets.values() for gene in genes}
    for index, gene in enumerate(used):
        if gene in pos_genes:
            values[index, :] *= weights[index]

    gene_index = {gene: index for index, gene in enumerate(used)}
    n_cells = values.shape[1]
    scores = np.zeros((len(gs_positive), n_cells), dtype=float)
    names = list(gs_positive.keys())
    for row_i, name in enumerate(names):
        pos_idx = [gene_index[gene] for gene in pos_sets[name]]
        neg_idx = [gene_index[gene] for gene in neg_sets[name]]
        if pos_idx:
            pos_sum = values[pos_idx, :].sum(axis=0) / np.sqrt(len(pos_idx))
        else:
            pos_sum = np.zeros(n_cells, dtype=float)
        if neg_idx:
            neg_sum = (-values[neg_idx, :]).sum(axis=0) / np.sqrt(len(neg_idx))
            neg_sum = np.nan_to_num(neg_sum, nan=0.0)
        else:
            neg_sum = 0.0
        scores[row_i, :] = pos_sum + neg_sum
    return pd.DataFrame(scores, index=names, columns=matrix.columns)


def _read_marker_db(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        table = pd.read_csv(path, sep="\t")
    elif suffix == ".csv":
        table = pd.read_csv(path)
    elif suffix in {".xlsx", ".xls"}:
        table = pd.read_excel(path)
    else:
        raise ValueError(f"ScType DB 支持 tsv/csv/xlsx，收到: {path}")
    missing = [column for column in _DB_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"ScType DB 缺少列: {', '.join(missing)}")
    return table


def _parse_genes(raw: Any) -> list[str]:
    text = "" if raw is None or (isinstance(raw, float) and np.isnan(raw)) else str(raw)
    text = text.replace(" ", "").replace("///", ",")
    genes: list[str] = []
    seen: set[str] = set()
    for part in text.split(","):
        gene = part.strip().upper()
        if not gene or gene in {"NA", "NAN", "NONE", "NULL"}:
            continue
        if gene not in seen:
            seen.add(gene)
            genes.append(gene)
    return genes


def _expression_matrix(adata: AnnData) -> pd.DataFrame:
    """Genes × cells from ``adata.X`` (log-normalized), unique gene symbols."""
    if "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    else:
        symbols = adata.var_names.astype(str).to_numpy()
    keep: list[int] = []
    names: list[str] = []
    seen: set[str] = set()
    for index, gene in enumerate(symbols):
        key = gene.upper()
        if key in seen or gene in {"", "nan", "None"}:
            continue
        seen.add(key)
        keep.append(index)
        names.append(gene)
    if not keep:
        raise ValueError("矩阵里没有基因符号。")
    block = adata.X[:, keep] if len(keep) != adata.n_vars else adata.X
    if sparse.issparse(block):
        values = np.asarray(block.toarray(), dtype=float)
    else:
        values = np.asarray(block, dtype=float)
    return pd.DataFrame(values.T, index=names, columns=adata.obs_names.astype(str))


def _cluster_order(clusters: np.ndarray) -> list[str]:
    return sorted(set(clusters), key=lambda value: int(value) if str(value).isdigit() else value)
