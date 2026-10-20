"""Read one DRAGEN CRISPR sample into an AnnData of gene counts."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.assign import (
    STATUS_NO_GUIDE,
    STATUS_SINGLE_NT,
    STATUS_SINGLE_TARGET,
    annotate_assignments,
    disambiguate_shared_guides,
    mito_fraction,
)


class InputMismatch(RuntimeError):
    """Matrix and assignment barcodes are not the same set."""

    def __init__(
        self,
        only_matrix: list[str],
        only_assignment: list[str],
        *,
        n_matrix: int = 0,
        n_assignment: int = 0,
    ) -> None:
        self.only_matrix = only_matrix
        self.only_assignment = only_assignment
        self.n_matrix = n_matrix
        self.n_assignment = n_assignment
        super().__init__(
            "矩阵与分配 barcode 集合不一致: "
            f"仅矩阵 {len(only_matrix)}，仅分配 {len(only_assignment)}。"
        )


def matrix_only_as_no_guide(only_matrix: list[str], only_assignment: list[str]) -> int:
    """How many matrix-only barcodes are kept as no_guide.

    Any assignment barcode missing from the matrix still fails closed. Matrix-only
    barcodes are kept only when that side is empty. This is not an intersection.
    """
    if only_assignment:
        return 0
    return len(only_matrix)


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def read_assignments(path: Path) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        dtype={
            "cell_barcode": str,
            "num_features": int,
            "feature_call": str,
            "num_transcripts": str,
        },
        keep_default_na=False,
    )
    expected = {"cell_barcode", "num_features", "feature_call", "num_transcripts"}
    missing = expected - set(table.columns)
    if missing:
        raise ValueError(f"分配文件缺少列: {sorted(missing)}")
    if table["cell_barcode"].duplicated().any():
        raise ValueError("分配文件中有重复的 cell_barcode。")
    return table


def read_reference(path: Path) -> pd.DataFrame:
    table = read_table(path)
    required = {"id", "name", "feature_type"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"参考表缺少列: {sorted(missing)}")
    return table


def read_features(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts.extend([""] * (3 - len(parts)))
            rows.append(parts[:3])
    return pd.DataFrame(rows, columns=["gene_id", "gene_symbol", "feature_type"])


def read_barcodes(path: Path) -> pd.Index:
    with gzip.open(path, "rt") as handle:
        values = [line.strip() for line in handle if line.strip()]
    if len(values) != len(set(values)):
        raise ValueError("表达矩阵 barcode 有重复。")
    return pd.Index(values, name="cell_barcode")


def matrix_header(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            rows, columns, entries = line.split()
            return int(rows), int(columns), int(entries)
    raise ValueError(f"无法读取 MatrixMarket 头: {path}")


def load_sample(
    paths: dict[str, Path],
    *,
    nt_control: str,
    feature_type_gex: str,
    feature_type_crispr: str,
    multi_guide_policy: str,
    guides_per_cell: int = 1,
    min_guide_umi: float | None,
    mito_prefix: str,
) -> tuple[AnnData, pd.DataFrame, dict[str, object]]:
    """Stream gene and CRISPR counts for cells that can enter QC.

    The MatrixMarket file is genes × cells. Only gene-expression rows of
    candidate cells are stored. CRISPR rows are kept as ``obsm['guide_umi']``.
    """
    reference = read_reference(paths["feature_reference"])
    if not (reference["feature_type"] == feature_type_crispr).all():
        unexpected = sorted(set(reference["feature_type"]) - {feature_type_crispr})
        raise ValueError(f"参考表含有未配置的 feature_type: {unexpected}")
    assignments = read_assignments(paths["assignments"])
    features = read_features(paths["features"])
    barcodes = read_barcodes(paths["barcodes"])
    n_rows, n_cols, n_entries = matrix_header(paths["matrix"])
    checks: dict[str, object] = {
        "matrix_rows": n_rows,
        "matrix_cols": n_cols,
        "matrix_entries": n_entries,
        "n_features": int(len(features)),
        "n_barcodes": int(len(barcodes)),
        "n_assignment_rows": int(len(assignments)),
        "feature_types": features["feature_type"].value_counts().to_dict(),
    }
    if n_rows != len(features) or n_cols != len(barcodes):
        raise ValueError(
            "矩阵维度与 feature / barcode 数量不一致: "
            f"matrix={n_rows}x{n_cols}, features={len(features)}, barcodes={len(barcodes)}"
        )
    matrix_set = set(barcodes.astype(str))
    assigned = set(assignments["cell_barcode"].astype(str))
    only_matrix = sorted(matrix_set - assigned)
    only_assignment = sorted(assigned - matrix_set)
    checks["n_barcodes_only_in_matrix"] = len(only_matrix)
    checks["n_barcodes_only_in_assignment"] = len(only_assignment)
    checks["n_matrix_only_as_no_guide"] = matrix_only_as_no_guide(only_matrix, only_assignment)
    checks["matrix_only_barcodes"] = [] if only_assignment else only_matrix
    if only_assignment:
        raise InputMismatch(
            only_matrix,
            only_assignment,
            n_matrix=len(barcodes),
            n_assignment=len(assignments),
        )
    gex = (features["feature_type"] == feature_type_gex).to_numpy()
    crispr = (features["feature_type"] == feature_type_crispr).to_numpy()
    mito_features = features["gene_symbol"].astype(str).str.startswith(mito_prefix).to_numpy() & gex
    if int(gex.sum()) == 0:
        raise ValueError("特征文件中没有 Gene Expression 行。")
    annotated = annotate_assignments(
        assignments,
        reference,
        nt_control=nt_control,
        multi_guide_policy=multi_guide_policy,
        guides_per_cell=guides_per_cell,
        min_guide_umi=min_guide_umi,
        matrix_barcodes=set(barcodes.astype(str)),
    )
    if only_matrix:
        annotated = pd.concat(
            [annotated, _no_guide_rows(only_matrix)], ignore_index=True
        )
    keep_status = annotated["assignment_status"].isin(
        [STATUS_SINGLE_TARGET, STATUS_SINGLE_NT]
    ) & (annotated["n_guides"].astype(int) == int(guides_per_cell))
    zero_guide = annotated["assignment_status"].eq(STATUS_NO_GUIDE) & annotated["n_guides"].astype(
        int
    ).eq(0)
    wanted = pd.Index(
        pd.concat(
            [
                annotated.loc[keep_status, "cell_barcode"],
                annotated.loc[zero_guide & ~keep_status, "cell_barcode"],
            ],
            ignore_index=True,
        ).astype(str)
    )
    checks["n_candidate_cells"] = int(len(wanted))

    gene_out = np.full(len(features), -1, dtype=np.int32)
    gene_out[np.flatnonzero(gex)] = np.arange(int(gex.sum()), dtype=np.int32)
    crispr_out = np.full(len(features), -1, dtype=np.int32)
    crispr_out[np.flatnonzero(crispr)] = np.arange(int(crispr.sum()), dtype=np.int32)
    in_matrix = np.asarray(barcodes.isin(wanted))
    original_positions = np.flatnonzero(in_matrix)
    kept_barcodes = barcodes[in_matrix]
    order = pd.Index(kept_barcodes.astype(str)).get_indexer(wanted)
    if np.any(order < 0):
        raise ValueError("候选细胞没有全部出现在矩阵 barcode 中。")
    cell_out = np.full(len(barcodes), -1, dtype=np.int32)
    cell_out[original_positions[order]] = np.arange(len(wanted), dtype=np.int32)

    print(
        f"  streaming {len(wanted):,} cells, {int(gex.sum()):,} genes, {int(crispr.sum()):,} guides",
        flush=True,
    )
    gene_matrix, guide_matrix, total_counts, n_genes, mito_counts = _stream_counts(
        paths["matrix"],
        gene_keep=gex,
        crispr_keep=crispr,
        cell_keep=in_matrix,
        gene_out=gene_out,
        crispr_out=crispr_out,
        cell_out=cell_out,
        n_genes=int(gex.sum()),
        n_guides=int(crispr.sum()),
        n_cells=len(wanted),
        n_barcodes=len(barcodes),
        mito_features=mito_features,
    )
    mito_frac = mito_fraction(mito_counts, total_counts)
    qc = pd.DataFrame(
        {
            "cell_barcode": barcodes.astype(str),
            "total_counts": total_counts,
            "n_genes_by_counts": n_genes,
            "mito_frac": mito_frac,
        }
    )
    annotated = annotated.merge(qc, on="cell_barcode", how="left")
    symbols = features.loc[gex, "gene_symbol"].astype(str).to_numpy()
    var = pd.DataFrame(
        {
            "gene_id": features.loc[gex, "gene_id"].astype(str).to_numpy(),
            "gene_symbol": symbols,
            "feature_type": feature_type_gex,
        },
        index=pd.Index(symbols),
    )
    crispr_features = features.loc[crispr, ["gene_id", "gene_symbol"]].reset_index(drop=True)
    crispr_ids = crispr_features["gene_id"].astype(str).tolist()
    crispr_names = crispr_features["gene_symbol"].astype(str).tolist()
    annotated = disambiguate_shared_guides(
        annotated,
        guide_matrix,
        wanted,
        crispr_ids,
        crispr_names,
    )
    obs = _obs_frame(annotated, wanted)
    adata = AnnData(X=gene_matrix, obs=obs, var=var)
    adata.var_names_make_unique(join="__dup")
    adata.var.index.name = None
    # Share the raw matrix until preprocessing copies it. A second in-memory copy
    # of this sample does not fit comfortably beside the streamed matrix.
    adata.layers["counts"] = adata.X
    adata.obsm["guide_umi"] = guide_matrix
    adata.uns["crispr_ids"] = crispr_ids
    adata.uns["crispr_names"] = crispr_names
    checks["guide_umi_match"] = _guide_umi_agreement(annotated, wanted)
    return adata, annotated, checks


def _no_guide_rows(barcodes: list[str]) -> pd.DataFrame:
    """Matrix cells absent from the positive assignment table."""
    rows = []
    for barcode in barcodes:
        rows.append(
            {
                "cell_barcode": barcode,
                "n_guides_called": 0,
                "n_guides": 0,
                "feature_call": "",
                "guide_list": "",
                "guide_list_json": "[]",
                "target_list_json": "[]",
                "guide_name": "",
                "guide_id": "",
                "guide_names": [],
                "guide_transcripts": [],
                "target_gene": pd.NA,
                "assignment_status": STATUS_NO_GUIDE,
                "multi_guide_subtype": pd.NA,
                "guide_umi_pass": True,
                "perturbation_status": STATUS_NO_GUIDE,
                "shared_guide_name": False,
                "in_matrix": True,
            }
        )
    return pd.DataFrame.from_records(rows)


def _obs_frame(annotated: pd.DataFrame, wanted: pd.Index) -> pd.DataFrame:
    columns = [
        "n_guides_called",
        "n_guides",
        "feature_call",
        "guide_list",
        "guide_name",
        "guide_id",
        "target_gene",
        "assignment_status",
        "perturbation_status",
        "shared_guide_name",
        "in_matrix",
        "guide_umi_observed",
    ]
    obs = annotated.set_index("cell_barcode").loc[wanted, columns].copy()
    obs.index = pd.Index(wanted, name="cell_barcode")
    obs["cell_barcode"] = pd.Index(wanted)
    obs["n_guides_called"] = obs["n_guides_called"].astype(int)
    obs["n_guides"] = obs["n_guides"].astype(int)
    obs["shared_guide_name"] = obs["shared_guide_name"].astype(bool)
    obs["in_matrix"] = obs["in_matrix"].astype(bool)
    return obs


def _stream_counts(
    path: Path,
    *,
    gene_keep: np.ndarray,
    crispr_keep: np.ndarray,
    cell_keep: np.ndarray,
    gene_out: np.ndarray,
    crispr_out: np.ndarray,
    cell_out: np.ndarray,
    n_genes: int,
    n_guides: int,
    n_cells: int,
    n_barcodes: int,
    mito_features: np.ndarray,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    total_counts = np.zeros(n_barcodes, dtype=np.float64)
    detected_genes = np.zeros(n_barcodes, dtype=np.int32)
    mito_counts = np.zeros(n_barcodes, dtype=np.float64)
    gene_rows: list[np.ndarray] = []
    gene_cols: list[np.ndarray] = []
    gene_data: list[np.ndarray] = []
    guide_rows: list[np.ndarray] = []
    guide_cols: list[np.ndarray] = []
    guide_data: list[np.ndarray] = []
    chunk_gene: list[tuple[int, int, float]] = []
    chunk_guide: list[tuple[int, int, float]] = []
    seen = 0

    def flush() -> None:
        if chunk_gene:
            gene_rows.append(np.asarray([item[0] for item in chunk_gene], dtype=np.int32))
            gene_cols.append(np.asarray([item[1] for item in chunk_gene], dtype=np.int32))
            gene_data.append(np.asarray([item[2] for item in chunk_gene], dtype=np.float32))
            chunk_gene.clear()
        if chunk_guide:
            guide_rows.append(np.asarray([item[0] for item in chunk_guide], dtype=np.int32))
            guide_cols.append(np.asarray([item[1] for item in chunk_guide], dtype=np.int32))
            guide_data.append(np.asarray([item[2] for item in chunk_guide], dtype=np.float32))
            chunk_guide.clear()

    with gzip.open(path, "rt") as handle:
        header_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            row_text, col_text, value_text = line.split()
            feature_i = int(row_text) - 1
            cell_i = int(col_text) - 1
            value = float(value_text)
            if gene_keep[feature_i]:
                total_counts[cell_i] += value
                detected_genes[cell_i] += 1
                if mito_features[feature_i]:
                    mito_counts[cell_i] += value
                if cell_keep[cell_i]:
                    chunk_gene.append((int(gene_out[feature_i]), int(cell_out[cell_i]), value))
            elif crispr_keep[feature_i] and cell_keep[cell_i]:
                chunk_guide.append((int(crispr_out[feature_i]), int(cell_out[cell_i]), value))
            if len(chunk_gene) >= 2_000_000 or len(chunk_guide) >= 2_000_000:
                flush()
            seen += 1
            if seen % 50_000_000 == 0:
                print(f"  matrix entries read: {seen:,}", flush=True)
    flush()
    gene_matrix = _assemble(gene_data, gene_rows, gene_cols, n_genes, n_cells).T.tocsr()
    guide_matrix = _assemble(guide_data, guide_rows, guide_cols, n_guides, n_cells).T.tocsr()
    return gene_matrix, guide_matrix, total_counts, detected_genes, mito_counts


def _assemble(
    data: list[np.ndarray],
    rows: list[np.ndarray],
    cols: list[np.ndarray],
    n_rows: int,
    n_cols: int,
) -> sparse.csr_matrix:
    if not data:
        return sparse.csr_matrix((n_rows, n_cols), dtype=np.float32)
    matrix = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_rows, n_cols),
        dtype=np.float32,
    )
    return matrix.tocsr()


def _guide_umi_agreement(annotated: pd.DataFrame, wanted: pd.Index) -> dict[str, float]:
    """Compare the called transcript count with the CRISPR matrix for resolved guides."""
    subset = annotated.set_index("cell_barcode").loc[wanted]
    observed = subset["guide_umi_observed"].to_numpy(dtype=float)
    called = []
    for raw in subset["guide_transcripts"]:
        if isinstance(raw, list) and len(raw) == 1:
            called.append(float(raw[0]))
        else:
            called.append(np.nan)
    called_array = np.asarray(called, dtype=float)
    ok = np.isfinite(observed) & np.isfinite(called_array) & (called_array > 0)
    if not np.any(ok):
        return {"n_compared": 0, "median_relative_difference": float("nan")}
    relative = np.abs(observed[ok] - called_array[ok]) / called_array[ok]
    return {
        "n_compared": int(ok.sum()),
        "median_relative_difference": float(np.median(relative)),
    }


def sum_by_group(
    matrix: sparse.spmatrix | np.ndarray, groups: pd.Series
) -> tuple[np.ndarray, pd.Index]:
    """Sum rows by group label. The grand total is unchanged."""
    labels = groups.astype(str)
    codes, uniques = pd.factorize(labels, sort=True)
    indicator = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float64), (codes, np.arange(len(codes)))),
        shape=(len(uniques), len(codes)),
    )
    if sparse.issparse(matrix):
        summed = indicator @ matrix.astype(np.float64)
        dense = np.asarray(summed.toarray() if sparse.issparse(summed) else summed)
    else:
        dense = np.asarray(indicator @ np.asarray(matrix, dtype=np.float64))
    return dense, pd.Index(uniques.astype(str))
