"""Load one DRAGEN sample without holding the unfiltered matrix."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

GEX_TYPE = "Gene Expression"


def sample_paths(input_dir: Path, sample: str) -> dict[str, Path]:
    prefix = input_dir / f"{sample}.scRNA"
    paths = {
        "matrix": Path(f"{prefix}.filtered.matrix.mtx.gz"),
        "barcodes": Path(f"{prefix}.filtered.barcodes.tsv.gz"),
        "features": Path(f"{prefix}.filtered.features.tsv.gz"),
        "reference": Path(f"{prefix}.feature_barcode_reference.csv"),
        "assignments": Path(f"{prefix}.positive_cell_guide_assignments.csv"),
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing DRAGEN outputs:\n" + "\n".join(missing))
    return paths


def read_reference(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"id", "name"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Feature reference is missing columns: {sorted(missing)}")
    return table


def read_assignments(path: Path) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        dtype={"cell_barcode": str, "num_features": int, "feature_call": str, "num_transcripts": str},
        keep_default_na=False,
    )
    expected = {"cell_barcode", "num_features", "feature_call", "num_transcripts"}
    missing = expected - set(table.columns)
    if missing:
        raise ValueError(f"Guide assignments are missing columns: {sorted(missing)}")
    if table["cell_barcode"].duplicated().any():
        raise ValueError("Guide assignments contain duplicate cell barcodes")
    return table


def read_features(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = parts + [""] * (3 - len(parts))
            rows.append(parts[:3])
    return pd.DataFrame(rows, columns=["gene_id", "gene_symbol", "feature_type"])


def read_barcodes(path: Path) -> pd.Index:
    with gzip.open(path, "rt") as handle:
        values = [line.strip() for line in handle if line.strip()]
    return pd.Index(values, name="cell_barcode")


def _stream_gex(
    matrix_path: Path,
    *,
    gene_keep: np.ndarray,
    cell_keep: np.ndarray,
    gene_out: np.ndarray,
    cell_out: np.ndarray,
) -> sparse.csr_matrix:
    """Stream a MatrixMarket file and keep selected gene rows and cell columns."""
    n_genes = int(gene_keep.sum())
    n_cells = int(cell_keep.sum())
    rows: list[int] = []
    cols: list[int] = []
    data: list[np.ndarray] = []
    chunk_rows: list[int] = []
    chunk_cols: list[int] = []
    chunk_data: list[float] = []
    seen = 0

    def flush() -> None:
        if not chunk_data:
            return
        rows.append(np.asarray(chunk_rows, dtype=np.int32))
        cols.append(np.asarray(chunk_cols, dtype=np.int32))
        data.append(np.asarray(chunk_data, dtype=np.float32))
        chunk_rows.clear()
        chunk_cols.clear()
        chunk_data.clear()

    with gzip.open(matrix_path, "rt") as handle:
        n_declared: int | None = None
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.split()
            if n_declared is None:
                n_declared = int(parts[2])
                continue
            gene_i = int(parts[0]) - 1
            cell_i = int(parts[1]) - 1
            if gene_keep[gene_i] and cell_keep[cell_i]:
                chunk_rows.append(int(gene_out[gene_i]))
                chunk_cols.append(int(cell_out[cell_i]))
                chunk_data.append(float(parts[2]))
                if len(chunk_data) >= 2_000_000:
                    flush()
            seen += 1
            if seen % 50_000_000 == 0:
                print(f"  matrix entries read: {seen:,}", flush=True)
    flush()
    if not data:
        return sparse.csr_matrix((n_genes, n_cells), dtype=np.float32)
    matrix = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_genes, n_cells),
        dtype=np.float32,
    )
    return matrix.tocsr()


def load_analysis_matrix(input_dir: Path, sample: str, assignments: pd.DataFrame) -> tuple[AnnData, pd.DataFrame]:
    """Expression matrix for Passing cells that also appear in the filtered matrix."""
    paths = sample_paths(input_dir, sample)
    features = read_features(paths["features"])
    barcodes = read_barcodes(paths["barcodes"])
    wanted = pd.Index(assignments["cell_barcode"])
    present = wanted.isin(barcodes)
    if not bool(present.all()):
        print(
            f"  assignment barcodes missing from the matrix: {int((~present).sum()):,}",
            flush=True,
        )
    gex = (features["feature_type"] == GEX_TYPE).to_numpy()
    gene_ids = features.loc[gex, "gene_id"].astype(str)
    symbols = features.loc[gex, "gene_symbol"].astype(str)
    # Map original row index -> compact GEX row.
    gene_out = np.full(features.shape[0], -1, dtype=np.int32)
    gene_out[np.flatnonzero(gex)] = np.arange(int(gex.sum()), dtype=np.int32)

    in_matrix = np.asarray(barcodes.isin(wanted))
    cell_out = np.full(len(barcodes), -1, dtype=np.int32)
    kept_barcodes = barcodes[in_matrix]
    # Preserve assignment order for the cells we keep.
    order = pd.Index(kept_barcodes).get_indexer(wanted[np.asarray(present)])
    # cell_out is indexed by original barcode position.
    original_positions = np.flatnonzero(in_matrix)
    # Rebuild in assignment order: position original_positions[order[k]] gets column k.
    cell_out[original_positions[order]] = np.arange(int(present.sum()), dtype=np.int32)

    print(
        f"  streaming gene expression for {int(present.sum()):,} cells and {int(gex.sum()):,} genes",
        flush=True,
    )
    matrix = _stream_gex(
        paths["matrix"],
        gene_keep=gex,
        cell_keep=in_matrix,
        gene_out=gene_out,
        cell_out=cell_out,
    ).T.tocsr()
    # Columns were written in assignment order of barcodes that exist.
    ordered_barcodes = wanted[present]
    var = pd.DataFrame({"gene_id": gene_ids.to_numpy(), "gene_symbol": symbols.to_numpy()})
    var.index = pd.Index(symbols.to_numpy())
    obs = pd.DataFrame(index=pd.Index(ordered_barcodes, name="cell_barcode"))
    adata = AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata, assignments.loc[np.asarray(present)].reset_index(drop=True)
