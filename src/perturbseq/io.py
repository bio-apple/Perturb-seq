from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse

SAMPLE_FILES = {
    "matrix": "{sample}.scRNA.filtered.matrix.mtx.gz",
    "barcodes": "{sample}.scRNA.filtered.barcodes.tsv.gz",
    "features": "{sample}.scRNA.filtered.features.tsv.gz",
    "feature_ref": "{sample}.scRNA.feature_barcode_reference.csv",
    "assignments": "{sample}.scRNA.positive_cell_guide_assignments.csv",
}

CRISPR_FEATURE_TOKENS = ("crispr", "guide")


def resolve_sample_files(input_dir: Path, sample_id: str) -> dict[str, Path]:
    files = {key: input_dir / name.format(sample=sample_id) for key, name in SAMPLE_FILES.items()}
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing DRAGEN files:\n" + "\n".join(missing))
    return files


def _read_table(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, **kwargs)


def read_barcodes(path: Path) -> pd.Index:
    table = _read_table(path, sep="\t", header=None, compression="infer")
    barcodes = table.iloc[:, 0].astype(str)
    if barcodes.duplicated().any():
        raise ValueError(f"Duplicate barcodes in {path}")
    return pd.Index(barcodes, name="cell_barcode")


def read_features(path: Path) -> pd.DataFrame:
    table = _read_table(path, sep="\t", header=None, compression="infer")
    n_cols = table.shape[1]
    if n_cols == 1:
        table.columns = ["gene_ids"]
        table["gene_symbols"] = table["gene_ids"]
        table["feature_types"] = "Gene Expression"
    elif n_cols == 2:
        table.columns = ["gene_ids", "gene_symbols"]
        table["feature_types"] = "Gene Expression"
    else:
        table = table.iloc[:, :3].copy()
        table.columns = ["gene_ids", "gene_symbols", "feature_types"]
    table["gene_ids"] = table["gene_ids"].astype(str)
    table["gene_symbols"] = table["gene_symbols"].astype(str)
    table["feature_types"] = table["feature_types"].astype(str)
    return table


def read_feature_reference(path: Path) -> pd.DataFrame:
    ref = _read_table(path)
    ref.columns = [str(col).strip() for col in ref.columns]
    required = {"id", "name"}
    missing = required.difference(ref.columns)
    if missing:
        raise ValueError(f"{path} is missing columns {sorted(missing)}")
    ref["id"] = ref["id"].astype(str)
    ref["name"] = ref["name"].astype(str)
    return ref


def read_guide_assignments(path: Path) -> pd.DataFrame:
    assignments = _read_table(path)
    assignments.columns = [str(col).strip() for col in assignments.columns]
    rename = {col: col.lower().replace(" ", "_") for col in assignments.columns}
    assignments = assignments.rename(columns=rename)
    required = {"cell_barcode", "num_features", "feature_call", "num_transcripts"}
    missing = required.difference(assignments.columns)
    if missing:
        raise ValueError(
            f"{path} is missing columns {sorted(missing)}; "
            "expected DRAGEN positive_cell_guide_assignments.csv"
        )
    assignments["cell_barcode"] = assignments["cell_barcode"].astype(str)
    assignments["num_features"] = pd.to_numeric(assignments["num_features"], errors="coerce").fillna(0).astype(int)
    assignments["feature_call"] = assignments["feature_call"].fillna("").astype(str)
    assignments["num_transcripts"] = assignments["num_transcripts"].fillna("").astype(str)
    if assignments["cell_barcode"].duplicated().any():
        raise ValueError(f"Duplicate cell_barcode values in {path}")
    return assignments.set_index("cell_barcode")


def is_crispr_feature(feature_type: str) -> bool:
    token = str(feature_type).lower()
    return any(part in token for part in CRISPR_FEATURE_TOKENS)


def load_filtered_mex(matrix_path: Path, barcodes_path: Path, features_path: Path) -> AnnData:
    """Load DRAGEN filtered MEX as cells x features.

    DRAGEN writes 10x-style sparse MTX (features x barcodes) plus prefixed filenames,
    so this does not use ``scanpy.read_10x_mtx``.
    """
    matrix = sc.read_mtx(matrix_path)
    x = matrix.X.tocsr() if sparse.issparse(matrix.X) else sparse.csr_matrix(matrix.X)
    barcodes = read_barcodes(barcodes_path)
    features = read_features(features_path)
    if x.shape[0] == len(features) and x.shape[1] == len(barcodes):
        x = x.T.tocsr()
    elif not (x.shape[0] == len(barcodes) and x.shape[1] == len(features)):
        raise ValueError(
            f"Matrix shape {tuple(x.shape)} does not match "
            f"{len(features)} features x {len(barcodes)} barcodes"
        )
    var = features.copy()
    var.index = pd.Index(var["gene_symbols"].astype(str), name=None)
    adata = AnnData(X=x, obs=pd.DataFrame(index=barcodes), var=var)
    adata.var_names_make_unique()
    return adata


def split_gene_and_crispr(adata: AnnData) -> tuple[AnnData, AnnData]:
    types = adata.var["feature_types"].astype(str)
    crispr_mask = types.map(is_crispr_feature).to_numpy()
    if crispr_mask.any() and (~crispr_mask).any():
        rna = adata[:, ~crispr_mask].copy()
        crispr = adata[:, crispr_mask].copy()
    elif crispr_mask.any():
        raise ValueError("Matrix contains only CRISPR features; gene expression is required")
    else:
        rna = adata.copy()
        crispr = adata[:, 0:0].copy()
    rna.var["feature_types"] = "Gene Expression"
    return rna, crispr


def load_dragen_sample(input_dir: Path, sample_id: str) -> tuple[AnnData, AnnData, pd.DataFrame, pd.DataFrame]:
    files = resolve_sample_files(input_dir, sample_id)
    combined = load_filtered_mex(files["matrix"], files["barcodes"], files["features"])
    rna, crispr = split_gene_and_crispr(combined)
    rna.obs["sample_id"] = sample_id
    if crispr.n_obs:
        crispr.obs["sample_id"] = sample_id
    assignments = read_guide_assignments(files["assignments"])
    feature_ref = read_feature_reference(files["feature_ref"])
    return rna, crispr, assignments, feature_ref


def write_h5ad(adata: AnnData, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if adata.var.index.name and adata.var.index.name in adata.var.columns:
        adata.var.index.name = None
    if adata.obs.index.name and adata.obs.index.name in adata.obs.columns:
        adata.obs.index.name = None
    adata.write_h5ad(path)
