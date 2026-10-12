"""One-sample tertiary analysis from the project glossary."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from sklearn.metrics import pairwise_distances
from statsmodels.stats.multitest import multipletests

from perturbseq.guides import (
    CONTROL_KEY,
    MIN_CONTROL_CELLS,
    MIN_PERTURBATION_CELLS,
    annotate_assignments,
    assignment_summary,
    guide_frequency,
    guide_target_mapping,
    has_target_columns,
    unmatched_records,
)
from perturbseq.io import load_analysis_matrix, read_assignments, read_barcodes, read_reference, sample_paths
from perturbseq.plots import (
    composition_bars,
    guide_charts,
    perturbation_composition_heatmap,
    qc_histograms,
    umaps,
)

N_PERMS = 1000
N_HVG = 2000
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 1.0
MIXSCAPE_NEIGHBORS = 20
SEED = 0
EVIDENCE = "within_sample_cell_level"


def run(input_dir: Path, sample: str, output_dir: Path, *, n_perms: int = N_PERMS) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = output_dir / "tables"
    figures = output_dir / "figures"
    tables.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    print(f"[1/6] Read {sample}", flush=True)
    sample_paths(input_dir, sample)
    reference = read_reference(input_dir / f"{sample}.scRNA.feature_barcode_reference.csv")
    assignments = read_assignments(input_dir / f"{sample}.scRNA.positive_cell_guide_assignments.csv")
    modeling = has_target_columns(reference)
    annotated = annotate_assignments(assignments, reference)
    n_features, frequency = guide_frequency(annotated)
    n_features.to_csv(tables / "guide_n_features.csv", index=False)
    frequency.to_csv(tables / "guide_frequency.csv", index=False)
    matrix_barcodes = read_barcodes(input_dir / f"{sample}.scRNA.filtered.barcodes.tsv.gz")
    mapping = guide_target_mapping(reference)
    unmatched = unmatched_records(annotated, reference, matrix_barcodes)
    summary = assignment_summary(annotated, reference, unmatched, matrix_barcodes)
    mapping.to_csv(tables / "guide_target_mapping.tsv", sep="\t", index=False)
    unmatched.to_csv(tables / "guide_unmatched_records.tsv", sep="\t", index=False)
    summary.to_csv(tables / "guide_assignment_summary.tsv", sep="\t", index=False)
    guide_charts(n_features, frequency, figures)

    print("[2/6] Gene expression for analysis cells", flush=True)
    adata, _kept = load_analysis_matrix(input_dir, sample, assignments)
    labels = annotated.set_index("cell_barcode").loc[adata.obs_names]
    for column in (
        "num_features",
        "guide",
        "target_gene_id",
        "target_gene_name",
        "icm_class",
        "perturbation",
    ):
        adata.obs[column] = labels[column].to_numpy()
    adata.obs["mixscape_class"] = ""
    adata.obs["cluster"] = ""
    _qc_metrics(adata)
    singlet = adata.obs["num_features"].to_numpy() == 1
    qc_histograms(adata, singlet, figures)

    print("[3/6] Expression space", flush=True)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sing = adata[singlet].copy()
    n_top = min(N_HVG, sing.n_vars)
    sc.pp.highly_variable_genes(sing, n_top_genes=n_top, flavor="seurat")
    n_comps = max(2, min(50, n_top - 1, sing.n_obs - 1))
    sc.pp.pca(sing, n_comps=n_comps, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(sing, n_neighbors=min(N_NEIGHBORS, sing.n_obs - 1), use_rep="X_pca", random_state=SEED)
    sc.tl.leiden(
        sing,
        resolution=LEIDEN_RESOLUTION,
        flavor="igraph",
        n_iterations=2,
        random_state=SEED,
        key_added="cluster",
    )
    sc.tl.umap(sing, random_state=SEED)
    _write_back(adata, sing, ["cluster"])
    _write_embedding(adata, sing, "X_pca")
    _write_embedding(adata, sing, "X_umap")
    markers = _marker_table(sing)
    markers.to_csv(tables / "markers.csv", index=False)
    composition = _composition(sing)
    composition.to_csv(tables / "composition.csv", index=False)
    composition_bars(composition, figures)
    pert_composition = _perturbation_composition(sing)
    pert_composition.to_csv(tables / "perturbation_composition.csv", index=False)
    perturbation_composition_heatmap(pert_composition, figures)

    summary: dict[str, object] = {
        "sample": sample,
        "n_passing_cells": int(len(assignments)),
        "n_analysis_cells": int(adata.n_obs),
        "n_single_guide_cells": int(singlet.sum()),
        "modeling": modeling,
        "evidence": EVIDENCE,
        "n_perms": n_perms,
    }
    n_control = int((sing.obs["perturbation"] == CONTROL_KEY).sum())
    summary["n_nontargeting"] = n_control
    run_mixscape = modeling and n_control >= MIN_CONTROL_CELLS
    run_distance = modeling and n_control > 0
    summary["mixscape"] = run_mixscape
    summary["edistance"] = run_distance

    eligible = _eligible_perturbations(sing.obs["perturbation"])
    summary["n_eligible_perturbations"] = int(len(eligible))

    if run_mixscape:
        print(f"[4/6] Mixscape on {len(eligible)} perturbations", flush=True)
        _run_mixscape(sing, eligible)
        _write_back(adata, sing, ["mixscape_class"])
    else:
        print("[4/6] Mixscape skipped", flush=True)

    counts = _mixscape_counts(sing)
    counts.to_csv(tables / "mixscape_counts.csv", index=False)
    umaps(adata, figures)

    if run_distance and eligible:
        print(f"[5/6] E-distance and E-test ({n_perms} permutations)", flush=True)
        effects = _etest(sing, eligible, n_perms=n_perms)
        effects.to_csv(tables / "edistance.csv", index=False)
    else:
        print("[5/6] E-distance skipped", flush=True)
        effects = pd.DataFrame()

    if run_mixscape:
        print("[6/6] Exploratory differential expression", flush=True)
        de_table = _exploratory_de(sing)
        de_table.to_csv(tables / "exploratory_de.csv.gz", index=False)
    else:
        print("[6/6] Exploratory differential expression skipped", flush=True)

    del sing
    gc.collect()
    adata.write_h5ad(output_dir / f"{sample}.h5ad")
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {output_dir}", flush=True)


def _qc_metrics(adata: AnnData) -> None:
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if sparse.issparse(counts):
        total = np.asarray(counts.sum(axis=1)).ravel()
        detected = np.asarray((counts > 0).sum(axis=1)).ravel()
    else:
        total = np.asarray(counts.sum(axis=1)).ravel()
        detected = np.asarray((counts > 0).sum(axis=1)).ravel()
    symbols = adata.var["gene_symbol"].astype(str)
    mito = symbols.str.startswith("MT-").to_numpy()
    ribo = symbols.str.startswith(("RPS", "RPL")).to_numpy()
    mito_sum = _column_sum(counts, mito)
    ribo_sum = _column_sum(counts, ribo)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mito = np.divide(mito_sum, total, out=np.zeros_like(total), where=total > 0) * 100
        pct_ribo = np.divide(ribo_sum, total, out=np.zeros_like(total), where=total > 0) * 100
    adata.obs["total_counts"] = total
    adata.obs["n_genes"] = detected
    adata.obs["pct_mito"] = pct_mito
    adata.obs["pct_ribo"] = pct_ribo


def _column_sum(matrix, columns: np.ndarray) -> np.ndarray:
    if not columns.any():
        n = matrix.shape[0]
        return np.zeros(n, dtype=np.float64)
    block = matrix[:, columns]
    if sparse.issparse(block):
        return np.asarray(block.sum(axis=1)).ravel()
    return np.asarray(block.sum(axis=1)).ravel()


def _write_back(adata: AnnData, sing: AnnData, columns: list[str]) -> None:
    for column in columns:
        adata.obs.loc[sing.obs_names, column] = sing.obs[column].astype(str).to_numpy()


def _write_embedding(adata: AnnData, sing: AnnData, key: str) -> None:
    values = np.asarray(sing.obsm[key])
    full = np.full((adata.n_obs, values.shape[1]), np.nan, dtype=np.float32)
    full[adata.obs_names.get_indexer(sing.obs_names)] = values
    adata.obsm[key] = full


def _marker_table(sing: AnnData) -> pd.DataFrame:
    expressed = _expressed(sing.X)
    view = sing[:, expressed].copy()
    if view.n_vars == 0 or view.obs["cluster"].nunique() < 2:
        return pd.DataFrame(columns=["cluster", "gene", "score", "logfoldchange", "pval"])
    sc.tl.rank_genes_groups(view, groupby="cluster", method="wilcoxon", use_raw=False, pts=False)
    return _rank_table(view, "cluster").rename(columns={"group": "cluster"})


def _composition(sing: AnnData) -> pd.DataFrame:
    counts = sing.obs["cluster"].value_counts().rename_axis("cluster").reset_index(name="n_cells")
    counts["proportion"] = counts["n_cells"] / counts["n_cells"].sum()
    return counts.sort_values("cluster")


def _perturbation_composition(sing: AnnData) -> pd.DataFrame:
    modeled = sing.obs["perturbation"].astype(str)
    keep = modeled != ""
    if not keep.any():
        return pd.DataFrame(columns=["perturbation"])
    counts = pd.crosstab(modeled[keep], sing.obs.loc[keep, "cluster"].astype(str))
    proportions = counts.div(counts.sum(axis=1), axis=0)
    proportions.index.name = "perturbation"
    return proportions.reset_index()


def _eligible_perturbations(perturbation: pd.Series) -> list[str]:
    counts = perturbation.value_counts()
    return sorted(
        gene
        for gene, n_cells in counts.items()
        if gene not in {"", CONTROL_KEY} and int(n_cells) >= MIN_PERTURBATION_CELLS
    )


def _run_mixscape(sing: AnnData, eligible: list[str]) -> None:
    import pertpy as pt

    keep = sing.obs["perturbation"].isin([CONTROL_KEY, *eligible]).to_numpy()
    hvg = sing.var["highly_variable"].to_numpy()
    model = sing[keep, hvg].copy()
    model.obsm["X_pca"] = np.asarray(sing.obsm["X_pca"][keep])
    mixscape = pt.tl.Mixscape()
    mixscape.perturbation_signature(
        model,
        pert_key="perturbation",
        control=CONTROL_KEY,
        n_neighbors=MIXSCAPE_NEIGHBORS,
        use_rep="X_pca",
    )
    mixscape.mixscape(
        model,
        pert_key="perturbation",
        control=CONTROL_KEY,
        layer="X_pert",
        perturbation_type="KO",
        random_state=SEED,
    )
    sing.obs["mixscape_class"] = ""
    sing.obs.loc[model.obs_names, "mixscape_class"] = model.obs["mixscape_class_global"].astype(str).to_numpy()


def _mixscape_counts(sing: AnnData) -> pd.DataFrame:
    frame = sing.obs
    genes = sorted(gene for gene in frame["perturbation"].unique() if gene not in {"", CONTROL_KEY})
    rows = []
    for gene in genes:
        cells = frame["perturbation"] == gene
        labels = frame.loc[cells, "mixscape_class"]
        rows.append(
            {
                "perturbation": gene,
                "n_cells": int(cells.sum()),
                "n_ko": int((labels == "KO").sum()),
                "n_np": int((labels == "NP").sum()),
            }
        )
    control = frame["perturbation"] == CONTROL_KEY
    rows.insert(
        0,
        {
            "perturbation": CONTROL_KEY,
            "n_cells": int(control.sum()),
            "n_ko": 0,
            "n_np": 0,
        },
    )
    return pd.DataFrame(rows)


def _etest(sing: AnnData, eligible: list[str], *, n_perms: int) -> pd.DataFrame:
    from pertpy.tools._distances._distances import Edistance

    metric = Edistance()
    embedding = np.asarray(sing.obsm["X_pca"])
    labels = sing.obs["perturbation"].astype(str).to_numpy()
    rng = np.random.default_rng(SEED)
    rows = []
    for index, gene in enumerate(eligible, start=1):
        mask = (labels == gene) | (labels == CONTROL_KEY)
        emb = embedding[mask]
        in_group = labels[mask] == gene
        pwd = pairwise_distances(emb, metric="euclidean")
        total = float(pwd.sum(dtype=np.float64))
        observed = float(metric.from_precomputed(pwd, in_group, total=total))
        n_greater = 0
        n_group = int(in_group.sum())
        for _ in range(n_perms):
            perm = np.zeros(in_group.size, dtype=bool)
            perm[rng.choice(in_group.size, size=n_group, replace=False)] = True
            if metric.from_precomputed(pwd, perm, total=total) > observed:
                n_greater += 1
        rows.append(
            {
                "perturbation": gene,
                "edistance": observed,
                "n_cells": n_group,
                "pvalue": max(n_greater, 1) / n_perms,
                "evidence": EVIDENCE,
            }
        )
        if index % 25 == 0 or index == len(eligible):
            print(f"  E-test {index}/{len(eligible)}", flush=True)
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    _reject, adjusted, _alphac, _alphacbonf = multipletests(table["pvalue"].to_numpy(), method="fdr_bh")
    table["bh_etest"] = adjusted
    return table.sort_values("edistance", ascending=False)


def _exploratory_de(sing: AnnData) -> pd.DataFrame:
    labels = sing.obs["mixscape_class"].astype(str)
    groups = pd.Series(index=sing.obs_names, dtype=object)
    groups.loc[labels == "KO"] = sing.obs.loc[labels == "KO", "perturbation"].to_numpy()
    groups.loc[labels == "NT"] = CONTROL_KEY
    keep = groups.notna()
    if keep.sum() == 0 or not (groups == CONTROL_KEY).any() or not (labels == "KO").any():
        return pd.DataFrame(columns=["perturbation", "gene", "score", "logfoldchange", "pval", "bh_exploratory", "evidence"])
    view = sing[keep].copy()
    view.obs["de_group"] = groups[keep].astype(str).to_numpy()
    expressed = _expressed(view.X)
    view = view[:, expressed].copy()
    present = [gene for gene, n_cells in view.obs["de_group"].value_counts().items() if gene != CONTROL_KEY and n_cells > 0]
    if not present:
        return pd.DataFrame(columns=["perturbation", "gene", "score", "logfoldchange", "pval", "bh_exploratory", "evidence"])
    sc.tl.rank_genes_groups(
        view,
        groupby="de_group",
        groups=present,
        reference=CONTROL_KEY,
        method="wilcoxon",
        use_raw=False,
        pts=False,
    )
    table = _rank_table(view, "de_group").rename(columns={"group": "perturbation"})
    table["bh_exploratory"] = table.groupby("perturbation", observed=True)["pval"].transform(_bh)
    table["evidence"] = EVIDENCE
    return table


def _expressed(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        totals = np.asarray(matrix.sum(axis=0)).ravel()
    else:
        totals = np.asarray(matrix.sum(axis=0)).ravel()
    return totals > 0


def _rank_table(adata: AnnData, groupby: str) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    frames = []
    for group in groups:
        frames.append(
            pd.DataFrame(
                {
                    "group": group,
                    "gene": result["names"][group],
                    "score": result["scores"][group],
                    "logfoldchange": result["logfoldchanges"][group],
                    "pval": result["pvals"][group],
                }
            )
        )
    table = pd.concat(frames, ignore_index=True)
    table = table.replace([np.inf, -np.inf], np.nan)
    table = table.loc[table["gene"].notna() & (table["gene"].astype(str) != "nan")]
    return table


def _bh(pvalues: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=pvalues.index, dtype=float)
    finite = np.isfinite(pvalues.to_numpy(dtype=float))
    if finite.sum() == 0:
        return adjusted
    _reject, values, _a, _b = multipletests(pvalues.to_numpy(dtype=float)[finite], method="fdr_bh")
    adjusted.iloc[np.flatnonzero(finite)] = values
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser(description="Tertiary analysis of one DRAGEN CRISPR sample")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.input, args.sample, args.output)


if __name__ == "__main__":
    main()
