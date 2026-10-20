"""Mixscape labels, descriptive CPM effects, and perturbation distances."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from pertpy.tools import Distance, Mixscale, Mixscape, PseudobulkSpace
from scipy import sparse
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from perturbseq.assign import STATUS_NT, STATUS_TARGET, control_cell_mask
from perturbseq.config import ZERO_GUIDE_LABEL, control_spec, expected_expression
from perturbseq.io import sum_by_group

PERTURBED = "perturbed"
NON_PERTURBED = "non_perturbed"
AMBIGUOUS = "ambiguous"
NOT_EVALUATED = "not_evaluated"
SKIPPED_NO_DE = "skipped_no_de_genes"
DESCRIPTIVE_UNCERTAINTY = "descriptive_estimate_uncertainty_not_significance_test"


def perturbation_type(config: dict[str, Any]) -> str:
    mode = config["experiment"]["crispr_mode"]
    return {"CRISPR_KO": "KO", "CRISPRi": "KD", "CRISPRa": "OE"}.get(mode, "PT")


def de_gate_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Thresholds the Mixscale tutorial passes to ``RunMixscale``.

    The tutorial calls ``min.de.genes = 5``, ``logfc.threshold = 0.2``,
    ``max.de.genes = 100``, and ``fine.mode = F``, with ``split.by = NULL``.
    ``pval.cutoff`` is omitted there; the function default is 0.05.
    Mixscape's own log fold-change is not read here.
    """
    mix = config.get("mixscale", {})
    return {
        "min_de_genes": int(mix.get("min_de_genes", 5)),
        "logfc_threshold": float(mix.get("logfc_threshold", 0.2)),
        "pval_cutoff": float(mix.get("pval_cutoff", 0.05)),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
        "fine_mode": False,
    }


# CalcPerturbSig in the Mixscape vignette uses ndims=40 and num.neighbors=20.
# RunMixscape's call omits logfc.threshold; the function default is 0.25.
# https://satijalab.org/seurat/articles/mixscape_vignette
# https://satijalab.org/seurat/reference/runmixscape
VIGNETTE_SIGNATURE_NDIMS = 40
VIGNETTE_CONTROL_NEIGHBORS = 20


def mixscape_gate_parameters(config: dict[str, Any]) -> dict[str, Any]:
    """Mixscape classification thresholds from the Seurat vignette.

    ``logfc_threshold`` defaults to RunMixscape / MixscapeHeatmap's 0.25.
    Mixscale's 0.2 cutoff is not read here. ``fine_mode`` stays off.
    """
    mix = config.get("mixscape", {})
    perturb = config.get("perturb", {})
    if not isinstance(mix, dict):
        mix = {}
    if not isinstance(perturb, dict):
        perturb = {}
    return {
        "min_de_genes": int(mix.get("min_de_genes", perturb.get("mixscape_min_de_genes", 5))),
        "logfc_threshold": float(mix.get("logfc_threshold", perturb.get("mixscape_logfc", 0.25))),
        "pval_cutoff": float(mix.get("pval_cutoff", perturb.get("mixscape_pval", 0.05))),
        "max_de_genes": int(mix.get("max_de_genes", 100)),
        "fine_mode": False,
    }


def signature_n_dims(n_obs: int, n_vars: int, requested: int = VIGNETTE_SIGNATURE_NDIMS) -> int:
    """CalcPerturbSig ndims. Cap only when the matrix is smaller than the request."""
    return max(1, min(int(requested), int(n_obs) - 1, int(n_vars)))


def differential_gene_gate(
    adata: AnnData,
    *,
    control: str,
    min_de_genes: int,
    logfc_threshold: float,
    pval_cutoff: float,
    max_de_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Wilcoxon each target against NT. The table stores counts, not p-values.

    A target passes only when at least ``min_de_genes`` genes clear the absolute
    log fold-change and adjusted p-value cutoffs. ``max_de_genes`` only caps the
    gene list used later for scoring.
    """
    labels = adata.obs["target_gene"].astype(str)
    targets = [str(gene) for gene in pd.unique(labels) if str(gene) != control]
    n_nt = int((labels == control).sum())
    rows: list[dict[str, Any]] = []
    gene_rows: list[dict[str, Any]] = []
    eligible: list[str] = []
    for gene in targets:
        n_cells = int((labels == gene).sum())
        if n_cells < 2 or n_nt < 2:
            rows.append(
                _gate_row(
                    gene,
                    0,
                    min_de_genes,
                    logfc_threshold,
                    pval_cutoff,
                    max_de_genes,
                    passed=False,
                )
            )
            continue
        eligible.append(gene)
    if eligible:
        work = AnnData(
            X=adata.X,
            obs=pd.DataFrame({"target_gene": labels.to_numpy()}, index=adata.obs_names),
            var=pd.DataFrame(index=adata.var_names),
        )
        if "log1p" in adata.uns:
            work.uns["log1p"] = adata.uns["log1p"]
        scored: list[dict[str, Any]] = []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                sc.tl.rank_genes_groups(
                    work,
                    groupby="target_gene",
                    groups=eligible,
                    reference=control,
                    method="wilcoxon",
                    use_raw=False,
                )
            result = work.uns["rank_genes_groups"]
            symbol_by_name = _symbol_by_var_name(adata)
            for gene in eligible:
                names = np.asarray(result["names"][gene]).astype(str)
                signed = np.asarray(result["logfoldchanges"][gene], dtype=float)
                adjusted = np.asarray(result["pvals_adj"][gene], dtype=float)
                raw_p = np.asarray(result["pvals"][gene], dtype=float)
                keep = (np.abs(signed) >= logfc_threshold) & (adjusted < pval_cutoff)
                kept_names = names[keep]
                kept_signed = signed[keep]
                kept_adjusted = adjusted[keep]
                kept_p = raw_p[keep]
                order = np.argsort(kept_p, kind="stable")
                kept_names = kept_names[order]
                kept_signed = kept_signed[order]
                kept_adjusted = kept_adjusted[order]
                n_de = int(kept_names.size)
                passed = n_de >= min_de_genes
                scored.append(
                    _gate_row(
                        gene,
                        n_de,
                        min_de_genes,
                        logfc_threshold,
                        pval_cutoff,
                        max_de_genes,
                        passed=passed,
                    )
                )
                if not passed:
                    continue
                capped = slice(0, max_de_genes)
                for rank, name, log_fc, adj_p in zip(
                    range(int(kept_names[capped].size)),
                    kept_names[capped],
                    kept_signed[capped],
                    kept_adjusted[capped],
                    strict=True,
                ):
                    gene_rows.append(
                        {
                            "target_gene": gene,
                            "gene": str(name),
                            "gene_symbol": symbol_by_name.get(str(name), str(name)),
                            "rank": rank,
                            "logfc": float(log_fc),
                            "adj_p": float(adj_p),
                        }
                    )
        except Exception:
            scored = [
                _gate_row(
                    gene,
                    0,
                    min_de_genes,
                    logfc_threshold,
                    pval_cutoff,
                    max_de_genes,
                    passed=False,
                )
                for gene in eligible
            ]
            gene_rows = []
        rows.extend(scored)
    table = pd.DataFrame(
        rows,
        columns=[
            "target_gene",
            "n_de_genes",
            "min_de_genes",
            "logfc_threshold",
            "pval_cutoff",
            "max_de_genes",
            "fine_mode",
            "passed",
            "reason",
        ],
    )
    genes = pd.DataFrame(
        gene_rows,
        columns=["target_gene", "gene", "gene_symbol", "rank", "logfc", "adj_p"],
    )
    return table, genes


def _symbol_by_var_name(adata: AnnData) -> dict[str, str]:
    """Map matrix names to gene symbols. Guide ids are not used."""
    if "gene_symbol" not in adata.var.columns:
        return {}
    mapping: dict[str, str] = {}
    for name, symbol in zip(
        adata.var_names.astype(str), adata.var["gene_symbol"].astype(str), strict=False
    ):
        mapping.setdefault(name, symbol)
    return mapping


def _gate_row(
    gene: str,
    n_de_genes: int,
    min_de_genes: int,
    logfc_threshold: float,
    pval_cutoff: float,
    max_de_genes: int,
    *,
    passed: bool,
) -> dict[str, Any]:
    return {
        "target_gene": gene,
        "n_de_genes": int(n_de_genes),
        "min_de_genes": int(min_de_genes),
        "logfc_threshold": float(logfc_threshold),
        "pval_cutoff": float(pval_cutoff),
        "max_de_genes": int(max_de_genes),
        "fine_mode": False,
        "passed": bool(passed),
        "reason": "" if passed else SKIPPED_NO_DE,
    }


def _gene_lists(frame: pd.DataFrame) -> dict[str, list[str]]:
    lists: dict[str, list[str]] = {}
    if frame.empty:
        return lists
    for gene, group in frame.groupby("target_gene", sort=False):
        ordered = group.sort_values("rank")
        lists[str(gene)] = ordered["gene"].astype(str).tolist()
    return lists


def run_mixscape(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Gate on DE genes, then perturbation_signature and Mixscape for passing targets.

    No replicate column is invented. A real cell-type column with more than one
    label is the only optional split. Targets below ``min_de_genes`` stay
    ``skipped_no_de_genes`` and are not forced to NP. Posterior probability
    above 0.5 is the software's perturbed class. The signature UMAP is stored
    separately from the expression UMAP.
    """
    spec = control_spec(config)
    control = spec["control_label"]
    control_cells = control_cell_mask(adata.obs, spec)
    if spec["temporary"] and int(control_cells.sum()) == 0:
        adata.obs["mixscape_class"] = NOT_EVALUATED
        adata.obs["mixscape_class_global"] = NOT_EVALUATED
        adata.obs["mixscape_class_p"] = np.nan
        adata.obs["mixscape_status"] = NOT_EVALUATED
        adata.uns["mixscape_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        _stamp_control_label(adata.obs, control_cells, "target_gene", control)
        adata.uns["control_source"] = spec["control_source"]
        adata.uns["control_label"] = control
        adata.uns["control_note"] = spec["note"]
    mixscape_settings = config.get("mixscape", {})
    seed = int(config["runtime"]["seed"])
    if "highly_variable" not in adata.var:
        raise ValueError("Mixscape 需要先完成高变基因选择。")
    hvg = adata.var["highly_variable"].to_numpy()
    if int(np.sum(hvg)) < 2:
        raise ValueError("高变基因不足，无法计算 perturbation signature。")
    params = mixscape_gate_parameters(config)
    gate, gene_table = differential_gene_gate(
        adata,
        control=control,
        min_de_genes=int(params["min_de_genes"]),
        logfc_threshold=float(params["logfc_threshold"]),
        pval_cutoff=float(params["pval_cutoff"]),
        max_de_genes=int(params["max_de_genes"]),
    )
    adata.uns["mixscape_de_gate"] = gate
    adata.uns["mixscape_de_genes"] = gene_table
    labels = adata.obs["target_gene"].astype(str)
    adata.obs["mixscape_class"] = NOT_EVALUATED
    adata.obs["mixscape_class_global"] = NOT_EVALUATED
    adata.obs["mixscape_class_p"] = np.nan
    adata.obs["mixscape_status"] = NOT_EVALUATED
    skipped_genes = set(gate.loc[~gate["passed"].astype(bool), "target_gene"].astype(str))
    skipped = labels.isin(skipped_genes)
    adata.obs.loc[skipped, "mixscape_class"] = SKIPPED_NO_DE
    adata.obs.loc[skipped, "mixscape_status"] = SKIPPED_NO_DE
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    label = perturbation_type(config)
    if not passed:
        nt = labels == control
        adata.obs.loc[nt, "mixscape_class"] = control
        adata.obs.loc[nt, "mixscape_status"] = (
            ZERO_GUIDE_LABEL if control == ZERO_GUIDE_LABEL else STATUS_NT
        )
        adata.uns["mixscape_scope"] = "unstratified"
        adata.uns["mixscape_umap_status"] = "no_targets_passed_de_gate"
        adata.uns["mixscape_parameters"] = {
            "pert_key": "target_gene",
            "control": control,
            "split_by": "",
            "layer": "X_pert",
            "posterior_cutoff": 0.5,
            "perturbation_type": label,
            "random_state": seed,
            "fine_mode": False,
            **params,
        }
        return adata
    keep = labels.isin(passed | {control}).to_numpy()
    work = adata[keep, hvg].copy()
    requested_dims = int(mixscape_settings.get("n_dims", VIGNETTE_SIGNATURE_NDIMS))
    n_dims = _prepare_signature_dims(work, requested_dims, seed)
    n_nt = int((work.obs["target_gene"].astype(str) == control).sum())
    if n_nt < 1:
        raise ValueError("没有 NT 细胞，无法计算 perturbation signature。")
    requested_neighbors = int(
        mixscape_settings.get("n_control_neighbors", VIGNETTE_CONTROL_NEIGHBORS)
    )
    n_neighbors = min(requested_neighbors, n_nt, max(work.n_obs - 1, 1))
    if n_neighbors < requested_neighbors:
        adata.uns["mixscape_neighbor_note"] = (
            f"n_control_neighbors reduced from {requested_neighbors} to {n_neighbors}"
        )
    split_by = _real_cell_type_split(work, config)
    signature_args: dict[str, Any] = {
        "pert_key": "target_gene",
        "control": control,
        "n_neighbors": n_neighbors,
        "n_dims": n_dims,
    }
    if split_by:
        signature_args["split_by"] = split_by
    model = Mixscape()
    model.perturbation_signature(work, **signature_args)
    embedding, umap_status = _signature_umap(work, seed, n_dims)
    mix_args: dict[str, Any] = {
        "pert_key": "target_gene",
        "control": control,
        "layer": "X_pert",
        "de_layer": None,
        "min_de_genes": params["min_de_genes"],
        "logfc_threshold": params["logfc_threshold"],
        "pval_cutoff": params["pval_cutoff"],
        "iter_num": int(mixscape_settings.get("iter_num", 10)),
        "test_method": "wilcoxon",
        "perturbation_type": label,
        "random_state": seed,
    }
    if split_by:
        mix_args["split_by"] = split_by
    model.mixscape(work, **mix_args)
    _write_signature_rows(adata, work)
    full_embedding = np.full((adata.n_obs, 2), np.nan)
    if embedding is not None:
        positions = adata.obs_names.get_indexer(work.obs_names)
        full_embedding[positions] = embedding
        adata.obsm["X_umap_mixscape"] = full_embedding
    adata.obs.loc[work.obs_names, "mixscape_class"] = work.obs["mixscape_class"].astype(str).to_numpy()
    adata.obs.loc[work.obs_names, "mixscape_class_global"] = (
        work.obs["mixscape_class_global"].astype(str).to_numpy()
    )
    probability_key = f"mixscape_class_p_{label.lower()}"
    if probability_key in work.obs:
        adata.obs.loc[work.obs_names, "mixscape_class_p"] = work.obs[probability_key].to_numpy()
    adata.obs.loc[work.obs_names, "mixscape_status"] = [
        _normalize_mixscape(value, control, label)
        for value in work.obs["mixscape_class"].astype(str)
    ]
    if "mixscape" in work.uns:
        adata.uns["mixscape"] = work.uns["mixscape"]
    adata.uns["mixscape_scope"] = "unstratified" if not split_by else str(split_by)
    adata.uns["mixscape_umap_status"] = umap_status if embedding is not None else umap_status
    adata.uns["mixscape_parameters"] = {
        "pert_key": "target_gene",
        "control": control,
        "n_neighbors": n_neighbors,
        "n_dims": n_dims,
        "split_by": split_by or "",
        "layer": "X_pert",
        "posterior_cutoff": 0.5,
        "perturbation_type": label,
        "random_state": seed,
        "umap_neighbors_metric": "cosine",
        "umap_key": "X_umap_mixscape",
        "iter_num": int(mixscape_settings.get("iter_num", 10)),
        "test_method": "wilcoxon",
        "fine_mode": False,
        **params,
    }
    return adata


def _real_cell_type_split(adata: AnnData, config: dict[str, Any]) -> str | None:
    """Return a real cell-type column. Never invent a replicate split."""
    key = None
    for section in ("preprocess", "perturb", "mixscape"):
        candidate = config.get(section, {}).get("cell_type_key")
        if candidate:
            key = str(candidate)
            break
    if key is None and "cell_type" in adata.obs:
        labels = set(adata.obs["cell_type"].astype(str).unique())
        if labels - {"ALL", "nan", "None", "", "<NA>"}:
            key = "cell_type"
    if key is None:
        return None
    if key not in adata.obs:
        raise ValueError(f"cell_type_key={key} 不在 obs 中。")
    labels = set(adata.obs[key].astype(str).unique())
    if len(labels - {"ALL", "nan", "None", "", "<NA>"}) < 2:
        return None
    return key


def _prepare_signature_dims(work: AnnData, requested: int, seed: int) -> int:
    """Width of the PCA used to find NT neighbors.

    Pertpy uses ``X_pca`` when the gene count exceeds Scanpy's PC threshold.
    A shorter preprocess PCA is rebuilt on a z-scored copy so the search can use
    the vignette width. ``work.X`` stays log-normalized: the signature itself is
    the difference of unscaled expression. The cap applies only when the matrix
    itself is smaller.
    """
    n_dims = signature_n_dims(work.n_obs, work.n_vars, requested)
    if work.n_vars <= int(sc.settings.N_PCS):
        return n_dims
    pca_limit = min(work.n_obs - 1, work.n_vars - 1)
    target = min(n_dims, pca_limit)
    if target < 2:
        return n_dims
    have = int(work.obsm["X_pca"].shape[1]) if "X_pca" in work.obsm else 0
    if have < target:
        scaled = work.copy()
        sc.pp.scale(scaled)
        sc.pp.pca(scaled, n_comps=target, random_state=seed)
        work.obsm["X_pca"] = np.asarray(scaled.obsm["X_pca"])
        have = int(work.obsm["X_pca"].shape[1])
    return min(n_dims, have)


def _signature_umap(work: AnnData, seed: int, n_dims: int) -> tuple[np.ndarray | None, str]:
    """PCA and UMAP on the perturbation signature. Vignette RunUMAP uses dims 1:40."""
    if "X_pert" not in work.layers or work.n_obs < 4:
        return None, "not_run_too_few_cells"
    vis = work.copy()
    vis.X = vis.layers["X_pert"].copy()
    try:
        n_comps = max(2, min(int(n_dims), vis.n_obs - 1, vis.n_vars - 1))
        sc.pp.pca(vis, n_comps=n_comps, random_state=seed)
        n_neighbors = min(15, vis.n_obs - 1)
        sc.pp.neighbors(vis, n_neighbors=n_neighbors, metric="cosine", random_state=seed)
        sc.tl.umap(vis, random_state=seed)
    except Exception as exc:
        return None, f"signature_umap_failed:{exc}"
    return np.asarray(vis.obsm["X_umap"], dtype=float), "signature_space"


def _write_signature_rows(adata: AnnData, work: AnnData) -> None:
    parent_genes = adata.var_names.get_indexer(work.var_names)
    parent_cells = adata.obs_names.get_indexer(work.obs_names)
    if np.any(parent_genes < 0) or np.any(parent_cells < 0):
        raise ValueError("X_pert 无法对齐回原矩阵。")
    pert = work.layers["X_pert"]
    matrix = pert.tocoo() if sparse.issparse(pert) else sparse.coo_matrix(np.asarray(pert))
    rows = parent_cells[np.asarray(matrix.row)]
    columns = parent_genes[np.asarray(matrix.col)]
    full = sparse.coo_matrix((matrix.data, (rows, columns)), shape=(adata.n_obs, adata.n_vars))
    adata.layers["X_pert"] = full.tocsr()


def _normalize_mixscape(value: str, control: str, label: str) -> str:
    if value == control or value.endswith(f" {control}"):
        if control == ZERO_GUIDE_LABEL:
            return ZERO_GUIDE_LABEL
        return STATUS_NT
    if value.endswith(f" {label}") or value == label:
        return PERTURBED
    if value.endswith(" NP") or value == "NP":
        return NON_PERTURBED
    if value == control:
        return STATUS_NT
    return AMBIGUOUS


def run_mixscale(adata: AnnData, config: dict[str, Any]) -> AnnData:
    """Score only targets that passed the Mixscale DE gate.

    Control cells are 0. Targets with fewer than ``min_de_genes`` keep a missing
    score rather than the software fallback of 1. No replicate split is invented.
    The Mixscape classification gate is left as Mixscape wrote it.
    These scores are not written into the descriptive effect tables.
    """
    spec = control_spec(config)
    control = spec["control_label"]
    control_cells = control_cell_mask(adata.obs, spec)
    if spec["temporary"] and int(control_cells.sum()) == 0:
        adata.uns["mixscale_skip_reason"] = "no_zero_guide_cells"
        return adata
    if spec["temporary"]:
        _stamp_control_label(adata.obs, control_cells, "target_gene", control)
    params = de_gate_parameters(config)
    gate, gene_table = differential_gene_gate(
        adata,
        control=control,
        min_de_genes=int(params["min_de_genes"]),
        logfc_threshold=float(params["logfc_threshold"]),
        pval_cutoff=float(params["pval_cutoff"]),
        max_de_genes=int(params["max_de_genes"]),
    )
    adata.uns["mixscale_de_gate"] = gate
    adata.uns["mixscale_de_genes"] = gene_table
    gene_lists = _gene_lists(gene_table)
    passed = set(gate.loc[gate["passed"].astype(bool), "target_gene"].astype(str))
    passed.discard(control)
    ready = set(passed)
    skipped = sorted(set(gate["target_gene"].astype(str)) - ready - {control})
    scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    labels = adata.obs["target_gene"].astype(str)
    min_de_genes = int(params["min_de_genes"])
    if ready and "X_pert" in adata.layers:
        keep = labels.isin(ready | {control}).to_numpy()
        if "highly_variable" in adata.var and int(adata.var["highly_variable"].sum()) >= 2:
            work = adata[keep, adata.var["highly_variable"].to_numpy()].copy()
        else:
            work = adata[keep].copy()
        available = set(work.var_names.astype(str))
        supplied = {
            gene: [name for name in gene_lists.get(gene, []) if name in available][
                : int(params["max_de_genes"])
            ]
            for gene in ready
        }
        mix_args: dict[str, Any] = {
            "pert_key": "target_gene",
            "control": control,
            "layer": "X_pert",
            "min_de_genes": min_de_genes,
            "max_de_genes": int(params["max_de_genes"]),
            "logfc_threshold": float(params["logfc_threshold"]),
            "pval_cutoff": float(params["pval_cutoff"]),
            "fine_mode": False,
            "random_state": int(config["runtime"]["seed"]),
        }
        if supplied and all(len(genes) >= min_de_genes for genes in supplied.values()):
            mix_args["de_genes_by_target"] = supplied
        Mixscale().mixscale(work, **mix_args)
        written = pd.Series(work.obs["mixscale_score"].to_numpy(dtype=float), index=work.obs_names)
        scored_targets = {str(gene) for gene in work.uns.get("mixscale", {})}
        used = work.uns.get("mixscale_de_genes", {})
        work_labels = labels.reindex(written.index).astype(str)
        for gene in ready:
            markers = list(used.get(gene, ()))
            if gene not in scored_targets or len(markers) < min_de_genes:
                written.loc[work_labels == gene] = np.nan
        scores.loc[written.index] = written.to_numpy()
    scores.loc[~labels.isin(ready | {control})] = np.nan
    adata.obs["mixscale_score"] = scores.to_numpy(dtype=float)
    adata.uns["mixscale_unscored_targets"] = skipped
    adata.uns["mixscale_parameters"] = {
        "layer": "X_pert",
        "control": control,
        "split_by": "",
        "meaning": "perturbation_strength_not_editing_success",
        **params,
    }
    return adata


def mixscape_fraction(adata: AnnData) -> pd.DataFrame:
    targeting = adata.obs[adata.obs["perturbation_status"] == STATUS_TARGET]
    rows = []
    for target, group in targeting.groupby("target_gene", sort=False):
        counts = group["mixscape_status"].value_counts()
        total = int(len(group))
        rows.append(
            {
                "target_gene": target,
                "n_cells": total,
                "n_perturbed": int(counts.get(PERTURBED, 0)),
                "n_non_perturbed": int(counts.get(NON_PERTURBED, 0)),
                "n_ambiguous": int(counts.get(AMBIGUOUS, 0)),
                "fraction_perturbed": float(counts.get(PERTURBED, 0) / total) if total else 0.0,
                "fraction_non_perturbed": float(counts.get(NON_PERTURBED, 0) / total)
                if total
                else 0.0,
                "fraction_ambiguous": float(counts.get(AMBIGUOUS, 0) / total) if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def guide_pseudobulk(
    adata: AnnData, config: dict[str, Any], *, branch: str = "main"
) -> tuple[AnnData, pd.DataFrame]:
    """Sum raw counts by real guide id.

    ``main`` uses every analysis-eligible cell. ``sensitivity`` keeps perturbed
    targeting cells and the same NT cells. A low guide count does not drop the
    cell from the target aggregate; it only marks the guide ineligible.
    """
    control = control_spec(config)["control_label"]
    identified = adata.obs["guide_id"].astype(str).ne("")
    if branch == "sensitivity":
        status = adata.obs["mixscape_status"].astype(str)
        use = identified & (
            (
                (adata.obs["perturbation_status"].astype(str) == STATUS_TARGET)
                & (status == PERTURBED)
            )
            | control_cell_mask(adata.obs, control_spec(config))
        )
    else:
        use = identified
    selected = adata[use.to_numpy()].copy()
    if selected.n_obs == 0:
        empty = AnnData(np.zeros((0, adata.n_vars)), var=adata.var.copy())
        return empty, pd.DataFrame(columns=["guide_id", "target_gene", "n_cells", "total_counts"])
    group_key = selected.obs["guide_id"].astype(str)
    matrix = selected.layers["counts"]
    summed, guide_ids = sum_by_group(matrix, group_key)
    meta = (
        selected.obs.groupby("guide_id", sort=False)
        .agg(
            target_gene=("target_gene", "first"),
            n_cells=("target_gene", "size"),
            guide_name=("guide_name", "first"),
        )
        .reindex(guide_ids)
    )
    meta["total_counts"] = summed.sum(axis=1)
    meta["is_control"] = meta["target_gene"].astype(str) == control
    pseudo = AnnData(X=np.rint(summed).astype(int), obs=meta, var=selected.var.copy())
    pseudo.layers["counts_exact"] = summed
    return pseudo, meta.reset_index().rename(columns={"index": "guide_id"})


def perturbation_distance(
    pseudobulk: AnnData,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Energy distance between targets, using guide pseudobulks as the sample.

    Guides are not biological replicates. Without an explicit
    ``perturb.min_guide_replicates`` this returns empty tables instead of
    inventing a replicate count. The cell-level E-test is a separate result.
    """
    perturb_cfg = config.get("perturb", {})
    if "min_guide_replicates" not in perturb_cfg:
        return pd.DataFrame(), pd.DataFrame()
    control = control_spec(config)["control_label"]
    min_reps = int(perturb_cfg["min_guide_replicates"])
    metric = str(config["perturb"].get("distance_metric", "edistance"))
    seed = int(config["runtime"]["seed"])
    counts = pseudobulk.obs["target_gene"].astype(str).value_counts()
    keep_targets = [
        target for target, count in counts.items() if target != control and int(count) >= min_reps
    ]
    if len(keep_targets) < 3:
        return pd.DataFrame(), pd.DataFrame()
    subset = pseudobulk[pseudobulk.obs["target_gene"].astype(str).isin(keep_targets)].copy()
    exact = subset.layers["counts_exact"] if "counts_exact" in subset.layers else subset.X
    normalized = sc.pp.normalize_total(AnnData(exact), target_sum=1e4, inplace=False)["X"]
    logged = np.log1p(np.asarray(normalized))
    n_comps = max(2, min(30, subset.n_obs - 1, subset.n_vars - 1))
    embedding = PCA(n_components=n_comps, random_state=seed).fit_transform(logged)
    embedded = AnnData(logged, obs=subset.obs.copy())
    embedded.obsm["X_pca"] = embedding
    distance = Distance(metric=metric, obsm_key="X_pca")
    matrix = distance.pairwise(embedded, groupby="target_gene", show_progressbar=False, n_jobs=1)
    frame = pd.DataFrame(matrix)
    frame = frame.loc[frame.index, frame.index]
    values = frame.to_numpy(dtype=float)
    values = (values + values.T) / 2
    np.fill_diagonal(values, 0.0)
    ordered = _cluster_order(values)
    names = frame.index.to_numpy()[ordered]
    clustered = pd.DataFrame(values[np.ix_(ordered, ordered)], index=names, columns=names)
    coordinates = _embed_distance(clustered.to_numpy(dtype=float), seed)
    embedding_table = pd.DataFrame(
        {"target_gene": names, "dimension_1": coordinates[:, 0], "dimension_2": coordinates[:, 1]}
    )
    return clustered, embedding_table


def _cluster_order(values: np.ndarray) -> np.ndarray:
    if values.shape[0] < 3:
        return np.arange(values.shape[0])
    condensed = squareform(values, checks=False)
    tree = linkage(condensed, method="average")
    return leaves_list(tree)


def _embed_distance(values: np.ndarray, seed: int) -> np.ndarray:
    model = MDS(
        n_components=2,
        dissimilarity="precomputed",
        random_state=seed,
        normalized_stress="auto",
    )
    return np.asarray(model.fit_transform(values))


def target_expression(adata: AnnData, config: dict[str, Any]) -> pd.DataFrame:
    """Guide-level mean of each target gene. Cells are not treated as replicates."""
    direction, sentence = expected_expression(config)
    symbols = pd.Series(adata.var["gene_symbol"].astype(str).to_numpy(), index=adata.var_names)
    symbol_to_column = {}
    for column, symbol in enumerate(symbols):
        symbol_to_column.setdefault(symbol, column)
    nt_cells = adata.obs["perturbation_status"].astype(str) == STATUS_NT
    rows = []
    expression = adata.X
    for target, group in adata.obs[adata.obs["perturbation_status"] == STATUS_TARGET].groupby(
        "target_gene"
    ):
        column = symbol_to_column.get(str(target))
        if column is None:
            rows.append(
                {
                    "target_gene": target,
                    "gene_found": False,
                    "n_guides": int(group["guide_id"].nunique()),
                    "mean_nt": np.nan,
                    "mean_targeting": np.nan,
                    "effect_log1p": np.nan,
                    "direction": direction,
                    "direction_text": sentence,
                }
            )
            continue
        nt_values = _column(expression, np.flatnonzero(nt_cells.to_numpy()), column)
        positions = adata.obs_names.get_indexer(group.index)
        target_values = _column(expression, positions, column)
        guide_means = (
            pd.DataFrame({"guide_id": group["guide_id"].to_numpy(), "value": target_values})
            .groupby("guide_id")["value"]
            .mean()
        )
        nt_mean = float(np.mean(nt_values)) if nt_values.size else np.nan
        target_mean = float(guide_means.mean()) if len(guide_means) else np.nan
        rows.append(
            {
                "target_gene": target,
                "gene_found": True,
                "n_guides": int(guide_means.size),
                "n_cells": int(len(group)),
                "mean_nt": nt_mean,
                "mean_targeting": target_mean,
                "effect_log1p": _difference(target_mean, nt_mean),
                "direction": direction,
                "direction_text": sentence,
            }
        )
    return pd.DataFrame(rows)


def guide_consistency(adata: AnnData) -> pd.DataFrame:
    """Descriptive per-guide log2 difference of the target gene versus NT cells."""
    symbols = adata.var["gene_symbol"].astype(str)
    symbol_to_column = {}
    for column, symbol in enumerate(symbols):
        symbol_to_column.setdefault(symbol, column)
    nt_positions = np.flatnonzero(
        adata.obs["perturbation_status"].astype(str).to_numpy() == STATUS_NT
    )
    rows = []
    for target, group in adata.obs[adata.obs["perturbation_status"] == STATUS_TARGET].groupby(
        "target_gene"
    ):
        column = symbol_to_column.get(str(target))
        if column is None:
            continue
        nt_mean = (
            float(np.mean(_column(adata.X, nt_positions, column))) if nt_positions.size else np.nan
        )
        positions = adata.obs_names.get_indexer(group.index)
        values = _column(adata.X, positions, column)
        frame = pd.DataFrame(
            {
                "guide_id": group["guide_id"].astype(str).to_numpy(),
                "guide_name": group["guide_name"].astype(str).to_numpy(),
                "value": values,
            }
        )
        for guide_id, guide_rows in frame.groupby("guide_id"):
            mean_value = float(guide_rows["value"].mean())
            rows.append(
                {
                    "target_gene": target,
                    "guide_id": guide_id,
                    "guide_name": guide_rows["guide_name"].iloc[0],
                    "n_cells": int(len(guide_rows)),
                    "mean_expression": mean_value,
                    "mean_nt": nt_mean,
                    "effect_log1p": _difference(mean_value, nt_mean),
                }
            )
    return pd.DataFrame(rows)


def _column(matrix: sparse.spmatrix | np.ndarray, rows: np.ndarray, column: int) -> np.ndarray:
    if rows.size == 0:
        return np.array([], dtype=float)
    block = matrix[rows, column]
    if sparse.issparse(block):
        return np.asarray(block.toarray()).ravel()
    return np.asarray(block).ravel()


def _difference(target_mean: float, nt_mean: float) -> float:
    if not np.isfinite(target_mean) or not np.isfinite(nt_mean):
        return float("nan")
    return float(target_mean - nt_mean)


TUTORIAL_DE_COLUMNS = [
    "variable",
    "baseMean",
    "log_fc",
    "lfcSE",
    "stat",
    "p_value",
    "adj_p_value",
    "contrast",
    "branch",
    "cell_type",
]

# pertpy DGE notebook order: edgeR, then PyDESeq2. Both are pseudobulk GLMs.
TUTORIAL_PSEUDOBULK_METHODS = ("edger", "pydeseq2")


def tutorial_differential_expression(
    adata: AnnData, config: dict[str, Any], *, branch: str
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Pseudobulk by biological sample, then the notebook's EdgeR and PyDESeq2 fits.

    The pertpy differential-expression notebook aggregates counts per biological
    sample and fits ``EdgeR`` and ``PyDESeq2`` (``fit``, then ``test_contrasts``
    or ``compare_groups``). It does not run a cell-level rank test. pertpy 1.4
    builds that pseudobulk with ``PseudobulkSpace.compute``. This pipeline has
    one biological sample, so neither model is estimated and no p-values are
    written. Guides are not replicates. Results stay out of the descriptive
    effect tables.
    """
    from perturbseq.preprocess import assign_cell_type

    spec = control_spec(config)
    control = spec["control_label"]
    sample_id = str(config["experiment"]["sample_id"])
    assign_cell_type(adata, config)
    statuses: list[dict[str, Any]] = []
    frames: dict[str, list[pd.DataFrame]] = {method: [] for method in TUTORIAL_PSEUDOBULK_METHODS}
    for cell_type, subset in _analysis_populations(adata):
        mask = _contrast_mask(subset, branch, spec)
        n_cells = int(mask.sum())
        guides = (
            subset.obs.loc[mask, "guide_id"].astype(str).nunique()
            if n_cells and "guide_id" in subset.obs
            else 0
        )
        n_target = 0
        if n_cells:
            n_target = int((subset.obs.loc[mask, "target_gene"].astype(str) != control).sum())
        if n_cells == 0 or n_target == 0:
            _skip_pseudobulk(
                statuses, sample_id, branch, cell_type, n_cells, int(guides), "insufficient_cells"
            )
            continue
        work = subset[mask].copy()
        work.obs["biological_sample"] = sample_id
        n_samples = int(work.obs["biological_sample"].astype(str).nunique())
        sample_count = int(config["experiment"].get("biological_sample_count", 1))
        if sample_count < 2 or n_samples < 2:
            _skip_pseudobulk(
                statuses,
                sample_id,
                branch,
                cell_type,
                n_cells,
                int(guides),
                "single_sample_de_not_estimable",
            )
            continue
        for method in TUTORIAL_PSEUDOBULK_METHODS:
            try:
                frames[method].append(
                    _pseudobulk_contrast(method, work, control, branch, cell_type)
                )
                statuses.append(
                    _comparison_row(
                        sample_id,
                        branch,
                        cell_type,
                        "",
                        method,
                        "estimated",
                        "",
                        n_cells,
                        int(guides),
                    )
                )
            except Exception as exc:
                statuses.append(
                    _comparison_row(
                        sample_id,
                        branch,
                        cell_type,
                        "",
                        method,
                        "skipped",
                        f"model_not_estimable:{exc}",
                        n_cells,
                        int(guides),
                    )
                )
    results = {
        method: (
            pd.concat(rows, ignore_index=True)
            if rows
            else pd.DataFrame(columns=TUTORIAL_DE_COLUMNS)
        )
        for method, rows in frames.items()
    }
    return results, pd.DataFrame(statuses)


def _skip_pseudobulk(
    statuses: list[dict[str, Any]],
    sample_id: str,
    branch: str,
    cell_type: str,
    n_cells: int,
    n_guides: int,
    reason: str,
) -> None:
    for method in TUTORIAL_PSEUDOBULK_METHODS:
        statuses.append(
            _comparison_row(
                sample_id,
                branch,
                cell_type,
                "",
                method,
                "skipped",
                reason,
                n_cells,
                n_guides,
            )
        )


def _analysis_populations(adata: AnnData) -> list[tuple[str, AnnData]]:
    if str(adata.uns.get("cell_type_scope", "unstratified")) != "annotated":
        return [("ALL", adata)]
    column = str(adata.uns.get("cell_type_col", "cell_type"))
    groups = []
    for name, index in adata.obs.groupby(column, observed=True).groups.items():
        label = str(name)
        if label in {"nan", "None", "<NA>", ""}:
            continue
        groups.append((label, adata[index]))
    return groups or [("ALL", adata)]


def _contrast_mask(adata: AnnData, branch: str, spec: dict[str, Any]) -> np.ndarray:
    control = str(spec["control_label"])
    target_gene = adata.obs["target_gene"].astype(str).to_numpy()
    if "perturbation_status" in adata.obs:
        nt = control_cell_mask(adata.obs, spec)
        target = adata.obs["perturbation_status"].astype(str).to_numpy() == STATUS_TARGET
    else:
        nt = target_gene == control
        target = ~nt
    if branch == "sensitivity":
        if "mixscape_status" not in adata.obs:
            return np.zeros(adata.n_obs, dtype=bool)
        target = target & (adata.obs["mixscape_status"].astype(str).to_numpy() == PERTURBED)
    return nt | target


def _pseudobulk_contrast(
    method: str, work: AnnData, control: str, branch: str, cell_type: str
) -> pd.DataFrame:
    """Notebook call on pertpy 1.4: PseudobulkSpace sum, then fit and compare_groups.

    The published notebook aggregates with decoupler and then calls
    ``EdgeR`` / ``PyDESeq2`` ``fit`` plus ``test_contrasts`` or ``compare_groups``.
    ``compare_groups`` is the multi-group form of that contrast and itself calls
    ``test_contrasts``. pertpy 1.4 builds the same sum with ``PseudobulkSpace``.
    """
    from pertpy.tools import EdgeR, PyDESeq2

    model_cls = {"edger": EdgeR, "pydeseq2": PyDESeq2}[method]
    pdata = PseudobulkSpace().compute(
        work,
        target_col="biological_sample",
        groups_col="target_gene",
        layer_key="counts" if "counts" in work.layers else None,
        mode="sum",
    )
    groups = [group for group in pdata.obs["target_gene"].astype(str).unique() if group != control]
    model = model_cls(adata=pdata, design="~target_gene")
    model.fit()
    result = model.compare_groups(
        pdata, column="target_gene", baseline=control, groups_to_compare=groups
    )
    result = result.copy()
    result["branch"] = branch
    result["cell_type"] = cell_type
    return result


def _bootstrap_parameters(config: dict[str, Any]) -> tuple[int, float, int]:
    statistics = config.get("statistics", {})
    if not isinstance(statistics, dict):
        statistics = {}
    runtime = config.get("runtime", {})
    if not isinstance(runtime, dict):
        runtime = {}
    return (
        int(statistics.get("bootstrap_n", 200)),
        float(statistics.get("bootstrap_ci", 0.95)),
        int(runtime.get("seed", 0)),
    )


def _bootstrap_rng(seed: int, branch: str, role: str, label: str) -> np.random.Generator:
    tokens = [int(seed), *role.encode(), 0, *branch.encode(), 0, *label.encode()]
    return np.random.default_rng(np.random.SeedSequence(tokens))


def _bootstrap_pool(
    counts: sparse.spmatrix | np.ndarray,
    obs: pd.DataFrame,
    positions: np.ndarray,
    *,
    cell_only: bool,
) -> tuple[str, int, sparse.spmatrix | np.ndarray]:
    """Raw-count libraries for one side of the contrast.

    At least two distinct guide ids become one row per guide. A single guide, and
    every zero-guide control, stays one row per cell.
    """
    n_genes = int(counts.shape[1])
    positions = np.asarray(positions, dtype=int)
    if positions.size == 0:
        return "cell", 0, np.zeros((0, n_genes), dtype=np.float64)
    guide_ids = obs["guide_id"].astype(str).to_numpy()[positions]
    distinct = pd.unique(guide_ids[guide_ids != ""])
    if cell_only or int(distinct.size) < 2:
        block = counts[positions]
        if sparse.issparse(block) and not sparse.isspmatrix_csr(block):
            block = block.tocsr()
        return "cell", int(positions.size), block
    keep = guide_ids != ""
    summed, _guide_ids = sum_by_group(counts[positions[keep]], pd.Series(guide_ids[keep]))
    return "guide", int(summed.shape[0]), summed


def _resample_count_sums(
    libraries: sparse.spmatrix | np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw the observed number of units with replacement and sum their counts."""
    n_units = int(libraries.shape[0])
    n_genes = int(libraries.shape[1])
    draws = rng.integers(0, n_units, size=(n_boot, n_units))
    weights = np.zeros((n_boot, n_units), dtype=np.float64)
    np.add.at(weights, (np.arange(n_boot)[:, None], draws), 1.0)
    if sparse.issparse(libraries):
        summed = weights @ libraries
    else:
        summed = weights @ np.asarray(libraries, dtype=np.float64)
    summed = np.reshape(np.asarray(summed, dtype=np.float64), (n_boot, n_genes))
    return summed, summed.sum(axis=1)


def _replicate_log2fc(
    target_counts: np.ndarray,
    target_library: np.ndarray,
    control_counts: np.ndarray,
    control_library: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """CPM log2FC of each replicate. Both-zero raw counts stay missing."""
    target_cpm = np.zeros(target_counts.shape, dtype=np.float64)
    control_cpm = np.zeros(control_counts.shape, dtype=np.float64)
    target_ok = target_library > 0
    control_ok = control_library > 0
    target_cpm[target_ok] = (1e6 * target_counts[target_ok]) / target_library[target_ok, None]
    control_cpm[control_ok] = (1e6 * control_counts[control_ok]) / control_library[control_ok, None]
    missing = (target_counts == 0) & (control_counts == 0)
    missing |= (~target_ok)[:, None] | (~control_ok)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        log2fc = np.log2((target_cpm + alpha) / (control_cpm + alpha))
    return np.where(missing, np.nan, log2fc)


def _missing_bootstrap(
    n_genes: int, status: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.full(n_genes, np.nan),
        np.full(n_genes, np.nan),
        np.full(n_genes, status, dtype=object),
    )


def _percentile_bounds(
    replicates: np.ndarray, point: np.ndarray, level: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tail = (1.0 - float(level)) / 2.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lower, upper = np.nanpercentile(replicates, [100.0 * tail, 100.0 * (1.0 - tail)], axis=0)
    lower = np.asarray(lower, dtype=np.float64).reshape(point.shape)
    upper = np.asarray(upper, dtype=np.float64).reshape(point.shape)
    point_missing = ~np.isfinite(np.asarray(point, dtype=np.float64))
    undefined = ~point_missing & ~(np.isfinite(lower) & np.isfinite(upper))
    lower = np.where(point_missing | undefined, np.nan, lower)
    upper = np.where(point_missing | undefined, np.nan, upper)
    status = np.full(point.shape, "computed", dtype=object)
    status[point_missing] = "point_estimate_missing"
    status[undefined] = "bootstrap_interval_undefined"
    return lower, upper, status


def _effect_bootstrap(
    counts: sparse.spmatrix | np.ndarray,
    obs: pd.DataFrame,
    positions: np.ndarray,
    point: np.ndarray,
    control_libraries: sparse.spmatrix | np.ndarray,
    control_n_units: int,
    control_cache: list[tuple[np.ndarray, np.ndarray] | None],
    *,
    alpha: float,
    n_boot: int,
    level: float,
    seed: int,
    branch: str,
    target_name: str,
) -> tuple[np.ndarray, np.ndarray, str, np.ndarray]:
    """Percentile interval for one target. Both sides are redrawn each replicate."""
    target_unit, target_n_units, target_libraries = _bootstrap_pool(
        counts, obs, positions, cell_only=False
    )
    n_genes = int(np.asarray(point).shape[0])
    if target_n_units < 2 and control_n_units < 2:
        lower, upper, status = _missing_bootstrap(
            n_genes, "fewer_than_2_target_and_control_sampling_units"
        )
        return lower, upper, target_unit, status
    if target_n_units < 2:
        lower, upper, status = _missing_bootstrap(n_genes, "fewer_than_2_target_sampling_units")
        return lower, upper, target_unit, status
    if control_n_units < 2:
        lower, upper, status = _missing_bootstrap(n_genes, "fewer_than_2_control_sampling_units")
        return lower, upper, target_unit, status
    if control_cache[0] is None:
        control_cache[0] = _resample_count_sums(
            control_libraries, n_boot, _bootstrap_rng(seed, branch, "control", "control")
        )
    control_counts, control_library = control_cache[0]
    target_counts, target_library = _resample_count_sums(
        target_libraries, n_boot, _bootstrap_rng(seed, branch, "target", target_name)
    )
    replicates = _replicate_log2fc(
        target_counts, target_library, control_counts, control_library, alpha
    )
    lower, upper, status = _percentile_bounds(replicates, np.asarray(point, dtype=np.float64), level)
    return lower, upper, target_unit, status


def descriptive_effects(
    adata: AnnData, config: dict[str, Any], *, branch: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """CPM log2 fold changes for one biological sample.

    p-values, FDR, and the inferential ``ci_lower`` / ``ci_upper`` stay missing.
    A within-sample bootstrap percentile interval is written beside them.
    Without a designed NT, the denominator is the pooled CPM of zero-guide cells.
    """
    spec = control_spec(config)
    control = spec["control_label"]
    sample_id = str(config["experiment"]["sample_id"])
    cell_type = "ALL"
    columns = _effect_columns()
    nt_mask = control_cell_mask(adata.obs, spec)
    if spec["temporary"] and int(nt_mask.sum()) == 0:
        status = _comparison_row(
            sample_id, branch, cell_type, "", "effect", "skipped", "no_zero_guide_cells", 0, 0
        )
        return pd.DataFrame(columns=columns), pd.DataFrame(), pd.DataFrame([status])
    alpha = float(config.get("statistics", {}).get("pseudocount_cpm", 1.0))
    min_cells = int(config["qc"]["min_cells_per_target"])
    min_guide_cells = int(config["qc"]["min_cells_per_guide"])
    min_nt = int(config.get("statistics", {}).get("min_nt_cells", min_cells))
    target_mask = adata.obs["perturbation_status"].astype(str).to_numpy() == STATUS_TARGET
    if branch == "sensitivity":
        if "mixscape_status" not in adata.obs:
            status = _comparison_row(
                sample_id, branch, cell_type, "", "effect", "skipped", "mixscape_not_run", 0, 0
            )
            return pd.DataFrame(columns=columns), pd.DataFrame(), pd.DataFrame([status])
        target_mask = target_mask & (
            adata.obs["mixscape_status"].astype(str).to_numpy() == PERTURBED
        )
    if int(nt_mask.sum()) < min_nt:
        status = _comparison_row(
            sample_id,
            branch,
            cell_type,
            control,
            "effect",
            "skipped",
            "insufficient_nt_cells",
            int(nt_mask.sum()),
            0,
        )
        return pd.DataFrame(columns=columns), pd.DataFrame(), pd.DataFrame([status])
    selected = adata[nt_mask | target_mask].copy()
    selected_control = control_cell_mask(selected.obs, spec)
    _stamp_control_label(selected.obs, selected_control, "target_gene", control)
    counts = selected.layers["counts"] if "counts" in selected.layers else selected.X
    summed, labels = sum_by_group(counts, selected.obs["target_gene"].astype(str))
    label_index = {str(label): index for index, label in enumerate(labels)}
    if control not in label_index:
        status = _comparison_row(
            sample_id,
            branch,
            cell_type,
            control,
            "effect",
            "skipped",
            "insufficient_nt_cells",
            0,
            0,
        )
        return pd.DataFrame(columns=columns), pd.DataFrame(), pd.DataFrame([status])
    nt_counts = np.asarray(summed[label_index[control]], dtype=float).ravel()
    nt_library = float(nt_counts.sum())
    if nt_library <= 0:
        status = _comparison_row(
            sample_id, branch, cell_type, control, "effect", "skipped", "zero_library_size", 0, 0
        )
        return pd.DataFrame(columns=columns), pd.DataFrame(), pd.DataFrame([status])
    nt_cpm = 1e6 * nt_counts / nt_library
    gene_ids = selected.var["gene_id"].astype(str).to_numpy()
    gene_symbols = selected.var["gene_symbol"].astype(str).to_numpy()
    nt_obs = selected.obs[selected.obs["target_gene"].astype(str) == control]
    n_nt_cells = int(len(nt_obs))
    nt_guide_sizes = nt_obs.groupby(nt_obs["guide_id"].astype(str)).size()
    n_nt_eligible = int((nt_guide_sizes[nt_guide_sizes.index != ""] >= min_guide_cells).sum())
    n_boot, ci_level, seed = _bootstrap_parameters(config)
    control_positions = np.flatnonzero(selected.obs["target_gene"].astype(str).to_numpy() == control)
    _unit, control_n_units, control_libraries = _bootstrap_pool(
        counts,
        selected.obs,
        control_positions,
        cell_only=bool(spec["temporary"]),
    )
    control_cache: list[tuple[np.ndarray, np.ndarray] | None] = [None]
    frames = []
    guide_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for target, group in selected.obs.groupby("target_gene", sort=False):
        target_name = str(target)
        if target_name == control:
            continue
        guide_sizes = group.groupby(group["guide_id"].astype(str)).size()
        n_guides = int((guide_sizes.index != "").sum())
        n_eligible = int((guide_sizes[guide_sizes.index != ""] >= min_guide_cells).sum())
        n_cells = int(len(group))
        if n_cells < min_cells:
            comparisons.append(
                _comparison_row(
                    sample_id,
                    branch,
                    cell_type,
                    target_name,
                    "effect",
                    "skipped",
                    "insufficient_cells",
                    n_cells,
                    n_guides,
                )
            )
            continue
        target_counts = np.asarray(summed[label_index[target_name]], dtype=float).ravel()
        target_library = float(target_counts.sum())
        if target_library <= 0:
            comparisons.append(
                _comparison_row(
                    sample_id,
                    branch,
                    cell_type,
                    target_name,
                    "effect",
                    "skipped",
                    "zero_library_size",
                    n_cells,
                    n_guides,
                )
            )
            continue
        target_cpm = 1e6 * target_counts / target_library
        both_zero = (target_counts == 0) & (nt_counts == 0)
        log2fc = np.log2((target_cpm + alpha) / (nt_cpm + alpha))
        log2fc = np.where(both_zero, np.nan, log2fc)
        positions = selected.obs_names.get_indexer(group.index)
        boot_lower, boot_upper, boot_unit, boot_status = _effect_bootstrap(
            counts,
            selected.obs,
            positions,
            log2fc,
            control_libraries,
            control_n_units,
            control_cache,
            alpha=alpha,
            n_boot=n_boot,
            level=ci_level,
            seed=seed,
            branch=branch,
            target_name=target_name,
        )
        frames.append(
            pd.DataFrame(
                {
                    "sample_id": sample_id,
                    "branch": branch,
                    "cell_type": cell_type,
                    "target_gene": target_name,
                    "gene_id": gene_ids,
                    "gene_symbol": gene_symbols,
                    "n_target_cells": n_cells,
                    "n_nt_cells": n_nt_cells,
                    "n_target_guides_total": n_guides,
                    "n_target_guides_eligible": n_eligible,
                    "n_nt_guides_eligible": n_nt_eligible,
                    "target_library_size": target_library,
                    "nt_library_size": nt_library,
                    "normalization": "cpm",
                    "pseudocount_cpm": alpha,
                    "log2fc_descriptive": log2fc,
                    "pvalue": pd.NA,
                    "fdr": pd.NA,
                    "ci_lower": pd.NA,
                    "ci_upper": pd.NA,
                    "log2fc_bootstrap_ci_lower": boot_lower,
                    "log2fc_bootstrap_ci_upper": boot_upper,
                    "bootstrap_unit": boot_unit,
                    "bootstrap_n": n_boot,
                    "bootstrap_ci_level": ci_level,
                    "uncertainty_meaning": DESCRIPTIVE_UNCERTAINTY,
                    "bootstrap_status": boot_status,
                    "inference_status": "not_run_single_sample_descriptive",
                    "reason": "",
                    "effect_status": np.where(both_zero, "not_detected_both", "descriptive"),
                }
            )
        )
        guide_rows.extend(
            _guide_summaries(
                selected,
                group,
                nt_cpm,
                nt_counts,
                alpha,
                min_guide_cells,
                sample_id,
                branch,
                cell_type,
                target_name,
            )
        )
        comparisons.append(
            _comparison_row(
                sample_id,
                branch,
                cell_type,
                target_name,
                "effect",
                "descriptive",
                "not_run_single_sample_descriptive",
                n_cells,
                n_guides,
            )
        )
    if not comparisons and int(target_mask.sum()) == 0:
        comparisons.append(
            _comparison_row(
                sample_id,
                branch,
                cell_type,
                "",
                "effect",
                "skipped",
                "insufficient_cells",
                0,
                0,
            )
        )
    effects = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
    return effects, pd.DataFrame(guide_rows), pd.DataFrame(comparisons)


def target_direction_table(adata: AnnData, config: dict[str, Any], *, branch: str) -> pd.DataFrame:
    """Interpret the target gene's own RNA change under the resolved CRISPR mode."""
    spec = control_spec(config)
    if "perturbation_status" not in adata.obs or "target_gene" not in adata.obs:
        return pd.DataFrame()
    settings = config.get("target_expression", {})
    min_fraction = float(settings.get("min_nt_detection_fraction", 0.05))
    min_cpm = float(settings.get("min_nt_mean_cpm", 1.0))
    tolerance = float(settings.get("direction_log2fc_tolerance", 0.1))
    alpha = float(config.get("statistics", {}).get("pseudocount_cpm", 1.0))
    resolved = str(config["experiment"]["resolved_expected_direction"])
    mode = str(config["experiment"]["crispr_mode"])
    symbols = adata.var["gene_symbol"].astype(str).to_numpy()
    symbol_to_column = {}
    for column, symbol in enumerate(symbols):
        symbol_to_column.setdefault(symbol, column)
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    library = np.asarray(counts.sum(axis=1)).ravel()
    nt_positions = np.flatnonzero(control_cell_mask(adata.obs, spec))
    rows = []
    targeting = adata.obs[adata.obs["perturbation_status"].astype(str) == STATUS_TARGET]
    if branch == "sensitivity" and "mixscape_status" in adata.obs:
        targeting = targeting[targeting["mixscape_status"].astype(str) == PERTURBED]
    for target, group in targeting.groupby("target_gene", sort=False):
        column = symbol_to_column.get(str(target))
        status = "evaluable"
        reason = ""
        if column is None:
            status, reason = "not_evaluable", "gene_not_in_rna"
            log2fc = np.nan
            nt_fraction = np.nan
            nt_mean = np.nan
            target_mean = np.nan
        else:
            nt_counts = _column(counts, nt_positions, column)
            nt_library = library[nt_positions]
            target_positions = adata.obs_names.get_indexer(group.index)
            target_counts = _column(counts, target_positions, column)
            target_library = library[target_positions]
            nt_fraction = float(np.mean(nt_counts > 0)) if nt_counts.size else 0.0
            nt_mean = _mean_cpm(nt_counts, nt_library)
            target_mean = _mean_cpm(target_counts, target_library)
            if nt_counts.size == 0 or len(group) == 0:
                status, reason = "not_evaluable", "insufficient_cells"
            elif nt_fraction < min_fraction or nt_mean < min_cpm:
                status, reason = "not_evaluable", "nt_expression_below_threshold"
            nt_sum = float(nt_counts.sum())
            target_sum = float(target_counts.sum())
            nt_lib = float(nt_library.sum())
            target_lib = float(target_library.sum())
            if status == "evaluable" and nt_lib > 0 and target_lib > 0:
                log2fc = float(
                    np.log2(
                        ((1e6 * target_sum / target_lib) + alpha)
                        / ((1e6 * nt_sum / nt_lib) + alpha)
                    )
                )
            else:
                log2fc = np.nan
        rows.append(
            {
                "sample_id": config["experiment"]["sample_id"],
                "branch": branch,
                "cell_type": "ALL",
                "crispr_mode": mode,
                "resolved_expected_direction": resolved,
                "target_gene": target,
                "target_expression_status": status,
                "reason": reason,
                "n_target_cells": int(len(group)),
                "n_nt_cells": int(nt_positions.size),
                "nt_detection_fraction": nt_fraction,
                "nt_mean_cpm": nt_mean,
                "target_mean_cpm": target_mean,
                "log2fc_descriptive": log2fc,
                "direction_label": (
                    "not_evaluable"
                    if status != "evaluable"
                    else direction_label(float(log2fc), resolved, tolerance)
                ),
            }
        )
    return pd.DataFrame(rows)


def direction_label(log2fc: float, resolved: str, tolerance: float) -> str:
    if resolved == "not_required":
        return "direction_not_required"
    if resolved == "unknown":
        return "direction_not_assessed"
    if not np.isfinite(log2fc):
        return "not_evaluable"
    if resolved == "bidirectional":
        if log2fc > tolerance:
            return "observed_up"
        if log2fc < -tolerance:
            return "observed_down"
        return "no_clear_change"
    expected_down = resolved == "down"
    if abs(log2fc) <= tolerance:
        return "no_clear_directional_change"
    increased = log2fc > tolerance
    agrees = increased != expected_down
    return "consistent_with_expected_direction" if agrees else "opposite_direction"


def _mean_cpm(counts: np.ndarray, library: np.ndarray) -> float:
    """Mean CPM. A cell with library size 0 contributes NA, not zero expression."""
    if counts.size == 0:
        return float("nan")
    observed = np.asarray(counts)
    library_size = np.asarray(library)
    cpm = np.full(observed.shape, np.nan, dtype=np.float64)
    np.divide(1e6 * observed, library_size, out=cpm, where=library_size > 0)
    if not np.isfinite(cpm).any():
        return float("nan")
    return float(np.nanmean(cpm))


def _guide_summaries(
    selected: AnnData,
    group: pd.DataFrame,
    nt_cpm: np.ndarray,
    nt_counts: np.ndarray,
    alpha: float,
    min_guide_cells: int,
    sample_id: str,
    branch: str,
    cell_type: str,
    target_name: str,
) -> list[dict[str, Any]]:
    counts = selected.layers["counts"] if "counts" in selected.layers else selected.X
    positions = selected.obs_names.get_indexer(group.index)
    block = counts[positions]
    summed, guide_ids = sum_by_group(block, group["guide_id"].astype(str).reset_index(drop=True))
    rows = []
    for index, guide_id in enumerate(guide_ids):
        guide_counts = np.asarray(summed[index], dtype=float).ravel()
        library = float(guide_counts.sum())
        n_cells = int((group["guide_id"].astype(str) == str(guide_id)).sum())
        symbols = selected.var["gene_symbol"].astype(str).to_numpy()
        symbol_hits = np.flatnonzero(symbols == target_name)
        if library <= 0:
            median = minimum = maximum = target_lfc = np.nan
            n_up = n_down = 0
        else:
            cpm = 1e6 * guide_counts / library
            both_zero = (guide_counts == 0) & (nt_counts == 0)
            full_lfc = np.log2((cpm + alpha) / (nt_cpm + alpha))
            target_lfc = float(full_lfc[symbol_hits[0]]) if symbol_hits.size else np.nan
            if symbol_hits.size and both_zero[symbol_hits[0]]:
                target_lfc = np.nan
            log2fc = full_lfc[~both_zero]
            finite = log2fc[np.isfinite(log2fc)]
            median = float(np.median(finite)) if finite.size else np.nan
            minimum = float(np.min(finite)) if finite.size else np.nan
            maximum = float(np.max(finite)) if finite.size else np.nan
            n_up = int((finite > 0).sum())
            n_down = int((finite < 0).sum())
        rows.append(
            {
                "sample_id": sample_id,
                "branch": branch,
                "cell_type": cell_type,
                "target_gene": target_name,
                "guide_id": guide_id,
                "n_cells": n_cells,
                "library_size": library,
                "guide_eligible": n_cells >= min_guide_cells,
                "effect_median": median,
                "effect_min": minimum,
                "effect_max": maximum,
                "n_positive": n_up,
                "n_negative": n_down,
                "target_gene_log2fc": target_lfc,
                "biological_sample_count": 1,
            }
        )
    return rows


def _stamp_control_label(obs: pd.DataFrame, mask: np.ndarray, column: str, label: str) -> None:
    """Write the control label without requiring it to already be a category."""
    values = obs[column].astype(str).to_numpy(copy=True)
    values[np.asarray(mask)] = label
    obs[column] = values


def _effect_columns() -> list[str]:
    return [
        "sample_id",
        "branch",
        "cell_type",
        "target_gene",
        "gene_id",
        "gene_symbol",
        "n_target_cells",
        "n_nt_cells",
        "n_target_guides_total",
        "n_target_guides_eligible",
        "n_nt_guides_eligible",
        "target_library_size",
        "nt_library_size",
        "normalization",
        "pseudocount_cpm",
        "log2fc_descriptive",
        "pvalue",
        "fdr",
        "ci_lower",
        "ci_upper",
        "log2fc_bootstrap_ci_lower",
        "log2fc_bootstrap_ci_upper",
        "bootstrap_unit",
        "bootstrap_n",
        "bootstrap_ci_level",
        "uncertainty_meaning",
        "bootstrap_status",
        "inference_status",
        "reason",
        "effect_status",
    ]


def _comparison_row(
    sample_id: str,
    branch: str,
    cell_type: str,
    target: str,
    module: str,
    status: str,
    reason: str,
    n_cells: int,
    n_guides: int,
) -> dict[str, Any]:
    return {
        "sample_id": sample_id,
        "branch": branch,
        "cell_type": cell_type,
        "target_gene": target,
        "module": module,
        "status": status,
        "reason": reason,
        "n_cells": n_cells,
        "n_guides": n_guides,
        "biological_sample_count": 1,
    }
