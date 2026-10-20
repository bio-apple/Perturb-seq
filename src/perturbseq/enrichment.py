"""Exploratory Enrichr on a descriptive log2FC list. Not a replicate-level GSEA."""

from __future__ import annotations

from typing import Any

import pandas as pd


def enrich_perturbations(
    effects: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, str]:
    """Exploratory Enrichr from |CPM log2FC|, not a replicate-level FDR.

    The gene list is the input a preranked GSEA could use. This function does
    not run that GSEA, and Enrichr's adjusted P is not a biological-replicate FDR.
    """
    settings = config.get("enrichment", {})
    if settings.get("enabled") is False:
        return pd.DataFrame(), "disabled"
    database = settings.get("database") or settings.get("gene_sets")
    if not database:
        return pd.DataFrame(columns=["target_gene", "direction"]), "database_not_configured"
    if effects.empty or "log2fc_descriptive" not in effects.columns:
        return pd.DataFrame(), "no_effect_table"
    lfc_limit = float(config["plotting"]["volcano_lfc"])
    min_genes = int(settings.get("min_genes", 10))
    max_targets = int(settings.get("max_targets", 40))
    ranked = []
    for target, table in effects.groupby("target_gene", sort=False):
        ranked_genes = _ranked(table, lfc_limit)
        ranked.append((int((ranked_genes["direction"] != "none").sum()), target, ranked_genes))
    ranked.sort(key=lambda item: item[0], reverse=True)
    chosen = [item for item in ranked if item[0] >= min_genes][:max_targets]
    if not chosen:
        return pd.DataFrame(), "no_effect_ranked_genes"
    try:
        import gseapy as gp
    except ImportError:
        return pd.DataFrame(), "gseapy_not_installed"
    organism = {"human": "human", "mouse": "mouse"}.get(config["experiment"]["organism"])
    if organism is None:
        return pd.DataFrame(), "no_enrichr_organism"
    library = str(database)
    frames = []
    for _count, target, ranked_genes in chosen:
        for direction in ("up", "down"):
            genes = (
                ranked_genes.loc[ranked_genes["direction"] == direction, "gene"].astype(str).tolist()
            )
            if len(genes) < min_genes:
                frames.append(
                    pd.DataFrame(
                        [
                            {
                                "target_gene": target,
                                "direction": direction,
                                "term": "",
                                "overlap": 0,
                                "pvalue": pd.NA,
                                "fdr": pd.NA,
                                "genes": "",
                                "reason": "insufficient_genes",
                            }
                        ]
                    )
                )
                continue
            try:
                result = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    organism=organism,
                    outdir=None,
                    no_plot=True,
                    cutoff=1,
                )
            except Exception as exc:
                frames.append(
                    pd.DataFrame(
                        [
                            {
                                "target_gene": target,
                                "direction": direction,
                                "term": "",
                                "overlap": 0,
                                "pvalue": pd.NA,
                                "fdr": pd.NA,
                                "genes": "",
                                "reason": f"enrichr_failed: {exc}",
                            }
                        ]
                    )
                )
                continue
            table = result.results.copy()
            if table.empty:
                continue
            overlap = table["Overlap"].astype(str).str.split("/", n=1).str[0]
            frames.append(
                pd.DataFrame(
                    {
                        "target_gene": target,
                        "direction": direction,
                        "term": table["Term"].astype(str),
                        "overlap": overlap.astype(int),
                        "enrichment_pvalue": table["P-value"],
                        "enrichment_fdr": table["Adjusted P-value"],
                        "genes": table["Genes"].astype(str),
                        "reason": "tested",
                    }
                )
            )
    if not frames:
        return pd.DataFrame(), "no_enrichment_results"
    return pd.concat(frames, ignore_index=True), "tested"


def _ranked(table: pd.DataFrame, lfc_limit: float) -> pd.DataFrame:
    if table.empty or "effect_status" not in table.columns:
        return pd.DataFrame(columns=["gene", "direction"])
    usable = table[table["effect_status"].astype(str) == "descriptive"].copy()
    direction = pd.Series("none", index=usable.index)
    lfc = pd.to_numeric(usable["log2fc_descriptive"], errors="coerce")
    direction = direction.mask(lfc >= lfc_limit, "up")
    direction = direction.mask(lfc <= -lfc_limit, "down")
    gene = usable["gene_symbol"] if "gene_symbol" in usable.columns else usable["gene"]
    return pd.DataFrame({"gene": gene.astype(str), "direction": direction})
