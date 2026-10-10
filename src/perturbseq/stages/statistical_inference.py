"""Stage 3b — statistical inference.

Owns: DE contrast orchestration and design checks for the pipeline stage.
Does NOT own: DE implementations (``statistics``) or Mixscape/E-distance (``perturbation``).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from anndata import AnnData

from perturbseq._deps import warn_if_pertpy_missing
from perturbseq.plots import plot_de_heatmap, plot_volcano
from perturbseq.report import completed_status, skipped_status
from perturbseq.stages._helpers import _ko_label, _mark_stage
from perturbseq.statistics import (
    check_experimental_design,
    resolve_de_covariates,
    run_de_contrasts,
)

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def stage_statistical_inference(
    rna: AnnData,
    config: PipelineConfig,
    report: dict,
    figures: Path,
    tables: Path,
    *,
    mixscape_ok: bool,
    edistances: pd.DataFrame,
) -> tuple[int, dict]:
    """Stage 3B: pseudobulk / DE / FDR (conceptually parallel to 3A)."""
    print("[3b/5] Statistical inference (pseudobulk DE / design checks)", flush=True)
    de_covs = resolve_de_covariates(rna, getattr(config, "de_covariates", True))
    de_prefer_pseudobulk = bool(getattr(config, "de_prefer_pseudobulk", True))
    n_pseudo_replicates = getattr(config, "n_pseudo_replicates", None)
    design = check_experimental_design(
        rna,
        replicate_col=config.replicate_col,
        groupby="gene_target",
        covariates=de_covs,
        n_pseudo_replicates=n_pseudo_replicates,
        de_prefer_pseudobulk=de_prefer_pseudobulk,
    )
    report["experimental_design"] = design
    needs_deseq2 = design.get("recommended_method", "").startswith("pydeseq2") and not config.skip_de
    if needs_deseq2:
        de_warn = warn_if_pertpy_missing(need_mixscape=False, need_distance=False, need_deseq2=True)
        if de_warn:
            warnings.warn(de_warn, UserWarning, stacklevel=2)
            report["pertpy_warning_de"] = de_warn
            print(f"WARNING: {de_warn}", flush=True)

    de_tables: list[pd.DataFrame] = []
    if config.skip_de:
        report["n_de_contrasts"] = 0
        report["de"] = skipped_status(
            "user_skip",
            detail="Differential expression skipped (--skip-de).",
            evidence_level=design.get("evidence_level"),
        )
    else:
        counts = rna.copy()
        if mixscape_ok:
            counts.obs["de_group"] = counts.obs["mixscape_class"].astype(str)
            groups = [
                g
                for g in counts.obs["de_group"].unique()
                if str(g).endswith(f" {config.perturbation_type}")
            ]
            reference = config.control
            report["de_scope"] = "mixscape_class (depends on Mixscape classification)"
            report.setdefault("matrix_provenance", {})["de"] = {
                "groupby": "mixscape_class",
                "depends_on_mixscape": True,
                "counts_layer": "layers['counts'] for PyDESeq2; log1p X for Wilcoxon",
                "covariates": de_covs,
            }
        else:
            counts.obs["de_group"] = counts.obs["gene_target"].astype(str)
            groups = [g for g in counts.obs["de_group"].unique() if g not in {config.control, "unassigned"}]
            reference = config.control
            report["de_scope"] = "gene_target (Mixscape not applied)"
            report.setdefault("matrix_provenance", {})["de"] = {
                "groupby": "gene_target",
                "depends_on_mixscape": False,
                "counts_layer": "layers['counts'] for PyDESeq2; log1p X for Wilcoxon",
                "covariates": de_covs,
            }
        if not edistances.empty:
            ranked = [g for g in edistances.index if g in set(map(str, groups))]
            if not ranked and mixscape_ok:
                gene_ranks = list(edistances.index)
                ranked = [
                    g
                    for gene in gene_ranks
                    for g in groups
                    if str(g) == _ko_label(str(gene), config.perturbation_type)
                ]
            groups = ranked[: config.de_top_n] if ranked else list(groups)[: config.de_top_n]
        else:
            groups = (
                counts.obs["de_group"].astype(str).value_counts().reindex(groups).sort_values(ascending=False).index.tolist()
            )[: config.de_top_n]
        report["de_groups"] = [str(g) for g in groups]
        report["de_covariates_used"] = de_covs
        contrast_tables, de_errors = run_de_contrasts(
            counts,
            groups,
            reference=reference,
            replicate_col=config.replicate_col,
            groupby="de_group",
            min_cells=config.min_cells_per_pert,
            n_jobs=config.n_jobs,
            covariates=de_covs,
            n_pseudo_replicates=n_pseudo_replicates,
            de_prefer_pseudobulk=de_prefer_pseudobulk,
            random_state=config.random_state,
        )
        if de_errors:
            report.setdefault("de_errors", {}).update(de_errors)
        contrast_evidence: dict[str, str] = {}
        for group, table in contrast_tables:
            safe = str(group).replace(" ", "_")
            table.to_csv(tables / f"de_{safe}.csv", index=False)
            plot_volcano(table, figures, safe)
            de_tables.append(table.head(50))
            if "evidence_level" in table.columns and len(table):
                contrast_evidence[str(group)] = str(table["evidence_level"].iloc[0])
        if de_tables:
            pd.concat(de_tables, ignore_index=True).to_csv(tables / "de_top50_concat.csv", index=False)
        if contrast_tables:
            if plot_de_heatmap(contrast_tables, figures):
                report["steps"].append("de_heatmap")
        report["n_de_contrasts"] = len(de_tables)
        report["steps"].append("de")
        de_detail = design.get("note")
        if de_detail:
            report["de_note"] = de_detail  # legacy alias for older consumers
        evidence_level = design.get("evidence_level")
        # If every contrast failed and none produced a table, treat as failed/skip-like.
        if not de_tables and de_errors:
            first_err = next(iter(de_errors.values()), "DE produced no contrasts")
            reason = "missing_pertpy" if "pertpy" in str(first_err).lower() else "failed"
            report["de"] = skipped_status(
                reason,
                detail=str(first_err),
                evidence_level=evidence_level,
            )
        else:
            report["de"] = completed_status(
                n_contrasts=len(de_tables),
                scope=report.get("de_scope"),
                groups=report.get("de_groups"),
                covariates=de_covs,
                design_formula=design.get("design_formula"),
                depends_on_mixscape=bool(mixscape_ok),
                errors=de_errors or None,
                detail=de_detail,
                evidence_level=evidence_level,
                contrast_evidence_levels=contrast_evidence or None,
            )

    _mark_stage(report, "3b_statistical_inference")
    return int(report.get("n_de_contrasts") or 0), design
