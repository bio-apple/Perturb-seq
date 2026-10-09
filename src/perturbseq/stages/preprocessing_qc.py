"""Stage 2 — preprocessing & QC.

Owns: wiring guide annotation, cell/sample QC, normalize/PCA/UMAP, optional cell annotation.
Does NOT own: the algorithms themselves (``guides``, ``qc``, ``preprocessing``, ``cell_annotation``).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from anndata import AnnData

from perturbseq.cell_annotation import annotate_cells, annotation_summary, write_annotation_tables
from perturbseq.composition import append_composition_audit
from perturbseq.guide_reassignment import run_guide_reassignment
from perturbseq.guides import annotate_guides, filter_singlets
from perturbseq.plots import plot_cell_annotation, plot_guide_composition, plot_qc, plot_umap
from perturbseq.preprocessing import preprocess_rna
from perturbseq.qc import add_qc_metrics, filter_cells, sample_qc_summary, write_qc_filter_by_guide
from perturbseq.report import completed_status, skipped_status
from perturbseq.stages._helpers import _embedding_provenance, _mark_stage

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def stage_preprocessing_qc(
    rna: AnnData,
    config: PipelineConfig,
    report: dict,
    figures: Path,
    composition_rows: list[pd.DataFrame],
) -> AnnData:
    """Stage 2: guide assignment validation, cell QC, sample QC, normalize/PCA/UMAP."""
    print("[2/5] Preprocessing & QC", flush=True)
    assignments = rna.uns.pop("_perturbseq_assignments")
    feature_ref = rna.uns.pop("_perturbseq_feature_ref")
    crispr = rna.uns.pop("_perturbseq_crispr", None)
    rna = annotate_guides(rna, assignments, feature_ref, config.control_patterns)
    rna = add_qc_metrics(rna)
    plot_qc(rna, figures)
    plot_guide_composition(rna, figures)
    report["guide_counts"] = rna.obs["num_features"].value_counts().sort_index().to_dict()
    report["gene_target_counts"] = rna.obs["gene_target"].value_counts().to_dict()
    report["sample_qc"] = sample_qc_summary(rna)
    append_composition_audit(composition_rows, rna, "loaded")

    aware = bool(getattr(config, "perturbation_aware_qc", False))
    rna, qc_log = filter_cells(
        rna,
        n_mads=config.n_mads,
        min_cells=config.min_cells,
        perturbation_aware=aware,
        control=config.control,
    )
    by_guide = qc_log.pop("qc_filter_by_guide", None)
    tables = config.output_dir / "tables"
    if isinstance(by_guide, pd.DataFrame) and not by_guide.empty:
        filter_path = tables / "qc_filter_by_guide.csv"
        write_qc_filter_by_guide(by_guide, filter_path)
        report["qc_filter_by_guide"] = str(filter_path)
        # Serializable summary for report.json (full table is on disk).
        qc_log["qc_filter_by_guide_path"] = str(filter_path)
        flagged = by_guide.loc[by_guide["cytotoxicity_qc_flag"]]
        if not flagged.empty:
            labels = (
                flagged.loc[flagged["category"] == "guide_id", "value"].astype(str).tolist()
                or flagged["value"].astype(str).tolist()
            )
            report.setdefault("composition_notes", []).append(
                f"QC cytotoxicity flags (≥90% cells removed) for: {', '.join(labels[:20])}"
                + ("…" if len(labels) > 20 else "")
            )
    report["qc"] = {k: v for k, v in qc_log.items() if k != "qc_filter_by_guide"}
    source = qc_log.get("qc_threshold_source", "all_cells")
    report.setdefault("composition_notes", []).append(
        f"QC removed {qc_log.get('n_cells_removed_qc', 0)} cells "
        f"(perturbation_aware_qc={aware}, thresholds={source}); "
        f"gene set also filtered (min_cells={config.min_cells})."
    )
    if aware and not source.startswith("control"):
        report["composition_notes"].append(
            "perturbation_aware_qc requested but fell back to global MAD "
            f"({qc_log.get('perturbation_aware_fallback', 'insufficient NT cells')})."
        )
    append_composition_audit(composition_rows, rna, "after_qc")
    n_before_singlet = int(rna.n_obs)
    rna = filter_singlets(rna, singlet_only=config.singlet_only)
    report["n_after_singlet"] = int(rna.n_obs)
    report["n_removed_nonsinglet"] = n_before_singlet - int(rna.n_obs)
    if config.singlet_only:
        report["composition_notes"].append(
            f"Singlet filter removed {report['n_removed_nonsinglet']} cells (num_features != 1)."
        )
    append_composition_audit(composition_rows, rna, "after_singlet")

    # Optional tertiary guide reassignment (compare / apply) after singlet filter.
    if getattr(config, "guide_reassign", "off") and config.guide_reassign != "off":
        if crispr is None or getattr(crispr, "n_vars", 0) == 0:
            report["guide_reassignment"] = {
                "skipped": True,
                "reason": "no_crispr_matrix",
                "mode": config.guide_reassign,
            }
        else:
            print(f"Guide reassignment ({config.guide_reassign})", flush=True)
            re_summary = run_guide_reassignment(
                rna,
                crispr,
                mode=config.guide_reassign,
                feature_ref=feature_ref,
                control_patterns=config.control_patterns,
                min_umi=float(getattr(config, "guide_reassign_min_umi", 1.0)),
                random_state=config.random_state,
                tables_dir=tables,
                figures_dir=figures,
            )
            report["guide_reassignment"] = re_summary
            report["steps"].append("guide_reassignment")
            if re_summary.get("applied"):
                report["guide_counts"] = rna.obs["num_features"].value_counts().sort_index().to_dict()
                report["gene_target_counts"] = rna.obs["gene_target"].value_counts().to_dict()
                append_composition_audit(composition_rows, rna, "after_guide_reassignment")
                report.setdefault("composition_notes", []).append(
                    f"Guide reassignment applied ({re_summary.get('applied_method')}); "
                    "DRAGEN originals kept in guide_id_dragen / gene_target_dragen."
                )
    else:
        report["guide_reassignment"] = {"skipped": True, "reason": "user_off", "mode": "off"}

    if rna.n_obs < 20:
        raise ValueError(
            f"Too few cells: only {rna.n_obs} remain after QC/singlet filters (need ≥20).\n"
            "Suggested commands / knobs:\n"
            f"  python -m perturbseq run ... --n-mads {max(config.n_mads, 8)}  # relax MAD QC\n"
            "  python -m perturbseq run ... --keep-multiplets  # keep multi-guide cells\n"
            f"  # or lower min_cells in YAML (currently {config.min_cells})\n"
            "  python -m perturbseq run ... --dry-run  # inspect resolved params"
        )

    rna = preprocess_rna(
        rna,
        n_top_genes=config.n_top_genes,
        n_pcs=config.n_pcs,
        leiden_resolution=config.leiden_resolution,
        random_state=config.random_state,
    )
    report.setdefault("matrix_provenance", {}).update(
        {
            "layers_counts": "raw_umi (preserved before normalize/log1p)",
            "X": "log1p_normalized",
            "X_pca_global": "PCA on HVG of log1p X; descriptive structure",
            "X_umap": "visualization_only; not perturbation-effect evidence",
            "leiden": "descriptive clusters; not perturbation-effect evidence",
            **_embedding_provenance(rna, "log1p_hvg", config),
        }
    )
    report["steps"].append("preprocess_hvg_pca_leiden")
    umap_color = ["leiden", "perturbation"]
    if rna.obs["gene_target"].nunique() <= 40:
        umap_color.insert(1, "gene_target")
    plot_umap(rna, figures, color=umap_color)
    report["steps"].append("umap")
    append_composition_audit(composition_rows, rna, "after_preprocess")

    if not config.skip_cell_annotation:
        print("Annotating cell cycle / cell states", flush=True)
        rna, markers = annotate_cells(rna)
        write_annotation_tables(rna, markers, config.output_dir / "tables")
        plot_cell_annotation(rna, figures)
        report["cell_annotation"] = completed_status(**annotation_summary(rna))
        report["steps"].append("cell_annotation")
        append_composition_audit(composition_rows, rna, "after_cell_annotation")
    else:
        report["cell_annotation"] = skipped_status(
            "user_skip",
            detail="Cell annotation skipped (--skip-cell-annotation).",
        )

    _mark_stage(report, "2_preprocessing_qc")
    return rna
