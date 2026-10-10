"""Stage 3a — perturbation modeling.

Owns: Mixscape / E-distance / guide-consistency orchestration for the pipeline stage.
Does NOT own: Mixscape/E-distance implementations (``perturbation``) or guide QC math (``guide_qc``).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from anndata import AnnData

from perturbseq._deps import PERTPY_MISSING_MSG, is_pertpy_import_error, warn_if_pertpy_missing
from perturbseq.composition import append_composition_audit
from perturbseq.guides import run_guide_qc, write_guide_qc_tables
from perturbseq.perturbation import (
    build_perturbation_space,
    estimate_mixscape_cost,
    filter_cells_for_mixscape_targets,
    merge_mixscape_annotations,
    run_edistance,
    run_mixscape,
    select_mixscape_targets_from_edistance,
)
from perturbseq.plots import plot_guide_qc, plot_mixscape_lda, plot_umap
from perturbseq.report import completed_status, mixscape_kd_caveat, skipped_status
from perturbseq.stages._helpers import (
    _compute_secondary_distances,
    _effective_mixscape_mode,
    _embedding_provenance,
    _mark_stage,
    _mixscape_subset,
    _pca_from_layer,
    _write_edistance_outputs,
)

if TYPE_CHECKING:
    from perturbseq.pipeline import PipelineConfig


def stage_perturbation_modeling(
    rna: AnnData,
    config: PipelineConfig,
    report: dict,
    figures: Path,
    tables: Path,
    composition_rows: list[pd.DataFrame],
) -> tuple[AnnData, AnnData, pd.DataFrame, bool, dict, pd.DataFrame]:
    """Stage 3A: Mixscape, E-distance, guide consistency (conceptually parallel to 3B)."""
    print("[3a/5] Perturbation modeling (Mixscape / E-distance / guide consistency)", flush=True)
    pertpy_warn = warn_if_pertpy_missing(
        need_mixscape=not config.skip_mixscape,
        need_distance=not config.skip_distance,
        need_deseq2=False,
    )
    if pertpy_warn:
        warnings.warn(pertpy_warn, UserWarning, stacklevel=2)
        report["pertpy_warning"] = pertpy_warn
        print(f"WARNING: {pertpy_warn}", flush=True)

    n_targets = int(
        rna.obs.loc[~rna.obs["gene_target"].isin([config.control, "unassigned"]), "gene_target"].nunique()
    )
    report["n_gene_targets"] = n_targets
    n_cells = int(rna.n_obs)
    mixscape_estimate = estimate_mixscape_cost(n_cells, n_targets, config.mixscape_max_targets)
    mixscape_mode = _effective_mixscape_mode(config)
    mixscape_ok = False
    skip_mixscape = mixscape_mode == "skip"
    mixscape_selected: list[str] | None = None
    pending_auto_skip = mixscape_mode == "auto" and n_targets > config.mixscape_max_targets

    if skip_mixscape:
        report["mixscape"] = skipped_status(
            "user_skip",
            detail=(
                "Mixscape was skipped (--skip-mixscape / mixscape_mode=skip). Downstream "
                "E-distance/DE use gene_target on all post-QC cells; no KO/NP classification filter."
            ),
            n_targets=n_targets,
            mixscape_max_targets=config.mixscape_max_targets,
            mixscape_mode=mixscape_mode,
            estimate=mixscape_estimate,
        )
    elif pending_auto_skip:
        skip_mixscape = True
        detail = (
            f"{n_targets} gene targets > mixscape_max_targets={config.mixscape_max_targets} "
            "(reason: n_targets>max). Use --force-mixscape / mixscape_mode=force for all targets, "
            "or mixscape_mode=subset with --mixscape-targets / --mixscape-top-n for a partial run. "
            "Genome-scale screens are quantified with E-distance by default."
        )
        report["mixscape"] = skipped_status(
            "too_many_targets",
            detail=detail,
            n_targets=n_targets,
            mixscape_max_targets=config.mixscape_max_targets,
            mixscape_mode=mixscape_mode,
            estimate=mixscape_estimate,
        )
        print(f"WARNING: Mixscape auto-skipped — {detail}", flush=True)

    kd_caveat = mixscape_kd_caveat(config.perturbation_type)
    if kd_caveat:
        caveats = list(report.get("statistical_caveats") or [])
        if kd_caveat not in caveats:
            caveats.append(kd_caveat)
        report["statistical_caveats"] = caveats

    edistances = pd.DataFrame()
    edist_pre = pd.DataFrame()
    if config.skip_distance:
        report["edistance"] = skipped_status(
            "user_skip",
            detail="E-distance / E-test skipped (--skip-distance).",
        )
    else:
        try:
            print("Computing E-distance (pre-Mixscape / log-norm PCA)", flush=True)
            pca_src = "log1p_hvg"
            emb_prov = _embedding_provenance(rna, pca_src, config)
            edist_pre, etest_pre = run_edistance(
                rna,
                groupby="gene_target",
                contrast=config.control,
                n_perms=config.n_perms,
                min_cells=config.min_cells_per_pert,
                random_state=config.random_state,
                pca_source=pca_src,
                n_jobs=config.n_jobs,
                power_min_cells=config.etest_power_min_cells,
                n_bootstrap=config.n_bootstrap,
            )
            sec_pre, sec_status_pre = _compute_secondary_distances(rna, config, pca_src)
            if not skip_mixscape:
                _write_edistance_outputs(
                    edist_pre,
                    etest_pre,
                    tables / "pre_mixscape",
                    figures / "pre_mixscape",
                    secondary=sec_pre,
                    secondary_status=sec_status_pre,
                )
                report["pre_mixscape"] = {
                    "n_cells": int(rna.n_obs),
                    "groupby": "gene_target",
                    "embedding": "X_pca",
                    **emb_prov,
                    "depends_on_mixscape": False,
                    "tables": str(tables / "pre_mixscape"),
                    "secondary_distances": sec_status_pre,
                }
                report["edistance"] = completed_status(
                    phase="pre_mixscape",
                    top=edist_pre.head(10).to_dict(),
                    depends_on_mixscape=False,
                    **emb_prov,
                    etest_power_min_cells=config.etest_power_min_cells,
                    n_bootstrap=config.n_bootstrap,
                    secondary_distances=sec_status_pre,
                )
            else:
                _write_edistance_outputs(
                    edist_pre,
                    etest_pre,
                    tables,
                    figures,
                    secondary=sec_pre,
                    secondary_status=sec_status_pre,
                )
                edistances = edist_pre
                report["edistance_top"] = edistances.head(10).to_dict()
                report.setdefault("matrix_provenance", {}).update(emb_prov)
                report.setdefault("matrix_provenance", {})["edistance"] = {
                    "embedding": "X_pca",
                    **emb_prov,
                    "cells": "all_post_qc_singlet",
                    "depends_on_mixscape": False,
                }
                report["edistance"] = completed_status(
                    phase="primary",
                    top=edistances.head(10).to_dict(),
                    depends_on_mixscape=False,
                    **emb_prov,
                    etest_power_min_cells=config.etest_power_min_cells,
                    n_bootstrap=config.n_bootstrap,
                    secondary_distances=sec_status_pre,
                )
            report["steps"].append("edistance_pre_or_only")
        except ImportError as exc:
            if is_pertpy_import_error(exc):
                msg = str(exc) if str(exc) else PERTPY_MISSING_MSG
                warnings.warn(msg, UserWarning, stacklevel=2)
                report["edistance"] = skipped_status("missing_pertpy", detail=msg)
                print(f"WARNING: E-distance skipped — {msg}", flush=True)
            else:
                raise
        except Exception as exc:  # noqa: BLE001
            report["edistance"] = skipped_status("failed", detail=str(exc))

    # Resolve subset target list (user list and/or top-N by pre-Mixscape E-distance).
    if not skip_mixscape and mixscape_mode == "subset":
        selected = [str(t) for t in config.mixscape_targets]
        top_n = config.mixscape_top_n
        if not selected:
            top_n = int(top_n) if top_n is not None else int(config.mixscape_max_targets)
        if top_n is not None and top_n > 0:
            ranked = select_mixscape_targets_from_edistance(
                edist_pre,
                top_n=int(top_n),
                control=config.control,
            )
            for gene in ranked:
                if gene not in selected:
                    selected.append(gene)
        if not selected:
            skip_mixscape = True
            detail = (
                "mixscape_mode=subset requires --mixscape-targets and/or --mixscape-top-n "
                "(needs a pre-Mixscape E-distance ranking; avoid --skip-distance)."
            )
            report["mixscape"] = skipped_status(
                "subset_no_targets",
                detail=detail,
                n_targets=n_targets,
                mixscape_max_targets=config.mixscape_max_targets,
                mixscape_mode=mixscape_mode,
                estimate=mixscape_estimate,
            )
            print(f"WARNING: Mixscape subset skipped — {detail}", flush=True)
        else:
            mixscape_selected = selected
            mixscape_estimate = estimate_mixscape_cost(
                n_cells=int(
                    (rna.obs["gene_target"].astype(str).isin([config.control, *selected])).sum()
                ),
                n_targets=len(selected),
                mixscape_max_targets=config.mixscape_max_targets,
            )
            print(
                f"Mixscape subset: {len(selected)} targets + {config.control} "
                f"(full library has {n_targets} targets)",
                flush=True,
            )

    analysis_source = rna
    if not skip_mixscape:
        try:
            mixscape_adata = rna
            if mixscape_selected is not None:
                mixscape_adata = filter_cells_for_mixscape_targets(
                    rna, config.control, mixscape_selected
                )
            run_mixscape(
                mixscape_adata,
                control=config.control,
                split_by=config.replicate_col,
                perturbation_type=config.perturbation_type,
            )
            if mixscape_selected is not None and mixscape_adata is not rna:
                merge_mixscape_annotations(rna, mixscape_adata)
            mixscape_ok = True
            analysis_source = mixscape_adata
            global_counts = mixscape_adata.obs["mixscape_class_global"].value_counts().to_dict()
            report["mixscape_global"] = global_counts
            mixscape_meta = dict(mixscape_adata.uns.get("_perturbseq_mixscape") or {})
            lda_ok = bool(mixscape_meta.get("lda")) and "mixscape_lda" in mixscape_adata.uns
            subset_detail = ""
            if mixscape_selected is not None:
                subset_detail = (
                    f" Subset Mixscape on {len(mixscape_selected)} targets + {config.control}; "
                    "labels/results do not cover the full gene library."
                )
            report["mixscape"] = completed_status(
                global_counts=global_counts,
                detail=(
                    "Mixscape succeeded. Primary E-distance/DE use KO+control cells and depend on "
                    "mixscape_class; pre_mixscape/ retains the unfiltered gene_target path."
                    + subset_detail
                ),
                n_targets=n_targets,
                mixscape_max_targets=config.mixscape_max_targets,
                mixscape_mode=mixscape_mode,
                subset=bool(mixscape_selected),
                selected_targets=list(mixscape_selected) if mixscape_selected else None,
                estimate=mixscape_estimate if mixscape_selected else None,
                lda=lda_ok,
                mixscale=bool(mixscape_meta.get("mixscale")),
                # KO / NP / control labels (Seurat mixscape_class.global analogue)
                class_global_labels=sorted(str(k) for k in global_counts),
            )
            plot_umap(
                mixscape_adata,
                figures / "mixscape",
                color=["mixscape_class_global", "perturbation"],
            )
            if lda_ok:
                plotted = plot_mixscape_lda(
                    mixscape_adata,
                    figures / "mixscape",
                    control=config.control,
                    perturbation_type=config.perturbation_type,
                )
                if not plotted:
                    warnings.warn(
                        "Mixscape LDA embedding present but lda_umap.png was not written",
                        UserWarning,
                        stacklevel=2,
                    )
            report["steps"].append("mixscape")
            append_composition_audit(composition_rows, mixscape_adata, "after_mixscape_classify")
        except ImportError as exc:
            if is_pertpy_import_error(exc):
                msg = str(exc) if str(exc) else PERTPY_MISSING_MSG
                warnings.warn(msg, UserWarning, stacklevel=2)
                detail = (
                    f"Mixscape skipped — pertpy missing. {msg} "
                    "Downstream uses gene_target on all post-QC cells; no KO/NP filter."
                )
                report["mixscape"] = skipped_status(
                    "missing_pertpy",
                    detail=detail,
                    n_targets=n_targets,
                    mixscape_max_targets=config.mixscape_max_targets,
                    mixscape_mode=mixscape_mode,
                    estimate=mixscape_estimate,
                )
                report["steps"].append("mixscape_skipped_missing_pertpy")
                print(f"WARNING: Mixscape skipped — {msg}", flush=True)
            else:
                raise
        except Exception as exc:  # noqa: BLE001
            detail = (
                f"Mixscape failed ({exc}). Downstream uses gene_target on all post-QC cells; "
                "no KO/NP filter."
            )
            report["mixscape"] = skipped_status(
                "failed",
                detail=detail,
                n_targets=n_targets,
                mixscape_max_targets=config.mixscape_max_targets,
                mixscape_mode=mixscape_mode,
                estimate=mixscape_estimate,
            )
            report["steps"].append("mixscape_failed")

    analysis_obj = _mixscape_subset(analysis_source, config.control) if mixscape_ok else rna
    pca_source = "log1p_hvg"
    if mixscape_ok:
        n_before_ko = int(analysis_source.n_obs)
        n_after_ko = int(analysis_obj.n_obs)
        before_counts = analysis_source.obs["mixscape_class_global"].astype(str).value_counts().to_dict()
        after_counts = analysis_obj.obs["mixscape_class_global"].astype(str).value_counts().to_dict()
        report["n_after_mixscape_ko_filter"] = n_after_ko
        report["n_removed_mixscape_non_ko"] = n_before_ko - n_after_ko
        mixscape_comp = {
            "n_before": n_before_ko,
            "n_after": n_after_ko,
            "n_removed_np_or_other": n_before_ko - n_after_ko,
            "mixscape_class_global_before": before_counts,
            "mixscape_class_global_after": after_counts,
            "inferential_subset": f"mixscape_class_global in {{{config.control}, KO}}",
        }
        report["mixscape_ko_filter"] = mixscape_comp
        if isinstance(report.get("mixscape"), dict) and report["mixscape"].get("skipped") is False:
            report["mixscape"]["ko_filter"] = mixscape_comp
        pd.DataFrame(
            [
                {"stage": "before_ko_filter", "label": k, "n_cells": int(v)}
                for k, v in before_counts.items()
            ]
            + [
                {"stage": "after_ko_filter", "label": k, "n_cells": int(v)}
                for k, v in after_counts.items()
            ]
        ).to_csv(tables / "mixscape_ko_filter_composition.csv", index=False)
        report.setdefault("composition_notes", []).append(
            f"Mixscape KO+{config.control} filter kept {n_after_ko}/{n_before_ko} cells "
            f"(removed {n_before_ko - n_after_ko} NP/other); "
            f"this KO-filtered set is the default inferential subset for E-distance."
        )
        if mixscape_selected is not None:
            report["composition_notes"].append(
                f"Mixscape subset mode: classified {len(mixscape_selected)} targets only "
                f"(library has {n_targets}); non-selected genes lack mixscape labels."
            )
        append_composition_audit(composition_rows, analysis_obj, "after_mixscape_ko_filter")
        if "X_pert" in analysis_obj.layers:
            _pca_from_layer(analysis_obj, "X_pert", config.n_pcs)
            pca_source = "X_pert"
            emb_post = _embedding_provenance(analysis_obj, pca_source, config)
            report.setdefault("matrix_provenance", {}).update(
                {
                    "X_pca_post_mixscape": "PCA recomputed on layers['X_pert']",
                    "X_pert": "Mixscape perturbation signature (pipeline perturbation embedding)",
                    "perturbation_embedding": "layers['X_pert']",
                    **emb_post,
                }
            )
            report["perturbation_embedding"] = {
                "key": "X_pert",
                "slot": "layers",
                "source": "Mixscape perturbation_signature",
                "used_for": ["mixscape_classification", "edistance_pca"],
                "note": (
                    "X_pert is the built-in low-dim perturbation signature used for "
                    "Mixscape and post-Mixscape E-distance PCA. External tools "
                    "(SCEPTRE / MIMOSCA / PerturbNet) are documented, not installed."
                ),
            }

    if mixscape_ok and not config.skip_distance:
        try:
            print("Computing E-distance (post-Mixscape / X_pert PCA)", flush=True)
            emb_prov = _embedding_provenance(analysis_obj, pca_source, config)
            edistances, etest = run_edistance(
                analysis_obj,
                groupby="gene_target",
                contrast=config.control,
                n_perms=config.n_perms,
                min_cells=config.min_cells_per_pert,
                random_state=config.random_state,
                pca_source=pca_source,
                n_jobs=config.n_jobs,
                power_min_cells=config.etest_power_min_cells,
                n_bootstrap=config.n_bootstrap,
            )
            sec_post, sec_status_post = _compute_secondary_distances(
                analysis_obj, config, pca_source
            )
            _write_edistance_outputs(
                edistances,
                etest,
                tables / "post_mixscape",
                figures / "post_mixscape",
                secondary=sec_post,
                secondary_status=sec_status_post,
            )
            _write_edistance_outputs(
                edistances,
                etest,
                tables,
                figures,
                secondary=sec_post,
                secondary_status=sec_status_post,
            )
            report["edistance_top"] = edistances.head(10).to_dict()
            report["post_mixscape"] = {
                "n_cells": int(analysis_obj.n_obs),
                "groupby": "gene_target",
                "embedding": "X_pca",
                **emb_prov,
                "depends_on_mixscape": True,
                "tables": str(tables / "post_mixscape"),
                "secondary_distances": sec_status_post,
            }
            report.setdefault("matrix_provenance", {}).update(emb_prov)
            report.setdefault("matrix_provenance", {})["edistance"] = {
                "embedding": "X_pca",
                **emb_prov,
                "cells": f"mixscape_class_global in {{{config.control}, KO}}",
                "depends_on_mixscape": True,
                "primary_tables": str(tables / "edistance.csv"),
            }
            report["edistance"] = completed_status(
                phase="post_mixscape",
                top=edistances.head(10).to_dict(),
                depends_on_mixscape=True,
                **emb_prov,
                etest_power_min_cells=config.etest_power_min_cells,
                n_bootstrap=config.n_bootstrap,
                secondary_distances=sec_status_post,
            )
            report["steps"].append("edistance_post_mixscape")
        except ImportError as exc:
            if is_pertpy_import_error(exc):
                msg = str(exc) if str(exc) else PERTPY_MISSING_MSG
                warnings.warn(msg, UserWarning, stacklevel=2)
                # Keep pre_mixscape results if any; mark primary path skipped.
                pre = report.get("edistance") if isinstance(report.get("edistance"), dict) else {}
                report["edistance"] = skipped_status(
                    "missing_pertpy",
                    detail=msg,
                )
                if pre.get("skipped") is False:
                    report["edistance"]["pre_phase"] = pre
                print(f"WARNING: post-Mixscape E-distance skipped — {msg}", flush=True)
            else:
                raise
        except Exception as exc:  # noqa: BLE001
            pre = report.get("edistance") if isinstance(report.get("edistance"), dict) else {}
            report["edistance"] = skipped_status("failed", detail=str(exc))
            if pre.get("skipped") is False:
                report["edistance"]["pre_phase"] = pre

    space_method = str(getattr(config, "perturbation_space", "pca_silhouette") or "pca_silhouette")
    try:
        space = build_perturbation_space(
            analysis_obj,
            method=space_method,
            groupby="gene_target",
            random_state=config.random_state,
        )
        space.clusters.to_csv(tables / "perturbation_clusters.csv")
        if space.embeddings is not None and not space.embeddings.empty:
            space.embeddings.to_csv(tables / "perturbation_space_embeddings.csv")
        analysis_obj.uns["perturbation_space"] = {
            "method": space.method,
            "backend": space.backend,
            **{
                k: v
                for k, v in space.metadata.items()
                if isinstance(v, (str, int, float, bool, type(None)))
            },
        }
        report.setdefault("matrix_provenance", {})["perturbation_clusters"] = {
            "method": space.method,
            "backend": space.backend,
            "embedding": space.metadata.get("embedding", "X_pca"),
            "pca_source": pca_source,
            "depends_on_mixscape": bool(mixscape_ok),
            **{k: v for k, v in space.metadata.items() if k != "embedding"},
        }
        n_clust = (
            int(space.clusters["pert_cluster"].nunique())
            if "pert_cluster" in space.clusters.columns
            else None
        )
        report["perturbation_clusters"] = completed_status(
            method=space.method,
            backend=space.backend,
            n_clusters=n_clust,
            depends_on_mixscape=bool(mixscape_ok),
            pca_source=pca_source,
            space_meta={
                k: v
                for k, v in space.metadata.items()
                if isinstance(v, (str, int, float, bool, type(None)))
            },
        )
    except ImportError as exc:
        if is_pertpy_import_error(exc):
            msg = (
                "pertpy is required for this perturbation-space path. "
                'Install with: pip install -e ".[de]" (or pip install "pertpy[de]>=1.3"), '
                "or use --perturbation-space kmeans / lr_classifier (local backends)."
            )
            report["perturbation_clusters"] = skipped_status(
                "missing_pertpy", detail=msg, method=space_method
            )
            if config.skip_mixscape and config.skip_distance:
                print(f"NOTE: {msg}", flush=True)
            else:
                warnings.warn(msg, UserWarning, stacklevel=2)
                print(f"WARNING: perturbation clustering skipped — {msg}", flush=True)
        else:
            raise
    except Exception as exc:  # noqa: BLE001
        report["perturbation_clusters"] = skipped_status(
            "failed", detail=str(exc), method=space_method
        )

    print("Running guide-level QC / consistency", flush=True)
    qc_filter_path = tables / "qc_filter_by_guide.csv"
    qc_filter_df = pd.read_csv(qc_filter_path) if qc_filter_path.is_file() else None
    guide_df, consistency_df, warnings_df, guide_summary = run_guide_qc(
        rna,
        control=config.control,
        guide_merge=config.guide_merge,
        perturbation_type=config.perturbation_type,
        on_target_lfc_cutoff=config.on_target_lfc_cutoff,
        on_target_min_fail_guides=config.on_target_min_fail_guides,
        qc_filter_df=qc_filter_df,
    )
    write_guide_qc_tables(guide_df, consistency_df, warnings_df, tables)
    plot_guide_qc(guide_df, consistency_df, figures)
    report["guide_qc"] = guide_summary
    report["steps"].append("guide_qc")

    _mark_stage(report, "3a_perturbation_modeling")
    return rna, analysis_obj, edistances, mixscape_ok, guide_summary, consistency_df
