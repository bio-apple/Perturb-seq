"""Thin orchestrator for the tertiary Perturb-seq pipeline.

Stage map (matches the engineering structure diagram)::

    1. Input & validation          → io
    2. Preprocessing & QC          → guides, qc, preprocessing
    3. Parallel tracks (sequential in-process):
         A. Perturbation modeling  → perturbation (+ guide consistency)
         B. Statistical inference  → statistics
    4. Robustness & evidence       → robustness
    5. Reproducible report         → report (H5AD / CSV / JSON / HTML / provenance)

Resume / audit artifacts are managed by ``perturbseq.repro.RunTracker``.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import pandas as pd
from anndata import AnnData, read_h5ad

from perturbseq.cell_annotation import annotate_cells, annotation_summary, write_annotation_tables
from perturbseq.composition import append_composition_audit, write_composition_audit
from perturbseq.guides import DEFAULT_CONTROL_PATTERNS, annotate_guides, filter_singlets, run_guide_qc, write_guide_qc_tables
from perturbseq.io import load_dragen_sample, validate_sample_inputs, write_h5ad
from perturbseq.perturbation import cluster_perturbations, run_edistance, run_mixscape
from perturbseq.plots import (
    plot_cell_annotation,
    plot_edistance,
    plot_guide_composition,
    plot_guide_qc,
    plot_qc,
    plot_umap,
    plot_volcano,
)
from perturbseq.preprocessing import preprocess_rna
from perturbseq.qc import add_qc_metrics, filter_cells, sample_qc_summary
from perturbseq.report import finalize_outputs, init_report
from perturbseq.repro import (
    STAGE_ORDER,
    RunTracker,
    config_params_dict,
    sha256_paths,
)
from perturbseq.robustness import integrate_evidence
from perturbseq.statistics import check_experimental_design, run_deseq2_or_wilcoxon

# Keep in sync with repro.STAGE_ORDER.
PIPELINE_STAGES = STAGE_ORDER


@dataclass
class PipelineConfig:
    input_dir: Path
    output_dir: Path
    sample_id: str = "sample1"
    singlet_only: bool = True
    n_mads: float = 5.0
    min_cells: int = 3
    n_top_genes: int = 2000
    n_pcs: int = 30
    leiden_resolution: float = 0.5
    control: str = "NT"
    replicate_col: str | None = None
    skip_mixscape: bool = False
    skip_cell_annotation: bool = False
    skip_distance: bool = False
    skip_de: bool = False
    n_perms: int = 200
    min_cells_per_pert: int = 10
    mixscape_max_targets: int = 40
    force_mixscape: bool = False
    de_top_n: int = 10
    perturbation_type: str = "KO"
    random_state: int = 0
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS
    extra: dict = field(default_factory=dict)


def config_params(config: PipelineConfig) -> dict[str, Any]:
    raw = asdict(config)
    raw.pop("extra", None)
    return config_params_dict(raw)


def pipeline_config_from_mapping(data: dict[str, Any]) -> PipelineConfig:
    """Build PipelineConfig from a flat mapping (YAML / CLI)."""
    known = {f.name for f in fields(PipelineConfig)}
    kwargs: dict[str, Any] = {}
    extra: dict[str, Any] = {}
    for key, value in data.items():
        if key not in known:
            extra[key] = value
            continue
        if key in {"input_dir", "output_dir"} and value is not None:
            kwargs[key] = Path(value)
        elif key == "control_patterns" and value is not None:
            if isinstance(value, str):
                kwargs[key] = tuple(p.strip() for p in value.split(",") if p.strip())
            else:
                kwargs[key] = tuple(value)
        else:
            kwargs[key] = value
    if "input_dir" not in kwargs or "output_dir" not in kwargs:
        raise ValueError("Config must provide input_dir and output_dir")
    if "control_patterns" not in kwargs:
        kwargs["control_patterns"] = DEFAULT_CONTROL_PATTERNS
    if extra:
        kwargs["extra"] = {**(kwargs.get("extra") or {}), **extra}
    return PipelineConfig(**kwargs)


def _mixscape_subset(adata: AnnData, control: str) -> AnnData:
    if "mixscape_class_global" not in adata.obs:
        return adata
    keep = adata.obs["mixscape_class_global"].astype(str).isin([control, "KO"])
    return adata[keep].copy()



def _ko_label(gene: str, perturbation_type: str) -> str:
    return f"{gene} {perturbation_type}"


def _write_edistance_outputs(
    edistances: pd.DataFrame,
    etest: pd.DataFrame | None,
    tables_dir: Path,
    figures_dir: Path,
) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    edistances.to_csv(tables_dir / "edistance.csv")
    plot_edistance(edistances, figures_dir)
    if etest is not None:
        etest.to_csv(tables_dir / "etest.csv")


def _pca_from_layer(adata: AnnData, layer: str, n_pcs: int) -> None:
    import scanpy as sc

    tmp = adata.copy()
    tmp.X = tmp.layers[layer]
    n_comps = min(n_pcs, max(2, tmp.n_obs - 1), max(2, tmp.n_vars - 1))
    sc.pp.pca(tmp, n_comps=n_comps)
    adata.obsm["X_pca"] = tmp.obsm["X_pca"]


def _mark_stage(report: dict, stage: str) -> None:
    stages = list(report.get("stages") or [])
    if stage not in stages:
        stages.append(stage)
    report["stages"] = stages


def _stage_params(config: PipelineConfig, stage: str) -> dict[str, Any]:
    common = {"sample_id": config.sample_id, "random_state": config.random_state}
    mapping = {
        "1_input_validation": {**common, "input_dir": str(config.input_dir)},
        "2_preprocessing_qc": {
            **common,
            "singlet_only": config.singlet_only,
            "n_mads": config.n_mads,
            "min_cells": config.min_cells,
            "n_top_genes": config.n_top_genes,
            "n_pcs": config.n_pcs,
            "leiden_resolution": config.leiden_resolution,
            "skip_cell_annotation": config.skip_cell_annotation,
            "control_patterns": list(config.control_patterns),
        },
        "3a_perturbation_modeling": {
            **common,
            "skip_mixscape": config.skip_mixscape,
            "skip_distance": config.skip_distance,
            "control": config.control,
            "replicate_col": config.replicate_col,
            "perturbation_type": config.perturbation_type,
            "mixscape_max_targets": config.mixscape_max_targets,
            "force_mixscape": config.force_mixscape,
            "n_perms": config.n_perms,
            "min_cells_per_pert": config.min_cells_per_pert,
            "n_pcs": config.n_pcs,
        },
        "3b_statistical_inference": {
            **common,
            "skip_de": config.skip_de,
            "control": config.control,
            "replicate_col": config.replicate_col,
            "perturbation_type": config.perturbation_type,
            "de_top_n": config.de_top_n,
            "min_cells_per_pert": config.min_cells_per_pert,
        },
        "4_robustness": {**common},
        "5_report": {**common, "sample_id": config.sample_id},
    }
    return mapping[stage]


def _save_checkpoint(tracker: RunTracker, stage: str, adata: AnnData) -> Path:
    path = tracker.checkpoint_path(stage)
    write_h5ad(adata, path)
    return path


def _load_checkpoint(tracker: RunTracker, stage: str) -> AnnData:
    path = tracker.checkpoint_path(stage)
    if not path.is_file():
        raise FileNotFoundError(f"Missing checkpoint for stage {stage}: {path}")
    return read_h5ad(path)


def _write_stage_state(tracker: RunTracker, stage: str, state: dict[str, Any]) -> Path:
    path = tracker.stage_dir(stage) / "state.json"
    path.write_text(json.dumps(state, indent=2, default=str))
    return path


def _read_stage_state(tracker: RunTracker, stage: str) -> dict[str, Any]:
    path = tracker.stage_dir(stage) / "state.json"
    if not path.is_file():
        return {}
    return json.loads(path.read_text())


def stage_input_validation(config: PipelineConfig, report: dict) -> tuple[AnnData, dict[str, Path]]:
    """Stage 1: DRAGEN MEX + guide reference + assignments."""
    print(f"[1/5] Input & validation: {config.sample_id} from {config.input_dir}", flush=True)
    files = validate_sample_inputs(config.input_dir, config.sample_id)
    rna, _crispr, assignments, feature_ref = load_dragen_sample(config.input_dir, config.sample_id)
    print(f"Loaded {rna.n_obs} cells × {rna.n_vars} genes", flush=True)
    report["n_cells_loaded"] = int(rna.n_obs)
    report["n_genes_loaded"] = int(rna.n_vars)
    report["input_files"] = {k: str(v) for k, v in files.items()}
    _mark_stage(report, "1_input_validation")
    # Carry assignments/ref via temporary attrs for stage 2 (avoid widening return API).
    rna.uns["_perturbseq_assignments"] = assignments
    rna.uns["_perturbseq_feature_ref"] = feature_ref
    return rna, files


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
    rna = annotate_guides(rna, assignments, feature_ref, config.control_patterns)
    rna = add_qc_metrics(rna)
    plot_qc(rna, figures)
    plot_guide_composition(rna, figures)
    report["guide_counts"] = rna.obs["num_features"].value_counts().sort_index().to_dict()
    report["gene_target_counts"] = rna.obs["gene_target"].value_counts().to_dict()
    report["sample_qc"] = sample_qc_summary(rna)
    append_composition_audit(composition_rows, rna, "loaded")

    rna, qc_log = filter_cells(rna, n_mads=config.n_mads, min_cells=config.min_cells)
    report["qc"] = qc_log
    report.setdefault("composition_notes", []).append(
        f"QC removed {qc_log.get('n_cells_removed_qc', 0)} cells; gene set also filtered (min_cells={config.min_cells})."
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
    if rna.n_obs < 20:
        raise ValueError(f"Only {rna.n_obs} cells remain after QC/singlet filters")

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
        report["cell_annotation"] = annotation_summary(rna)
        report["steps"].append("cell_annotation")
        append_composition_audit(composition_rows, rna, "after_cell_annotation")

    _mark_stage(report, "2_preprocessing_qc")
    return rna


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
    n_targets = int(
        rna.obs.loc[~rna.obs["gene_target"].isin([config.control, "unassigned"]), "gene_target"].nunique()
    )
    report["n_gene_targets"] = n_targets
    mixscape_ok = False
    skip_mixscape = config.skip_mixscape
    if skip_mixscape:
        report["mixscape_status"] = "skipped_by_config"
        report["mixscape_note"] = (
            "Mixscape was skipped (--skip-mixscape). Downstream E-distance/DE use gene_target "
            "on all post-QC cells; no KO/NP classification filter."
        )
    elif n_targets > config.mixscape_max_targets and not config.force_mixscape:
        skip_mixscape = True
        report["mixscape_status"] = "skipped_too_many_targets"
        report["mixscape_skipped"] = (
            f"{n_targets} gene targets > mixscape_max_targets={config.mixscape_max_targets}; "
            "use --force-mixscape to override. Genome-scale screens are quantified with E-distance."
        )
        report["mixscape_note"] = report["mixscape_skipped"]
    else:
        report["mixscape_status"] = "attempted"

    edistances = pd.DataFrame()
    if not config.skip_distance:
        try:
            print("Computing E-distance (pre-Mixscape / log-norm PCA)", flush=True)
            edist_pre, etest_pre = run_edistance(
                rna,
                groupby="gene_target",
                contrast=config.control,
                n_perms=config.n_perms,
                min_cells=config.min_cells_per_pert,
                random_state=config.random_state,
                pca_source="log1p_hvg",
            )
            if not skip_mixscape:
                _write_edistance_outputs(
                    edist_pre, etest_pre, tables / "pre_mixscape", figures / "pre_mixscape"
                )
                report["pre_mixscape"] = {
                    "n_cells": int(rna.n_obs),
                    "groupby": "gene_target",
                    "embedding": "X_pca",
                    "pca_source": "log1p_hvg",
                    "depends_on_mixscape": False,
                    "tables": str(tables / "pre_mixscape"),
                }
            else:
                _write_edistance_outputs(edist_pre, etest_pre, tables, figures)
                edistances = edist_pre
                report["edistance_top"] = edistances.head(10).to_dict()
                report.setdefault("matrix_provenance", {})["edistance"] = {
                    "embedding": "X_pca",
                    "pca_source": "log1p_hvg",
                    "cells": "all_post_qc_singlet",
                    "depends_on_mixscape": False,
                }
            report["steps"].append("edistance_pre_or_only")
        except Exception as exc:  # noqa: BLE001
            report["distance_error_pre"] = str(exc)

    if not skip_mixscape:
        try:
            run_mixscape(
                rna,
                control=config.control,
                split_by=config.replicate_col,
                perturbation_type=config.perturbation_type,
            )
            mixscape_ok = True
            report["mixscape_status"] = "ok"
            report["mixscape_global"] = rna.obs["mixscape_class_global"].value_counts().to_dict()
            report["mixscape_note"] = (
                "Mixscape succeeded. Primary E-distance/DE use KO+control cells and depend on "
                "mixscape_class; pre_mixscape/ retains the unfiltered gene_target path."
            )
            plot_umap(rna, figures / "mixscape", color=["mixscape_class_global", "perturbation"])
            report["steps"].append("mixscape")
            append_composition_audit(composition_rows, rna, "after_mixscape_classify")
        except Exception as exc:  # noqa: BLE001
            report["mixscape_status"] = "failed"
            report["mixscape_error"] = str(exc)
            report["mixscape_note"] = (
                f"Mixscape failed ({exc}). Downstream uses gene_target on all post-QC cells; "
                "no KO/NP filter."
            )
            report["steps"].append("mixscape_failed")

    analysis_obj = _mixscape_subset(rna, config.control) if mixscape_ok else rna
    pca_source = "log1p_hvg"
    if mixscape_ok:
        n_before_ko = int(rna.n_obs)
        n_after_ko = int(analysis_obj.n_obs)
        report["n_after_mixscape_ko_filter"] = n_after_ko
        report["n_removed_mixscape_non_ko"] = n_before_ko - n_after_ko
        report.setdefault("composition_notes", []).append(
            f"Mixscape KO+{config.control} filter kept {n_after_ko}/{n_before_ko} cells "
            f"(removed {n_before_ko - n_after_ko} NP/other)."
        )
        append_composition_audit(composition_rows, analysis_obj, "after_mixscape_ko_filter")
        if "X_pert" in analysis_obj.layers:
            _pca_from_layer(analysis_obj, "X_pert", config.n_pcs)
            pca_source = "X_pert"
            report.setdefault("matrix_provenance", {}).update(
                {
                    "X_pca_post_mixscape": "PCA recomputed on layers['X_pert']",
                    "X_pert": "Mixscape perturbation signature",
                }
            )

    if mixscape_ok and not config.skip_distance:
        try:
            print("Computing E-distance (post-Mixscape / X_pert PCA)", flush=True)
            edistances, etest = run_edistance(
                analysis_obj,
                groupby="gene_target",
                contrast=config.control,
                n_perms=config.n_perms,
                min_cells=config.min_cells_per_pert,
                random_state=config.random_state,
                pca_source=pca_source,
            )
            _write_edistance_outputs(edistances, etest, tables / "post_mixscape", figures / "post_mixscape")
            _write_edistance_outputs(edistances, etest, tables, figures)
            report["edistance_top"] = edistances.head(10).to_dict()
            report["post_mixscape"] = {
                "n_cells": int(analysis_obj.n_obs),
                "groupby": "gene_target",
                "embedding": "X_pca",
                "pca_source": pca_source,
                "depends_on_mixscape": True,
                "tables": str(tables / "post_mixscape"),
            }
            report.setdefault("matrix_provenance", {})["edistance"] = {
                "embedding": "X_pca",
                "pca_source": pca_source,
                "cells": f"mixscape_class_global in {{{config.control}, KO}}",
                "depends_on_mixscape": True,
                "primary_tables": str(tables / "edistance.csv"),
            }
            report["steps"].append("edistance_post_mixscape")
        except Exception as exc:  # noqa: BLE001
            report["distance_error"] = str(exc)

    try:
        pert_clusters = cluster_perturbations(analysis_obj, groupby="gene_target")
        pert_clusters.to_csv(tables / "perturbation_clusters.csv")
        report.setdefault("matrix_provenance", {})["perturbation_clusters"] = {
            "embedding": "X_pca",
            "pca_source": pca_source,
            "depends_on_mixscape": bool(mixscape_ok),
        }
    except Exception as exc:  # noqa: BLE001
        report["cluster_error"] = str(exc)

    print("Running guide-level QC / consistency", flush=True)
    guide_df, consistency_df, warnings_df, guide_summary = run_guide_qc(rna, control=config.control)
    write_guide_qc_tables(guide_df, consistency_df, warnings_df, tables)
    plot_guide_qc(guide_df, consistency_df, figures)
    report["guide_qc"] = guide_summary
    report["steps"].append("guide_qc")

    _mark_stage(report, "3a_perturbation_modeling")
    return rna, analysis_obj, edistances, mixscape_ok, guide_summary, consistency_df


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
    design = check_experimental_design(rna, replicate_col=config.replicate_col, groupby="gene_target")
    report["experimental_design"] = design

    de_tables: list[pd.DataFrame] = []
    if not config.skip_de:
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
        for group in groups:
            n_group = int((counts.obs["de_group"] == group).sum())
            if n_group < config.min_cells_per_pert:
                continue
            try:
                table = run_deseq2_or_wilcoxon(
                    counts,
                    group=group,
                    reference=reference,
                    replicate_col=config.replicate_col,
                    groupby="de_group",
                )
            except Exception as exc:  # noqa: BLE001
                report.setdefault("de_errors", {})[str(group)] = str(exc)
                continue
            if table.empty:
                continue
            safe = str(group).replace(" ", "_")
            table.to_csv(tables / f"de_{safe}.csv", index=False)
            plot_volcano(table, figures, safe)
            de_tables.append(table.head(50))
        if de_tables:
            pd.concat(de_tables, ignore_index=True).to_csv(tables / "de_top50_concat.csv", index=False)
        report["n_de_contrasts"] = len(de_tables)
        report["steps"].append("de")
        if design.get("note"):
            report["de_note"] = design["note"]
    else:
        report["n_de_contrasts"] = 0

    _mark_stage(report, "3b_statistical_inference")
    return int(report.get("n_de_contrasts") or 0), design


def stage_robustness(
    report: dict,
    *,
    guide_summary: dict,
    consistency_df: pd.DataFrame,
    edistances: pd.DataFrame,
    mixscape_ok: bool,
    n_de_contrasts: int,
    design: dict,
) -> None:
    """Stage 4: sensitivity / effect consistency / confidence flags."""
    print("[4/5] Robustness & evidence integration", flush=True)
    integrate_evidence(
        report,
        guide_summary=guide_summary,
        consistency_df=consistency_df,
        edistances=edistances,
        mixscape_ok=mixscape_ok,
        n_de_contrasts=n_de_contrasts,
        design=design,
    )
    _mark_stage(report, "4_robustness")


def stage_report(
    rna: AnnData,
    config: PipelineConfig,
    report: dict,
    *,
    input_files: dict[str, Path],
) -> dict:
    """Stage 5: H5AD, JSON/HTML reports, provenance manifest."""
    print("[5/5] Reproducible report", flush=True)
    report = finalize_outputs(
        rna,
        report,
        config.output_dir,
        config.sample_id,
        inputs={k: str(v) for k, v in input_files.items()},
    )
    _mark_stage(report, "5_report")
    return report


def run_pipeline(
    config: PipelineConfig,
    *,
    resume: bool = False,
    config_path: str | Path | None = None,
) -> dict:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    figures = config.output_dir / "figures"
    tables = config.output_dir / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)

    tracker = RunTracker(config.output_dir, resume=resume, seed=config.random_state)
    input_files = validate_sample_inputs(config.input_dir, config.sample_id)
    input_checksums = sha256_paths(input_files.values())
    tracker.set_run_metadata(
        inputs=input_checksums,
        params=config_params(config),
        config_path=str(config_path) if config_path else None,
    )

    report = init_report(config.sample_id, config)
    report["random_state"] = config.random_state
    report["versions"] = tracker.manifest["versions"]
    report["environment"] = tracker.manifest["environment"]
    report["step_kinds"] = {
        "qc": "filtering",
        "singlet": "filtering",
        "preprocess_hvg_pca_leiden": "descriptive",
        "umap": "descriptive_visualization",
        "cell_annotation": "descriptive",
        "mixscape": "filtering_classification",
        "edistance": "inferential",
        "de": "inferential_or_exploratory",
        "guide_qc": "descriptive",
    }
    report.setdefault("matrix_provenance", {})
    report.setdefault("composition_notes", [])
    report["analysis_dependencies_doc"] = "docs/ANALYSIS_DEPENDENCIES.md"
    composition_rows: list[pd.DataFrame] = []

    rna: AnnData | None = None
    edistances = pd.DataFrame()
    mixscape_ok = False
    guide_summary: dict = {}
    consistency_df = pd.DataFrame()
    n_de_contrasts = 0
    design: dict = {}

    try:
        # ---- 1_input_validation ----
        stage = "1_input_validation"
        params = _stage_params(config, stage)
        if tracker.should_skip(stage, params=params, input_checksums=input_checksums):
            rna = _load_checkpoint(tracker, stage)
            report["n_cells_loaded"] = int(rna.n_obs)
            report["n_genes_loaded"] = int(rna.n_vars)
            report["input_files"] = {k: str(v) for k, v in input_files.items()}
            _mark_stage(report, stage)
            tracker.logger.info("Resumed checkpoint for %s", stage)
        else:
            payload = tracker.begin_stage(stage, inputs=input_checksums, params=params)
            try:
                rna, input_files = stage_input_validation(config, report)
                ckpt = _save_checkpoint(tracker, stage, rna)
                tracker.finish_stage(
                    stage,
                    payload,
                    status="success",
                    outputs={"adata": str(ckpt), **{k: str(v) for k, v in input_files.items()}},
                )
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        assert rna is not None

        # ---- 2_preprocessing_qc ----
        stage = "2_preprocessing_qc"
        params = _stage_params(config, stage)
        stage_inputs = sha256_paths([tracker.checkpoint_path("1_input_validation")])
        if tracker.should_skip(stage, params=params, input_checksums=stage_inputs):
            rna = _load_checkpoint(tracker, stage)
            st = _read_stage_state(tracker, stage)
            report.update({k: st[k] for k in ("qc", "n_after_singlet", "n_removed_nonsinglet", "cell_annotation") if k in st})
            _mark_stage(report, stage)
        else:
            payload = tracker.begin_stage(stage, inputs=stage_inputs, params=params)
            try:
                if "_perturbseq_assignments" not in rna.uns:
                    rna = _load_checkpoint(tracker, "1_input_validation")
                rna = stage_preprocessing_qc(rna, config, report, figures, composition_rows)
                ckpt = _save_checkpoint(tracker, stage, rna)
                state = {
                    "qc": report.get("qc"),
                    "n_after_singlet": report.get("n_after_singlet"),
                    "n_removed_nonsinglet": report.get("n_removed_nonsinglet"),
                    "cell_annotation": report.get("cell_annotation"),
                }
                state_path = _write_stage_state(tracker, stage, state)
                tracker.finish_stage(
                    stage,
                    payload,
                    status="success",
                    outputs={"adata": str(ckpt), "state": str(state_path)},
                )
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        # ---- 3a_perturbation_modeling ----
        stage = "3a_perturbation_modeling"
        params = _stage_params(config, stage)
        stage_inputs = sha256_paths([tracker.checkpoint_path("2_preprocessing_qc")])
        if tracker.should_skip(stage, params=params, input_checksums=stage_inputs):
            rna = _load_checkpoint(tracker, stage)
            st = _read_stage_state(tracker, stage)
            mixscape_ok = bool(st.get("mixscape_ok"))
            guide_summary = st.get("guide_summary") or {}
            edist_path = tables / "edistance.csv"
            if edist_path.is_file():
                edistances = pd.read_csv(edist_path, index_col=0)
            cons_path = tables / "gene_guide_consistency.csv"
            if cons_path.is_file():
                consistency_df = pd.read_csv(cons_path)
            report["guide_qc"] = guide_summary
            report["n_gene_targets"] = st.get("n_gene_targets")
            _mark_stage(report, stage)
        else:
            payload = tracker.begin_stage(stage, inputs=stage_inputs, params=params)
            try:
                rna = _load_checkpoint(tracker, "2_preprocessing_qc")
                rna, _analysis_obj, edistances, mixscape_ok, guide_summary, consistency_df = (
                    stage_perturbation_modeling(rna, config, report, figures, tables, composition_rows)
                )
                ckpt = _save_checkpoint(tracker, stage, rna)
                state = {
                    "mixscape_ok": mixscape_ok,
                    "guide_summary": guide_summary,
                    "n_gene_targets": report.get("n_gene_targets"),
                }
                state_path = _write_stage_state(tracker, stage, state)
                outputs = {"adata": str(ckpt), "state": str(state_path)}
                for name in ("edistance.csv", "etest.csv", "gene_guide_consistency.csv", "guide_qc.csv"):
                    p = tables / name
                    if p.is_file():
                        outputs[name] = str(p)
                tracker.finish_stage(stage, payload, status="success", outputs=outputs)
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        # ---- 3b_statistical_inference ----
        stage = "3b_statistical_inference"
        params = _stage_params(config, stage)
        stage_inputs = sha256_paths([tracker.checkpoint_path("3a_perturbation_modeling")])
        edist_csv = tables / "edistance.csv"
        if edist_csv.is_file():
            stage_inputs.update(sha256_paths([edist_csv]))
        if tracker.should_skip(stage, params=params, input_checksums=stage_inputs):
            rna = _load_checkpoint(tracker, "3a_perturbation_modeling")
            st = _read_stage_state(tracker, stage)
            n_de_contrasts = int(st.get("n_de_contrasts") or 0)
            design = st.get("design") or {}
            report["n_de_contrasts"] = n_de_contrasts
            report["experimental_design"] = design
            if "de_note" in st:
                report["de_note"] = st["de_note"]
            _mark_stage(report, stage)
        else:
            payload = tracker.begin_stage(stage, inputs=stage_inputs, params=params)
            try:
                rna = _load_checkpoint(tracker, "3a_perturbation_modeling")
                n_de_contrasts, design = stage_statistical_inference(
                    rna,
                    config,
                    report,
                    figures,
                    tables,
                    mixscape_ok=mixscape_ok,
                    edistances=edistances,
                )
                # AnnData unchanged; reuse 3a checkpoint path as this stage's artifact pointer.
                ckpt = tracker.checkpoint_path("3a_perturbation_modeling")
                state = {
                    "n_de_contrasts": n_de_contrasts,
                    "design": design,
                    "de_note": report.get("de_note"),
                    "de_groups": report.get("de_groups"),
                }
                state_path = _write_stage_state(tracker, stage, state)
                outputs = {"adata": str(ckpt), "state": str(state_path)}
                concat = tables / "de_top50_concat.csv"
                if concat.is_file():
                    outputs["de_top50_concat"] = str(concat)
                tracker.finish_stage(stage, payload, status="success", outputs=outputs)
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        # ---- 4_robustness ----
        stage = "4_robustness"
        params = _stage_params(config, stage)
        stage_inputs = sha256_paths(
            [
                tracker.stage_dir("3a_perturbation_modeling") / "state.json",
                tracker.stage_dir("3b_statistical_inference") / "state.json",
            ]
        )
        if tracker.should_skip(stage, params=params, input_checksums=stage_inputs):
            st = _read_stage_state(tracker, stage)
            if "confidence_flags" in st:
                report["confidence_flags"] = st["confidence_flags"]
            _mark_stage(report, stage)
        else:
            payload = tracker.begin_stage(stage, inputs=stage_inputs, params=params)
            try:
                stage_robustness(
                    report,
                    guide_summary=guide_summary,
                    consistency_df=consistency_df,
                    edistances=edistances,
                    mixscape_ok=mixscape_ok,
                    n_de_contrasts=n_de_contrasts,
                    design=design,
                )
                state = {"confidence_flags": report.get("confidence_flags")}
                state_path = _write_stage_state(tracker, stage, state)
                tracker.finish_stage(
                    stage,
                    payload,
                    status="success",
                    outputs={"state": str(state_path)},
                )
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        audit_path = tables / "composition_audit.csv"
        if composition_rows:
            composition_table = write_composition_audit(composition_rows, audit_path)
        elif audit_path.is_file():
            composition_table = pd.read_csv(audit_path)
        else:
            composition_table = write_composition_audit([], audit_path)
        report["composition_audit"] = str(audit_path)
        report["composition_stages"] = (
            composition_table["stage"].drop_duplicates().tolist() if not composition_table.empty else []
        )
        report["steps"].append("composition_audit")

        # ---- 5_report ----
        stage = "5_report"
        params = _stage_params(config, stage)
        stage_inputs = sha256_paths([tracker.checkpoint_path("3a_perturbation_modeling")])
        out_h5ad = config.output_dir / f"{config.sample_id}.tertiary.h5ad"
        if tracker.should_skip(stage, params=params, input_checksums=stage_inputs) and out_h5ad.is_file():
            _mark_stage(report, stage)
            report["output_h5ad"] = str(out_h5ad)
        else:
            payload = tracker.begin_stage(stage, inputs=stage_inputs, params=params)
            try:
                if rna is None or "gene_target" not in rna.obs:
                    rna = _load_checkpoint(tracker, "3a_perturbation_modeling")
                report = stage_report(rna, config, report, input_files=input_files)
                outputs = {
                    "tertiary_h5ad": str(out_h5ad),
                    "report_json": str(config.output_dir / "report.json"),
                }
                html_path = config.output_dir / "report.html"
                if html_path.is_file():
                    outputs["report_html"] = str(html_path)
                tracker.finish_stage(stage, payload, status="success", outputs=outputs)
            except Exception as exc:  # noqa: BLE001
                tracker.finish_stage(stage, payload, status="failed", error=str(exc))
                raise

        report_path = config.output_dir / "report.json"
        if report_path.is_file():
            try:
                report = json.loads(report_path.read_text())
            except json.JSONDecodeError:
                pass
        report["manifest"] = str(tracker.manifest_path)
        report["run_log"] = str(tracker.log_path)
        report_path.write_text(json.dumps(report, indent=2, default=str))
        tracker.manifest["output_h5ad"] = str(out_h5ad)
        tracker.finalize(status="success", report_path=str(report_path))
        return report
    except Exception:
        tracker.finalize(status="failed")
        raise


__all__ = [
    "PIPELINE_STAGES",
    "PipelineConfig",
    "config_params",
    "pipeline_config_from_mapping",
    "run_pipeline",
]
