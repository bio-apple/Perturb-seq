"""Thin orchestrator for the tertiary Perturb-seq pipeline.

Stage map (matches the engineering structure diagram)::

    1. Input & validation          → stages.input_validation → io
    2. Preprocessing & QC          → stages.preprocessing_qc → guides, qc, preprocessing
    3. Parallel tracks (sequential in-process):
         A. Perturbation modeling  → stages.perturbation_modeling → perturbation
         B. Statistical inference  → stages.statistical_inference → statistics
    4. Robustness & evidence       → stages.robustness → robustness
    5. Reproducible report         → stages.report → report

Owns: ``PipelineConfig``, dry-run planning, checkpoint/resume orchestration.
Does NOT own: stage science bodies (``stages/``) or domain algorithms.

Resume / audit artifacts are managed by ``perturbseq.repro.RunTracker``.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import pandas as pd
from anndata import AnnData, read_h5ad

from perturbseq.annotation_policy import SAMPLE_TYPES, resolve_cell_annotation_decision
from perturbseq.guide_reassignment import GUIDE_REASSIGN_MODES
from perturbseq.guides import DEFAULT_CONTROL_PATTERNS
from perturbseq.io import validate_sample_inputs, write_h5ad
from perturbseq.perturbation import (
    DEFAULT_ETEST_POWER_MIN_CELLS,
    DEFAULT_N_BOOTSTRAP,
    DEFAULT_SECONDARY_DISTANCE_METRICS,
    PERTURBATION_SPACE_METHODS,
)
from perturbseq.report import init_report
from perturbseq.repro import (
    STAGE_ORDER,
    RunTracker,
    config_params_dict,
    sha256_paths,
)
from perturbseq.composition import write_composition_audit
from perturbseq.stages import (
    stage_input_validation,
    stage_perturbation_modeling,
    stage_preprocessing_qc,
    stage_report,
    stage_robustness,
    stage_statistical_inference,
)
from perturbseq.stages._helpers import (
    _effective_mixscape_mode,
    _embedding_provenance,
    _mark_stage,
)

# Keep in sync with repro.STAGE_ORDER.
PIPELINE_STAGES = STAGE_ORDER


@dataclass
class PipelineConfig:
    input_dir: Path
    output_dir: Path
    sample_id: str = "sample1"
    singlet_only: bool = True
    n_mads: float = 5.0
    # False = global MAD (default, backward-compatible). True = fit MAD on NT/control
    # only then apply to all cells — preferred for Perturb-seq (avoids dropping strong phenotypes).
    perturbation_aware_qc: bool = False
    min_cells: int = 3
    n_top_genes: int = 2000
    n_pcs: int = 30
    leiden_resolution: float = 0.5
    control: str = "NT"
    replicate_col: str | None = None
    skip_mixscape: bool = False
    # None = decide from sample_type; True/False = explicit --skip/--run-cell-annotation.
    skip_cell_annotation: bool | None = None
    # cell_line (default) → skip annotation; primary/mixed → run; unknown → skip + note.
    sample_type: str = "cell_line"
    skip_distance: bool = False
    skip_de: bool = False
    n_perms: int = 200
    min_cells_per_pert: int = 10
    # Peidli et al.: prefer ≥50–100 cells/pert for trusted E-test; ~200 more stable.
    etest_power_min_cells: int = DEFAULT_ETEST_POWER_MIN_CELLS
    secondary_distance_metrics: tuple[str, ...] = DEFAULT_SECONDARY_DISTANCE_METRICS
    # Cell bootstrap replicates for E-distance CI (0 = skip). Modest default for speed.
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP
    mixscape_max_targets: int = 40
    force_mixscape: bool = False
    # auto | skip | force | subset — CLI --skip/--force still work; subset = NT+selected targets.
    mixscape_mode: str = "auto"
    mixscape_targets: tuple[str, ...] = ()
    mixscape_top_n: int | None = None
    de_top_n: int = 10
    # Cap biologist report plots (guide consistency + target validation) for speed.
    report_plot_top_n: int = 15
    n_jobs: int = 1
    perturbation_type: str = "KO"
    # Perturbation-space construction: pca_silhouette (default) | kmeans | lr_classifier
    perturbation_space: str = "pca_silhouette"
    random_state: int = 0
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS
    # none = flag inconsistency only (default); equal|umi|confidence|umi_confidence
    # add weighted gene-level summaries without dropping per-guide rows.
    guide_merge: str = "none"
    # On-target efficiency proxy (guide QC): |log2FC| vs control in expected direction.
    on_target_lfc_cutoff: float = 0.25
    # Gene flagged potential_low_efficiency when ≥ this many adequate guides fail on-target.
    on_target_min_fail_guides: int = 2
    # Optional tertiary guide re-call vs DRAGEN: off | compare | apply_max | apply_gmm
    guide_reassign: str = "off"
    guide_reassign_min_umi: float = 1.0
    # Include cell-level covariates in PyDESeq2 design when available (Wilcoxon: documented ignore).
    # True = default candidates (phase, pct_counts_mt, log_n_counts); False/[] = off; list = explicit.
    de_covariates: bool | tuple[str, ...] = True
    # Prefer sample_id(+perturbation) sum pseudobulk + PyDESeq2 even without bio reps.
    # Wilcoxon remains a loud exploratory fallback when <2 pseudobulk units exist.
    de_prefer_pseudobulk: bool = True
    # Optional technical pseudo-replicates (N≥2) for single-sample exploratory PyDESeq2.
    # None/0 = off. Not a substitute for biological replicates.
    n_pseudo_replicates: int | None = None
    extra: dict = field(default_factory=dict)


def config_params(config: PipelineConfig) -> dict[str, Any]:
    raw = asdict(config)
    raw.pop("extra", None)
    return config_params_dict(raw)


def dry_run_plan(config: PipelineConfig) -> dict[str, Any]:
    """Validate inputs and return planned stages + resolved params (no analysis)."""
    files = validate_sample_inputs(config.input_dir, config.sample_id)
    stages: list[dict[str, Any]] = []
    ann_decision = resolve_cell_annotation_decision(
        config.sample_type,
        skip_cell_annotation=config.skip_cell_annotation,
    )
    for name in PIPELINE_STAGES:
        notes: list[str] = []
        if name == "2_preprocessing_qc" and ann_decision["skip"]:
            notes.append(f"skip cell_annotation ({ann_decision['reason']})")
        if name == "3a_perturbation_modeling":
            if config.skip_mixscape:
                notes.append("skip mixscape")
            if config.skip_distance:
                notes.append("skip edistance")
            notes.append(f"perturbation_space={config.perturbation_space}")
        if name == "3b_statistical_inference" and config.skip_de:
            notes.append("skip de")
        stages.append({"stage": name, "action": "run", "notes": notes})
    return {
        "sample_id": config.sample_id,
        "input_dir": str(config.input_dir),
        "output_dir": str(config.output_dir),
        "input_files": {k: str(v) for k, v in files.items()},
        "stages": stages,
        "params": config_params(config),
    }


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
        elif key in {"control_patterns", "secondary_distance_metrics", "mixscape_targets"} and value is not None:
            if isinstance(value, str):
                kwargs[key] = tuple(p.strip() for p in value.split(",") if p.strip())
            else:
                kwargs[key] = tuple(str(p).strip() for p in value if str(p).strip())
        elif key == "de_covariates" and value is not None:
            if isinstance(value, bool):
                kwargs[key] = value
            elif isinstance(value, str):
                low = value.strip().lower()
                if low in {"true", "yes", "1"}:
                    kwargs[key] = True
                elif low in {"false", "no", "0", "none", ""}:
                    kwargs[key] = False
                else:
                    kwargs[key] = tuple(p.strip() for p in value.split(",") if p.strip())
            else:
                kwargs[key] = tuple(str(p).strip() for p in value if str(p).strip())
        elif key == "de_prefer_pseudobulk" and value is not None:
            if isinstance(value, bool):
                kwargs[key] = value
            else:
                kwargs[key] = str(value).strip().lower() in {"true", "yes", "1"}
        elif key == "n_pseudo_replicates":
            if value is None or value == "" or value is False:
                kwargs[key] = None
            else:
                n = int(value)
                kwargs[key] = n if n >= 2 else None
        elif key == "mixscape_mode" and value is not None:
            kwargs[key] = str(value).strip().lower()
        elif key == "sample_type" and value is not None:
            kwargs[key] = str(value).strip().lower()
        elif key == "skip_cell_annotation":
            if value is None or value == "":
                kwargs[key] = None
            elif isinstance(value, bool):
                kwargs[key] = value
            else:
                low = str(value).strip().lower()
                if low in {"true", "yes", "1"}:
                    kwargs[key] = True
                elif low in {"false", "no", "0"}:
                    kwargs[key] = False
                elif low in {"null", "none", "auto", ""}:
                    kwargs[key] = None
                else:
                    raise ValueError(
                        f"skip_cell_annotation must be true/false/null, got {value!r}"
                    )
        elif key == "perturbation_space" and value is not None:
            kwargs[key] = str(value).strip().lower()
        elif key == "guide_reassign" and value is not None:
            # YAML 1.1 treats bare `off`/`on` as bool; coerce to mode strings.
            if value is False:
                kwargs[key] = "off"
            elif value is True:
                raise ValueError(
                    "guide_reassign got boolean true; use off|compare|apply_max|apply_gmm "
                    '(quote "off" in YAML to avoid YAML 1.1 bool parsing)'
                )
            else:
                kwargs[key] = str(value).strip().lower()
        else:
            kwargs[key] = value
    if "input_dir" not in kwargs or "output_dir" not in kwargs:
        raise ValueError(
            "Config must provide input_dir and output_dir.\n"
            "Suggested commands:\n"
            "  python -m perturbseq config --generate -o configs/pipeline.yaml\n"
            "  python -m perturbseq run --input-dir <dir> --output-dir <out> --dry-run"
        )
    if "control_patterns" not in kwargs:
        kwargs["control_patterns"] = DEFAULT_CONTROL_PATTERNS
    if "secondary_distance_metrics" not in kwargs:
        kwargs["secondary_distance_metrics"] = DEFAULT_SECONDARY_DISTANCE_METRICS
    if "mixscape_mode" in kwargs and kwargs["mixscape_mode"] is not None:
        mode = str(kwargs["mixscape_mode"]).lower()
        allowed = {"auto", "skip", "force", "subset"}
        if mode not in allowed:
            raise ValueError(f"mixscape_mode must be one of {sorted(allowed)}, got {kwargs['mixscape_mode']!r}")
        kwargs["mixscape_mode"] = mode
    if "sample_type" in kwargs and kwargs["sample_type"] is not None:
        st = str(kwargs["sample_type"]).lower()
        if st not in SAMPLE_TYPES:
            raise ValueError(f"sample_type must be one of {SAMPLE_TYPES}, got {kwargs['sample_type']!r}")
        kwargs["sample_type"] = st
    if "perturbation_space" in kwargs and kwargs["perturbation_space"] is not None:
        space = str(kwargs["perturbation_space"]).lower()
        if space not in PERTURBATION_SPACE_METHODS:
            raise ValueError(
                f"perturbation_space must be one of {PERTURBATION_SPACE_METHODS}, got {kwargs['perturbation_space']!r}"
            )
        kwargs["perturbation_space"] = space
    if "guide_merge" in kwargs and kwargs["guide_merge"] is not None:
        from perturbseq.guide_qc import GUIDE_MERGE_MODES

        mode = str(kwargs["guide_merge"]).lower()
        if mode not in GUIDE_MERGE_MODES:
            raise ValueError(f"guide_merge must be one of {GUIDE_MERGE_MODES}, got {kwargs['guide_merge']!r}")
        kwargs["guide_merge"] = mode
    if "guide_reassign" in kwargs and kwargs["guide_reassign"] is not None:
        mode = str(kwargs["guide_reassign"]).lower()
        if mode not in GUIDE_REASSIGN_MODES:
            raise ValueError(
                f"guide_reassign must be one of {GUIDE_REASSIGN_MODES}, got {kwargs['guide_reassign']!r}"
            )
        kwargs["guide_reassign"] = mode
    if extra:
        kwargs["extra"] = {**(kwargs.get("extra") or {}), **extra}
    return PipelineConfig(**kwargs)


def _stage_params(config: PipelineConfig, stage: str) -> dict[str, Any]:
    common = {"sample_id": config.sample_id, "random_state": config.random_state}
    mapping = {
        "1_input_validation": {**common, "input_dir": str(config.input_dir)},
        "2_preprocessing_qc": {
            **common,
            "singlet_only": config.singlet_only,
            "n_mads": config.n_mads,
            "perturbation_aware_qc": config.perturbation_aware_qc,
            "min_cells": config.min_cells,
            "n_top_genes": config.n_top_genes,
            "n_pcs": config.n_pcs,
            "leiden_resolution": config.leiden_resolution,
            "sample_type": config.sample_type,
            "skip_cell_annotation": config.skip_cell_annotation,
            "control_patterns": list(config.control_patterns),
            "guide_reassign": config.guide_reassign,
            "guide_reassign_min_umi": config.guide_reassign_min_umi,
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
            "mixscape_mode": config.mixscape_mode,
            "mixscape_targets": list(config.mixscape_targets),
            "mixscape_top_n": config.mixscape_top_n,
            "n_perms": config.n_perms,
            "min_cells_per_pert": config.min_cells_per_pert,
            "etest_power_min_cells": config.etest_power_min_cells,
            "secondary_distance_metrics": list(config.secondary_distance_metrics),
            "n_bootstrap": config.n_bootstrap,
            "n_pcs": config.n_pcs,
            "n_jobs": config.n_jobs,
            "guide_merge": config.guide_merge,
            "on_target_lfc_cutoff": config.on_target_lfc_cutoff,
            "on_target_min_fail_guides": config.on_target_min_fail_guides,
            "perturbation_space": config.perturbation_space,
        },
        "3b_statistical_inference": {
            **common,
            "skip_de": config.skip_de,
            "control": config.control,
            "replicate_col": config.replicate_col,
            "perturbation_type": config.perturbation_type,
            "de_top_n": config.de_top_n,
            "de_covariates": config.de_covariates
            if isinstance(config.de_covariates, bool)
            else list(config.de_covariates),
            "de_prefer_pseudobulk": config.de_prefer_pseudobulk,
            "n_pseudo_replicates": config.n_pseudo_replicates,
            "min_cells_per_pert": config.min_cells_per_pert,
            "n_jobs": config.n_jobs,
        },
        "4_robustness": {**common},
        "5_report": {
            **common,
            "sample_id": config.sample_id,
            "report_plot_top_n": config.report_plot_top_n,
        },
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
            report.update(
                {
                    k: st[k]
                    for k in ("qc", "n_after_singlet", "n_removed_nonsinglet", "cell_annotation")
                    if k in st
                }
            )
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
            for key in ("mixscape", "edistance", "perturbation_clusters", "mixscape_global", "edistance_top"):
                if key in st:
                    report[key] = st[key]
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
                    "mixscape": report.get("mixscape"),
                    "edistance": report.get("edistance"),
                    "perturbation_clusters": report.get("perturbation_clusters"),
                    "mixscape_global": report.get("mixscape_global"),
                    "edistance_top": report.get("edistance_top"),
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
            if "de" in st:
                report["de"] = st["de"]
            if "de_scope" in st:
                report["de_scope"] = st["de_scope"]
            if "de_groups" in st:
                report["de_groups"] = st["de_groups"]
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
                    "de_scope": report.get("de_scope"),
                    "de": report.get("de"),
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
    "dry_run_plan",
    "pipeline_config_from_mapping",
    "run_pipeline",
    "stage_input_validation",
    "stage_perturbation_modeling",
    "stage_preprocessing_qc",
    "stage_report",
    "stage_robustness",
    "stage_statistical_inference",
    # Compat for tests / internal callers
    "_effective_mixscape_mode",
    "_embedding_provenance",
]
