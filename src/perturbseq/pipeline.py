"""Run the single-sample tertiary analysis one step at a time."""

from __future__ import annotations

import argparse
import json
import platform
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from perturbseq import __version__
from perturbseq.assign import (
    LOW_LIBRARY_NO_GUIDE,
    STATUS_NO_GUIDE,
    STATUS_NT,
    STATUS_RETAINED,
    STATUS_TARGET,
    GuideReferenceError,
    apply_eligibility,
    capture_library_correlation,
    cells_per_guide,
    cells_per_target,
    derive_qc_thresholds,
    filter_stage_counts,
    low_library_no_guide_rule,
    multi_guide_composition,
    validate_guide_targets,
)
from perturbseq.config import (
    TEMPORARY_CONTROL_WARNING,
    ConfigError,
    apply_customer_arguments,
    control_spec,
    mito_prefix,
    output_dir,
    read_config,
    resolve_input_paths,
    ribo_prefixes,
    software_config_path,
    validate_config,
    verify_species_consistency,
)
from perturbseq.distance_tests import write_edistance_results
from perturbseq.io import (
    InputMismatch,
    load_sample,
    matrix_only_as_no_guide,
    read_features,
    read_reference,
)
from perturbseq.mixscale import run_mixscale
from perturbseq.mixscape import mixscape_fraction, run_mixscape
from perturbseq.perturb_effect import (
    descriptive_effects,
    guide_pseudobulk,
    perturbation_distance,
    target_direction_table,
    tutorial_differential_expression,
)
from perturbseq.preprocess import preprocess
from perturbseq.report import write_report

STEPS = ["io", "assign", "preprocess", "mixscape", "mixscale", "perturb", "report"]
PACKAGES = [
    "perturbseq",
    "scanpy",
    "anndata",
    "pertpy",
    "pandas",
    "numpy",
    "scipy",
    "pydeseq2",
    "scikit-learn",
    "matplotlib",
    "seaborn",
    "statsmodels",
    "PyYAML",
]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Single-sample Perturb-seq tertiary analysis")
    parser.add_argument(
        "--config",
        default=None,
        help="分析软件默认参数 YAML。缺省为仓库里的 config/config.yaml。样本和路径以命令行为准。",
    )
    parser.add_argument("--sample-id", required=True, help="样本名，用来拼接 DRAGEN 文件名")
    parser.add_argument("--organism", required=True, help="human 或 mouse")
    parser.add_argument("--input-dir", required=True, help="五个 DRAGEN 输入文件所在目录")
    parser.add_argument("--output-dir", required=True, help="结果输出目录")
    parser.add_argument("--guide-separator", required=True, help="多条 guide 之间的分隔符，例如 |")
    parser.add_argument(
        "--nt-control", default="", help="设计好的非靶向对照名。不填则用 zero-guide 细胞"
    )
    parser.add_argument(
        "--crispr-mode", default="CRISPR_KO", help="KO、i、a 或 other，缺省 CRISPR_KO"
    )
    parser.add_argument(
        "--crispr-description", default=None, help="crispr-mode 为 other 时的文字说明"
    )
    parser.add_argument(
        "--expected-direction",
        default=None,
        help="other 模式的预期方向：down、up、bidirectional 或 unknown",
    )
    parser.add_argument("--batch-key", default=None, help="obs 里的批次列。不填则不做 Harmony")
    parser.add_argument(
        "--cell-type-key", default=None, help="obs 里已有的细胞类型列。不填则全程 ALL"
    )
    parser.add_argument("--min-umi", type=float, default=None)
    parser.add_argument("--max-umi", type=float, default=None)
    parser.add_argument(
        "--max-mito-frac",
        type=float,
        default=None,
        help="线粒体 reads 比例上限，0 到 1。细胞 ≥ 该值时剔除。默认读配置（0.2 = 20%）。",
    )
    parser.add_argument(
        "--min-genes",
        type=float,
        default=None,
        help="表达基因数下限。细胞 n_genes ≤ 该值时剔除。默认读配置（200）。",
    )
    parser.add_argument(
        "--max-genes",
        type=float,
        default=None,
        help="表达基因数上限。细胞 n_genes ≥ 该值时剔除。默认读配置（6000）。",
    )
    parser.add_argument(
        "--max-ribo-frac",
        type=float,
        default=None,
        help="核糖体蛋白基因比例上限，0 到 1。不填则用 median + n_mad × MAD，不并入 UMI 绝对阈值。",
    )
    parser.add_argument("--steps", default="all", help="all, or a comma-separated subset")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    config_path = args.config or software_config_path()
    config = read_config(config_path)
    apply_customer_arguments(
        config,
        {
            "sample_id": args.sample_id,
            "organism": args.organism,
            "input_dir": args.input_dir,
            "output_dir": args.output_dir,
            "guide_separator": args.guide_separator,
            "nt_control": args.nt_control,
            "crispr_mode": args.crispr_mode,
            "crispr_mode_description": args.crispr_description,
            "expected_target_direction": args.expected_direction,
            "batch_key": args.batch_key,
            "cell_type_key": args.cell_type_key,
            "min_umi": args.min_umi,
            "max_umi": args.max_umi,
            "max_mito_frac": args.max_mito_frac,
            "min_genes": args.min_genes,
            "max_genes": args.max_genes,
            "max_ribo_frac": args.max_ribo_frac,
        },
    )
    validate_config(config)
    if args.validate_only:
        folder = output_dir(config)
        folder.mkdir(parents=True, exist_ok=True)
        _record_species(folder, config)
        _validate_guide_reference(folder, config)
        _validate_barcodes(folder, config)
        print(f"Validated {folder}", flush=True)
        return
    steps = (
        STEPS
        if args.steps.strip() == "all"
        else [step.strip() for step in args.steps.split(",") if step.strip()]
    )
    unknown = [step for step in steps if step not in STEPS]
    if unknown:
        raise SystemExit(f"未知步骤: {', '.join(unknown)}")
    run(config, steps, resume=args.resume)


def run(config: dict[str, Any], steps: list[str], *, resume: bool) -> None:
    folder = output_dir(config)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "tables").mkdir(exist_ok=True)
    (folder / "checkpoints").mkdir(exist_ok=True)
    (folder / "pseudobulk").mkdir(exist_ok=True)
    _record_species(folder, config)
    seed = int(config["runtime"]["seed"])
    np.random.seed(seed)
    info = _load_info(folder)
    info.update(
        {
            "version": __version__,
            "sample_id": config["experiment"]["sample_id"],
            "python": platform.python_version(),
            "packages": _package_versions(),
            "config": config,
            "seed": seed,
            "biological_sample_count": 1,
            "analysis_scope": "single_sample_exploratory",
            "document_version": config.get("document_version"),
            "resolved_organism": config["experiment"].get("resolved_organism", ""),
            "resolved_mito_prefix": config["experiment"].get("resolved_mito_prefix", ""),
            "species_check": config["experiment"].get("species_check", ""),
            "organism_evidence": config["experiment"].get("organism_evidence", []),
            "other_organelle_prefixes": config["experiment"].get("other_organelle_prefixes", []),
            **_control_summary(config),
        }
    )
    state: dict[str, Any] = {"adata": None}
    if resume:
        state["adata"] = _load_resume(folder, steps[0])
    started = time.time()
    info["temporary_control_warning"] = _publish_temporary_control_warning(folder, config)
    _write_info(folder, info, started)
    for step in steps:
        print(f"[{step}]", flush=True)
        step_started = time.time()
        if step == "io":
            state["adata"] = _step_io(folder, config)
            validation = json.loads((folder / "input_validation.json").read_text())
            info["input_validation"] = validation
            info["n_matrix_only_as_no_guide"] = validation.get("n_matrix_only_as_no_guide", 0)
        elif step == "assign":
            state["adata"] = _require(state, folder, "io")
            state["adata"] = _step_assign(folder, config, state["adata"])
            info.update(_control_summary(config, state["adata"]))
        elif step == "preprocess":
            state["adata"] = _require(state, folder, "assign")
            state["adata"] = preprocess(state["adata"], config)
            _score_cell_cycle(state["adata"], config)
            info["cell_cycle_status"] = state["adata"].uns.get("cell_cycle_status", "")
            info["cell_cycle_note"] = state["adata"].uns.get("cell_cycle_note", "")
            _write_h5ad(state["adata"], folder / "checkpoints" / "preprocess.h5ad")
        elif step == "mixscape":
            state["adata"] = _require(state, folder, "preprocess")
            state["adata"] = _step_mixscape(folder, config, state["adata"])
        elif step == "mixscale":
            state["adata"] = _require(state, folder, "mixscape")
            state["adata"] = _step_mixscale(folder, config, state["adata"])
        elif step == "perturb":
            state["adata"] = _require(state, folder, "mixscale")
            state["adata"] = _step_perturb(folder, config, state["adata"])
        elif step == "report":
            state["adata"] = (
                state["adata"] if state["adata"] is not None else _load_resume(folder, "report")
            )
            write_report(folder, config, state["adata"])
        info.setdefault("steps", {})[step] = {"seconds": round(time.time() - step_started, 3)}
        _write_info(folder, info, started)
    print(f"Wrote {folder}", flush=True)


def _step_io(folder: Path, config: dict[str, Any]):
    paths = resolve_input_paths(config)
    _validate_guide_reference(folder, config)
    min_guide_umi = config["qc"].get("min_guide_umi")
    try:
        adata, annotated, checks = load_sample(
            paths,
            nt_control=config["labels"]["nt_control"],
            feature_type_gex=config["labels"].get("feature_type_gex", "Gene Expression"),
            feature_type_crispr=config["labels"]["feature_type_crispr"],
            multi_guide_policy=config["labels"]["multi_guide_policy"],
            targets_per_cell=int(config["labels"].get("targets_per_cell", 1)),
            min_guide_umi=None if min_guide_umi is None else float(min_guide_umi),
            mito_prefix=mito_prefix(config),
            ribo_prefixes=ribo_prefixes(config),
        )
    except InputMismatch as exc:
        _write_barcode_differences(folder, exc.only_matrix, exc.only_assignment)
        _write_input_validation(
            folder,
            config,
            n_matrix=exc.n_matrix,
            n_assignment=exc.n_assignment,
            only_matrix=exc.only_matrix,
            only_assignment=exc.only_assignment,
        )
        raise SystemExit(str(exc)) from exc
    only_matrix = [str(barcode) for barcode in checks.pop("matrix_only_barcodes", [])]
    _write_barcode_differences(folder, only_matrix, [])
    _write_input_validation(
        folder,
        config,
        n_matrix=int(checks["n_barcodes"]),
        n_assignment=int(checks["n_assignment_rows"]),
        only_matrix=only_matrix,
        only_assignment=[],
    )
    if only_matrix:
        print(
            f"  matrix-only barcodes recorded as no_guide: {len(only_matrix)}",
            flush=True,
        )
    _calls(annotated).to_csv(folder / "tables" / "cell_calls.tsv", sep="\t", index=False)
    (folder / "tables" / "io_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    reference = read_reference(paths["feature_reference"])
    composition = multi_guide_composition(
        annotated,
        nt_control=config["labels"]["nt_control"],
        reference=reference,
    )
    composition.to_csv(folder / "tables" / "multi_guide_composition.tsv", sep="\t", index=False)
    _write_h5ad(adata, folder / "checkpoints" / "io.h5ad")
    return adata


def _step_assign(folder: Path, config: dict[str, Any], adata):
    calls = pd.read_csv(folder / "tables" / "cell_calls.tsv", sep="\t", keep_default_na=False)
    calls["guide_umi_pass"] = _as_bool(calls["guide_umi_pass"])
    for column in (
        "total_counts",
        "n_genes_by_counts",
        "mito_frac",
        "ribo_frac",
        "guide_umi_total",
        "grna_capture_efficiency",
    ):
        if column in calls.columns:
            calls[column] = pd.to_numeric(calls[column], errors="coerce")
    metrics = pd.DataFrame(
        {
            "total_counts": calls["total_counts"],
            "n_genes_by_counts": calls["n_genes_by_counts"],
            "pct_counts_mito": calls["mito_frac"],
            "pct_counts_ribo": calls["ribo_frac"],
        }
    )
    thresholds = derive_qc_thresholds(metrics, config["qc"])
    n_mad = float(config["qc"].get("n_mad", 3))
    umi_pair = ("min_umi", "max_umi")
    configured = all(config["qc"].get(key) is not None for key in umi_pair)
    library_rule = low_library_no_guide_rule(calls["total_counts"], calls["n_guides"], n_mad)
    correlation = capture_library_correlation(calls)
    correlation["low_library_no_guide_applied"] = bool(library_rule["applied"])
    correlation["low_library_no_guide_reason"] = library_rule["reason"]
    correlation["low_library_min_total_counts"] = library_rule["min_total_counts"]
    recorded = {
        **thresholds,
        "n_mad": n_mad,
        "source": "command_line" if configured else "median_mad",
        "gene_bounds_source": "configured",
        "gene_filter": "drop_n_genes_le_min_or_ge_max",
        "mito_frac_source": "configured",
        "mito_filter": "drop_mito_frac_ge_max",
        "ribo_frac_source": (
            "command_line" if config["qc"].get("max_ribo_frac") is not None else "median_mad"
        ),
        "ribo_prefixes": list(ribo_prefixes(config)),
        "ribo_match": "gene_symbol_prefix_case_insensitive",
        "low_library_no_guide": library_rule,
    }
    (folder / "tables" / "qc_thresholds.json").write_text(
        json.dumps(recorded, indent=2, ensure_ascii=False) + "\n"
    )
    correlation.to_csv(
        _qc_tables(folder) / "grna_capture_library_correlation.tsv", sep="\t", index=False
    )
    if library_rule["applied"]:
        print(
            f"  low_library_no_guide threshold: {library_rule['min_total_counts']}",
            flush=True,
        )
    else:
        print(f"  low_library_no_guide skipped: {library_rule['reason']}", flush=True)
    updated = apply_eligibility(
        calls,
        thresholds,
        min_cells_per_target=int(config["qc"]["min_cells_per_target"]),
        targets_per_cell=int(config["labels"].get("targets_per_cell", 1)),
        low_library_rule=library_rule,
    )
    updated.insert(0, "sample_id", config["experiment"]["sample_id"])
    spec = control_spec(config)
    temporary = pd.Series(False, index=updated.index)
    if spec["temporary"]:
        temporary = (
            updated["assignment_status"].astype(str).eq(STATUS_NO_GUIDE)
            & updated["n_guides"].astype(int).eq(0)
            & updated["qc_pass"].astype(bool)
            & updated["exclusion_reason_primary"].astype(str).ne(LOW_LIBRARY_NO_GUIDE)
        )
    updated["temporary_control"] = temporary.to_numpy()
    updated["control_source"] = spec["control_source"]
    updated["control_label"] = spec["control_label"]
    stages = filter_stage_counts(updated)
    retained = int((updated["exclusion_reason_primary"] == STATUS_RETAINED).sum())
    analysis_n = int(updated["analysis_eligible"].astype(bool).sum())
    targeting_n = int(
        (
            updated["analysis_eligible"].astype(bool)
            & (updated["perturbation_status"] == STATUS_TARGET)
        ).sum()
    )
    nt_n = int(
        (
            updated["analysis_eligible"].astype(bool)
            & (updated["perturbation_status"] == STATUS_NT)
        ).sum()
    )
    if analysis_n != retained or retained != targeting_n + nt_n:
        raise RuntimeError("细胞账本不守恒。")
    if len(updated) != int(adata.uns.get("n_input_barcodes", len(updated))):
        pass
    ledger = _cell_status(updated)
    ledger.to_csv(folder / "cell_status.tsv", sep="\t", index=False)
    metric_columns = [
        "cell_barcode",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mito",
        "pct_counts_ribo",
        "guide_umi_total",
        "grna_capture_efficiency",
        "n_guides",
        "n_targets",
        "exclusion_reason_primary",
        "perturbation_status",
    ]
    ledger[metric_columns].to_csv(_qc_tables(folder) / "qc_metrics.tsv", sep="\t", index=False)
    stages.to_csv(folder / "filter_stage_counts.tsv", sep="\t", index=False)
    assignment = (
        updated["assignment_status"]
        .value_counts()
        .rename_axis("assignment_status")
        .reset_index(name="n_cells")
    )
    assignment.to_csv(_qc_tables(folder) / "assignment_summary.tsv", sep="\t", index=False)
    cells_per_guide(updated).to_csv(
        _qc_tables(folder) / "cells_per_guide.tsv", sep="\t", index=False
    )
    cells_per_target(updated, int(config["qc"]["min_cells_per_guide"])).to_csv(
        _qc_tables(folder) / "cells_per_target.tsv", sep="\t", index=False
    )
    eligible = set(
        updated.loc[updated["analysis_eligible"].astype(bool), "cell_barcode"].astype(str)
    )
    control_ids = set(
        updated.loc[updated["temporary_control"].astype(bool), "cell_barcode"].astype(str)
    )
    adata = adata[adata.obs_names.astype(str).isin(eligible | control_ids)].copy()
    if control_ids:
        hit = np.isin(adata.obs_names.astype(str), list(control_ids))
        for column in ("target_gene", "perturbation_status"):
            values = adata.obs[column].astype(str).to_numpy(copy=True)
            values[hit] = spec["control_label"]
            adata.obs[column] = values
    adata.uns["analysis_scope"] = "single_sample_exploratory"
    adata.uns["config"] = config
    adata.uns["control_source"] = spec["control_source"]
    adata.uns["control_label"] = spec["control_label"]
    adata.uns["control_note"] = spec["note"]
    adata.uns["n_temporary_control_cells"] = int(len(control_ids))
    if spec["temporary"]:
        print(
            f"  temporary control {spec['control_label']}: {len(control_ids):,} zero-guide cells",
            flush=True,
        )
    _write_h5ad(adata, folder / "checkpoints" / "assign.h5ad")
    print(f"  analysis cells: {adata.n_obs:,}", flush=True)
    return adata


def _record_species(folder: Path, config: dict[str, Any]) -> None:
    """Check the declared species, then write it even when the check fails."""
    try:
        _resolve_species(config)
    except ConfigError:
        _write_resolved(folder, config)
        raise
    _write_resolved(folder, config)


def _resolve_species(config: dict[str, Any]) -> None:
    """Read gene ids and symbols and check them against the input organism."""
    paths = resolve_input_paths(config)
    features = read_features(paths["features"])
    gex_name = str(config["labels"].get("feature_type_gex", "Gene Expression"))
    if "feature_type" in features.columns:
        gex = features["feature_type"].astype(str).eq(gex_name)
        if bool(gex.any()):
            features = features.loc[gex]
    verify_species_consistency(
        config,
        features["gene_id"].astype(str),
        features["gene_symbol"].astype(str),
    )


def _step_mixscape(folder: Path, config: dict[str, Any], adata):
    """Independent Mixscape step: DE gate, X_pert signature, classification."""
    mixscape_on = config.get("mixscape", {}).get("enabled", True)
    spec = control_spec(config)
    if spec["temporary"]:
        print(f"  temporary control warning: {spec['note']}", flush=True)
    if not mixscape_on:
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
        adata.obs["mixscape_class_raw"] = "not_evaluated"
        _write_de_gate(folder, adata)
        fractions = mixscape_fraction(adata)
        fractions.to_csv(_mixscape_tables(folder) / "mixscape_summary.tsv", sep="\t", index=False)
        _write_h5ad(adata, folder / "checkpoints" / "mixscape.h5ad")
        return adata
    if not config.get("mixscape", {}).get("cell_type_key") and not config["perturb"].get(
        "cell_type_key"
    ):
        print("  Mixscape scope: ALL cells, not stratified", flush=True)
    try:
        adata = run_mixscape(adata, config)
        adata.obs["mixscape_class_raw"] = adata.obs["mixscape_class"].astype(str)
        if "mixscape_class_p" in adata.obs:
            adata.obs["mixscape_probability"] = adata.obs["mixscape_class_p"]
        _write_de_gate(folder, adata)
    except Exception as exc:
        print(f"  Mixscape failed: {exc}", flush=True)
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
        adata.obs["mixscape_class_raw"] = "not_evaluated"
        adata.uns["mixscape_skip_reason"] = f"mixscape_failed:{exc}"
        _write_de_gate(folder, adata)
    fractions = mixscape_fraction(adata)
    summary_path = _mixscape_tables(folder) / "mixscape_summary.tsv"
    fractions.to_csv(summary_path, sep="\t", index=False)
    for stale in (
        folder / "DGE" / "tables" / "mixscape_summary.tsv",
        folder / "DGE" / "mixscape" / "tables" / "mixscape_summary.tsv",
        folder / "tables" / "mixscape_fraction.tsv",
    ):
        if stale.is_file() and stale != summary_path:
            stale.unlink()
    _refresh_mixscape_ledger(folder, adata)
    _write_h5ad(adata, folder / "checkpoints" / "mixscape.h5ad")
    return adata


def _step_mixscale(folder: Path, config: dict[str, Any], adata):
    """Independent Mixscale step: own DE gate and mixscale_score on X_pert."""
    if "X_pert" not in adata.layers:
        print("  Mixscale skipped: missing X_pert (run mixscape first)", flush=True)
        adata.uns["mixscale_status"] = "skipped:missing_X_pert"
        _write_mixscale_gate(folder, adata)
        _write_h5ad(adata, folder / "checkpoints" / "mixscale.h5ad")
        return adata
    try:
        adata = run_mixscale(adata, config)
    except Exception as mixscale_exc:
        print(f"  Mixscale failed: {mixscale_exc}", flush=True)
        adata.uns["mixscale_status"] = f"failed:{mixscale_exc}"
    _write_mixscale_gate(folder, adata)
    _write_h5ad(adata, folder / "checkpoints" / "mixscale.h5ad")
    return adata


def _step_perturb(folder: Path, config: dict[str, Any], adata):
    """Descriptive effects, pseudobulk distances, and E-test (after Mixscape/Mixscale)."""
    comparisons = []
    spec = control_spec(config)
    if spec["temporary"]:
        print(f"  temporary control warning: {spec['note']}", flush=True)
    if "mixscape_status" not in adata.obs:
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
    main_effects, main_guides, main_status = descriptive_effects(adata, config, branch="main")
    sens_effects, sens_guides, sens_status = descriptive_effects(
        adata, config, branch="sensitivity"
    )
    descriptive = _descriptive_tables(folder)
    main_effects.to_csv(descriptive / "effect_main.tsv", sep="\t", index=False)
    sens_effects.to_csv(descriptive / "effect_sensitivity.tsv", sep="\t", index=False)
    guides = pd.concat([main_guides, sens_guides], ignore_index=True)
    guides.to_csv(descriptive / "guide_effects.tsv", sep="\t", index=False)
    direction = pd.concat(
        [
            target_direction_table(adata, config, branch="main"),
            target_direction_table(adata, config, branch="sensitivity"),
        ],
        ignore_index=True,
    )
    direction.to_csv(descriptive / "target_expression_effect.tsv", sep="\t", index=False)
    de_main, de_main_status = tutorial_differential_expression(adata, config, branch="main")
    de_sens, de_sens_status = tutorial_differential_expression(adata, config, branch="sensitivity")
    for branch_name, tables in (("main", de_main), ("sensitivity", de_sens)):
        for method, frame in tables.items():
            name = f"de_{method}_{branch_name}.tsv"
            frame.to_csv(descriptive / name, sep="\t", index=False)
            legacy = folder / "DGE" / "tables" / name
            if legacy.is_file():
                legacy.unlink()
    for name in (
        "effect_main.tsv",
        "effect_sensitivity.tsv",
        "guide_effects.tsv",
        "target_expression_effect.tsv",
    ):
        legacy = folder / name
        if legacy.is_file():
            legacy.unlink()
    comparisons_table = pd.concat(
        [main_status, sens_status, de_main_status, de_sens_status, pd.DataFrame(comparisons)],
        ignore_index=True,
    )
    comparisons_table.to_csv(folder / "comparison_status.tsv", sep="\t", index=False)
    pseudobulk, guide_table = guide_pseudobulk(adata, config, branch="main")
    guide_table.to_csv(folder / "pseudobulk" / "guide_metadata.tsv", sep="\t", index=False)
    _write_h5ad(pseudobulk, folder / "pseudobulk" / "guide_counts.h5ad")
    _write_guide_distance(folder, pseudobulk, config)
    write_edistance_results(folder, adata, config)
    _refresh_mixscape_ledger(folder, adata)
    _write_h5ad(adata, folder / "checkpoints" / "perturb.h5ad")
    return adata


def _write_guide_distance(folder: Path, pseudobulk, config: dict[str, Any]) -> None:
    """Skip the guide-pseudobulk distance when no replicate minimum was configured.

    Guides are not biological replicates, so a missing minimum is not filled in.
    The cell-level E-test is written later and is not this table.
    """
    tables = _edistance_tables(folder)
    skip_path = tables / "perturbation_distance_heatmap.skipped.txt"
    if "min_guide_replicates" not in config.get("perturb", {}):
        skip_path.write_text("min_guide_replicates_not_set\n")
        _retire_legacy_distance_skip(folder)
        return
    distance, embedding = perturbation_distance(pseudobulk, config)
    if distance.empty:
        skip_path.write_text("insufficient_cells\n")
        _retire_legacy_distance_skip(folder)
        return
    if skip_path.exists():
        skip_path.unlink()
    distance.to_csv(tables / "perturbation_distance.tsv", sep="\t")
    embedding.to_csv(tables / "perturbation_embedding.tsv", sep="\t", index=False)
    for stale in (
        folder / "perturbation_distance.tsv",
        folder / "perturbation_embedding.tsv",
        folder / "plot_data" / "perturbation_distance_heatmap.skipped.txt",
    ):
        if stale.is_file():
            stale.unlink()


def _score_cell_cycle(adata, config: dict[str, Any]) -> None:
    """Score human or mouse cell cycle on the log-normalized matrix.

    Uses scib gene lists. Does not alter ``X`` or regress out the cycle.
    """
    from perturbseq.plotting import score_cell_cycle_after_log1p

    organism = str(config["experiment"].get("organism", ""))
    try:
        status = score_cell_cycle_after_log1p(
            adata,
            organism=organism,
            random_state=int(config["runtime"]["seed"]),
        )
    except ValueError as exc:
        print(f"  cell cycle skipped: {exc}", flush=True)
        adata.uns["cell_cycle_status"] = "skipped"
        adata.uns["cell_cycle_note"] = str(exc)
        return
    adata.uns["cell_cycle_status"] = status
    notes = {
        "scored_after_log1p": "scib Tirosh/Regev 基因集；写入 phase、S_score、G2M_score。",
        "already_scored": "obs 已有 phase，跳过重复打分。",
        "skipped_unsupported_organism": "organism 不是 human 或 mouse。",
        "skipped_missing_log1p": "尚未 log1p，跳过打分。",
        "skipped_missing_cell_cycle_genes": "矩阵里几乎找不到该物种的细胞周期基因。",
        "skipped_scib_unavailable": "未安装 scib，无法读取细胞周期基因集。",
    }
    adata.uns["cell_cycle_note"] = notes.get(status, "")
    print(f"  cell cycle: {status}", flush=True)


def _qc_tables(folder: Path) -> Path:
    path = folder / "qc" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _mixscape_tables(folder: Path) -> Path:
    path = folder / "DGE" / "Mixscape" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _mixscale_tables(folder: Path) -> Path:
    path = folder / "DGE" / "mixscale" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _descriptive_tables(folder: Path) -> Path:
    path = folder / "DGE" / "descriptive" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _edistance_tables(folder: Path) -> Path:
    path = folder / "E-distance" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _retire_legacy_distance_skip(folder: Path) -> None:
    legacy = folder / "plot_data" / "perturbation_distance_heatmap.skipped.txt"
    if legacy.is_file():
        legacy.unlink()


def _write_de_gate(folder: Path, adata) -> None:
    gate = adata.uns.get("mixscape_de_gate")
    if not isinstance(gate, pd.DataFrame):
        return
    path = _mixscape_tables(folder) / "mixscape_de_gate.tsv"
    gate.to_csv(path, sep="\t", index=False)
    for stale in (
        folder / "DGE" / "tables" / "mixscape_de_gate.tsv",
        folder / "tables" / "mixscape_de_gate.tsv",
        folder / "DGE" / "mixscape" / "tables" / "mixscape_de_gate.tsv",
    ):
        if stale.is_file() and stale != path:
            stale.unlink()


def _write_mixscale_gate(folder: Path, adata) -> None:
    gate = adata.uns.get("mixscale_de_gate")
    if not isinstance(gate, pd.DataFrame):
        return
    gate.to_csv(_mixscale_tables(folder) / "mixscale_de_gate.tsv", sep="\t", index=False)


def _write_resolved(folder: Path, config: dict[str, Any]) -> None:
    public = {key: value for key, value in config.items() if not str(key).startswith("_")}
    (folder / "resolved_config.yaml").write_text(
        yaml.safe_dump(public, sort_keys=False, allow_unicode=True)
    )


def _validate_guide_reference(folder: Path, config: dict[str, Any]) -> None:
    paths = resolve_input_paths(config)
    reference = read_reference(paths["feature_reference"])
    try:
        validate_guide_targets(reference, nt_control=config["labels"]["nt_control"])
    except GuideReferenceError as exc:
        _write_guide_conflicts(folder, exc)
        raise SystemExit(str(exc)) from exc


def _write_guide_conflicts(folder: Path, exc: GuideReferenceError) -> None:
    (folder / "tables").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"guide_id": exc.guide_ids, "reason": exc.reason}).to_csv(
        folder / "tables" / "guide_reference_conflicts.tsv", sep="\t", index=False
    )


def _validate_barcodes(folder: Path, config: dict[str, Any]) -> None:
    from perturbseq.io import read_assignments, read_barcodes

    paths = resolve_input_paths(config)
    matrix_ids = set(read_barcodes(paths["barcodes"]).astype(str))
    assigned = set(read_assignments(paths["assignments"])["cell_barcode"].astype(str))
    only_matrix = sorted(matrix_ids - assigned)
    only_assignment = sorted(assigned - matrix_ids)
    _write_barcode_differences(folder, only_matrix, only_assignment)
    payload = _write_input_validation(
        folder,
        config,
        n_matrix=len(matrix_ids),
        n_assignment=len(assigned),
        only_matrix=only_matrix,
        only_assignment=only_assignment,
    )
    if payload["status"] == "failed":
        raise SystemExit(
            "矩阵与分配 barcode 集合不一致: "
            f"仅矩阵 {len(only_matrix)}，仅分配 {len(only_assignment)}。"
        )


def _write_input_validation(
    folder: Path,
    config: dict[str, Any],
    *,
    n_matrix: int,
    n_assignment: int,
    only_matrix: list[str],
    only_assignment: list[str],
) -> dict[str, Any]:
    failed = len(only_assignment) > 0
    payload = {
        "barcode_policy": "strict",
        "status": "failed" if failed else "passed",
        "n_matrix": n_matrix,
        "n_assignment": n_assignment,
        "n_only_matrix": len(only_matrix),
        "n_only_assignment": len(only_assignment),
        "n_matrix_only_as_no_guide": matrix_only_as_no_guide(only_matrix, only_assignment),
        "guide_separator": config["input"]["guide_separator"],
        "biological_sample_count": 1,
        **_control_summary(config),
    }
    (folder / "input_validation.json").write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def _control_summary(config: dict[str, Any], adata=None) -> dict[str, Any]:
    spec = control_spec(config)
    summary = {
        "control_source": spec["control_source"],
        "control_label": spec["control_label"],
        "control_note": spec["note"],
    }
    if adata is not None:
        summary["n_temporary_control_cells"] = int(adata.uns.get("n_temporary_control_cells", 0))
    return summary


def _publish_temporary_control_warning(folder: Path, config: dict[str, Any]) -> str:
    """Record the zero-guide warning in the log and on disk.

    A designed NT control removes that warning and the temporary-control QC figure.
    """
    path = folder / "tables" / "temporary_control_warning.txt"
    if not control_spec(config)["temporary"]:
        if path.is_file():
            path.unlink()
        from perturbseq.plotting import remove_temporary_control_qc

        remove_temporary_control_qc(folder)
        return ""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TEMPORARY_CONTROL_WARNING.strip() + "\n")
    print(f"  temporary control warning: {TEMPORARY_CONTROL_WARNING}", flush=True)
    return TEMPORARY_CONTROL_WARNING


def _write_barcode_differences(
    folder: Path, only_matrix: list[str], only_assignment: list[str]
) -> None:
    rows = [{"cell_barcode": barcode, "present_in": "matrix_only"} for barcode in only_matrix]
    rows.extend(
        {"cell_barcode": barcode, "present_in": "assignment_only"} for barcode in only_assignment
    )
    pd.DataFrame(rows, columns=["cell_barcode", "present_in"]).to_csv(
        folder / "barcode_differences.tsv", sep="\t", index=False
    )


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().isin(["true", "1"])


def _cell_status(table: pd.DataFrame) -> pd.DataFrame:
    ledger = table.copy()
    ledger["pct_counts_mito"] = ledger["mito_frac"]
    if "ribo_frac" in ledger.columns:
        ledger["pct_counts_ribo"] = ledger["ribo_frac"]
    for column in (
        "mixscape_class_raw",
        "mixscape_class",
        "mixscape_probability",
        "mixscape_status",
        "mixscape_reason",
        "sensitivity_eligible",
    ):
        if column not in ledger.columns:
            ledger[column] = pd.NA
    return ledger


def _refresh_mixscape_ledger(folder: Path, adata) -> None:
    path = folder / "cell_status.tsv"
    if not path.exists():
        return
    ledger = pd.read_csv(path, sep="\t", keep_default_na=False)
    labels = adata.obs.copy()
    labels["cell_barcode"] = adata.obs_names.astype(str)
    mapped = labels.set_index("cell_barcode")
    barcodes = ledger["cell_barcode"].astype(str)
    for column in ("mixscape_class", "mixscape_status", "mixscape_probability"):
        if column in mapped.columns:
            ledger[column] = barcodes.map(mapped[column])
    if "mixscape_class_raw" in mapped.columns:
        ledger["mixscape_class_raw"] = barcodes.map(mapped["mixscape_class_raw"])
    in_analysis = barcodes.isin(set(mapped.index.astype(str)))
    perturbed = ledger["mixscape_class"].astype(str).eq("perturbed") & in_analysis
    temporary = (
        ledger["temporary_control"].astype(str).str.lower().isin(["true", "1"])
        if "temporary_control" in ledger.columns
        else pd.Series(False, index=ledger.index)
    )
    nt = (
        ledger["perturbation_status"].astype(str).eq(STATUS_NT) | temporary.to_numpy()
    ) & in_analysis
    ledger["sensitivity_eligible"] = perturbed | nt
    ledger.loc[~in_analysis, "mixscape_status"] = ledger.loc[
        ~in_analysis, "mixscape_status"
    ].replace("", "not_evaluated")
    ledger.to_csv(path, sep="\t", index=False)


def _note_untested_targets(notes: pd.DataFrame, adata, control: str) -> pd.DataFrame:
    """Targets with too few perturbed guide replicates never enter the model."""
    present = (
        set(notes["target_gene"].astype(str))
        if not notes.empty and "target_gene" in notes
        else set()
    )
    targeting = adata.obs.loc[
        adata.obs["perturbation_status"].astype(str) == STATUS_TARGET, "target_gene"
    ]
    rows = []
    for target, count in targeting.astype(str).value_counts().items():
        if target in present or target == control:
            continue
        rows.append(
            {
                "target_gene": target,
                "reason": "insufficient_replicates",
                "n_guides": 0,
                "n_cells": int(count),
            }
        )
    if not rows:
        return notes
    return pd.concat([notes, pd.DataFrame(rows)], ignore_index=True)


def _require(state: dict[str, Any], folder: Path, previous: str):
    if state["adata"] is not None:
        return state["adata"]
    path = folder / "checkpoints" / f"{previous}.h5ad"
    if not path.exists():
        raise SystemExit(f"缺少上一步结果 {path}。加上 --resume 前需要先跑完 {previous}。")
    print(f"  resume {path}", flush=True)
    return sc.read_h5ad(path)


def _load_resume(folder: Path, step: str):
    previous = {
        "assign": "io",
        "preprocess": "assign",
        "mixscape": "preprocess",
        "mixscale": "mixscape",
        "perturb": "mixscale",
        "report": "perturb",
    }
    if step not in previous:
        return None
    path = folder / "checkpoints" / f"{previous[step]}.h5ad"
    fallback = {
        "mixscale": ["mixscape.h5ad", "preprocess.h5ad"],
        "perturb": ["mixscale.h5ad", "mixscape.h5ad", "preprocess.h5ad"],
        "report": [
            "perturb.h5ad",
            "mixscale.h5ad",
            "mixscape.h5ad",
            "preprocess.h5ad",
            "assign.h5ad",
        ],
    }
    candidates = (
        [path]
        if step not in fallback
        else [folder / "checkpoints" / name for name in fallback[step]]
    )
    for candidate in candidates:
        if candidate.exists():
            print(f"  resume {candidate}", flush=True)
            return sc.read_h5ad(candidate)
    return None


def _calls(annotated: pd.DataFrame) -> pd.DataFrame:
    drop = [
        column for column in ("guide_names", "guide_transcripts") if column in annotated.columns
    ]
    return annotated.drop(columns=drop)


def _write_h5ad(adata, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="lzf")


def _write_info(folder: Path, info: dict[str, Any], started: float) -> None:
    info["elapsed_seconds"] = round(time.time() - started, 3)
    (folder / "run_info.json").write_text(json.dumps(info, indent=2, default=str) + "\n")


def _load_info(folder: Path) -> dict[str, Any]:
    path = folder / "run_info.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _package_versions() -> dict[str, str]:
    found = {}
    for name in PACKAGES:
        if name == "perturbseq":
            found[name] = __version__
            continue
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            found[name] = "not_installed"
    return found


def _safe(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._+-" else "_" for character in value
    )
