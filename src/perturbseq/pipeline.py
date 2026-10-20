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
    STATUS_NO_GUIDE,
    STATUS_NT,
    STATUS_RETAINED,
    STATUS_TARGET,
    GuideReferenceError,
    apply_eligibility,
    cells_per_guide,
    cells_per_target,
    filter_stage_counts,
    multi_guide_composition,
    validate_guide_targets,
)
from perturbseq.config import (
    ConfigError,
    control_spec,
    load_config,
    mito_prefix,
    output_dir,
    resolve_input_paths,
    verify_species_consistency,
)
from perturbseq.distance_tests import write_edistance_results
from perturbseq.enrichment import enrich_perturbations
from perturbseq.io import (
    InputMismatch,
    load_sample,
    matrix_only_as_no_guide,
    read_features,
    read_reference,
)
from perturbseq.perturb_effect import (
    descriptive_effects,
    guide_pseudobulk,
    mixscape_fraction,
    perturbation_distance,
    run_mixscale,
    run_mixscape,
    target_direction_table,
    tutorial_differential_expression,
)
from perturbseq.preprocess import preprocess
from perturbseq.report import write_report
from perturbseq.responsiveness import run_augur

STEPS = ["io", "assign", "preprocess", "perturb", "augur", "enrich", "report"]
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
    "gseapy",
]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Single-sample Perturb-seq tertiary analysis")
    parser.add_argument("--config", required=True)
    parser.add_argument("--steps", default="all", help="all, or a comma-separated subset")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    config = load_config(args.config)
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
        elif step == "perturb":
            state["adata"] = _require(state, folder, "preprocess")
            state["adata"] = _step_perturb(folder, config, state["adata"])
        elif step == "augur":
            state["adata"] = _require(state, folder, "perturb")
            scores = run_augur(state["adata"], config)
            scores.to_csv(folder / "augur_scores.tsv", sep="\t", index=False)
            _write_h5ad(state["adata"], folder / "checkpoints" / "augur.h5ad")
        elif step == "enrich":
            effects = _read_table(folder / "effect_main.tsv")
            enrichment, reason = enrich_perturbations(effects, config)
            enrichment.to_csv(folder / "enrichment_results.tsv", sep="\t", index=False)
            (folder / "tables" / "enrichment_status.txt").write_text(reason + "\n")
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
            guides_per_cell=int(config["labels"].get("guides_per_cell", 1)),
            min_guide_umi=None if min_guide_umi is None else float(min_guide_umi),
            mito_prefix=mito_prefix(config),
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
    for column in ("total_counts", "n_genes_by_counts", "mito_frac"):
        calls[column] = pd.to_numeric(calls[column], errors="coerce")
    thresholds = {
        "min_umi": float(config["qc"]["min_umi"]),
        "max_umi": float(config["qc"]["max_umi"]),
        "min_genes": float(config["qc"]["min_genes"]),
        "max_mito_frac": float(config["qc"]["max_mito_frac"]),
    }
    (folder / "tables" / "qc_thresholds.json").write_text(json.dumps(thresholds, indent=2) + "\n")
    updated = apply_eligibility(
        calls,
        thresholds,
        min_cells_per_target=int(config["qc"]["min_cells_per_target"]),
        guides_per_cell=int(config["labels"].get("guides_per_cell", 1)),
    )
    updated.insert(0, "sample_id", config["experiment"]["sample_id"])
    spec = control_spec(config)
    temporary = pd.Series(False, index=updated.index)
    if spec["temporary"]:
        temporary = (
            updated["assignment_status"].astype(str).eq(STATUS_NO_GUIDE)
            & updated["n_guides"].astype(int).eq(0)
            & updated["qc_pass"].astype(bool)
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
    ledger[
        [
            "cell_barcode",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mito",
            "exclusion_reason_primary",
            "perturbation_status",
        ]
    ].to_csv(_qc_tables(folder) / "qc_metrics.tsv", sep="\t", index=False)
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
    control_ids = set(updated.loc[updated["temporary_control"].astype(bool), "cell_barcode"].astype(str))
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


def _step_perturb(folder: Path, config: dict[str, Any], adata):
    mixscape_on = config.get("mixscape", {}).get("enabled", True)
    comparisons = []
    spec = control_spec(config)
    if spec["temporary"]:
        print(f"  control: {spec['note']}", flush=True)
    if mixscape_on:
        if not config.get("mixscape", {}).get("cell_type_key") and not config["perturb"].get(
            "cell_type_key"
        ):
            print("  Mixscape scope: ALL cells, not stratified", flush=True)
        try:
            adata = run_mixscape(adata, config)
            adata.obs["mixscape_class_raw"] = adata.obs["mixscape_class"].astype(str)
            if "mixscape_class_p" in adata.obs:
                adata.obs["mixscape_probability"] = adata.obs["mixscape_class_p"]
            try:
                adata = run_mixscale(adata, config)
            except Exception as mixscale_exc:
                print(f"  Mixscale failed: {mixscale_exc}", flush=True)
                adata.uns["mixscale_status"] = f"failed:{mixscale_exc}"
            _write_de_gate(folder, adata)
        except Exception as exc:
            print(f"  Mixscape failed: {exc}", flush=True)
            adata.obs["mixscape_status"] = "not_evaluated"
            adata.obs["mixscape_class"] = "not_evaluated"
            comparisons.append(
                {
                    "branch": "sensitivity",
                    "cell_type": "ALL",
                    "target_gene": "",
                    "module": "mixscape",
                    "status": "skipped",
                    "reason": "mixscape_failed",
                    "n_cells": 0,
                    "n_guides": 0,
                    "biological_sample_count": 1,
                }
            )
            _write_de_gate(folder, adata)
    else:
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
    fractions = mixscape_fraction(adata)
    fractions.to_csv(_dge_tables(folder) / "mixscape_summary.tsv", sep="\t", index=False)
    main_effects, main_guides, main_status = descriptive_effects(adata, config, branch="main")
    sens_effects, sens_guides, sens_status = descriptive_effects(
        adata, config, branch="sensitivity"
    )
    main_effects.to_csv(folder / "effect_main.tsv", sep="\t", index=False)
    sens_effects.to_csv(folder / "effect_sensitivity.tsv", sep="\t", index=False)
    guides = pd.concat([main_guides, sens_guides], ignore_index=True)
    guides.to_csv(folder / "guide_effects.tsv", sep="\t", index=False)
    direction = pd.concat(
        [
            target_direction_table(adata, config, branch="main"),
            target_direction_table(adata, config, branch="sensitivity"),
        ],
        ignore_index=True,
    )
    direction.to_csv(folder / "target_expression_effect.tsv", sep="\t", index=False)
    de_main, de_main_status = tutorial_differential_expression(adata, config, branch="main")
    de_sens, de_sens_status = tutorial_differential_expression(adata, config, branch="sensitivity")
    for branch_name, tables in (("main", de_main), ("sensitivity", de_sens)):
        for method, frame in tables.items():
            frame.to_csv(
                _dge_tables(folder) / f"de_{method}_{branch_name}.tsv", sep="\t", index=False
            )
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
    plot_data = folder / "plot_data"
    if "min_guide_replicates" not in config.get("perturb", {}):
        plot_data.mkdir(exist_ok=True)
        (plot_data / "perturbation_distance_heatmap.skipped.txt").write_text(
            "min_guide_replicates_not_set\n"
        )
        return
    distance, embedding = perturbation_distance(pseudobulk, config)
    if distance.empty:
        plot_data.mkdir(exist_ok=True)
        (plot_data / "perturbation_distance_heatmap.skipped.txt").write_text(
            "insufficient_cells\n"
        )
        return
    distance.to_csv(folder / "perturbation_distance.tsv", sep="\t")
    embedding.to_csv(folder / "perturbation_embedding.tsv", sep="\t", index=False)


def _score_cell_cycle(adata, config: dict[str, Any]) -> None:
    """Score human cell cycle on the log-normalized matrix. Do not alter ``X``.

    Mouse has no gene set in this code. Do not reuse the human Seurat lists.
    """
    from perturbseq.plotting import score_human_cell_cycle_after_log1p

    organism = str(config["experiment"].get("organism", ""))
    if organism != "human":
        note = "代码只有人类 Seurat 细胞周期基因，没有小鼠基因集，因此跳过，不编造打分。"
        adata.uns["cell_cycle_status"] = "skipped_no_mouse_cell_cycle_genes"
        adata.uns["cell_cycle_note"] = note
        print(f"  cell cycle skipped: {note}", flush=True)
        return
    try:
        status = score_human_cell_cycle_after_log1p(
            adata,
            organism=organism,
            random_state=int(config["runtime"]["seed"]),
        )
    except ValueError as exc:
        print(f"  cell cycle skipped: {exc}", flush=True)
        adata.uns["cell_cycle_status"] = "skipped"
        return
    adata.uns["cell_cycle_status"] = status
    print(f"  cell cycle: {status}", flush=True)


def _qc_tables(folder: Path) -> Path:
    path = folder / "qc" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _dge_tables(folder: Path) -> Path:
    path = folder / "DGE" / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_de_gate(folder: Path, adata) -> None:
    gate = adata.uns.get("mixscape_de_gate")
    if not isinstance(gate, pd.DataFrame):
        return
    gate.to_csv(_dge_tables(folder) / "mixscape_de_gate.tsv", sep="\t", index=False)


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


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


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
        "perturb": "preprocess",
        "augur": "perturb",
        "enrich": "perturb",
        "report": "augur",
    }
    if step not in previous:
        return None
    path = folder / "checkpoints" / f"{previous[step]}.h5ad"
    fallback = {
        "report": ["augur.h5ad", "perturb.h5ad", "preprocess.h5ad", "assign.h5ad"],
        "enrich": ["perturb.h5ad"],
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
