"""Run the single-sample tertiary analysis one step at a time."""

from __future__ import annotations

import argparse
import json
import os
import platform
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from perturbseq import __version__
from perturbseq.assign import (
    DEFAULT_NT_CONTROL_ARG,
    LOW_LIBRARY_NO_GUIDE,
    STATUS_NO_GUIDE,
    STATUS_NT,
    STATUS_RETAINED,
    STATUS_TARGET,
    GuideReferenceError,
    apply_eligibility,
    bind_nt_control_from_reference,
    capture_library_correlation,
    cells_per_guide,
    cells_per_guide_name,
    cells_per_target,
    derive_qc_thresholds,
    filter_stage_counts,
    multi_guide_composition,
    reference_guide_names,
    reference_target_genes,
    validate_guide_targets,
    write_qc_figures,
)
from perturbseq.config import (
    TEMPORARY_CONTROL_WARNING,
    ConfigError,
    apply_customer_arguments,
    apply_thread_budget,
    control_spec,
    mito_prefix,
    nt_labels_from_config,
    output_dir,
    read_config,
    resolve_control_spec,
    resolve_input_paths,
    resolve_n_threads,
    restore_control_spec,
    software_config_path,
    validate_config,
    verify_species_consistency,
)
from perturbseq.io import (
    InputMismatch,
    load_sample,
    matrix_only_as_no_guide,
    read_features,
    read_reference,
)
from perturbseq.mixscale import run_mixscale, write_mixscale_figures
from perturbseq.mixscape import run_mixscape, write_mixscape_figures
from perturbseq.outputs import cleanup_lean_outputs
from perturbseq.preprocess import preprocess, write_preprocess_figures
from perturbseq.sctype import sctype_tissue_choices

STEPS = [
    "io",
    "assign",
    "preprocess",
    "mixscape",
    "mixscale",
]
# After preprocess these only read preprocess AnnData and write disjoint PNG dirs.
PARALLEL_ANALYSIS_STEPS = ("mixscape", "mixscale")
# Fixed run order: always Mixscape + Mixscale after preprocess (KO / i / a).
DEFAULT_STEPS = ["io", "assign", "preprocess", "mixscape", "mixscale"]
PACKAGES = [
    "perturbseq",
    "scanpy",
    "anndata",
    "pertpy",
    "pandas",
    "numpy",
    "scipy",
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
    parser.add_argument(
        "--organism", required=True, help="human 或 mouse（写入 experiment.organism）"
    )
    parser.add_argument("--input-dir", required=True, help="五个 DRAGEN 输入文件所在目录")
    parser.add_argument("--output-dir", required=True, help="结果输出目录")
    parser.add_argument("--guide-separator", required=True, help="多条 guide 之间的分隔符，例如 |")
    parser.add_argument(
        "--nt-control",
        default=DEFAULT_NT_CONTROL_ARG,
        help=(
            "在 *.scRNA.feature_barcode_reference.csv 的 "
            "id/name/target_gene_name/target_gene_id 上按正则 re.search 匹配设计 NT（忽略大小写）"
            f"（默认 {DEFAULT_NT_CONTROL_ARG}，忽略大小写；可用 / 或 , 自定义）。"
            "默认或自定义关键词均未匹配到时自动用 zero_guide；也可显式填 zero_guide"
        ),
    )
    parser.add_argument(
        "--crispr-mode",
        default="KO",
        choices=("KO", "i", "a"),
        help="KO、i 或 a，缺省 KO（映射 CRISPR_KO / CRISPRi / CRISPRa；预期方向 not_required / down / up；图标签 KO / KD / OE）",
    )
    parser.add_argument(
        "--cell-type-key", default=None, help="obs 里已有的细胞类型列。不填则全程 ALL"
    )
    parser.add_argument(
        "--max-mito-frac",
        type=float,
        default=None,
        help="线粒体 reads 比例上限，0 到 1。细胞 ≥ 该值时剔除。默认读配置（0.15 = 15%%）。",
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
    tissue_choices = sctype_tissue_choices()
    parser.add_argument(
        "--sctype-tissue",
        default=None,
        choices=tissue_choices,
        metavar="TISSUE",
        help=(
            "启用 ScType 并指定组织（choices 来自离线 ScTypeDB）。"
            "未传则跳过 ScType。含空格时请加引号，例如 --sctype-tissue 'Immune system'。"
            f"可选: {', '.join(tissue_choices)}"
        ),
    )
    parser.add_argument(
        "--sctype-db-path",
        default="data/sctype/ScTypeDB_full.tsv",
        help="ScType marker DB 路径（tsv/csv/xlsx）。默认 data/sctype/ScTypeDB_full.tsv。",
    )
    parser.add_argument(
        "--sctype-low-score-fraction",
        type=float,
        default=None,
        help=(
            "ScType 低分簇阈值：顶分 < n_cells * fraction → Unknown。"
            "仅在指定 --sctype-tissue 时生效；默认 0.25。"
        ),
    )
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
            "crispr_mode_description": None,
            "expected_target_direction": None,
            "cell_type_key": args.cell_type_key,
            "max_mito_frac": args.max_mito_frac,
            "min_genes": args.min_genes,
            "max_genes": args.max_genes,
            "sctype_tissue": args.sctype_tissue,
            "sctype_db_path": args.sctype_db_path,
            "sctype_low_score_fraction": args.sctype_low_score_fraction,
        },
    )
    validate_config(config)
    steps = list(DEFAULT_STEPS)
    mode = str(config["experiment"]["crispr_mode"])
    print(
        f"  steps → {', '.join(steps)} (crispr_mode={mode}; always mixscape ‖ mixscale)",
        flush=True,
    )
    run(config, steps)


def run(config: dict[str, Any], steps: list[str]) -> None:
    folder = output_dir(config)
    folder.mkdir(parents=True, exist_ok=True)
    if steps and steps[0] != "io":
        raise SystemExit(
            "不支持 checkpoints / --resume："
            f"请从 io 起连续运行（当前首步为 {steps[0]}）。"
        )
    n_threads = apply_thread_budget(config)
    print(f"  runtime.n_threads={n_threads} (host/3 default; BLAS/OpenMP/scanpy)", flush=True)
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
    state: dict[str, Any] = {"adata": None, "preprocess_adata": None}
    started = time.time()
    # Provisional until assign resolves from the ledger (may fall back to zero_guide).
    info["temporary_control_warning"] = _publish_temporary_control_warning(folder, config)
    _write_info(folder, info, started)
    remaining = list(steps)
    while remaining:
        step = remaining[0]
        if step in PARALLEL_ANALYSIS_STEPS:
            batch = [name for name in PARALLEL_ANALYSIS_STEPS if name in remaining]
            for name in batch:
                remaining.remove(name)
            print(
                f"[parallel] {', '.join(batch)} (in-memory preprocess; temp IPC only)",
                flush=True,
            )
            step_started = time.time()
            preprocess_adata = _ensure_preprocess_adata(state, config)
            timings, state["adata"] = _run_parallel_analysis(
                folder, config, batch, preprocess_adata
            )
            wall = round(time.time() - step_started, 3)
            for name in batch:
                info.setdefault("steps", {})[name] = {
                    "seconds": timings.get(name, wall),
                    "parallel_wall_seconds": wall,
                    "parallel_batch": list(batch),
                }
            info.setdefault("steps", {})["_parallel_analysis"] = {
                "steps": list(batch),
                "wall_seconds": wall,
            }
            _write_info(folder, info, started)
            continue
        remaining.pop(0)
        print(f"[{step}]", flush=True)
        step_started = time.time()
        if step == "io":
            state["adata"] = _step_io(folder, config)
            validation = json.loads((folder / "input_validation.json").read_text())
            info["input_validation"] = validation
            info["n_matrix_only_as_no_guide"] = validation.get("n_matrix_only_as_no_guide", 0)
        elif step == "assign":
            state["adata"] = _require(state, "io")
            state["adata"] = _step_assign(folder, config, state["adata"])
            info.update(_control_summary(config, state["adata"]))
            info["temporary_control_warning"] = _publish_temporary_control_warning(folder, config)
        elif step == "preprocess":
            state["adata"] = _require(state, "assign")
            state["adata"] = preprocess(state["adata"], config)
            _score_cell_cycle(state["adata"], config)
            info["cell_cycle_status"] = state["adata"].uns.get("cell_cycle_status", "")
            info["cell_cycle_note"] = state["adata"].uns.get("cell_cycle_note", "")
            write_preprocess_figures(folder, config, state["adata"])
            state["preprocess_adata"] = state["adata"]
        else:
            raise SystemExit(f"未知步骤: {step}")
        info.setdefault("steps", {})[step] = {"seconds": round(time.time() - step_started, 3)}
        _write_info(folder, info, started)
    removed = cleanup_lean_outputs(folder, protect=_lean_protect_paths(config, folder))
    print(f"Wrote {folder} (lean cleanup removed {removed} intermediate files)", flush=True)


def _lean_protect_paths(config: dict[str, Any], folder: Path) -> list[Path]:
    """Keep DRAGEN inputs when the output dir is the sample input dir."""
    paths = list(resolve_input_paths(config).values())
    sample_id = str(config["experiment"]["sample_id"])
    input_dir = Path(config["input"].get("input_dir", "."))
    if input_dir.resolve() == Path(folder).resolve():
        paths.extend(input_dir.glob(f"{sample_id}.scRNA.*"))
    return paths


def _ensure_preprocess_adata(state: dict[str, Any], config: dict[str, Any]):
    """Return preprocess AnnData from the in-memory run state."""
    cached = state.get("preprocess_adata")
    if cached is not None:
        restore_control_spec(config, cached)
        return cached
    adata = _require(state, "preprocess")
    restore_control_spec(config, adata)
    state["preprocess_adata"] = adata
    return adata


def _run_parallel_analysis(
    folder: Path,
    config: dict[str, Any],
    batch: list[str],
    preprocess_adata,
) -> tuple[dict[str, float], Any]:
    """Run mixscape / mixscale; IPC via a temporary directory (not output checkpoints)."""
    public_config = {key: value for key, value in config.items() if not str(key).startswith("_")}
    with tempfile.TemporaryDirectory(prefix="perturbseq_parallel_") as tmp:
        tmp_dir = Path(tmp)
        preprocess_path = tmp_dir / "preprocess.h5ad"
        print("  writing temporary preprocess for workers…", flush=True)
        _write_h5ad(preprocess_adata, preprocess_path)
        product_paths = {step: tmp_dir / f"{step}.h5ad" for step in batch}
        if len(batch) == 1:
            result = _parallel_step_worker(
                (
                    batch[0],
                    str(folder),
                    public_config,
                    str(preprocess_path),
                    str(product_paths[batch[0]]),
                )
            )
            timings = {result["step"]: float(result["seconds"])}
        else:
            timings = {}
            errors: list[str] = []
            max_workers = min(len(batch), resolve_n_threads(config))
            with ProcessPoolExecutor(max_workers=max_workers) as pool:
                futures = {
                    pool.submit(
                        _parallel_step_worker,
                        (
                            step,
                            str(folder),
                            public_config,
                            str(preprocess_path),
                            str(product_paths[step]),
                        ),
                    ): step
                    for step in batch
                }
                for future in as_completed(futures):
                    step = futures[future]
                    try:
                        payload = future.result()
                    except Exception as exc:  # noqa: BLE001
                        errors.append(f"{step}: {exc}")
                        continue
                    timings[str(payload["step"])] = float(payload["seconds"])
            if errors:
                raise SystemExit("并行分析失败: " + "; ".join(errors))
        merged = _merge_post_preprocess_products(preprocess_adata, product_paths)
        # Mixscale worker runs without Mixscape labels; redraw when both are present.
        if (
            "mixscale" in batch
            and "mixscape_status" in merged.obs
            and "mixscale_score" in merged.obs
            and "mixscale_de_genes" in merged.uns
        ):
            print("  redrawing Mixscale figures with Mixscape NT/NP/KO labels…", flush=True)
            write_mixscale_figures(folder, config, merged)
    return timings, merged


def _parallel_step_worker(
    payload: tuple[str, str, dict[str, Any], str, str],
) -> dict[str, Any]:
    """Process entrypoint: one analysis step from a temporary preprocess.h5ad."""
    step, folder_str, config, preprocess_path, product_path = payload
    os.environ.setdefault("MPLBACKEND", "Agg")
    apply_thread_budget(config)
    folder = Path(folder_str)
    seed = int(config["runtime"]["seed"])
    np.random.seed(seed)
    print(f"[{step}] start (worker)", flush=True)
    started = time.time()
    adata = sc.read_h5ad(preprocess_path)
    restore_control_spec(config, adata)
    if step == "mixscape":
        adata = _step_mixscape(folder, config, adata)
        _write_h5ad(adata, Path(product_path))
    elif step == "mixscale":
        adata = _step_mixscale(folder, config, adata)
        _write_h5ad(adata, Path(product_path))
    else:
        raise ValueError(f"unsupported parallel step: {step}")
    seconds = round(time.time() - started, 3)
    print(f"[{step}] done ({seconds}s)", flush=True)
    return {"step": step, "seconds": seconds}


def _merge_post_preprocess_products(preprocess_adata, product_paths: dict[str, Path]):
    """Merge Mixscape / Mixscale products onto the preprocess AnnData."""
    import gc

    adata = preprocess_adata.copy()
    mixscape_path = Path(product_paths["mixscape"]) if "mixscape" in product_paths else None
    if mixscape_path is not None and mixscape_path.exists():
        mix = sc.read_h5ad(mixscape_path, backed="r")
        for column in (
            "mixscape_status",
            "mixscape_class",
            "mixscape_class_raw",
            "mixscape_class_global",
            "mixscape_class_p",
            "mixscape_probability",
            "mixscape_reason",
        ):
            if column in mix.obs.columns:
                adata.obs[column] = mix.obs[column].to_numpy()
        if "X_pert" in mix.layers:
            adata.layers["X_pert"] = mix.layers["X_pert"][...]
        if "X_umap_mixscape" in mix.obsm:
            adata.obsm["X_umap_mixscape"] = np.asarray(mix.obsm["X_umap_mixscape"])
        for key in list(mix.uns.keys()):
            if str(key).startswith("mixscape"):
                adata.uns[key] = mix.uns[key]
        mix.file.close()
        del mix
        gc.collect()
    mixscale_path = Path(product_paths["mixscale"]) if "mixscale" in product_paths else None
    if mixscale_path is not None and mixscale_path.exists():
        from perturbseq.mixscale import MIXSCALE_LAYER

        scale = sc.read_h5ad(mixscale_path, backed="r")
        if "mixscale_score" in scale.obs.columns:
            adata.obs["mixscale_score"] = scale.obs["mixscale_score"].to_numpy()
        if MIXSCALE_LAYER in scale.layers:
            adata.layers[MIXSCALE_LAYER] = scale.layers[MIXSCALE_LAYER][...]
        for key in list(scale.uns.keys()):
            if str(key).startswith("mixscale"):
                adata.uns[key] = scale.uns[key]
        scale.file.close()
        del scale
        gc.collect()
    return adata


def _step_io(folder: Path, config: dict[str, Any]):
    paths = resolve_input_paths(config)
    reference = read_reference(paths["feature_reference"])
    bind_nt_control_from_reference(config, reference)
    if str(config["labels"].get("nt_mode", "")) == "zero_guide":
        if config["labels"].get("nt_unmatched_fallback"):
            keys = config["labels"].get("nt_keywords") or []
            print(
                "  NT control: zero_guide (temporary; no reference match for "
                f"keywords {list(keys)})",
                flush=True,
            )
        else:
            print("  NT control: zero_guide (temporary)", flush=True)
    else:
        print(
            "  NT control: designed "
            f"{config['labels']['nt_control']!r} "
            f"({len(config['labels'].get('nt_labels') or [])} matched labels)",
            flush=True,
        )
    _validate_guide_reference(folder, config)
    try:
        adata, annotated, checks = load_sample(
            paths,
            nt_control=config["labels"]["nt_control"],
            nt_labels=nt_labels_from_config(config),
            feature_type_gex=config["labels"].get("feature_type_gex", "Gene Expression"),
            feature_type_crispr=config["labels"]["feature_type_crispr"],
            multi_guide_policy=config["labels"]["multi_guide_policy"],
            targets_per_cell=int(config["labels"].get("targets_per_cell", 1)),
            mito_prefix=mito_prefix(config),
            guide_separator=str(config["input"]["guide_separator"]),
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
    calls = _calls(annotated)
    # Keep the annotated ledger on the in-memory AnnData for assign in the same run.
    adata.uns["cell_calls"] = calls
    adata.uns["io_checks"] = checks
    adata.uns["multi_guide_composition"] = multi_guide_composition(
        annotated,
        nt_control=config["labels"]["nt_control"],
        nt_labels=nt_labels_from_config(config),
        reference=reference,
        guide_separator=str(config["input"]["guide_separator"]),
    )
    return adata


def _load_cell_calls(adata) -> pd.DataFrame:
    stored = adata.uns.get("cell_calls") if adata is not None else None
    if not isinstance(stored, pd.DataFrame) or stored.empty:
        raise SystemExit("缺少 cell_calls：请从 io 起连续运行（账本在内存 AnnData.uns 中传递）。")
    calls = stored.copy()
    for column in (
        "total_counts",
        "n_genes_by_counts",
        "mito_frac",
        "guide_umi_total",
        "grna_capture_efficiency",
    ):
        if column in calls.columns:
            calls[column] = pd.to_numeric(calls[column], errors="coerce")
    return calls


def _step_assign(folder: Path, config: dict[str, Any], adata):
    calls = _load_cell_calls(adata)
    thresholds = derive_qc_thresholds(config["qc"])
    # low_library_no_guide stage retained in ledger; exclusion count stays 0.
    library_rule = {
        "applied": False,
        "reason": "disabled_mad_qc_removed",
        "min_total_counts": None,
        "source": "disabled",
        "note": "low_library_no_guide 不应用排除；账本阶段保留，新增排除数为 0。",
    }
    correlation = capture_library_correlation(calls)
    correlation["low_library_no_guide_applied"] = False
    correlation["low_library_no_guide_reason"] = library_rule["reason"]
    # Keep a float column (np.nan), not None — nested None/DataFrame break h5ad/JSON.
    correlation["low_library_min_total_counts"] = np.nan
    recorded = {
        **thresholds,
        "library_umi_filter": "none",
        "gene_bounds_source": "configured",
        "gene_filter": "drop_n_genes_le_min_or_ge_max",
        "mito_frac_source": "configured",
        "mito_filter": "drop_mito_frac_ge_max",
        "low_library_no_guide": {
            **library_rule,
            "min_total_counts": (
                np.nan
                if library_rule.get("min_total_counts") is None
                else library_rule["min_total_counts"]
            ),
        },
        # Nested DataFrame is not JSON-serializable (figure_context) and is fragile in h5ad.
        "grna_capture_library_correlation": _uns_scalar_row(correlation),
    }
    print(f"  low_library_no_guide skipped: {library_rule['reason']}", flush=True)
    updated = apply_eligibility(
        calls,
        thresholds,
        min_cells_per_target=int(config["qc"]["min_cells_per_target"]),
        targets_per_cell=int(config["labels"].get("targets_per_cell", 1)),
        low_library_rule=library_rule,
    )
    updated.insert(0, "sample_id", config["experiment"]["sample_id"])
    spec = resolve_control_spec(config, updated)
    if spec["temporary"] and spec.get("fallback_from"):
        print(
            f"  control fallback: --nt-control={spec['fallback_from']!r} "
            "has no usable NT cells → temporary zero_guide",
            flush=True,
        )
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
    reference = read_reference(resolve_input_paths(config)["feature_reference"])
    targets = cells_per_target(
        updated,
        int(config["qc"]["min_cells_per_guide"]),
        all_targets=reference_target_genes(
            reference,
            nt_control=str(config["labels"].get("nt_control", "") or ""),
            nt_labels=nt_labels_from_config(config),
        ),
    )
    guide_name_counts = cells_per_guide_name(
        updated,
        all_names=reference_guide_names(reference),
    )
    control_ids = set(
        updated.loc[updated["temporary_control"].astype(bool), "cell_barcode"].astype(str)
    )
    # Preprocess / Leiden / UMAP / SingleR: every QC-pass cell already in the
    # AnnData (io candidates = single assignment ∪ zero-guide). Mixscape /
    # Mixscale still restrict to analysis_eligible ∪ contrast via obs flags.
    qc_ids = set(updated.loc[updated["qc_pass"].astype(bool), "cell_barcode"].astype(str))
    adata = adata[adata.obs_names.astype(str).isin(qc_ids)].copy()
    if "cell_calls" in adata.uns:
        del adata.uns["cell_calls"]
    by_barcode = updated.set_index(updated["cell_barcode"].astype(str))
    aligned = by_barcode.reindex(adata.obs_names.astype(str))
    adata.obs["qc_pass"] = aligned["qc_pass"].astype(bool).to_numpy()
    adata.obs["analysis_eligible"] = aligned["analysis_eligible"].astype(bool).to_numpy()
    adata.obs["exclusion_reason_primary"] = (
        aligned["exclusion_reason_primary"].astype(str).to_numpy()
    )
    adata.obs["temporary_control"] = aligned["temporary_control"].astype(bool).to_numpy()
    if control_ids:
        hit = np.isin(adata.obs_names.astype(str), list(control_ids))
        for column in ("target_gene", "perturbation_status"):
            values = adata.obs[column].astype(str).to_numpy(copy=True)
            values[hit] = spec["control_label"]
            adata.obs[column] = values
    adata.uns["analysis_scope"] = "single_sample_exploratory"
    adata.uns["embedding_cell_scope"] = "qc_pass"
    adata.uns["n_embedding_cells"] = int(adata.n_obs)
    adata.uns["n_analysis_eligible"] = analysis_n
    adata.uns["config"] = config
    adata.uns["control_source"] = spec["control_source"]
    adata.uns["control_label"] = spec["control_label"]
    adata.uns["control_temporary"] = bool(spec["temporary"])
    adata.uns["control_phrase"] = spec["control_phrase"]
    adata.uns["control_note"] = spec["note"]
    adata.uns["control_fallback_from"] = str(spec.get("fallback_from", ""))
    adata.uns["n_temporary_control_cells"] = int(len(control_ids)) if spec["temporary"] else 0
    adata.uns["qc_thresholds"] = recorded
    adata.uns["cells_per_target"] = targets
    adata.uns["cells_per_guide_name"] = guide_name_counts
    adata.uns["filter_stage_counts"] = filter_stage_counts(updated)
    adata.uns["cells_per_guide"] = cells_per_guide(updated)
    if spec["temporary"]:
        print(
            f"  temporary control {spec['control_label']}: {len(control_ids):,} zero-guide cells",
            flush=True,
        )
    else:
        print(
            f"  designed NT control {spec['control_label']}: "
            f"{int((updated['analysis_eligible'].astype(bool) & (updated['perturbation_status'] == STATUS_NT)).sum()):,} cells",
            flush=True,
        )
    _write_resolved(folder, config)
    write_qc_figures(
        folder,
        config,
        adata,
        cell_status=ledger,
        targets=targets,
        guide_names=guide_name_counts,
    )
    print(
        f"  embedding cells (QC pass): {adata.n_obs:,}; "
        f"analysis_eligible (perturbation): {analysis_n:,}",
        flush=True,
    )
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
    """Independent Mixscape step: signature on unscaled log1p, then pertpy mixscape."""
    mixscape_on = config.get("mixscape", {}).get("enabled", True)
    spec = control_spec(config)
    if spec["temporary"]:
        print(f"  temporary control warning: {spec['note']}", flush=True)
    if not mixscape_on:
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
        adata.obs["mixscape_class_raw"] = "not_evaluated"
        write_mixscape_figures(folder, config, adata)
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
    except Exception as exc:
        print(f"  Mixscape failed: {exc}", flush=True)
        adata.obs["mixscape_status"] = "not_evaluated"
        adata.obs["mixscape_class"] = "not_evaluated"
        adata.obs["mixscape_class_raw"] = "not_evaluated"
        adata.uns["mixscape_skip_reason"] = f"mixscape_failed:{exc}"
    write_mixscape_figures(folder, config, adata)
    return adata


def _step_mixscale(folder: Path, config: dict[str, Any], adata):
    """Independent Mixscale step: own signature and mixscale_score. Does not read Mixscape."""
    spec = control_spec(config)
    if spec["temporary"]:
        print(f"  temporary control warning: {spec['note']}", flush=True)
    try:
        adata = run_mixscale(adata, config)
    except Exception as mixscale_exc:
        print(f"  Mixscale failed: {mixscale_exc}", flush=True)
        adata.uns["mixscale_status"] = f"failed:{mixscale_exc}"
    write_mixscale_figures(folder, config, adata)
    return adata


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
        "skipped_unsupported_organism": "experiment.organism 不是 human 或 mouse。",
        "skipped_missing_log1p": "尚未 log1p，跳过打分。",
        "skipped_missing_cell_cycle_genes": "矩阵里几乎找不到该物种的细胞周期基因。",
        "skipped_scib_unavailable": "未安装 scib，无法读取细胞周期基因集。",
    }
    adata.uns["cell_cycle_note"] = notes.get(status, "")
    print(f"  cell cycle: {status}", flush=True)


def _write_resolved(folder: Path, config: dict[str, Any]) -> None:
    public = {key: value for key, value in config.items() if not str(key).startswith("_")}
    (folder / "resolved_config.yaml").write_text(
        yaml.safe_dump(public, sort_keys=False, allow_unicode=True)
    )


def _validate_guide_reference(folder: Path, config: dict[str, Any]) -> None:
    paths = resolve_input_paths(config)
    reference = read_reference(paths["feature_reference"])
    try:
        validate_guide_targets(
            reference,
            nt_control=config["labels"]["nt_control"],
            nt_labels=nt_labels_from_config(config),
        )
    except GuideReferenceError as exc:
        _write_guide_conflicts(folder, exc)
        raise SystemExit(str(exc)) from exc


def _write_guide_conflicts(folder: Path, exc: GuideReferenceError) -> None:
    (folder / "tables").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"guide_id": exc.guide_ids, "reason": exc.reason}).to_csv(
        folder / "tables" / "guide_reference_conflicts.tsv", sep="\t", index=False
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
        "control_temporary": bool(spec["temporary"]),
        "control_phrase": spec["control_phrase"],
        "control_note": spec["note"],
        "control_fallback_from": str(spec.get("fallback_from", "")),
    }
    if adata is not None:
        summary["n_temporary_control_cells"] = int(adata.uns.get("n_temporary_control_cells", 0))
    return summary


def _publish_temporary_control_warning(folder: Path, config: dict[str, Any]) -> str:
    """Record the zero-guide warning in the log and on disk.

    A designed NT control removes that warning and the temporary-control QC figure.
    After assign, the note may include a data-level fallback from a missing NT label.
    """
    path = folder / "tables" / "temporary_control_warning.txt"
    spec = control_spec(config)
    if not spec["temporary"]:
        if path.is_file():
            path.unlink()
        from perturbseq.plotting import remove_temporary_control_qc

        remove_temporary_control_qc(folder)
        return ""
    text = str(spec.get("note") or TEMPORARY_CONTROL_WARNING).strip()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n")
    print(f"  temporary control warning: {text}", flush=True)
    return text


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


def _require(state: dict[str, Any], previous: str):
    if state["adata"] is not None:
        return state["adata"]
    raise SystemExit(
        f"缺少上一步「{previous}」的内存结果：请从 io 起连续运行（不支持 checkpoints / --resume）。"
    )


def _calls(annotated: pd.DataFrame) -> pd.DataFrame:
    drop = [
        column for column in ("guide_names", "guide_transcripts") if column in annotated.columns
    ]
    return annotated.drop(columns=drop)


def _uns_scalar_row(frame: pd.DataFrame) -> dict[str, Any]:
    """One-row QC table as a plain dict safe for nested uns / JSON / h5ad."""
    if frame.empty:
        return {}
    out: dict[str, Any] = {}
    for key, value in frame.iloc[0].items():
        name = str(key)
        if isinstance(value, (bool, np.bool_)):
            out[name] = bool(value)
        elif isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
            out[name] = int(value)
        elif isinstance(value, (float, np.floating)):
            out[name] = float(value)
        elif isinstance(value, np.generic):
            out[name] = value.item()
        else:
            try:
                missing = bool(pd.isna(value))
            except (TypeError, ValueError):
                missing = value is None
            if missing:
                out[name] = np.nan
            elif isinstance(value, str):
                out[name] = value
            else:
                out[name] = str(value)
    return out


def _write_h5ad(adata, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # h5py vlen strings reject pd.NA / None in object columns (e.g. target_gene).
    adata.obs = _h5_safe_strings(adata.obs)
    adata.var = _h5_safe_strings(adata.var)
    adata.uns.update({key: _h5_safe_uns(value) for key, value in list(adata.uns.items())})
    adata.write_h5ad(path, compression="lzf")


def _h5_safe_uns(value: Any) -> Any:
    """Recursively make uns values writable (nested DataFrames / object NA)."""
    if isinstance(value, pd.DataFrame):
        return _h5_safe_strings(value)
    if isinstance(value, dict):
        return {str(key): _h5_safe_uns(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_h5_safe_uns(item) for item in value]
    if isinstance(value, tuple):
        return [_h5_safe_uns(item) for item in value]
    return value


def _h5_safe_strings(frame: pd.DataFrame) -> pd.DataFrame:
    """Fill missing values in string/object columns so AnnData can write h5ad."""
    out = frame.copy()
    for column in out.columns:
        series = out[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            if series.isna().any():
                if "" not in series.cat.categories:
                    series = series.cat.add_categories([""])
                out[column] = series.fillna("")
            continue
        if series.dtype != object and not pd.api.types.is_string_dtype(series):
            continue
        non_null = series[series.notna()]
        if len(non_null) and isinstance(non_null.iloc[0], (list, tuple, dict, set, np.ndarray)):
            continue
        out[column] = series.astype("string").fillna("").astype(str)
    return out


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
