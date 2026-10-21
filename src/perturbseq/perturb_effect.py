"""Descriptive CPM effects.

Mixscape and Mixscale live in ``perturbseq.mixscape`` and ``perturbseq.mixscale``.
Cell-level E-test lives in ``perturbseq.distance_tests``.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from perturbseq.assign import STATUS_NT, STATUS_TARGET, control_cell_mask
from perturbseq.config import control_spec, expected_expression
from perturbseq.de_gate import stamp_control_label
from perturbseq.io import sum_by_group
from perturbseq.mixscape import PERTURBED

DESCRIPTIVE_UNCERTAINTY = "descriptive_estimate_uncertainty_not_significance_test"


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


def _missing_bootstrap(n_genes: int, status: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    lower, upper, status = _percentile_bounds(
        replicates, np.asarray(point, dtype=np.float64), level
    )
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
    stamp_control_label(selected.obs, selected_control, "target_gene", control)
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
    control_positions = np.flatnonzero(
        selected.obs["target_gene"].astype(str).to_numpy() == control
    )
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
                    "control_source": spec["control_source"],
                    "control_label": control,
                    "temporary_control": bool(spec["temporary"]),
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
        "control_source",
        "control_label",
        "temporary_control",
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


def write_descriptive_figures(
    output_dir: Path,
    config: dict[str, Any],
    *,
    effects: pd.DataFrame | None = None,
    direction: pd.DataFrame | None = None,
    guides: pd.DataFrame | None = None,
) -> None:
    """Write descriptive DE PNGs after the perturb step (E-test is separate)."""
    from perturbseq.plotting import (
        draw_descriptive_heatmap,
        draw_guide_consistency,
        draw_target_expression,
        open_figure_store,
    )

    output_dir = Path(output_dir)
    store = open_figure_store(output_dir, config, None)
    # Cell-level Wilcoxon is computed in the dge step (no separate volcano PNG).
    store.skip("volcano", "wilcoxon_in_dge_step")
    store.skip("target_effect_forest", "wilcoxon_in_dge_step")
    if effects is None:
        effects = _read_first(
            output_dir / "DGE" / "tables" / "effect_main.tsv",
            output_dir / "DGE" / "descriptive" / "tables" / "effect_main.tsv",
            output_dir / "effect_main.tsv",
        )
    draw_descriptive_heatmap(store, effects)
    if direction is None:
        direction = _read_first(
            output_dir / "DGE" / "tables" / "target_expression_effect.tsv",
            output_dir / "DGE" / "descriptive" / "tables" / "target_expression_effect.tsv",
            output_dir / "target_expression_effect.tsv",
            output_dir / "tables" / "target_expression.tsv",
        )
    draw_target_expression(store, direction)
    if guides is None:
        guides = _read_first(
            output_dir / "DGE" / "tables" / "guide_effects.tsv",
            output_dir / "DGE" / "descriptive" / "tables" / "guide_effects.tsv",
            output_dir / "guide_effects.tsv",
        )
    draw_guide_consistency(store, guides)


def _read_first(*paths: Path) -> pd.DataFrame:
    for path in paths:
        if path.exists():
            return pd.read_csv(path, sep="\t")
    return pd.DataFrame()
