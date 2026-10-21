"""Cell-level differential expression via pertpy WilcoxonTest.

Independent of Mixscape / Mixscale. Outputs under ``wilcoxon/``. Aligns with
pertpy Simple statistical tests:

https://pertpy.scverse.org/en/stable/tutorials/notebooks/differential_gene_expression.html

Single-sample cell-level tests treat cells as observations; p-values are not
biological-replicate inference.

pertpy ``log_fc`` is ``log2(mean_group) - log2(mean_control)`` on the current
matrix (here log1p), with no pseudocount. Zero means yield ±inf / NaN; those
rows are set to NA after the call (p-values unchanged).
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbseq.assign import STATUS_TARGET, control_cell_mask
from perturbseq.config import control_spec
from perturbseq.de_gate import stamp_control_label

WILCOXON_COLUMNS = [
    "variable",
    "log_fc",
    "p_value",
    "statistic",
    "comparison",
    "adj_p_value",
    "target_gene",
    "cell_type",
    "control",
    "control_source",
    "temporary_control",
    "n_target_cells",
    "n_control_cells",
    "method",
    "inference_note",
]

INFERENCE_NOTE = (
    "pertpy WilcoxonTest (Mann-Whitney U) on log1p matrix; "
    "log_fc = log2(mean_log1p_group) - log2(mean_log1p_control) without "
    "pseudocount; non-finite log_fc set to NA; single-sample cell-level test, "
    "not biological-replicate DE"
)

LOG_FC_NOTE = (
    "pertpy formula on adata.X (log1p): log2(mean_g1)-log2(mean_g0); "
    "no pseudocount; non-finite values replaced with NA"
)


def run_wilcoxon_de(adata: AnnData, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare each target gene to the current control with WilcoxonTest.compare_groups.

    Uses ``adata.X`` (log1p after preprocess). Does not read Mixscape labels or
    ``X_pert``. Returns ``(results, status)``.
    """
    from pertpy.tools import WilcoxonTest

    from perturbseq.preprocess import assign_cell_type

    assign_cell_type(adata, config)
    spec = control_spec(config)
    control = str(spec["control_label"])
    control_source = str(spec["control_source"])
    temporary = bool(spec["temporary"])
    if temporary:
        stamp_control_label(
            adata.obs, control_cell_mask(adata.obs, spec), "target_gene", control
        )
    sample_id = str(config["experiment"]["sample_id"])
    min_cells = int(config.get("qc", {}).get("min_cells_per_target", 20))
    min_nt = int(config.get("statistics", {}).get("min_nt_cells", min_cells))
    results: list[pd.DataFrame] = []
    statuses: list[dict[str, Any]] = []

    for cell_type, subset in _analysis_populations(adata):
        control_mask = control_cell_mask(subset.obs, spec)
        n_control = int(control_mask.sum())
        if temporary and n_control == 0:
            statuses.append(
                _status_row(
                    sample_id,
                    cell_type,
                    "",
                    "skipped",
                    "no_zero_guide_cells",
                    0,
                    n_control,
                )
            )
            continue
        if n_control < min_nt:
            statuses.append(
                _status_row(
                    sample_id,
                    cell_type,
                    "",
                    "skipped",
                    "insufficient_control_cells",
                    0,
                    n_control,
                )
            )
            continue
        if "perturbation_status" in subset.obs:
            targeting = subset.obs["perturbation_status"].astype(str).eq(STATUS_TARGET)
        else:
            targeting = subset.obs["target_gene"].astype(str).ne(control)
        targets = sorted(
            {
                str(gene)
                for gene in subset.obs.loc[targeting, "target_gene"].astype(str)
                if gene and gene != control
            }
        )
        if not targets:
            statuses.append(
                _status_row(
                    sample_id,
                    cell_type,
                    "",
                    "skipped",
                    "no_target_cells",
                    0,
                    n_control,
                )
            )
            continue
        for gene in targets:
            gene_mask = subset.obs["target_gene"].astype(str).eq(gene)
            n_target = int(gene_mask.sum())
            if n_target < min_cells:
                statuses.append(
                    _status_row(
                        sample_id,
                        cell_type,
                        gene,
                        "skipped",
                        "insufficient_target_cells",
                        n_target,
                        n_control,
                    )
                )
                continue
            work = subset[control_mask | gene_mask].copy()
            if temporary:
                stamp_control_label(
                    work.obs, control_cell_mask(work.obs, spec), "target_gene", control
                )
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message="divide by zero encountered in log2",
                        category=RuntimeWarning,
                    )
                    warnings.filterwarnings(
                        "ignore",
                        message="invalid value encountered in subtract",
                        category=RuntimeWarning,
                    )
                    frame = WilcoxonTest.compare_groups(
                        work,
                        column="target_gene",
                        baseline=control,
                        groups_to_compare=[gene],
                    )
            except Exception as exc:
                statuses.append(
                    _status_row(
                        sample_id,
                        cell_type,
                        gene,
                        "skipped",
                        f"wilcoxon_failed:{exc}",
                        n_target,
                        n_control,
                    )
                )
                continue
            frame = _sanitize_wilcoxon_log_fc(frame)
            frame["target_gene"] = gene
            frame["cell_type"] = cell_type
            frame["control"] = control
            frame["control_source"] = control_source
            frame["temporary_control"] = temporary
            frame["n_target_cells"] = n_target
            frame["n_control_cells"] = n_control
            frame["method"] = "WilcoxonTest"
            frame["inference_note"] = INFERENCE_NOTE
            results.append(frame)
            statuses.append(
                _status_row(
                    sample_id,
                    cell_type,
                    gene,
                    "estimated",
                    "",
                    n_target,
                    n_control,
                )
            )

    table = (
        pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=WILCOXON_COLUMNS)
    )
    status = pd.DataFrame(statuses)
    n_log_fc = int(len(table)) if not table.empty else 0
    n_log_fc_na = (
        int(pd.to_numeric(table["log_fc"], errors="coerce").isna().sum()) if n_log_fc else 0
    )
    adata.uns["wilcoxon_de"] = {
        "method": "pertpy.tools.WilcoxonTest",
        "layer": "X",
        "control": control,
        "control_source": control_source,
        "temporary_control": temporary,
        "control_note": str(spec.get("note", "")),
        "n_contrasts": int(status["status"].eq("estimated").sum()) if not status.empty else 0,
        "n_log_fc_rows": n_log_fc,
        "n_log_fc_nonfinite_set_na": n_log_fc_na,
        "log_fc_note": LOG_FC_NOTE,
        "inference_note": INFERENCE_NOTE,
        "reference": (
            "https://pertpy.scverse.org/en/stable/tutorials/notebooks/"
            "differential_gene_expression.html"
        ),
    }
    return table, status


def _sanitize_wilcoxon_log_fc(frame: pd.DataFrame) -> pd.DataFrame:
    """Replace non-finite pertpy log_fc (±inf / NaN from zero means) with NA."""
    frame = frame.copy()
    if "log_fc" not in frame.columns:
        return frame
    values = pd.to_numeric(frame["log_fc"], errors="coerce").to_numpy(dtype=float, copy=True)
    values[~np.isfinite(values)] = np.nan
    frame["log_fc"] = values
    return frame


def _analysis_populations(adata: AnnData) -> list[tuple[str, AnnData]]:
    if str(adata.uns.get("cell_type_scope", "unstratified")) != "annotated":
        return [("ALL", adata)]
    column = str(adata.uns.get("cell_type_col", "cell_type"))
    groups = []
    for name, index in adata.obs.groupby(column, observed=True).groups.items():
        label = str(name)
        if label in {"nan", "None", "<NA>", ""}:
            continue
        groups.append((label, adata[index]))
    return groups or [("ALL", adata)]


def _status_row(
    sample_id: str,
    cell_type: str,
    target: str,
    status: str,
    reason: str,
    n_target: int,
    n_control: int,
) -> dict[str, Any]:
    return {
        "sample_id": sample_id,
        "branch": "main",
        "cell_type": cell_type,
        "target_gene": target,
        "module": "wilcoxon",
        "status": status,
        "reason": reason,
        "n_cells": int(n_target + n_control),
        "n_target_cells": int(n_target),
        "n_control_cells": int(n_control),
        "n_guides": 0,
        "biological_sample_count": 1,
    }
