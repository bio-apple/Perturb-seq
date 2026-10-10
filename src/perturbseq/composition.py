from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

DEFAULT_COMPOSITION_COLUMNS = (
    "gene_target",
    "guide_id",
    "perturbation",
    "num_features",
    "sample_id",
    "leiden",
    "cell_state",
    "phase",
    "mixscape_class_global",
    "mixscape_class",
)


def composition_snapshot(
    adata: AnnData,
    stage: str,
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Long-form cell counts by category at a pipeline stage.

    Used to surface how QC / singlet / Mixscape (and similar) filters change
    sample composition before interpreting perturbation effects.
    """
    n = int(adata.n_obs)
    rows: list[dict] = [
        {
            "stage": stage,
            "category": "n_cells",
            "value": "all",
            "n_cells": n,
            "fraction": 1.0,
        }
    ]
    wanted = list(columns) if columns is not None else list(DEFAULT_COMPOSITION_COLUMNS)
    for col in wanted:
        if col not in adata.obs:
            continue
        counts = adata.obs[col].astype(str).value_counts(dropna=False)
        for value, count in counts.items():
            rows.append(
                {
                    "stage": stage,
                    "category": col,
                    "value": str(value),
                    "n_cells": int(count),
                    "fraction": float(count) / n if n else 0.0,
                }
            )
    return pd.DataFrame(rows)


def append_composition_audit(
    rows: list[pd.DataFrame],
    adata: AnnData,
    stage: str,
    columns: Sequence[str] | None = None,
) -> None:
    rows.append(composition_snapshot(adata, stage, columns=columns))


def write_composition_audit(rows: list[pd.DataFrame], path: Path) -> pd.DataFrame:
    if not rows:
        table = pd.DataFrame(columns=["stage", "category", "value", "n_cells", "fraction"])
    else:
        table = pd.concat(rows, ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table


def composition_by_perturbation(
    adata: AnnData,
    *,
    cluster_key: str = "leiden",
    perturbation_key: str = "gene_target",
) -> pd.DataFrame:
    """Per-perturbation counts and fractions within each cluster (推荐 cell-composition module)."""
    if cluster_key not in adata.obs or perturbation_key not in adata.obs:
        return pd.DataFrame(columns=["cluster", "perturbation", "n_cells", "fraction_in_cluster"])
    frame = adata.obs[[cluster_key, perturbation_key]].copy()
    frame.columns = ["cluster", "perturbation"]
    frame["cluster"] = frame["cluster"].astype(str)
    frame["perturbation"] = frame["perturbation"].astype(str)
    counts = (
        frame.groupby(["cluster", "perturbation"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = counts.groupby("cluster", observed=True)["n_cells"].transform("sum")
    counts["fraction_in_cluster"] = counts["n_cells"] / totals.replace(0, np.nan)
    counts["fraction_in_cluster"] = counts["fraction_in_cluster"].fillna(0.0)
    return counts


def write_composition_by_perturbation(table: pd.DataFrame, path: Path) -> pd.DataFrame:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table


def guide_assignment_summary(adata: AnnData, *, control: str = "NT") -> pd.DataFrame:
    """Compact assignment summary for gRNA QC (singlet rate, NTC, top guides)."""
    n = int(adata.n_obs)
    rows: list[dict] = [{"metric": "n_cells", "value": str(n), "n_cells": n, "fraction": 1.0}]
    if "num_features" in adata.obs:
        for val, count in adata.obs["num_features"].value_counts().sort_index().items():
            rows.append(
                {
                    "metric": "num_features",
                    "value": str(val),
                    "n_cells": int(count),
                    "fraction": float(count) / n if n else 0.0,
                }
            )
        singlets = int((adata.obs["num_features"] == 1).sum())
        rows.append(
            {
                "metric": "singlet_rate",
                "value": "num_features==1",
                "n_cells": singlets,
                "fraction": float(singlets) / n if n else 0.0,
            }
        )
    if "gene_target" in adata.obs:
        n_ntc = int(adata.obs["gene_target"].astype(str).eq(str(control)).sum())
        rows.append(
            {
                "metric": "ntc_cells",
                "value": str(control),
                "n_cells": n_ntc,
                "fraction": float(n_ntc) / n if n else 0.0,
            }
        )
        n_targets = int(adata.obs["gene_target"].astype(str).nunique())
        rows.append(
            {
                "metric": "n_gene_targets",
                "value": str(n_targets),
                "n_cells": n,
                "fraction": 1.0,
            }
        )
    if "guide_id" in adata.obs:
        rows.append(
            {
                "metric": "n_guides",
                "value": str(int(adata.obs["guide_id"].astype(str).nunique())),
                "n_cells": n,
                "fraction": 1.0,
            }
        )
    return pd.DataFrame(rows)


def write_guide_assignment_summary(table: pd.DataFrame, path: Path) -> pd.DataFrame:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    return table
