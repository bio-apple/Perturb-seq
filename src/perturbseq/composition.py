from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

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
