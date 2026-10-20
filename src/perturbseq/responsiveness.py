"""Score how separable each perturbation is from non-targeting cells."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from pertpy.tools import Augur

from perturbseq.assign import control_cell_mask
from perturbseq.config import control_spec
from perturbseq.preprocess import assign_cell_type

CLASSIFIERS = {
    "random_forest": "random_forest_classifier",
    "logistic_regression": "logistic_regression_classifier",
}


def run_augur(adata: AnnData, config: dict[str, Any]) -> pd.DataFrame:
    """Augur AUC of perturbation versus NT, stratified by a real cell-type column.

    Without a cell-type annotation every cell is ALL and the population is
    recorded as unstratified. AUC measures label separability in the current
    data, not editing efficiency. The pertpy tutorial call is
    ``Augur.load(..., cell_type_col=...)`` then ``predict``.
    """
    settings = config["augur"]
    if str(settings.get("branch", "main")) != "main":
        raise ValueError("本版 Augur 只定义了 main 分支。")
    assign_cell_type(adata, config)
    scope = str(adata.uns.get("cell_type_scope", "unstratified"))
    cell_type_col = str(adata.uns.get("cell_type_col", "cell_type"))
    spec = control_spec(config)
    control = spec["control_label"]
    if spec["temporary"] and int(control_cell_mask(adata.obs, spec).sum()) == 0:
        return pd.DataFrame(
            columns=[
                "target_gene",
                "cell_type",
                "cell_type_scope",
                "cell_type_col",
                "auc",
                "auc_sd",
                "n_cells",
                "n_target_cells",
                "n_nt_cells",
                "n_folds",
                "seed",
                "reason",
                "auc_meaning",
            ]
        )
    seed = int(config["runtime"]["seed"])
    subsample_size = int(settings["subsample_size"])
    folds = int(settings["n_folds"])
    n_subsamples = int(settings["n_subsamples"])
    feature_perc = float(settings.get("feature_perc", 0.5))
    n_threads = int(settings.get("n_threads", 1))
    select_variance = bool(settings.get("select_variance_features", True))
    classifier_name = str(settings["classifier"])
    model = Augur(
        CLASSIFIERS[classifier_name],
        n_estimators=int(settings.get("n_estimators", 100)),
        max_features="sqrt",
        random_state=seed,
    )
    parameters = {
        "seed": seed,
        "n_subsamples": n_subsamples,
        "subsample_size": subsample_size,
        "n_folds": folds,
        "feature_perc": feature_perc,
        "n_threads": n_threads,
        "select_variance_features": select_variance,
        "classifier": classifier_name,
        "n_estimators": int(settings.get("n_estimators", 100)),
        "augur_mode": "default",
        "cell_type_col": cell_type_col,
        "cell_type_scope": scope,
        "cell_type_note": str(adata.uns.get("cell_type_note", "")),
        "label": (
            "perturbation_vs_zero_guide_cells" if spec["temporary"] else "perturbation_vs_designed_nt"
        ),
        "auc_meaning": "label_separability_not_editing_efficiency",
        "branch": "main",
    }
    if spec["temporary"]:
        values = adata.obs["target_gene"].astype(str).to_numpy(copy=True)
        values[control_cell_mask(adata.obs, spec)] = control
        adata.obs["target_gene"] = values
    nt_mask = control_cell_mask(adata.obs, spec)
    targets = (
        adata.obs.loc[adata.obs["perturbation_status"].astype(str) == "targeting", "target_gene"]
        .astype(str)
        .unique()
    )
    rows: list[dict[str, object]] = []
    for index, target in enumerate(targets, start=1):
        target_mask = adata.obs["target_gene"].astype(str).to_numpy() == target
        subset = adata[target_mask | nt_mask].copy()
        subset.obs["augur_label"] = np.where(
            subset.obs["target_gene"].astype(str).to_numpy() == control, control, target
        )
        eligible: list[str] = []
        for cell_type, group in subset.obs.groupby(cell_type_col, observed=True):
            label = str(cell_type)
            group_target = group["target_gene"].astype(str) != control
            n_target = int(group_target.sum())
            n_nt = int((~group_target).sum())
            base = {
                **parameters,
                "target_gene": target,
                "cell_type": label if label not in {"nan", "None", "<NA>"} else "unannotated",
                "n_cells": n_target + n_nt,
                "n_target_cells": n_target,
                "n_nt_cells": n_nt,
            }
            if label in {"nan", "None", "<NA>", ""}:
                rows.append({**base, "auc": np.nan, "auc_sd": np.nan, "reason": "not_annotated"})
                continue
            if n_target < subsample_size or n_nt < subsample_size:
                rows.append(
                    {**base, "auc": np.nan, "auc_sd": np.nan, "reason": "insufficient_cells"}
                )
                continue
            eligible.append(label)
        if not eligible:
            continue
        if index == 1 or index % 25 == 0:
            print(f"  Augur {index}/{len(targets)}: {target}", flush=True)
        piece = subset[subset.obs[cell_type_col].astype(str).isin(eligible)].copy()
        try:
            loaded = model.load(piece, label_col="augur_label", cell_type_col=cell_type_col)
            _loaded, scores = model.predict(
                loaded,
                n_subsamples=n_subsamples,
                subsample_size=subsample_size,
                folds=folds,
                feature_perc=feature_perc,
                n_threads=n_threads,
                select_variance_features=select_variance,
                random_state=seed,
                augur_mode="default",
            )
        except Exception as exc:
            for label in eligible:
                rows.append(
                    _failed_row(parameters, piece, cell_type_col, target, control, label, exc)
                )
            continue
        summary = scores.get("summary_metrics")
        for label in eligible:
            auc = np.nan
            if (
                isinstance(summary, pd.DataFrame)
                and label in summary.columns
                and "mean_augur_score" in summary.index
            ):
                auc = float(summary.loc["mean_augur_score", label])
            group = piece.obs[cell_type_col].astype(str) == label
            n_target = int((piece.obs.loc[group, "target_gene"].astype(str) != control).sum())
            n_nt = int(group.sum()) - n_target
            fold_values = _fold_scores(scores, label)
            rows.append(
                {
                    **parameters,
                    "target_gene": target,
                    "cell_type": label,
                    "n_cells": n_target + n_nt,
                    "n_target_cells": n_target,
                    "n_nt_cells": n_nt,
                    "auc": auc,
                    "auc_sd": float(np.std(fold_values, ddof=1))
                    if len(fold_values) > 1
                    else np.nan,
                    "reason": "tested" if np.isfinite(auc) else "augur_skipped",
                }
            )
    columns = [
        "target_gene",
        "cell_type",
        "cell_type_scope",
        "cell_type_col",
        "auc",
        "auc_sd",
        "n_cells",
        "n_target_cells",
        "n_nt_cells",
        "n_folds",
        "seed",
        "reason",
        "auc_meaning",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
    table = pd.DataFrame(rows).sort_values("auc", ascending=False, na_position="last")
    return table.reset_index(drop=True)


def _failed_row(
    parameters: dict[str, object],
    piece: AnnData,
    cell_type_col: str,
    target: str,
    control: str,
    label: str,
    exc: Exception,
) -> dict[str, object]:
    group = piece.obs[cell_type_col].astype(str) == label
    n_target = int((piece.obs.loc[group, "target_gene"].astype(str) != control).sum())
    n_nt = int(group.sum()) - n_target
    return {
        **parameters,
        "target_gene": target,
        "cell_type": label,
        "n_cells": n_target + n_nt,
        "n_target_cells": n_target,
        "n_nt_cells": n_nt,
        "auc": np.nan,
        "auc_sd": np.nan,
        "reason": f"augur_failed:{exc}",
    }


def _fold_scores(scores: dict[str, Any], cell_type: str) -> np.ndarray:
    full = scores.get("full_results")
    if not isinstance(full, pd.DataFrame) or "augur_score" not in full.columns:
        return np.array([], dtype=float)
    values = full
    if "cell_type" in full.columns:
        values = full[full["cell_type"].astype(str) == cell_type]
    scores_only = values["augur_score"].to_numpy(dtype=float)
    return scores_only[np.isfinite(scores_only)]
