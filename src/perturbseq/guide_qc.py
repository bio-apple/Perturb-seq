"""Guide-level QC and perturbation efficacy (separate from gene-level evidence)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

DEFAULT_MIN_CELLS = 10
DEFAULT_MIN_MEDIAN_UMI = 5.0
DEFAULT_MIN_DETECTION_RATE = 0.5
DEFAULT_LFC_NOISE = 0.25
DEFAULT_ON_TARGET_LFC_CUTOFF = 0.25
DEFAULT_ON_TARGET_MIN_FAIL_GUIDES = 2
DEFAULT_PERTURBATION_TYPE = "KO"
# none = flag only (default, conservative); other modes add weighted summaries
# without discarding per-guide rows or clearing inconsistency flags.
GUIDE_MERGE_MODES = ("none", "equal", "umi", "confidence", "umi_confidence")
DEFAULT_GUIDE_MERGE = "none"
_CONFIDENCE_WEIGHT = {"high": 3.0, "medium": 2.0, "low": 1.0}
# Tokens whose expected on-target mRNA change is upregulation.
_ACTIVATION_TYPES = frozenset({"A", "CRISPRA", "ACTIVATION", "OE", "OVEREXPRESSION"})


def _to_1d(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.todense()).ravel()
    arr = np.asarray(matrix)
    return arr.ravel() if arr.ndim > 1 else arr


def _gene_index(adata: AnnData, gene: str) -> int | None:
    if gene not in adata.var_names:
        return None
    loc = adata.var_names.get_loc(gene)
    return int(loc) if isinstance(loc, (int, np.integer)) else int(loc.start)


def _mean_expr(adata: AnnData, cell_mask: np.ndarray, gene_idx: int | None) -> float:
    if gene_idx is None or not cell_mask.any():
        return float("nan")
    return float(_to_1d(adata.X[cell_mask, gene_idx]).mean())


def _control_baselines(obs: pd.DataFrame, control: str) -> dict[str, float]:
    ctrl = obs["gene_target"].astype(str) == control
    if not ctrl.any():
        return {}
    out: dict[str, float] = {"n_cells": float(ctrl.sum())}
    for col in ("n_counts", "n_genes", "pct_counts_mt", "guide_umi"):
        if col in obs.columns:
            out[col] = float(obs.loc[ctrl, col].astype(float).mean())
    return out


def _effect_direction(log2fc: float, noise: float = DEFAULT_LFC_NOISE) -> str:
    if not np.isfinite(log2fc):
        return "unknown"
    if log2fc <= -noise:
        return "down"
    if log2fc >= noise:
        return "up"
    return "none"


def expected_on_target_direction(perturbation_type: str = DEFAULT_PERTURBATION_TYPE) -> str:
    """Expected target-gene mRNA direction for this perturbation modality."""
    token = str(perturbation_type or DEFAULT_PERTURBATION_TYPE).strip().upper()
    if token in _ACTIVATION_TYPES:
        return "up"
    return "down"


def on_target_effect_score(
    log2fc: float,
    perturbation_type: str = DEFAULT_PERTURBATION_TYPE,
) -> float:
    """Higher = stronger on-target change in the expected direction (KO/KD/i ↓, CRISPRa ↑)."""
    if not np.isfinite(log2fc):
        return float("nan")
    expected = expected_on_target_direction(perturbation_type)
    return float(-log2fc if expected == "down" else log2fc)


def on_target_pass(
    log2fc: float,
    perturbation_type: str = DEFAULT_PERTURBATION_TYPE,
    lfc_cutoff: float = DEFAULT_ON_TARGET_LFC_CUTOFF,
) -> bool:
    """Pragmatic on-target check: signed log2FC vs control crosses cutoff in expected direction."""
    if not np.isfinite(log2fc) or not np.isfinite(lfc_cutoff):
        return False
    expected = expected_on_target_direction(perturbation_type)
    if expected == "down":
        return bool(log2fc <= -abs(lfc_cutoff))
    return bool(log2fc >= abs(lfc_cutoff))


def _guide_assignment_adequate(row: pd.Series) -> bool:
    n_ok = int(row["n_cells"]) >= int(row["min_cells_threshold"])
    umi = float(row["median_guide_umi"])
    umi_ok = (not np.isfinite(umi)) or umi >= float(row["min_median_umi_threshold"])
    det = float(row["detection_rate"])
    det_ok = (not np.isfinite(det)) or det >= float(row["min_detection_rate_threshold"])
    return bool(n_ok and umi_ok and det_ok)


def _interpretation(
    *,
    n_cells: int,
    min_cells: int,
    median_umi: float,
    min_umi: float,
    detection_rate: float,
    min_detection: float,
    target_in_matrix: bool,
    direction: str,
    gene_consistent: bool | None,
    n_guides_for_gene: int,
    low_efficiency_guide: bool = False,
) -> str:
    """Distinguish low-efficiency / technical failure from unclear biology."""
    if n_cells < min_cells:
        return "too_few_cells"
    if median_umi < min_umi:
        return "low_guide_umi"
    if detection_rate < min_detection:
        return "low_detection"
    # Prefer low-efficiency over "no phenotype" when this guide fails on-target.
    if low_efficiency_guide:
        return "potential_low_efficiency_guide"
    # This guide's own on-target signal beats gene-level inconsistency labeling.
    if target_in_matrix and direction == "down":
        return "target_knockdown_detected"
    if target_in_matrix and direction == "up":
        return "target_upregulated"
    if gene_consistent is False and n_guides_for_gene >= 2:
        return "inconsistent_guides"
    if target_in_matrix and direction == "none":
        return "no_target_effect_with_adequate_guides"
    if not target_in_matrix:
        return "target_not_in_expression_matrix"
    return "adequate_assignment_unclear_effect"


def _gene_on_target_interpretation(
    *,
    potential_low_efficiency: bool,
    guides_consistent: bool | None,
    n_on_target_pass: int,
    n_adequate_guides: int,
    majority_direction: str,
    expected_direction: str,
) -> str:
    """Gene-level label: low-efficiency ≠ no phenotype."""
    if potential_low_efficiency:
        return "potential_low_efficiency_guides"
    if guides_consistent is False:
        return "inconsistent_guides"
    if n_on_target_pass >= 1 and majority_direction == expected_direction:
        return "on_target_effect_detected"
    if n_adequate_guides >= 1 and majority_direction == "none":
        return "no_clear_on_target_effect"
    if n_adequate_guides == 0:
        return "no_adequate_guides"
    return "unclear_on_target"


def compute_guide_metrics(
    adata: AnnData,
    control: str = "NT",
    sample_col: str = "sample_id",
    min_cells: int = DEFAULT_MIN_CELLS,
    min_median_umi: float = DEFAULT_MIN_MEDIAN_UMI,
    min_detection_rate: float = DEFAULT_MIN_DETECTION_RATE,
    lfc_noise: float = DEFAULT_LFC_NOISE,
    perturbation_type: str = DEFAULT_PERTURBATION_TYPE,
    on_target_lfc_cutoff: float = DEFAULT_ON_TARGET_LFC_CUTOFF,
) -> pd.DataFrame:
    """Per-guide QC / efficacy metrics (does not merge guides into gene-level cells)."""
    obs = adata.obs
    required = {"guide_id", "gene_target"}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(
            f"Guide column mismatch in AnnData.obs: missing {sorted(missing)}. "
            f"Found: {list(obs.columns)}. "
            "Expected tertiary h5ad from `python -m perturbseq run` "
            "(obs must include guide_id, gene_target).\n"
            "Suggested commands:\n"
            "  python -m perturbseq run --input-dir <dir> --output-dir <out>\n"
            "  python -m perturbseq guide-qc --h5ad <sample>.tertiary.h5ad --output-dir <out>"
        )

    baselines = _control_baselines(obs, control)
    ctrl_mask = obs["gene_target"].astype(str).to_numpy() == control
    has_sample = sample_col in obs.columns
    n_samples = int(obs[sample_col].astype(str).nunique()) if has_sample else 1
    samples = obs[sample_col].astype(str) if has_sample else pd.Series(["all"] * adata.n_obs, index=obs.index)

    gene_targets = sorted(
        {
            str(g)
            for g in obs["gene_target"].astype(str).unique()
            if g not in {control, "unassigned", ""}
        }
    )
    gene_idx = {g: _gene_index(adata, g) for g in gene_targets}
    ctrl_expr = {
        g: _mean_expr(adata, ctrl_mask, idx) if idx is not None else float("nan")
        for g, idx in gene_idx.items()
    }

    has_mixscape = "mixscape_class_global" in obs.columns
    has_phase = "phase" in obs.columns
    has_state = "cell_state" in obs.columns

    rows: list[dict] = []
    guide_ids = [g for g in obs["guide_id"].astype(str).unique() if g and g.lower() != "nan"]
    for guide in sorted(guide_ids):
        mask = obs["guide_id"].astype(str).to_numpy() == guide
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        sub = obs.loc[mask]
        gene = str(sub["gene_target"].astype(str).mode().iloc[0])
        is_control = gene == control or guide.lower().startswith("non-targeting")
        median_umi = float(sub["guide_umi"].astype(float).median()) if "guide_umi" in sub else float("nan")
        mean_umi = float(sub["guide_umi"].astype(float).mean()) if "guide_umi" in sub else float("nan")
        mean_n_features = float(sub["num_features"].astype(float).mean()) if "num_features" in sub else float("nan")
        frac_singlet = float((sub["num_features"].astype(int) == 1).mean()) if "num_features" in sub else float("nan")

        detected_samples = int(samples.loc[mask].nunique()) if n_cells else 0
        detection_rate = detected_samples / n_samples if n_samples else float("nan")

        target_in_matrix = (not is_control) and gene_idx.get(gene) is not None
        target_mean = float("nan")
        target_log2fc = float("nan")
        if target_in_matrix:
            target_mean = _mean_expr(adata, mask, gene_idx[gene])
            base = ctrl_expr.get(gene, float("nan"))
            if np.isfinite(target_mean) and np.isfinite(base) and base > 0:
                target_log2fc = float(np.log2((target_mean + 1e-6) / (base + 1e-6)))
            elif np.isfinite(target_mean) and np.isfinite(base):
                target_log2fc = float(target_mean - base)

        direction = "control" if is_control else _effect_direction(target_log2fc, lfc_noise)
        score = (
            float("nan")
            if is_control
            else on_target_effect_score(target_log2fc, perturbation_type)
        )
        passed = (
            False
            if is_control or not target_in_matrix
            else on_target_pass(target_log2fc, perturbation_type, on_target_lfc_cutoff)
        )

        row: dict = {
            "guide_id": guide,
            "gene_target": gene,
            "is_control": is_control,
            "n_cells": n_cells,
            "n_samples_detected": detected_samples,
            "n_samples_total": n_samples,
            "detection_rate": detection_rate,
            "median_guide_umi": median_umi,
            "mean_guide_umi": mean_umi,
            "mean_num_features": mean_n_features,
            "frac_singlet": frac_singlet,
            "assignment_confidence": (
                "high"
                if n_cells >= min_cells and median_umi >= min_median_umi and frac_singlet >= 0.9
                else "low"
                if n_cells < min_cells or median_umi < min_median_umi
                else "medium"
            ),
            "target_in_matrix": target_in_matrix,
            "target_expr_mean": target_mean,
            "target_expr_log2fc_vs_control": target_log2fc,
            "target_effect_direction": direction,
            "perturbation_type": str(perturbation_type),
            "expected_on_target_direction": (
                "control" if is_control else expected_on_target_direction(perturbation_type)
            ),
            "on_target_score": score,
            "on_target_lfc_cutoff": float(on_target_lfc_cutoff),
            "on_target_pass": passed,
        }

        for col, key in (
            ("n_counts", "mean_n_counts"),
            ("n_genes", "mean_n_genes"),
            ("pct_counts_mt", "mean_pct_mt"),
        ):
            if col in sub.columns:
                val = float(sub[col].astype(float).mean())
                row[key] = val
                base = baselines.get(col)
                row[f"delta_{key}_vs_control"] = float(val - base) if base is not None else float("nan")

        if has_mixscape and not is_control:
            mg = sub["mixscape_class_global"].astype(str)
            row["frac_mixscape_ko"] = float((mg == "KO").mean())
            row["frac_mixscape_np"] = float((mg == "NP").mean())
        if has_phase:
            row["top_phase"] = str(sub["phase"].astype(str).value_counts().index[0])
            row["frac_top_phase"] = float(sub["phase"].astype(str).value_counts(normalize=True).iloc[0])
        if has_state:
            row["top_cell_state"] = str(sub["cell_state"].astype(str).value_counts().index[0])
            row["frac_top_cell_state"] = float(sub["cell_state"].astype(str).value_counts(normalize=True).iloc[0])

        row["min_cells_threshold"] = min_cells
        row["min_median_umi_threshold"] = min_median_umi
        row["min_detection_rate_threshold"] = min_detection_rate
        rows.append(row)

    guide_df = pd.DataFrame(rows)
    if guide_df.empty:
        return guide_df

    # Adequate assignment + target in matrix + fail on-target → low-efficiency candidate.
    adequate = guide_df.apply(_guide_assignment_adequate, axis=1)
    guide_df["low_efficiency_guide"] = (
        (~guide_df["is_control"])
        & guide_df["target_in_matrix"].astype(bool)
        & adequate
        & (~guide_df["on_target_pass"].astype(bool))
    )

    # Fill interpretation after gene-level consistency is known (placeholder; updated by run_guide_qc).
    guide_df["interpretation"] = "pending"
    return guide_df


def compute_gene_guide_consistency(
    guide_df: pd.DataFrame,
    lfc_noise: float = DEFAULT_LFC_NOISE,
    on_target_min_fail_guides: int = DEFAULT_ON_TARGET_MIN_FAIL_GUIDES,
) -> pd.DataFrame:
    """Compare effect directions across guides targeting the same gene (no cell merging)."""
    if guide_df.empty:
        return pd.DataFrame()

    rows: list[dict] = []
    working = guide_df.loc[~guide_df["is_control"]].copy()
    for gene, sub in working.groupby("gene_target", sort=True):
        directions = sub["target_effect_direction"].tolist()
        scored = [d for d in directions if d in {"down", "up", "none"}]
        n_guides = len(sub)
        n_with_cells = int((sub["n_cells"] > 0).sum())
        if not scored:
            consistent = None
            majority = "unknown"
        else:
            majority = max(set(scored), key=scored.count)
            # "none" does not count as disagreeing with another "none"; opposing up/down is inconsistent.
            signed = [d for d in scored if d in {"up", "down"}]
            if len(set(signed)) > 1:
                consistent = False
            elif signed and any(d == "none" for d in scored):
                # Some guides show effect, others do not → inconsistent efficacy
                consistent = False
            else:
                consistent = True

        lfcs = sub["target_expr_log2fc_vs_control"].astype(float)
        finite = lfcs[np.isfinite(lfcs)]
        lfc_spread = float(finite.max() - finite.min()) if len(finite) else float("nan")
        outlier_guides: list[str] = []
        if len(finite) >= 3:
            med = float(finite.median())
            mad = float(np.median(np.abs(finite - med))) * 1.4826
            if mad > 0:
                for guide_id, lfc in zip(sub["guide_id"], lfcs, strict=False):
                    if np.isfinite(lfc) and abs(lfc - med) > max(2.5 * mad, lfc_noise * 2):
                        outlier_guides.append(str(guide_id))

        n_adequate = int(
            (
                (sub["n_cells"] >= sub["min_cells_threshold"])
                & (sub["median_guide_umi"] >= sub["min_median_umi_threshold"])
            ).sum()
        )
        n_pass = int(sub["on_target_pass"].astype(bool).sum()) if "on_target_pass" in sub.columns else 0
        n_low_eff = (
            int(sub["low_efficiency_guide"].astype(bool).sum())
            if "low_efficiency_guide" in sub.columns
            else 0
        )
        potential_low_eff = bool(n_low_eff >= int(on_target_min_fail_guides))
        ptype = (
            str(sub["perturbation_type"].iloc[0])
            if "perturbation_type" in sub.columns
            else DEFAULT_PERTURBATION_TYPE
        )
        expected = expected_on_target_direction(ptype)
        low_eff_ids = (
            sub.loc[sub["low_efficiency_guide"].astype(bool), "guide_id"].astype(str).tolist()
            if "low_efficiency_guide" in sub.columns
            else []
        )
        gene_interp = _gene_on_target_interpretation(
            potential_low_efficiency=potential_low_eff,
            guides_consistent=consistent,
            n_on_target_pass=n_pass,
            n_adequate_guides=n_adequate,
            majority_direction=str(majority),
            expected_direction=expected,
        )

        rows.append(
            {
                "gene_target": gene,
                "n_guides": n_guides,
                "n_guides_with_cells": n_with_cells,
                "total_cells": int(sub["n_cells"].sum()),
                "directions": "|".join(f"{g}:{d}" for g, d in zip(sub["guide_id"], directions, strict=False)),
                "majority_direction": majority,
                "guides_consistent": consistent,
                "target_lfc_spread": lfc_spread,
                "median_target_log2fc": float(finite.median()) if len(finite) else float("nan"),
                "outlier_guides": "|".join(outlier_guides),
                "n_adequate_guides": n_adequate,
                "n_on_target_pass": n_pass,
                "n_low_efficiency_guides": n_low_eff,
                "low_efficiency_guides": "|".join(low_eff_ids),
                "potential_low_efficiency": potential_low_eff,
                "expected_on_target_direction": expected,
                "on_target_min_fail_guides": int(on_target_min_fail_guides),
                "interpretation": gene_interp,
            }
        )
    return pd.DataFrame(rows)


def _confidence_weight(label: object) -> float:
    return float(_CONFIDENCE_WEIGHT.get(str(label).lower(), 1.0))


def _guide_weights(sub: pd.DataFrame, mode: str) -> np.ndarray:
    n = len(sub)
    if mode == "equal":
        return np.ones(n, dtype=float)
    umi = (
        sub["median_guide_umi"].astype(float).to_numpy()
        if "median_guide_umi" in sub.columns
        else np.ones(n, dtype=float)
    )
    umi = np.where(np.isfinite(umi) & (umi > 0), umi, 0.0)
    conf = (
        sub["assignment_confidence"].map(_confidence_weight).to_numpy(dtype=float)
        if "assignment_confidence" in sub.columns
        else np.ones(n, dtype=float)
    )
    if mode == "umi":
        return umi
    if mode == "confidence":
        return conf
    if mode == "umi_confidence":
        return umi * conf
    raise ValueError(f"Unknown guide_merge mode for weights: {mode}")


def compute_weighted_guide_summary(
    guide_df: pd.DataFrame,
    guide_merge: str = DEFAULT_GUIDE_MERGE,
    lfc_noise: float = DEFAULT_LFC_NOISE,
) -> pd.DataFrame:
    """Per-gene weighted effect summary (does not merge cells or drop guide rows).

    ``guide_merge=none`` returns an empty frame (conservative default: flag only).
    Other modes compute weighted mean target log2FC and a weight-voted direction.
    """
    mode = str(guide_merge or DEFAULT_GUIDE_MERGE).lower()
    if mode not in GUIDE_MERGE_MODES:
        raise ValueError(f"guide_merge must be one of {GUIDE_MERGE_MODES}, got {guide_merge!r}")
    if mode == "none" or guide_df.empty:
        return pd.DataFrame(
            columns=[
                "gene_target",
                "guide_merge",
                "weighted_target_log2fc",
                "weighted_direction",
                "weighted_direction_vote",
                "weight_sum",
                "n_guides_weighted",
                "guide_weights",
            ]
        )

    rows: list[dict] = []
    working = guide_df.loc[~guide_df["is_control"]].copy()
    for gene, sub in working.groupby("gene_target", sort=True):
        weights = _guide_weights(sub, mode)
        lfcs = sub["target_expr_log2fc_vs_control"].astype(float).to_numpy()
        finite = np.isfinite(lfcs) & (weights > 0)
        if finite.any():
            w = weights[finite]
            wsum = float(w.sum())
            wlfc = float(np.dot(lfcs[finite], w) / wsum) if wsum > 0 else float("nan")
        else:
            wsum = float(weights.sum())
            wlfc = float("nan")
        wdir = _effect_direction(wlfc, lfc_noise)

        # Weight-voted direction among scored guides
        vote: dict[str, float] = {}
        for d, w in zip(sub["target_effect_direction"].tolist(), weights, strict=False):
            if d in {"up", "down", "none"} and w > 0:
                vote[d] = vote.get(d, 0.0) + float(w)
        voted = max(vote, key=vote.get) if vote else "unknown"

        weight_labels = [
            f"{g}:{w:.4g}" for g, w in zip(sub["guide_id"].astype(str), weights, strict=False)
        ]
        rows.append(
            {
                "gene_target": gene,
                "guide_merge": mode,
                "weighted_target_log2fc": wlfc,
                "weighted_direction": wdir,
                "weighted_direction_vote": voted,
                "weight_sum": wsum if finite.any() else float(weights.sum()),
                "n_guides_weighted": int((weights > 0).sum()),
                "guide_weights": "|".join(weight_labels),
            }
        )
    return pd.DataFrame(rows)


def attach_weighted_summary(
    consistency_df: pd.DataFrame,
    weighted_df: pd.DataFrame,
    guide_merge: str = DEFAULT_GUIDE_MERGE,
) -> pd.DataFrame:
    """Join weighted columns onto gene_guide_consistency; keep inconsistency flags."""
    mode = str(guide_merge or DEFAULT_GUIDE_MERGE).lower()
    if consistency_df.empty:
        return weighted_df.copy() if not weighted_df.empty else consistency_df.copy()
    out = consistency_df.copy()
    out["guide_merge"] = mode
    if weighted_df.empty or mode == "none":
        return out
    cols = [
        "weighted_target_log2fc",
        "weighted_direction",
        "weighted_direction_vote",
        "weight_sum",
        "n_guides_weighted",
        "guide_weights",
    ]
    merge_cols = ["gene_target"] + [c for c in cols if c in weighted_df.columns]
    out = out.drop(columns=[c for c in cols if c in out.columns], errors="ignore")
    return out.merge(weighted_df[merge_cols], on="gene_target", how="left")


def build_qc_warnings(
    guide_df: pd.DataFrame,
    consistency_df: pd.DataFrame,
) -> pd.DataFrame:
    """Flag insufficient cells, low UMI/detection, outlier / inconsistent guide effects."""
    warnings: list[dict] = []
    if guide_df.empty:
        return pd.DataFrame(columns=["level", "entity", "warning", "detail"])

    for _, row in guide_df.iterrows():
        if row["is_control"]:
            continue
        guide = row["guide_id"]
        if row["n_cells"] < row["min_cells_threshold"]:
            warnings.append(
                {
                    "level": "guide",
                    "entity": guide,
                    "gene_target": row["gene_target"],
                    "warning": "insufficient_cells",
                    "detail": f"n_cells={row['n_cells']} < {row['min_cells_threshold']}",
                }
            )
        if np.isfinite(row["median_guide_umi"]) and row["median_guide_umi"] < row["min_median_umi_threshold"]:
            warnings.append(
                {
                    "level": "guide",
                    "entity": guide,
                    "gene_target": row["gene_target"],
                    "warning": "low_guide_umi",
                    "detail": f"median_guide_umi={row['median_guide_umi']:.2f} < {row['min_median_umi_threshold']}",
                }
            )
        if np.isfinite(row["detection_rate"]) and row["detection_rate"] < row["min_detection_rate_threshold"]:
            warnings.append(
                {
                    "level": "guide",
                    "entity": guide,
                    "gene_target": row["gene_target"],
                    "warning": "low_detection",
                    "detail": (
                        f"detection_rate={row['detection_rate']:.2f} "
                        f"({row['n_samples_detected']}/{row['n_samples_total']} samples)"
                    ),
                }
            )
        if "delta_mean_n_counts_vs_control" in row and "mean_n_counts" in row:
            delta = row["delta_mean_n_counts_vs_control"]
            mean_counts = row["mean_n_counts"]
            if np.isfinite(delta) and np.isfinite(mean_counts):
                ctrl_mean = mean_counts - delta
                if ctrl_mean > 0 and mean_counts < 0.5 * ctrl_mean:
                    warnings.append(
                        {
                            "level": "guide",
                            "entity": guide,
                            "gene_target": row["gene_target"],
                            "warning": "cytotoxicity_outlier",
                            "detail": f"mean_n_counts={mean_counts:.0f} < 50% of control ({ctrl_mean:.0f})",
                        }
                    )
        if bool(row.get("cytotoxicity_qc_flag")):
            frac = row.get("frac_removed_qc", float("nan"))
            n_before = row.get("n_cells_before_qc", float("nan"))
            n_removed = row.get("n_cells_removed_qc", float("nan"))
            warnings.append(
                {
                    "level": "guide",
                    "entity": guide,
                    "gene_target": row["gene_target"],
                    "warning": "cytotoxicity_qc_depletion",
                    "detail": (
                        f"QC removed {frac:.0%} of cells "
                        f"({n_removed:.0f}/{n_before:.0f}) — strong cytotoxicity / phenotype signal; "
                        f"do not treat as silent technical discard"
                    ),
                }
            )
        if bool(row.get("low_efficiency_guide")):
            score = row.get("on_target_score", float("nan"))
            lfc = row.get("target_expr_log2fc_vs_control", float("nan"))
            expected = row.get("expected_on_target_direction", "down")
            cutoff = row.get("on_target_lfc_cutoff", DEFAULT_ON_TARGET_LFC_CUTOFF)
            warnings.append(
                {
                    "level": "guide",
                    "entity": guide,
                    "gene_target": row["gene_target"],
                    "warning": "low_efficiency_guide",
                    "detail": (
                        f"on_target_pass=False (expected {expected}, "
                        f"log2FC={lfc if np.isfinite(lfc) else 'nan'}, "
                        f"score={score if np.isfinite(score) else 'nan'}, "
                        f"cutoff={cutoff}); potential low-efficiency guide — "
                        f"do not conclude no phenotype for the gene"
                    ),
                }
            )

    if not consistency_df.empty:
        for _, row in consistency_df.iterrows():
            if row.get("potential_low_efficiency"):
                warnings.append(
                    {
                        "level": "gene",
                        "entity": row["gene_target"],
                        "gene_target": row["gene_target"],
                        "warning": "potential_low_efficiency_guides",
                        "detail": (
                            f"{row.get('n_low_efficiency_guides', 0)}/{row['n_guides']} guides "
                            f"fail on-target (expected {row.get('expected_on_target_direction', 'down')}); "
                            f"low_efficiency_guides={row.get('low_efficiency_guides', '')} — "
                            f"do not conclude no phenotype"
                        ),
                    }
                )
            if row["guides_consistent"] is False:
                warnings.append(
                    {
                        "level": "gene",
                        "entity": row["gene_target"],
                        "gene_target": row["gene_target"],
                        "warning": "inconsistent_guides",
                        "detail": row["directions"],
                    }
                )
            if row.get("outlier_guides"):
                warnings.append(
                    {
                        "level": "gene",
                        "entity": row["gene_target"],
                        "gene_target": row["gene_target"],
                        "warning": "outlier_guide_effect",
                        "detail": f"outliers={row['outlier_guides']}; spread={row['target_lfc_spread']}",
                    }
                )
            if row["n_adequate_guides"] == 0 and row["n_guides"] > 0:
                warnings.append(
                    {
                        "level": "gene",
                        "entity": row["gene_target"],
                        "gene_target": row["gene_target"],
                        "warning": "no_adequate_guides",
                        "detail": "all guides fail cell-count or UMI thresholds",
                    }
                )
            elif (
                not bool(row.get("potential_low_efficiency"))
                and row["n_adequate_guides"] >= 1
                and row["majority_direction"] == "none"
                and row["guides_consistent"] is True
            ):
                warnings.append(
                    {
                        "level": "gene",
                        "entity": row["gene_target"],
                        "gene_target": row["gene_target"],
                        "warning": "no_clear_effect_adequate_guides",
                        "detail": (
                            "adequate guides but no clear on-target change — "
                            "review assay limits; with ≥2 failing guides see potential_low_efficiency instead"
                        ),
                    }
                )

    return pd.DataFrame(warnings)


def _apply_interpretations(guide_df: pd.DataFrame, consistency_df: pd.DataFrame) -> pd.DataFrame:
    consistency_map = (
        consistency_df.set_index("gene_target")["guides_consistent"].to_dict() if not consistency_df.empty else {}
    )
    n_guides_map = consistency_df.set_index("gene_target")["n_guides"].to_dict() if not consistency_df.empty else {}
    low_eff_map = (
        consistency_df.set_index("gene_target")["potential_low_efficiency"].to_dict()
        if not consistency_df.empty and "potential_low_efficiency" in consistency_df.columns
        else {}
    )
    interpretations = []
    for _, row in guide_df.iterrows():
        if row["is_control"]:
            interpretations.append("control")
            continue
        gene = row["gene_target"]
        interpretations.append(
            _interpretation(
                n_cells=int(row["n_cells"]),
                min_cells=int(row["min_cells_threshold"]),
                median_umi=float(row["median_guide_umi"]),
                min_umi=float(row["min_median_umi_threshold"]),
                detection_rate=float(row["detection_rate"]),
                min_detection=float(row["min_detection_rate_threshold"]),
                target_in_matrix=bool(row["target_in_matrix"]),
                direction=str(row["target_effect_direction"]),
                gene_consistent=consistency_map.get(gene),
                n_guides_for_gene=int(n_guides_map.get(gene, 1)),
                low_efficiency_guide=bool(row.get("low_efficiency_guide", False)),
            )
        )
    guide_df = guide_df.copy()
    guide_df["interpretation"] = interpretations
    guide_df["guides_consistent_for_gene"] = guide_df["gene_target"].map(consistency_map)
    guide_df["potential_low_efficiency_for_gene"] = guide_df["gene_target"].map(low_eff_map)
    return guide_df


def summarize_guide_qc(
    guide_df: pd.DataFrame,
    consistency_df: pd.DataFrame,
    warnings_df: pd.DataFrame,
    *,
    guide_merge: str = DEFAULT_GUIDE_MERGE,
) -> dict:
    warning_counts = (
        warnings_df["warning"].value_counts().to_dict() if not warnings_df.empty else {}
    )
    n_consistent = int(consistency_df["guides_consistent"].eq(True).sum()) if not consistency_df.empty else 0
    n_inconsistent = int(consistency_df["guides_consistent"].eq(False).sum()) if not consistency_df.empty else 0
    n_unknown = int(consistency_df["guides_consistent"].isna().sum()) if not consistency_df.empty else 0
    n_low_eff_genes = (
        int(consistency_df["potential_low_efficiency"].eq(True).sum())
        if not consistency_df.empty and "potential_low_efficiency" in consistency_df.columns
        else 0
    )
    n_low_eff_guides = (
        int(guide_df["low_efficiency_guide"].eq(True).sum())
        if not guide_df.empty and "low_efficiency_guide" in guide_df.columns
        else 0
    )
    return {
        "n_guides": int(len(guide_df)),
        "n_guides_noncontrol": int((~guide_df["is_control"]).sum()) if not guide_df.empty else 0,
        "n_genes_with_multi_guides": int((consistency_df["n_guides"] >= 2).sum()) if not consistency_df.empty else 0,
        "n_genes_consistent_guides": n_consistent,
        "n_genes_inconsistent_guides": n_inconsistent,
        "n_genes_consistency_unknown": n_unknown,
        "n_low_efficiency_guides": n_low_eff_guides,
        "n_genes_potential_low_efficiency": n_low_eff_genes,
        "n_warnings": int(len(warnings_df)),
        "warning_counts": {str(k): int(v) for k, v in warning_counts.items()},
        "interpretation_counts": (
            guide_df.loc[~guide_df["is_control"], "interpretation"].value_counts().to_dict()
            if not guide_df.empty
            else {}
        ),
        "guide_merge": str(guide_merge or DEFAULT_GUIDE_MERGE).lower(),
    }


def merge_qc_filter_stats(
    guide_df: pd.DataFrame,
    filter_df: pd.DataFrame | None,
    *,
    control: str = "NT",
) -> pd.DataFrame:
    """Join per-guide QC removal stats from ``qc_filter_by_guide.csv`` onto guide_qc.

    Guides fully removed by QC (absent from post-QC AnnData) are appended as stub rows
    so cytotoxicity is never a silent discard.
    """
    if filter_df is None or filter_df.empty:
        return guide_df
    if "category" not in filter_df.columns or "value" not in filter_df.columns:
        return guide_df
    guide_rows = filter_df.loc[filter_df["category"].astype(str) == "guide_id"].copy()
    if guide_rows.empty:
        return guide_df
    cols = [
        c
        for c in (
            "n_cells_before_qc",
            "n_cells_after_qc",
            "n_cells_removed_qc",
            "frac_removed_qc",
            "cytotoxicity_qc_flag",
        )
        if c in guide_rows.columns
    ]
    if not cols:
        return guide_df

    out = guide_df.copy() if not guide_df.empty else pd.DataFrame()
    if not out.empty:
        out = out.drop(columns=[c for c in cols if c in out.columns], errors="ignore")
        merge_src = guide_rows.rename(columns={"value": "guide_id"})[["guide_id", *cols]]
        out = out.merge(merge_src, on="guide_id", how="left")
    else:
        out = pd.DataFrame(columns=["guide_id", "gene_target", "is_control", *cols])

    present = set(out["guide_id"].astype(str)) if "guide_id" in out.columns else set()
    stubs: list[dict] = []
    for _, row in guide_rows.iterrows():
        gid = str(row["value"])
        if gid in present:
            continue
        gene = str(row["gene_target"]) if "gene_target" in row.index and pd.notna(row["gene_target"]) else ""
        is_control = gene == control or gid.lower().startswith("non-targeting")
        stub = {
            "guide_id": gid,
            "gene_target": gene or "unknown",
            "is_control": is_control,
            "n_cells": 0,
            "median_guide_umi": float("nan"),
            "detection_rate": float("nan"),
            "min_cells_threshold": DEFAULT_MIN_CELLS,
            "min_median_umi_threshold": DEFAULT_MIN_MEDIAN_UMI,
            "min_detection_rate_threshold": DEFAULT_MIN_DETECTION_RATE,
            "low_efficiency_guide": False,
            "interpretation": "cytotoxicity_qc_depletion" if row.get("cytotoxicity_qc_flag") else "qc_filtered_out",
        }
        for c in cols:
            stub[c] = row[c]
        stubs.append(stub)
    if stubs:
        out = pd.concat([out, pd.DataFrame(stubs)], ignore_index=True)

    if "cytotoxicity_qc_flag" in out.columns:
        out["cytotoxicity_qc_flag"] = out["cytotoxicity_qc_flag"].fillna(False).astype(bool)
        if "is_control" in out.columns:
            cyto = out["cytotoxicity_qc_flag"].astype(bool) & (~out["is_control"].astype(bool))
            out.loc[cyto, "interpretation"] = "cytotoxicity_qc_depletion"
    return out


def run_guide_qc(
    adata: AnnData,
    control: str = "NT",
    sample_col: str = "sample_id",
    min_cells: int = DEFAULT_MIN_CELLS,
    min_median_umi: float = DEFAULT_MIN_MEDIAN_UMI,
    min_detection_rate: float = DEFAULT_MIN_DETECTION_RATE,
    lfc_noise: float = DEFAULT_LFC_NOISE,
    guide_merge: str = DEFAULT_GUIDE_MERGE,
    perturbation_type: str = DEFAULT_PERTURBATION_TYPE,
    on_target_lfc_cutoff: float = DEFAULT_ON_TARGET_LFC_CUTOFF,
    on_target_min_fail_guides: int = DEFAULT_ON_TARGET_MIN_FAIL_GUIDES,
    qc_filter_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    mode = str(guide_merge or DEFAULT_GUIDE_MERGE).lower()
    if mode not in GUIDE_MERGE_MODES:
        raise ValueError(f"guide_merge must be one of {GUIDE_MERGE_MODES}, got {guide_merge!r}")
    guide_df = compute_guide_metrics(
        adata,
        control=control,
        sample_col=sample_col,
        min_cells=min_cells,
        min_median_umi=min_median_umi,
        min_detection_rate=min_detection_rate,
        lfc_noise=lfc_noise,
        perturbation_type=perturbation_type,
        on_target_lfc_cutoff=on_target_lfc_cutoff,
    )
    consistency_df = compute_gene_guide_consistency(
        guide_df,
        lfc_noise=lfc_noise,
        on_target_min_fail_guides=on_target_min_fail_guides,
    )
    weighted_df = compute_weighted_guide_summary(guide_df, guide_merge=mode, lfc_noise=lfc_noise)
    consistency_df = attach_weighted_summary(consistency_df, weighted_df, guide_merge=mode)
    guide_df = _apply_interpretations(guide_df, consistency_df)
    guide_df = merge_qc_filter_stats(guide_df, qc_filter_df, control=control)
    warnings_df = build_qc_warnings(guide_df, consistency_df)
    summary = summarize_guide_qc(guide_df, consistency_df, warnings_df, guide_merge=mode)
    summary["weighted"] = {
        "enabled": mode != "none",
        "n_genes": int(len(weighted_df)),
    }
    summary["on_target"] = {
        "perturbation_type": str(perturbation_type),
        "expected_direction": expected_on_target_direction(perturbation_type),
        "lfc_cutoff": float(on_target_lfc_cutoff),
        "min_fail_guides": int(on_target_min_fail_guides),
    }
    return guide_df, consistency_df, warnings_df, summary


def write_guide_qc_tables(
    guide_df: pd.DataFrame,
    consistency_df: pd.DataFrame,
    warnings_df: pd.DataFrame,
    tables: Path,
    *,
    weighted_df: pd.DataFrame | None = None,
) -> None:
    tables.mkdir(parents=True, exist_ok=True)
    guide_df.to_csv(tables / "guide_qc.csv", index=False)
    consistency_df.to_csv(tables / "gene_guide_consistency.csv", index=False)
    if warnings_df.empty and len(warnings_df.columns) == 0:
        warnings_df = pd.DataFrame(columns=["level", "entity", "gene_target", "warning", "detail"])
    warnings_df.to_csv(tables / "qc_warnings.csv", index=False)
    # Companion weighted table when merge mode produced rows (mode != none).
    if weighted_df is not None and not weighted_df.empty:
        weighted_df.to_csv(tables / "gene_guide_weighted.csv", index=False)
    elif (
        not consistency_df.empty
        and "guide_merge" in consistency_df.columns
        and consistency_df["guide_merge"].astype(str).ne("none").any()
        and "weighted_target_log2fc" in consistency_df.columns
    ):
        cols = [
            c
            for c in (
                "gene_target",
                "guide_merge",
                "weighted_target_log2fc",
                "weighted_direction",
                "weighted_direction_vote",
                "weight_sum",
                "n_guides_weighted",
                "guide_weights",
            )
            if c in consistency_df.columns
        ]
        consistency_df[cols].to_csv(tables / "gene_guide_weighted.csv", index=False)
