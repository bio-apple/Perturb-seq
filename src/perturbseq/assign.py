"""Parse guide calls and keep a complete, exclusive cell ledger."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

STATUS_NO_GUIDE = "no_guide"
STATUS_SINGLE_TARGET = "single_targeting"
STATUS_SINGLE_NT = "single_nt"
STATUS_MULTI = "multi_guide"
STATUS_TARGET = "targeting"
STATUS_NT = "nt"
STATUS_RETAINED = "retained"
LOW_LIBRARY_NO_GUIDE = "low_library_no_guide"

EXCLUSION_PRIORITY = (
    LOW_LIBRARY_NO_GUIDE,
    "no_guide",
    "multi_guide",
    "low_guide_umi",
    "qc_failed",
    "low_cells_target",
)

MULTI_ALL_NT = "all_nt"
MULTI_SAME_TARGET = "same_target"
MULTI_MIXED_NT = "targeting_nt_mixed"
MULTI_DIFFERENT_TARGETS = "multiple_targets"


class GuideReferenceError(ValueError):
    """A guide id does not have exactly one target_gene_name."""

    def __init__(self, guide_ids: list[str], reason: str) -> None:
        self.guide_ids = list(guide_ids)
        self.reason = reason
        super().__init__(reason + ": " + ", ".join(self.guide_ids))


def name_token_width(reference: pd.DataFrame) -> int:
    widths = reference["name"].astype(str).map(lambda value: value.count("|") + 1)
    return int(widths.max()) if len(widths) else 1


def segment_feature_call(feature_call: str, names: set[str], max_tokens: int) -> list[str]:
    """Split a DRAGEN feature_call on the guide separator.

    Guide names may themselves contain '|'. The longest reference name that
    consumes the next tokens is taken, so a name is not cut at its internal pipe.
    """
    if feature_call == "":
        return []
    tokens = feature_call.split("|")
    parts: list[str] = []
    index = 0
    while index < len(tokens):
        matched: str | None = None
        upper = min(max_tokens, len(tokens) - index)
        for width in range(upper, 0, -1):
            candidate = "|".join(tokens[index : index + width])
            if candidate in names:
                matched = candidate
                index += width
                break
        if matched is None:
            parts.append("|".join(tokens[index:]))
            break
        parts.append(matched)
    return parts


def parse_transcripts(raw: str) -> list[float]:
    if raw == "":
        return []
    values = []
    for field in raw.split("|"):
        if field == "":
            raise ValueError("num_transcripts 含空字段。")
        try:
            values.append(float(field))
        except ValueError as exc:
            raise ValueError(f"num_transcripts 含非数值 UMI: {field}") from exc
    return values


def build_guide_catalog(reference: pd.DataFrame) -> pd.DataFrame:
    """One row per reference guide, with the configured NT label applied later."""
    catalog = pd.DataFrame(
        {
            "guide_id": reference["id"].astype(str),
            "guide_name": reference["name"].astype(str),
            "target_gene_name": reference.get(
                "target_gene_name", pd.Series("", index=reference.index)
            )
            .astype(str)
            .str.strip(),
            "target_gene_id": reference.get("target_gene_id", pd.Series("", index=reference.index))
            .astype(str)
            .str.strip(),
            "feature_type": reference.get(
                "feature_type", pd.Series("", index=reference.index)
            ).astype(str),
        }
    )
    return catalog.reset_index(drop=True)


def validate_guide_targets(reference: pd.DataFrame, *, nt_control: str) -> pd.DataFrame:
    """Use target_gene_name as the only target when that column exists.

    A blank target is an error. The one exception is a guide whose name equals
    the configured NT label: that guide is the designed control even when
    target_gene_name is blank. A blank is never treated as non-targeting for
    any other guide.
    References without the column keep that name-based control path.
    """
    catalog = build_guide_catalog(reference)
    if "target_gene_name" not in reference.columns:
        return catalog
    conflicts = sorted(
        {
            str(guide_id)
            for guide_id, group in catalog.groupby("guide_id", sort=False)
            if group["target_gene_name"].nunique(dropna=False) > 1
        }
    )
    if conflicts:
        raise GuideReferenceError(conflicts, "同一 guide id 对应多个 target_gene_name")
    empty_ids = sorted(
        {
            str(row.guide_id)
            for row in catalog.itertuples(index=False)
            if row.target_gene_name == "" and row.guide_name != nt_control
        }
    )
    if empty_ids:
        raise GuideReferenceError(empty_ids, "target_gene_name 为空")
    return catalog


def target_by_guide_name(catalog: pd.DataFrame) -> dict[str, str]:
    """Map each guide name to the column target. Do not parse the name or id."""
    mapping: dict[str, str] = {}
    conflicts: list[str] = []
    for name, group in catalog.groupby("guide_name", sort=False):
        values = set(group["target_gene_name"].astype(str))
        if len(values) != 1:
            conflicts.extend(str(guide_id) for guide_id in group["guide_id"])
            continue
        mapping[str(name)] = next(iter(values))
    if conflicts:
        raise GuideReferenceError(
            sorted(set(conflicts)), "同一 guide name 对应多个 target_gene_name"
        )
    return mapping


def annotate_assignments(
    assignments: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    nt_control: str,
    multi_guide_policy: str,
    min_guide_umi: float | None,
    matrix_barcodes: set[str],
    targets_per_cell: int = 1,
) -> pd.DataFrame:
    """Classify every assignment row into four exclusive categories.

    Classification uses the number of unique target genes, not the number of
    guides. Two guides aimed at the same gene are ``single_targeting``.
    ``targets_per_cell`` is the unique-target count kept for later analysis; it
    does not add a fifth assignment class. Cells with more unique targets stay
    ``multi_guide``.
    """
    if multi_guide_policy != "drop":
        raise ValueError("本版只支持 multi_guide_policy=drop。")
    if targets_per_cell < 1:
        raise ValueError("targets_per_cell 必须是正整数。")
    catalog = validate_guide_targets(reference, nt_control=nt_control)
    names = set(catalog["guide_name"])
    max_tokens = name_token_width(reference)
    shared_names = set(catalog.loc[catalog["guide_name"].duplicated(keep=False), "guide_name"])
    target_by_name = target_by_guide_name(catalog)
    id_by_unique_name = {}
    for name, group in catalog.groupby("guide_name", sort=False):
        if len(group) == 1:
            id_by_unique_name[str(name)] = str(group["guide_id"].iloc[0])

    rows: list[dict[str, Any]] = []
    for record in assignments.itertuples(index=False):
        call = str(record.feature_call)
        parts = segment_feature_call(call, names, max_tokens)
        transcripts = parse_transcripts(str(record.num_transcripts))
        known = all(part in names for part in parts)
        aligned = len(parts) == int(record.num_features) and len(transcripts) == len(parts)
        barcode = str(record.cell_barcode)
        if barcode not in matrix_barcodes:
            raise ValueError(f"分配 barcode 不在矩阵中: {barcode}")
        if not known or not aligned:
            raise ValueError(f"guide 解析与 num_features 不一致: {barcode}")
        status, guide_name, target, targets = _assignment_for_call(
            parts=parts,
            nt_control=nt_control,
            target_by_name=target_by_name,
        )
        guide_id = id_by_unique_name.get(guide_name, "")
        if guide_name in shared_names:
            guide_id = ""
        guide_umi_pass = min_guide_umi is None or all(
            weight >= float(min_guide_umi) for weight in transcripts
        )
        subtype = (
            _multi_subtype(parts, nt_control, target_by_name) if status == STATUS_MULTI else pd.NA
        )
        analysis_status = {
            STATUS_SINGLE_TARGET: STATUS_TARGET,
            STATUS_SINGLE_NT: STATUS_NT,
        }.get(status, status)
        unique_targets = list(dict.fromkeys(target for target in targets if target))
        rows.append(
            {
                "cell_barcode": barcode,
                "n_guides_called": int(record.num_features),
                "n_guides": len(parts),
                "n_targets": len(unique_targets),
                "feature_call": call,
                "guide_list": "|".join(parts),
                "guide_list_json": json.dumps(parts),
                "target_list_json": json.dumps(targets),
                "guide_name": guide_name,
                "guide_id": guide_id,
                "guide_names": parts,
                "guide_transcripts": transcripts,
                "target_gene": target
                if status in {STATUS_SINGLE_TARGET, STATUS_SINGLE_NT}
                else pd.NA,
                "assignment_status": status,
                "multi_guide_subtype": subtype,
                "guide_umi_pass": guide_umi_pass or not parts,
                "perturbation_status": analysis_status,
                "shared_guide_name": guide_name in shared_names,
                "in_matrix": True,
            }
        )
    return pd.DataFrame.from_records(rows)


def _assignment_for_call(
    *,
    parts: list[str],
    nt_control: str,
    target_by_name: dict[str, str],
) -> tuple[str, str, str, list[str]]:
    """Assign by unique target genes. Same-gene dual guides stay single-target."""
    if not parts:
        return STATUS_NO_GUIDE, "", "", []
    targets = [_target_of(part, nt_control, target_by_name) for part in parts]
    for guide_name, target in zip(parts, targets, strict=True):
        if target == "":
            raise ValueError(f"guide 没有靶基因: {guide_name}")
    unique = list(dict.fromkeys(targets))
    if len(unique) > 1:
        return STATUS_MULTI, "", "", targets
    target = unique[0]
    guide_name = parts[0]
    label = nt_control.strip() if isinstance(nt_control, str) else ""
    if label and target == label:
        return STATUS_SINGLE_NT, guide_name, label, targets
    return STATUS_SINGLE_TARGET, guide_name, target, targets


def _target_of(guide_name: str, nt_control: str, target_by_name: dict[str, str]) -> str:
    target = target_by_name.get(guide_name, "").strip()
    if _is_nt(guide_name, target, nt_control):
        return nt_control
    return target


def _multi_subtype(parts: list[str], nt_control: str, target_by_name: dict[str, str]) -> str:
    nt_flags = []
    targets = []
    for part in parts:
        target = target_by_name.get(part, "")
        is_nt = _is_nt(part, target, nt_control)
        nt_flags.append(is_nt)
        if not is_nt and target:
            targets.append(target)
    if all(nt_flags):
        return MULTI_ALL_NT
    if any(nt_flags) and targets:
        return MULTI_MIXED_NT
    if len(set(targets)) == 1 and targets:
        return MULTI_SAME_TARGET
    return MULTI_DIFFERENT_TARGETS


def control_cell_mask(obs: pd.DataFrame, spec: dict[str, Any]) -> np.ndarray:
    """Cells used as the contrast. Designed NT and zero-guide cells stay separate.

    Temporary control cells are assignment ``no_guide`` with zero guides and
    ``qc_pass``. That covers an empty assignment call and a matrix barcode
    missing from the positive assignment table. ``multi_guide`` and
    ``single_targeting`` are not controls. A configured NT label does not
    recruit zero-guide cells. Cells excluded as ``low_library_no_guide`` stay
    out of both contrasts. A real ``single_nt`` cell is not removed by that rule.
    """
    n = len(obs)
    if n == 0:
        return np.zeros(0, dtype=bool)
    shallow = _low_library_no_guide_rows(obs)
    if spec.get("temporary"):
        if "assignment_status" in obs.columns and "n_guides" in obs.columns:
            assignment = obs["assignment_status"].astype(str)
            n_guides = pd.to_numeric(obs["n_guides"], errors="coerce").fillna(0)
            mask = _writable_bool(
                (assignment.eq(STATUS_NO_GUIDE) & n_guides.eq(0)).to_numpy()
            )
        else:
            label = str(spec.get("control_label", ""))
            mask = _writable_bool(
                _obs_eq(obs, "target_gene", label) | _obs_eq(obs, "perturbation_status", label)
            )
        qc_pass = _flag_column(obs, "qc_pass")
        if qc_pass is not None:
            mask &= qc_pass
        return mask & ~shallow
    label = str(spec.get("control_label", ""))
    designed = _writable_bool(
        _obs_eq(obs, "perturbation_status", STATUS_NT) | _obs_eq(obs, "target_gene", label)
    )
    if "assignment_status" in obs.columns:
        designed &= ~_obs_eq(obs, "assignment_status", STATUS_NO_GUIDE)
    return designed & ~(shallow & ~_single_nt_rows(obs))


def _writable_bool(values: np.ndarray) -> np.ndarray:
    """Copy boolean masks so in-place ``&=`` / ``|=`` works under pandas CoW."""
    return np.array(values, dtype=bool, copy=True)


def _flag_column(obs: pd.DataFrame, column: str) -> np.ndarray | None:
    if column not in obs.columns:
        return None
    values = obs[column]
    if pd.api.types.is_bool_dtype(values):
        return _writable_bool(values.fillna(False).to_numpy(dtype=bool))
    return _writable_bool(
        values.astype(str).str.strip().str.lower().isin(["true", "1"]).to_numpy()
    )


def _low_library_no_guide_rows(obs: pd.DataFrame) -> np.ndarray:
    if "exclusion_reason_primary" not in obs.columns:
        return np.zeros(len(obs), dtype=bool)
    return _writable_bool(
        obs["exclusion_reason_primary"].astype(str).eq(LOW_LIBRARY_NO_GUIDE).to_numpy()
    )


def _single_nt_rows(obs: pd.DataFrame) -> np.ndarray:
    status = _obs_eq(obs, "assignment_status", STATUS_SINGLE_NT)
    return status | _obs_eq(obs, "perturbation_status", STATUS_NT)


def _obs_eq(obs: pd.DataFrame, column: str, value: str) -> np.ndarray:
    if column not in obs.columns:
        return np.zeros(len(obs), dtype=bool)
    return _writable_bool(obs[column].astype(str).eq(value).to_numpy())


def _is_nt(guide_name: str, target_gene_name: str, nt_control: str) -> bool:
    label = nt_control.strip() if isinstance(nt_control, str) else ""
    if not label:
        return False
    return guide_name == label or target_gene_name == label


def multi_guide_composition(
    annotated: pd.DataFrame,
    *,
    nt_control: str,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    """Classify true multi-target calls before they are dropped.

    Same-gene dual guides are ``single_targeting`` and are not counted here.
    """
    catalog = validate_guide_targets(reference, nt_control=nt_control)
    target_by_name = target_by_guide_name(catalog)
    names = set(catalog["guide_name"])
    max_tokens = name_token_width(reference)
    counts = {MULTI_ALL_NT: 0, MULTI_SAME_TARGET: 0, MULTI_MIXED_NT: 0, MULTI_DIFFERENT_TARGETS: 0}
    multi = annotated[annotated["assignment_status"].astype(str).eq(STATUS_MULTI)]
    for record in multi.itertuples(index=False):
        parts = segment_feature_call(str(record.feature_call), names, max_tokens)
        if len(parts) < 2 or any(part not in names for part in parts):
            continue
        targets = []
        nt_flags = []
        for part in parts:
            target = target_by_name.get(part, "")
            is_nt = _is_nt(part, target, nt_control)
            nt_flags.append(is_nt)
            if not is_nt and target:
                targets.append(target)
        if all(nt_flags):
            counts[MULTI_ALL_NT] += 1
        elif any(nt_flags) and targets:
            counts[MULTI_MIXED_NT] += 1
        elif len(set(targets)) == 1 and targets:
            counts[MULTI_SAME_TARGET] += 1
        elif len(set(targets)) > 1:
            counts[MULTI_DIFFERENT_TARGETS] += 1
    return pd.DataFrame({"composition": list(counts), "n_cells": list(counts.values())})


def mito_fraction(mito_counts: np.ndarray, total_counts: np.ndarray) -> np.ndarray:
    """Share of library size on a 0–1 scale.

    A library size of 0 has no denominator, so the value is missing rather than 0.
    """
    part = np.asarray(mito_counts, dtype=np.float64).ravel()
    total = np.asarray(total_counts, dtype=np.float64).ravel()
    fraction = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(part, total, out=fraction, where=total > 0)
    return fraction


def capture_efficiency(guide_umi: np.ndarray, total_counts: np.ndarray) -> np.ndarray:
    """CRISPR guide UMI sum divided by Gene Expression library size.

    A guide-UMI sum of 0 is efficiency 0, including when the library size is 0.
    A positive guide-UMI sum with library size 0 has no denominator and is missing.
    """
    guide = np.asarray(guide_umi, dtype=np.float64).ravel()
    total = np.asarray(total_counts, dtype=np.float64).ravel()
    efficiency = np.full(guide.shape, np.nan, dtype=np.float64)
    finite_guide = np.isfinite(guide)
    efficiency[finite_guide & (guide == 0)] = 0.0
    positive = finite_guide & (guide > 0) & np.isfinite(total) & (total > 0)
    np.divide(guide, total, out=efficiency, where=positive)
    return efficiency


def capture_library_correlation(table: pd.DataFrame) -> pd.DataFrame:
    """Spearman correlation of capture efficiency and library size before QC pass.

    The p-value stays in this QC table. It is not an effect-table p-value.
    """
    row: dict[str, Any] = {
        "x": "total_counts",
        "y": "grna_capture_efficiency",
        "population": "before_qc_pass",
        "n_cells_input": int(len(table)),
        "n_cells": 0,
        "spearman_r": np.nan,
        "spearman_pvalue": np.nan,
        "status": "skipped",
        "reason": "",
    }
    needed = {"total_counts", "grna_capture_efficiency"}
    if not needed.issubset(table.columns):
        row["reason"] = "missing_capture_columns"
        return pd.DataFrame([row])
    library = pd.to_numeric(table["total_counts"], errors="coerce").to_numpy(dtype=float)
    efficiency = pd.to_numeric(table["grna_capture_efficiency"], errors="coerce").to_numpy(
        dtype=float
    )
    ok = np.isfinite(library) & np.isfinite(efficiency)
    row["n_cells"] = int(ok.sum())
    if int(ok.sum()) < 2:
        row["reason"] = "fewer_than_2_cells"
        return pd.DataFrame([row])
    if float(np.var(library[ok])) == 0.0 or float(np.var(efficiency[ok])) == 0.0:
        row["reason"] = "zero_variance"
        return pd.DataFrame([row])
    from scipy.stats import spearmanr

    result = spearmanr(library[ok], efficiency[ok])
    statistic = getattr(result, "statistic", None)
    if statistic is None:
        statistic = result.correlation
    row["spearman_r"] = float(statistic)
    row["spearman_pvalue"] = float(result.pvalue)
    row["status"] = "computed"
    return pd.DataFrame([row])


def derive_qc_thresholds(qc: dict[str, Any]) -> dict[str, float]:
    """Resolve QC cutoffs from config / CLI: gene bounds and mito fraction only.

    Library UMI / nCount (``total_counts``) is computed for metrics and plots
    but is not used as an absolute cell filter.
    """
    if qc.get("min_genes") is None or qc.get("max_genes") is None:
        raise ValueError("缺少 qc.min_genes 或 qc.max_genes。")
    if qc.get("max_mito_frac") is None:
        raise ValueError("缺少 qc.max_mito_frac。")
    return {
        "min_genes": float(qc["min_genes"]),
        "max_genes": float(qc["max_genes"]),
        "max_mito_frac": float(qc["max_mito_frac"]),
    }


def _low_library_no_guide_mask(
    table: pd.DataFrame,
    total_counts: np.ndarray,
    rule: dict[str, Any] | None,
) -> np.ndarray:
    """Shallow cells with no detected guide. Cells that have a guide are left alone."""
    empty = np.zeros(len(table), dtype=bool)
    if not rule or not rule.get("applied") or rule.get("min_total_counts") is None:
        return empty
    if "assignment_status" not in table.columns or "n_guides" not in table.columns:
        return empty
    status = table["assignment_status"].astype(str).to_numpy()
    n_guides = pd.to_numeric(table["n_guides"], errors="coerce").fillna(0).to_numpy(dtype=float)
    no_guide = (status == STATUS_NO_GUIDE) & (n_guides == 0)
    floor = float(rule["min_total_counts"])
    return no_guide & np.isfinite(total_counts) & (total_counts < floor)


def apply_eligibility(
    annotated: pd.DataFrame,
    thresholds: dict[str, float],
    *,
    min_cells_per_target: int,
    cell_type: str = "ALL",
    targets_per_cell: int = 1,
    low_library_rule: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Attach QC, every exclusion reason, and one primary reason per cell.

    ``low_library_rule`` is kept for the filter-stage ledger; pass ``None`` or
    ``applied=False`` so that stage adds zero exclusions. Retention uses unique
    target-gene count (``n_targets``), not guide count.
    """
    table = annotated.copy()
    for column in ("total_counts", "n_genes_by_counts", "mito_frac"):
        if column not in table.columns:
            raise ValueError(f"缺少 {column}，不能把未计算的 QC 默认为通过。")
        table[column] = pd.to_numeric(table[column], errors="coerce")
    total = table["total_counts"].to_numpy(dtype=float)
    genes = table["n_genes_by_counts"].to_numpy(dtype=float)
    mito = table["mito_frac"].to_numpy(dtype=float)
    finite = np.isfinite(total) & np.isfinite(genes) & np.isfinite(mito)
    if "max_genes" not in thresholds:
        raise ValueError("缺少 max_genes，不能把未计算的基因数 QC 默认为通过。")
    # Drop cells with ≤ min_genes or ≥ max_genes expressed genes.
    low_genes = finite & (genes <= float(thresholds["min_genes"]))
    high_genes = finite & (genes >= float(thresholds["max_genes"]))
    # Drop cells with mitochondrial read fraction ≥ max_mito_frac (default 0.2).
    high_mito = finite & (mito >= float(thresholds["max_mito_frac"]))
    low_library = _low_library_no_guide_mask(table, total, low_library_rule)
    qc_pass = finite & ~low_genes & ~high_genes & ~high_mito & ~low_library
    qc_reasons = []
    for index in range(len(table)):
        reasons = []
        if not finite[index]:
            reasons.append("qc_not_computable")
        if low_genes[index]:
            reasons.append("low_genes")
        if high_genes[index]:
            reasons.append("high_genes")
        if high_mito[index]:
            reasons.append("high_mito")
        if low_library[index]:
            reasons.append(LOW_LIBRARY_NO_GUIDE)
        qc_reasons.append(json.dumps(reasons))
    table["qc_pass"] = qc_pass
    table["qc_fail_reasons_json"] = qc_reasons
    table["cell_type"] = cell_type
    if "n_targets" not in table.columns:
        raise ValueError("缺少 n_targets，不能按靶基因数判定资格。")
    single = table["assignment_status"].isin([STATUS_SINGLE_TARGET, STATUS_SINGLE_NT])
    allowed_count = table["n_targets"].astype(int) == int(targets_per_cell)
    base = single & allowed_count & table["guide_umi_pass"].astype(bool) & table["qc_pass"]
    table["base_eligible"] = base.to_numpy()
    targeting_base = base & (table["assignment_status"] == STATUS_SINGLE_TARGET)
    grouped = table.loc[targeting_base].groupby(["target_gene", "cell_type"], dropna=False).size()
    low_keys = {key for key, count in grouped.items() if int(count) < min_cells_per_target}
    pair = pd.MultiIndex.from_arrays(
        [table["target_gene"].astype(str), table["cell_type"].astype(str)]
    )
    low_cells = targeting_base.to_numpy() & pair.isin(low_keys)
    all_reasons = []
    primary = []
    for index, record in enumerate(table.itertuples(index=False)):
        reasons = []
        if low_library[index]:
            reasons.append(LOW_LIBRARY_NO_GUIDE)
        if record.assignment_status == STATUS_NO_GUIDE:
            reasons.append("no_guide")
        if record.assignment_status == STATUS_MULTI or int(record.n_targets) > int(
            targets_per_cell
        ):
            reasons.append("multi_guide")
        if not bool(record.guide_umi_pass):
            reasons.append("low_guide_umi")
        if not bool(record.qc_pass):
            reasons.append("qc_failed")
        if low_cells[index]:
            reasons.append("low_cells_target")
        chosen = STATUS_RETAINED
        for reason in EXCLUSION_PRIORITY:
            if reason in reasons:
                chosen = reason
                break
        all_reasons.append(json.dumps(reasons))
        primary.append(chosen)
    table["exclusion_reasons_json"] = all_reasons
    table["exclusion_reason_primary"] = primary
    table["analysis_eligible"] = table["exclusion_reason_primary"] == STATUS_RETAINED
    return table


def filter_stage_counts(table: pd.DataFrame) -> pd.DataFrame:
    """Newly excluded cells at each priority stage. Earlier stages are already gone."""
    remaining = len(table)
    rows = [
        {
            "stage": "input",
            "n_before": remaining,
            "n_newly_excluded": 0,
            "n_after": remaining,
            "primary_reason": "",
        }
    ]
    counts = table["exclusion_reason_primary"].value_counts()
    for reason in EXCLUSION_PRIORITY:
        removed = int(counts.get(reason, 0))
        after = remaining - removed
        rows.append(
            {
                "stage": reason,
                "n_before": remaining,
                "n_newly_excluded": removed,
                "n_after": after,
                "primary_reason": reason,
            }
        )
        remaining = after
    retained = int((table["exclusion_reason_primary"] == STATUS_RETAINED).sum())
    if remaining != retained:
        raise RuntimeError(f"阶段计数 {remaining} 与 retained {retained} 不一致。")
    return pd.DataFrame(rows)


def disambiguate_shared_guides(
    annotated: pd.DataFrame,
    guide_umi: sparse.spmatrix,
    barcodes: pd.Index,
    crispr_ids: list[str],
    crispr_names: list[str],
) -> pd.DataFrame:
    """Assign a reference id when several guides share one name.

    The shared name comes from the reference table. The CRISPR count matrix
    picks the id with the most UMIs. There is no fixed non-targeting name.
    """
    table = annotated.copy()
    if guide_umi.shape[1] == 0:
        table["guide_umi_observed"] = np.nan
        return table
    positions = pd.Index(barcodes.astype(str)).get_indexer(table["cell_barcode"].astype(str))
    names = np.asarray(crispr_names)
    ids = np.asarray(crispr_ids)
    columns_by_name = {
        name: np.flatnonzero(names == name)
        for name in set(table.loc[table["shared_guide_name"], "guide_name"].astype(str))
    }
    id_to_column = {guide_id: index for index, guide_id in enumerate(crispr_ids)}
    name_to_columns = {name: np.flatnonzero(names == name) for name in set(names.tolist()) if name}
    resolved: list[str] = []
    umis: list[float] = []
    for row_index, record in enumerate(table.itertuples(index=False)):
        matrix_row = int(positions[row_index])
        status = str(record.perturbation_status)
        if matrix_row < 0 or status not in {STATUS_NT, STATUS_TARGET}:
            resolved.append(str(record.guide_id))
            umis.append(np.nan)
            continue
        parts = _record_guide_names(record)
        if len(parts) > 1:
            columns = np.concatenate(
                [name_to_columns.get(part, np.array([], dtype=int)) for part in parts]
            )
            if columns.size == 0:
                resolved.append("")
                umis.append(0.0)
                continue
            values = _row_values(guide_umi, matrix_row, columns)
            if values.size == 0 or float(np.nanmax(values)) <= 0:
                resolved.append("")
                umis.append(0.0)
                continue
            winner = int(np.argmax(values))
            resolved.append(str(ids[columns[winner]]))
            umis.append(float(values[winner]))
            continue
        if bool(record.shared_guide_name):
            columns = columns_by_name.get(str(record.guide_name), np.array([], dtype=int))
            values = _row_values(guide_umi, matrix_row, columns)
            if values.size == 0 or float(np.nanmax(values)) <= 0:
                resolved.append("")
                umis.append(0.0)
                continue
            winner = int(np.argmax(values))
            resolved.append(str(ids[columns[winner]]))
            umis.append(float(values[winner]))
            continue
        column = id_to_column.get(str(record.guide_id))
        resolved.append(str(record.guide_id))
        if column is None:
            umis.append(np.nan)
            continue
        values = _row_values(guide_umi, matrix_row, np.asarray([column]))
        umis.append(float(values[0]) if values.size else np.nan)
    table["guide_id"] = resolved
    table["guide_umi_observed"] = umis
    return table


def _record_guide_names(record: Any) -> list[str]:
    names = getattr(record, "guide_names", None)
    if names is None:
        name = str(getattr(record, "guide_name", "") or "")
        return [name] if name else []
    if isinstance(names, str):
        return [names] if names else []
    return [str(name) for name in list(names) if str(name)]


def _row_values(matrix: sparse.spmatrix, row: int, columns: np.ndarray) -> np.ndarray:
    if columns.size == 0:
        return np.array([], dtype=float)
    block = matrix[row, columns]
    if sparse.issparse(block):
        return np.asarray(block.toarray()).ravel()
    return np.asarray(block).ravel()


def cells_per_guide(table: pd.DataFrame) -> pd.DataFrame:
    used = _analysis_rows(table)
    if used.empty:
        return pd.DataFrame(columns=["guide_id", "guide_name", "target_gene", "n_cells"])
    grouped = (
        used.groupby(["guide_id", "guide_name", "target_gene"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    return grouped


def _analysis_rows(table: pd.DataFrame) -> pd.DataFrame:
    if "analysis_eligible" in table.columns:
        return table[table["analysis_eligible"].astype(bool)].copy()
    return table[table["perturbation_status"].isin([STATUS_TARGET, STATUS_NT])].copy()


def cells_per_target(table: pd.DataFrame, min_cells_per_guide: int) -> pd.DataFrame:
    used = _analysis_rows(table)
    if used.empty:
        return pd.DataFrame(
            columns=[
                "target_gene",
                "n_cells",
                "n_guides",
                "n_guides_passing",
                "perturbation_status",
            ]
        )
    rows = []
    for target, group in used.groupby("target_gene", sort=False):
        guide_counts = group.groupby("guide_id").size()
        status = group["perturbation_status"].iloc[0]
        rows.append(
            {
                "target_gene": target,
                "n_cells": int(len(group)),
                "n_guides": int((guide_counts.index != "").sum()),
                "n_guides_passing": int(
                    (guide_counts[guide_counts.index != ""] >= min_cells_per_guide).sum()
                ),
                "perturbation_status": status,
            }
        )
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def write_qc_figures(
    output_dir: Path,
    config: dict[str, Any],
    adata: AnnData | None = None,
    *,
    cell_status: pd.DataFrame | None = None,
    targets: pd.DataFrame | None = None,
) -> None:
    """Write QC / assignment PNGs at the end of the assign step."""
    from perturbseq.config import control_spec
    from perturbseq.plotting import (
        draw_cells_per_target,
        draw_guides_per_cell,
        draw_qc_metrics_panel,
        draw_qc_scatter,
        draw_scatter_ncount_nfeature_gex,
        draw_temporary_control_qc,
        open_figure_store,
        remove_qc_violin_products,
        remove_temporary_control_qc,
    )

    output_dir = Path(output_dir)
    remove_qc_violin_products(output_dir)
    store = open_figure_store(output_dir, config, adata)
    if cell_status is None:
        ledger_path = output_dir / "cell_status.tsv"
        cell_status = (
            pd.read_csv(ledger_path, sep="\t") if ledger_path.exists() else pd.DataFrame()
        )
    if not cell_status.empty and "n_guides" in cell_status.columns:
        draw_guides_per_cell(store, cell_status["n_guides"])
    else:
        store.skip("guides_per_cell", "missing_fields")
    if targets is None:
        targets_path = output_dir / "qc" / "tables" / "cells_per_target.tsv"
        targets = (
            pd.read_csv(targets_path, sep="\t") if targets_path.exists() else pd.DataFrame()
        )
    draw_cells_per_target(store, targets, int(config["qc"]["min_cells_per_target"]))
    qc_table, qc_note = _qc_population(cell_status, adata)
    if qc_table.empty:
        for name in (
            "qc_metrics_panel",
            "qc_scatter",
            "scatter_ncount_nfeature_gex",
        ):
            store.skip(name, "missing_fields")
    else:
        draw_qc_metrics_panel(store, qc_table, config, population=qc_note)
        draw_qc_scatter(store, qc_table, int(config["runtime"]["seed"]), population=qc_note)
        draw_scatter_ncount_nfeature_gex(
            store, qc_table, int(config["runtime"]["seed"]), population=qc_note
        )
    if control_spec(config)["temporary"]:
        draw_temporary_control_qc(store, cell_status)
    else:
        remove_temporary_control_qc(output_dir)
    # plot_manifest is ephemeral; lean cleanup removes it after the run.


def _qc_population(cell_status: pd.DataFrame, adata: AnnData | None) -> tuple[pd.DataFrame, str]:
    """Prefer the full annotated ledger, which is recorded before analysis filtering."""
    mito = {"mito_frac", "pct_counts_mito"}
    needed = {"total_counts", "n_genes_by_counts"}
    if (
        not cell_status.empty
        and needed.issubset(cell_status.columns)
        and mito & set(cell_status.columns)
    ):
        return cell_status, "before filtering; full annotated cells"
    if adata is not None:
        frame = adata.obs.reset_index(drop=False)
        if needed.issubset(frame.columns) and mito & set(frame.columns):
            return frame, "analysis object; full annotated QC columns were absent"
    return pd.DataFrame(), ""
