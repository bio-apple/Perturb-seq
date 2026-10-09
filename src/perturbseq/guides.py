"""Guide parsing, NT/control identification, and guide-level QC.

Guide assignment / parsing lives here. Guide-level efficacy QC is implemented in
``guide_qc`` and re-exported for a single import surface.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbseq.guide_qc import (
    compute_gene_guide_consistency,
    compute_guide_metrics,
    run_guide_qc,
    summarize_guide_qc,
    write_guide_qc_tables,
)

DEFAULT_CONTROL_PATTERNS = (
    r"^nt$",
    r"^ntc$",
    r"^nt_",
    r"^negctrl",
    r"^neg_ctrl",
    r"non[-_ ]?targeting",
    r"^control$",
    r"^ctrl$",
)

_GUIDE_INDEX_SUFFIX = re.compile(r"^([A-Za-z0-9.-]+)_\d+$")


def _matches_control(value: str, patterns: tuple[str, ...]) -> bool:
    text = str(value).strip().lower()
    if not text:
        return False
    return any(re.search(pattern, text) for pattern in patterns)


def parse_gene_target(feature_name: str, control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS) -> str:
    name = str(feature_name).strip()
    if not name or name.lower() in {"nan", "none"}:
        return "unassigned"
    if _matches_control(name, control_patterns):
        return "NT"
    if re.search(r"\|design_\d+$", name, flags=re.IGNORECASE):
        return name.split("|", 1)[0]
    if "/" in name:
        return name.split("/", 1)[0]
    match = _GUIDE_INDEX_SUFFIX.match(name)
    if match:
        gene = match.group(1)
        return "NT" if _matches_control(gene, control_patterns) else gene
    return name


def _guide_umi_total(num_transcripts: str) -> int:
    if not num_transcripts:
        return 0
    values = []
    for token in str(num_transcripts).split("|"):
        token = token.strip()
        if not token:
            continue
        values.append(int(float(token)))
    return int(sum(values))


def split_feature_calls(feature_call: str, known_names: set[str] | None = None) -> list[str]:
    """Split DRAGEN feature_call values.

    Guide *names* may themselves contain ``|`` (e.g. ``NEK4|design_3``), so a naive
    split is wrong. Prefer longest-prefix matches against the feature reference.
    """
    call = str(feature_call).strip()
    if not call:
        return []
    known = {str(name) for name in (known_names or set()) if str(name)}
    if call in known:
        return [call]
    if not known:
        if re.search(r"\|design_\d+(\||$)", call, flags=re.IGNORECASE):
            return re.findall(r"(?:Non-Targeting_Control|[A-Za-z0-9.-]+\|design_\d+)", call)
        return [token for token in call.split("|") if token]
    names = sorted(known, key=len, reverse=True)
    remaining = call
    parts: list[str] = []
    while remaining:
        if remaining.startswith("|"):
            remaining = remaining[1:]
            continue
        match = next((name for name in names if remaining == name or remaining.startswith(name + "|")), None)
        if match is None:
            parts.append(remaining)
            break
        parts.append(match)
        remaining = remaining[len(match) :]
    return parts


def build_feature_target_map(
    feature_ref: pd.DataFrame,
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS,
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    gene_col = next(
        (col for col in ("target_gene", "gene_target", "gene", "target_gene_name") if col in feature_ref.columns),
        None,
    )
    for _, row in feature_ref.iterrows():
        names = {str(row["id"]), str(row["name"])}
        if gene_col and pd.notna(row[gene_col]) and str(row[gene_col]).strip():
            target = str(row[gene_col]).strip()
            if _matches_control(target, control_patterns):
                target = "NT"
        else:
            target = parse_gene_target(row["name"], control_patterns)
        for name in names:
            mapping[name] = target
    return mapping


def annotate_guides(
    adata: AnnData,
    assignments: pd.DataFrame,
    feature_ref: pd.DataFrame,
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS,
) -> AnnData:
    """Join DRAGEN guide calls onto gene-expression cells."""
    mapping = build_feature_target_map(feature_ref, control_patterns)
    known_names = set(mapping)
    joined = adata.obs.join(assignments, how="left")
    joined["num_features"] = joined["num_features"].fillna(0).astype(int)
    joined["feature_call"] = joined["feature_call"].fillna("")
    joined["num_transcripts"] = joined["num_transcripts"].fillna("")
    split_calls = [split_feature_calls(call, known_names) for call in joined["feature_call"]]
    joined["guide_id"] = [guides[0] if guides else "" for guides in split_calls]
    joined["guide_umi"] = joined["num_transcripts"].map(_guide_umi_total)
    targets = []
    for guides in split_calls:
        if not guides:
            targets.append("unassigned")
            continue
        gene_targets = [mapping.get(guide, parse_gene_target(guide, control_patterns)) for guide in guides]
        unique = list(dict.fromkeys(gene_targets))
        targets.append("|".join(unique) if unique else "unassigned")
    joined["gene_target"] = targets
    joined.loc[joined["num_features"] == 0, "gene_target"] = "unassigned"
    nt_mask = joined["gene_target"].map(lambda value: _matches_control(value, control_patterns) or value == "NT")
    joined["perturbation"] = np.where(joined["gene_target"].eq("unassigned"), "unassigned", "perturbed")
    joined.loc[nt_mask, ["gene_target", "perturbation"]] = ["NT", "NT"]
    adata.obs = joined
    return adata


def filter_singlets(adata: AnnData, singlet_only: bool = True) -> AnnData:
    if not singlet_only:
        return adata
    return adata[adata.obs["num_features"] == 1].copy()


__all__ = [
    "DEFAULT_CONTROL_PATTERNS",
    "annotate_guides",
    "build_feature_target_map",
    "compute_gene_guide_consistency",
    "compute_guide_metrics",
    "filter_singlets",
    "parse_gene_target",
    "run_guide_qc",
    "split_feature_calls",
    "summarize_guide_qc",
    "write_guide_qc_tables",
]
