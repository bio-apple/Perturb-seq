"""Guide call parsing. Counts follow num_features, not the number of '|' characters."""

from __future__ import annotations

import numpy as np
import pandas as pd

NT_PATTERNS = ("negctrl", "non-targeting", "nontargeting", "ntc")
INTRON_PATTERN = "intronctrl"
NT_CLASS = "非靶向对照"
PERT_CLASS = "受扰动细胞"
INTRON_CLASS = "内含子对照"
CONTROL_KEY = "NT"
TARGET_COLUMNS = ("target_gene_id", "target_gene_name")
MIN_CONTROL_CELLS = 20
MIN_PERTURBATION_CELLS = 10


def has_target_columns(reference: pd.DataFrame) -> bool:
    return set(TARGET_COLUMNS).issubset(reference.columns)


def _hit(value: str, pattern: str) -> bool:
    return pattern in value.lower()


def is_nontargeting(feature_call: str) -> bool:
    return any(_hit(feature_call, pattern) for pattern in NT_PATTERNS)


def is_intron(feature_call: str) -> bool:
    return _hit(feature_call, INTRON_PATTERN)


def segment_feature_call(feature_call: str, names_longest_first: list[str]) -> list[str]:
    """Split a DRAGEN feature_call into reference `name`s.

    Guide names themselves may contain '|', and '|' is also the separator
    between guides. Match the longest reference name first.
    """
    if feature_call == "":
        return []
    parts: list[str] = []
    remaining = feature_call
    while remaining:
        if remaining.startswith("|"):
            remaining = remaining[1:]
            continue
        match = next(
            (name for name in names_longest_first if remaining == name or remaining.startswith(name + "|")),
            None,
        )
        if match is None:
            parts.append(remaining)
            break
        parts.append(match)
        remaining = remaining[len(match) :]
    return parts


def _lookup_tables(reference: pd.DataFrame) -> tuple[dict[str, int], dict[str, int], list[str]]:
    id_row = {str(row.id): i for i, row in enumerate(reference.itertuples(index=False))}
    name_row: dict[str, int] = {}
    for i, name in enumerate(reference["name"].astype(str)):
        # Shared names (the non-targeting guides) carry the same empty target.
        name_row.setdefault(name, i)
    names = sorted(set(reference["name"].astype(str)), key=len, reverse=True)
    return id_row, name_row, names


def _target_from_row(reference: pd.DataFrame, row_index: int, modeling: bool) -> tuple[str, str]:
    if not modeling:
        return "", ""
    gene_id = str(reference.iloc[row_index]["target_gene_id"]).strip()
    gene_name = str(reference.iloc[row_index]["target_gene_name"]).strip()
    if gene_id == "" or gene_name == "":
        return "", ""
    return gene_id, gene_name


def annotate_assignments(assignments: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    """Attach guide name, class, and target gene to every assignment row."""
    modeling = has_target_columns(reference)
    id_row, name_row, names = _lookup_tables(reference)
    records: list[dict[str, object]] = []
    for row in assignments.itertuples(index=False):
        call = str(row.feature_call)
        parts = segment_feature_call(call, names)
        transcripts = _transcripts(str(row.num_transcripts), len(parts))
        guide = parts[0] if int(row.num_features) == 1 and parts else ""
        row_index = None
        if guide:
            row_index = id_row.get(guide)
            if row_index is None:
                row_index = name_row.get(guide)
        gene_id, gene_name = ("", "")
        if row_index is not None:
            gene_id, gene_name = _target_from_row(reference, row_index, modeling)
        icm = ""
        perturbation = ""
        if int(row.num_features) == 1 and guide:
            if is_nontargeting(call):
                icm = NT_CLASS
                perturbation = CONTROL_KEY
                gene_id, gene_name = "", ""
            elif is_intron(call):
                icm = INTRON_CLASS
                gene_id, gene_name = "", ""
            else:
                icm = PERT_CLASS
                perturbation = gene_name
        records.append(
            {
                "cell_barcode": row.cell_barcode,
                "num_features": int(row.num_features),
                "feature_call": call,
                "guide": guide,
                "guide_parts": parts,
                "guide_transcripts": transcripts,
                "target_gene_id": gene_id,
                "target_gene_name": gene_name,
                "icm_class": icm,
                "perturbation": perturbation,
            }
        )
    return pd.DataFrame.from_records(records)


def _transcripts(raw: str, n_parts: int) -> list[float]:
    if raw == "" or n_parts == 0:
        return []
    fields = raw.split("|")
    values = [float(field) if field else 0.0 for field in fields]
    if len(values) != n_parts:
        return values
    return values


def guide_target_mapping(reference: pd.DataFrame) -> pd.DataFrame:
    """One reference row per guide. A target gene requires both columns to be non-empty."""
    modeling = has_target_columns(reference)
    rows = []
    for record in reference.itertuples(index=False):
        gene_id = ""
        gene_name = ""
        if modeling:
            gene_id = str(record.target_gene_id).strip()
            gene_name = str(record.target_gene_name).strip()
            if gene_id == "" or gene_name == "":
                gene_id, gene_name = "", ""
        rows.append(
            {
                "guide_id": str(record.id),
                "guide": str(record.name),
                "target_gene_id": gene_id,
                "target_gene_name": gene_name,
            }
        )
    return pd.DataFrame(rows, columns=["guide_id", "guide", "target_gene_id", "target_gene_name"])


def unmatched_records(
    annotated: pd.DataFrame,
    reference: pd.DataFrame,
    matrix_barcodes: pd.Index,
) -> pd.DataFrame:
    """Barcodes and guides that do not meet on both sides of the sample."""
    columns = ["record_type", "cell_barcode", "guide_id", "guide", "detail"]
    rows: list[dict[str, str]] = []
    assigned = set(annotated["cell_barcode"].astype(str))
    matrix = [str(barcode) for barcode in matrix_barcodes]
    matrix_set = set(matrix)
    for barcode in matrix:
        if barcode not in assigned:
            rows.append(
                {
                    "record_type": "barcode_not_in_assignment",
                    "cell_barcode": barcode,
                    "guide_id": "",
                    "guide": "",
                    "detail": "filtered matrix barcode has no assignment row",
                }
            )
    for barcode in annotated["cell_barcode"].astype(str):
        if barcode not in matrix_set:
            rows.append(
                {
                    "record_type": "barcode_not_in_matrix",
                    "cell_barcode": barcode,
                    "guide_id": "",
                    "guide": "",
                    "detail": "assignment barcode is absent from the filtered matrix",
                }
            )
    known = set(reference["id"].astype(str)) | set(reference["name"].astype(str))
    seen: set[str] = set()
    for record in annotated.itertuples(index=False):
        parts = [str(part) for part in record.guide_parts]
        seen.update(parts)
        for part in parts:
            if part not in known:
                rows.append(
                    {
                        "record_type": "guide_not_in_reference",
                        "cell_barcode": str(record.cell_barcode),
                        "guide_id": "",
                        "guide": part,
                        "detail": "feature_call fragment matched neither id nor name",
                    }
                )
        n_called = int(record.num_features)
        if (n_called == 0 and str(record.feature_call)) or (n_called > 0 and len(parts) != n_called):
            rows.append(
                {
                    "record_type": "guide_not_in_reference",
                    "cell_barcode": str(record.cell_barcode),
                    "guide_id": "",
                    "guide": str(record.feature_call),
                    "detail": f"segmented {len(parts)} names but num_features is {n_called}",
                }
            )
    for record in reference.itertuples(index=False):
        guide_id = str(record.id)
        guide = str(record.name)
        if guide_id not in seen and guide not in seen:
            rows.append(
                {
                    "record_type": "guide_not_in_assignment",
                    "cell_barcode": "",
                    "guide_id": guide_id,
                    "guide": guide,
                    "detail": "reference guide was not called in any passing cell",
                }
            )
    return pd.DataFrame(rows, columns=columns)


def assignment_summary(
    annotated: pd.DataFrame,
    reference: pd.DataFrame,
    unmatched: pd.DataFrame,
    matrix_barcodes: pd.Index,
) -> pd.DataFrame:
    """Cell classes from num_features, plus counts of unmatched records."""
    known = set(reference["id"].astype(str)) | set(reference["name"].astype(str))
    classes = [_assignment_class(record, known) for record in annotated.itertuples(index=False)]
    counts = pd.Series(classes).value_counts()
    assigned = set(annotated["cell_barcode"].astype(str))
    matrix_set = {str(barcode) for barcode in matrix_barcodes}
    unmatched_counts = unmatched["record_type"].value_counts() if len(unmatched) else pd.Series(dtype=int)
    items = [
        ("no_guide", int(counts.get("no_guide", 0)), "cells"),
        ("single_guide", int(counts.get("single_guide", 0)), "cells"),
        ("multiple_guides", int(counts.get("multiple_guides", 0)), "cells"),
        ("ambiguous", int(counts.get("ambiguous", 0)), "cells"),
        ("passing_cells", int(len(annotated)), "cells"),
        ("analysis_cells", int(len(assigned & matrix_set)), "cells"),
        ("barcode_not_in_assignment", int(unmatched_counts.get("barcode_not_in_assignment", 0)), "records"),
        ("barcode_not_in_matrix", int(unmatched_counts.get("barcode_not_in_matrix", 0)), "records"),
        ("guide_not_in_reference", int(unmatched_counts.get("guide_not_in_reference", 0)), "records"),
        ("guide_not_in_assignment", int(unmatched_counts.get("guide_not_in_assignment", 0)), "records"),
    ]
    return pd.DataFrame(items, columns=["item", "n", "unit"])


def _assignment_class(record, known: set[str]) -> str:
    parts = [str(part) for part in record.guide_parts]
    n_called = int(record.num_features)
    unmatched = any(part not in known for part in parts) or len(parts) != n_called
    if n_called == 0:
        return "ambiguous" if str(record.feature_call) else "no_guide"
    if unmatched:
        return "ambiguous"
    if n_called == 1:
        return "single_guide"
    return "multiple_guides"


def guide_frequency(annotated: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """num_features histogram, and per-guide cell counts plus transcript sums.

    The same guide name contributes at most one cell. Transcript counts of
    repeated names in that cell are still added.
    """
    n_features = (
        annotated["num_features"]
        .value_counts()
        .rename_axis("num_features")
        .reset_index(name="n_cells")
        .sort_values("num_features")
    )
    cells: dict[str, int] = {}
    transcripts: dict[str, float] = {}
    for parts, weights in zip(annotated["guide_parts"], annotated["guide_transcripts"], strict=True):
        paired = list(zip(parts, weights, strict=False))
        for name in dict.fromkeys(name for name, _weight in paired):
            cells[name] = cells.get(name, 0) + 1
            transcripts[name] = transcripts.get(name, 0.0) + sum(weight for guide, weight in paired if guide == name)
    frequency = pd.DataFrame(
        {
            "guide": list(cells),
            "n_cells": [cells[name] for name in cells],
            "transcript_sum": [transcripts.get(name, 0.0) for name in cells],
        }
    ).sort_values(["transcript_sum", "n_cells"], ascending=False)
    return n_features, frequency


def singlet_mask(annotated: pd.DataFrame) -> np.ndarray:
    return annotated["num_features"].to_numpy() == 1
