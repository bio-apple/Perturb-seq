"""Unit tests for sample_type → cell-annotation policy (no scanpy runtime needed)."""

from __future__ import annotations

import pytest

from perturbseq.annotation_policy import SAMPLE_TYPES, resolve_cell_annotation_decision


@pytest.mark.parametrize(
    "sample_type,expect_skip,reason_prefix",
    [
        ("cell_line", True, "sample_type=cell_line"),
        ("primary", False, "sample_type=primary"),
        ("mixed", False, "sample_type=mixed"),
        ("unknown", True, "sample_type=unknown"),
    ],
)
def test_sample_type_policy(sample_type: str, expect_skip: bool, reason_prefix: str) -> None:
    d = resolve_cell_annotation_decision(sample_type)
    assert d["skip"] is expect_skip
    assert d["reason"] == reason_prefix
    assert d["source"] == "sample_type"
    assert d["sample_type"] == sample_type


def test_skip_override_wins() -> None:
    d = resolve_cell_annotation_decision("primary", skip_cell_annotation=True)
    assert d["skip"] is True
    assert d["reason"] == "user_skip"
    assert d["source"] == "override"


def test_run_override_wins() -> None:
    d = resolve_cell_annotation_decision("cell_line", skip_cell_annotation=False)
    assert d["skip"] is False
    assert d["reason"] == "user_run"
    assert d["source"] == "override"


def test_invalid_sample_type() -> None:
    with pytest.raises(ValueError, match="sample_type"):
        resolve_cell_annotation_decision("organoid")


def test_sample_types_constant() -> None:
    assert set(SAMPLE_TYPES) == {"cell_line", "primary", "mixed", "unknown"}
