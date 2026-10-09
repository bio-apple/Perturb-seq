"""Machine-safe analysis status blocks in report.json."""

from __future__ import annotations

import json

from perturbseq.cli import main
from perturbseq.report import (
    ANALYSIS_STATUS_KEYS,
    completed_status,
    skipped_status,
    status_note,
)


def test_skipped_status_shape():
    block = skipped_status("user_skip", detail="Mixscape skipped (--skip-mixscape).")
    assert block == {
        "skipped": True,
        "reason": "user_skip",
        "detail": "Mixscape skipped (--skip-mixscape).",
    }
    assert "user_skip" in (status_note(block) or "")


def test_completed_status_shape():
    block = completed_status(n_contrasts=3, scope="gene_target")
    assert block["skipped"] is False
    assert block["reason"] is None
    assert block["n_contrasts"] == 3


def test_pipeline_report_has_explicit_skip_blocks(tmp_path):
    demo = tmp_path / "demo"
    out = tmp_path / "results"
    assert main(["write-demo", "--output-dir", str(demo)]) == 0
    rc = main(
        [
            "run",
            "--input-dir",
            str(demo),
            "--output-dir",
            str(out),
            "--sample-id",
            "sample1",
            "--skip-mixscape",
            "--skip-distance",
            "--skip-de",
            "--skip-cell-annotation",
            "--n-perms",
            "10",
        ]
    )
    assert rc == 0
    report = json.loads((out / "report.json").read_text())
    for key in ANALYSIS_STATUS_KEYS:
        assert key in report, f"missing analysis status key: {key}"
        block = report[key]
        assert isinstance(block, dict), f"{key} must be an object"
        assert "skipped" in block and isinstance(block["skipped"], bool)
        if block["skipped"]:
            assert isinstance(block.get("reason"), str) and block["reason"]
            assert block.get("detail") is None or isinstance(block["detail"], str)

    assert report["mixscape"]["skipped"] is True
    assert report["mixscape"]["reason"] == "user_skip"
    assert report["edistance"]["skipped"] is True
    assert report["edistance"]["reason"] == "user_skip"
    assert report["de"]["skipped"] is True
    assert report["de"]["reason"] == "user_skip"
    assert report["cell_annotation"]["skipped"] is True
    assert report["cell_annotation"]["reason"] == "user_skip"
    # Default pca_silhouette now has a local backend (no pertpy required).
    if report["perturbation_clusters"]["skipped"]:
        assert report["perturbation_clusters"]["reason"] in {"missing_pertpy", "failed", "user_skip"}
    else:
        assert report["perturbation_clusters"].get("method") == "pca_silhouette"

    html = (out / "report.html").read_text()
    assert "mixscape" in html
    assert "user_skip" in html or "skipped" in html.lower()
