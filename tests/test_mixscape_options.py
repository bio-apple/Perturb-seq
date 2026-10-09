"""Mixscape skip payload shape, cost estimate, and subset filtering helpers."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from perturbseq.cli import main
from perturbseq.perturbation import (
    estimate_mixscape_cost,
    filter_cells_for_mixscape_targets,
    select_mixscape_targets_from_edistance,
)
from perturbseq.pipeline import PipelineConfig, _effective_mixscape_mode, pipeline_config_from_mapping
from perturbseq.report import mixscape_kd_caveat, skipped_status


def test_estimate_mixscape_cost_shape():
    est = estimate_mixscape_cost(n_cells=10_000, n_targets=785, mixscape_max_targets=40)
    assert est["n_cells"] == 10_000
    assert est["n_targets"] == 785
    assert est["mixscape_max_targets"] == 40
    assert est["approx_work_units"] == 10_000 * 785
    assert est["approx_relative_to_threshold"] == round(785 / 40, 2)
    assert "O(n_cells" in est["note"]
    assert est["approx_memory_hint_gb"] >= 0.5


def test_skipped_status_accepts_estimate_extra():
    est = estimate_mixscape_cost(100, 50, 40)
    block = skipped_status(
        "too_many_targets",
        detail="n_targets>max",
        n_targets=50,
        mixscape_max_targets=40,
        estimate=est,
    )
    assert block["skipped"] is True
    assert block["reason"] == "too_many_targets"
    assert block["n_targets"] == 50
    assert block["mixscape_max_targets"] == 40
    assert block["estimate"]["approx_work_units"] == 5000


def test_filter_cells_for_mixscape_targets():
    obs = pd.DataFrame(
        {
            "gene_target": ["NT", "NT", "A", "B", "C", "unassigned"],
            "perturbation": ["NT", "NT", "A KO", "B KO", "C KO", "unassigned"],
        },
        index=[f"c{i}" for i in range(6)],
    )
    adata = ad.AnnData(X=np.ones((6, 3)), obs=obs)
    sub = filter_cells_for_mixscape_targets(adata, "NT", ["A", "C"])
    assert set(sub.obs["gene_target"].astype(str)) == {"NT", "A", "C"}
    assert sub.n_obs == 4


def test_select_mixscape_targets_from_edistance():
    edist = pd.DataFrame({"edistance": [3.0, 2.0, 1.0]}, index=["GENE_A", "GENE_B", "GENE_C"])
    assert select_mixscape_targets_from_edistance(edist, top_n=2, control="NT") == [
        "GENE_A",
        "GENE_B",
    ]
    assert select_mixscape_targets_from_edistance(edist, top_n=0, control="NT") == []


def test_effective_mixscape_mode_flags(tmp_path):
    base = dict(input_dir=tmp_path / "in", output_dir=tmp_path / "out")
    assert _effective_mixscape_mode(PipelineConfig(**base, skip_mixscape=True)) == "skip"
    assert _effective_mixscape_mode(PipelineConfig(**base, force_mixscape=True)) == "force"
    assert (
        _effective_mixscape_mode(PipelineConfig(**base, mixscape_targets=("TP53",))) == "subset"
    )
    assert _effective_mixscape_mode(PipelineConfig(**base, mixscape_top_n=5)) == "subset"
    assert _effective_mixscape_mode(PipelineConfig(**base, mixscape_mode="subset")) == "subset"
    assert _effective_mixscape_mode(PipelineConfig(**base)) == "auto"


def test_mixscape_kd_caveat():
    assert mixscape_kd_caveat("KO") is None
    note = mixscape_kd_caveat("KD")
    assert note is not None
    assert "knockdown" in note.lower() or "CRISPRi" in note
    assert mixscape_kd_caveat("CRISPRi") is not None


def test_pipeline_config_mixscape_subset_fields(tmp_path):
    cfg = pipeline_config_from_mapping(
        {
            "input_dir": tmp_path / "in",
            "output_dir": tmp_path / "out",
            "mixscape_mode": "subset",
            "mixscape_targets": "TP53,MDM2",
            "mixscape_top_n": 10,
        }
    )
    assert cfg.mixscape_mode == "subset"
    assert cfg.mixscape_targets == ("TP53", "MDM2")
    assert cfg.mixscape_top_n == 10


def test_auto_skip_report_includes_estimate(tmp_path):
    """Smoke: user skip still has explicit status; too_many_targets payload unit-tested above."""
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
            "--perturbation-type",
            "KD",
        ]
    )
    assert rc == 0
    report = json.loads((out / "report.json").read_text())
    mix = report["mixscape"]
    assert mix["skipped"] is True
    assert mix["reason"] == "user_skip"
    assert "estimate" in mix
    assert mix["estimate"]["approx_work_units"] == mix["estimate"]["n_cells"] * mix["estimate"]["n_targets"]
    assert report.get("statistical_caveats")
    assert any("knockdown" in str(c).lower() or "CRISPRi" in str(c) for c in report["statistical_caveats"])
    html = (out / "report.html").read_text()
    assert "Statistical caveats" in html or "knockdown" in html.lower()
