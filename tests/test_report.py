"""Tests for per-perturbation report aggregation from tiny fake tables."""

from __future__ import annotations

import json

import pandas as pd

from perturbseq.report import (
    build_perturbation_summaries,
    load_result_tables,
    summarize_perturbation,
    write_analysis_report,
)


def _write_fake_results(tmp_path):
    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    (figures / "volcano_GENEA.png").write_bytes(b"fake")

    pd.DataFrame(
        {
            "gene_target": ["GENEA", "GENEB", "GENEC"],
            "edistance": [3.0, 1.5, 0.2],
            "n_cells": [40, 25, 12],
        }
    ).to_csv(tables / "edistance.csv", index=False)

    pd.DataFrame(
        {
            "distance": [3.0, 1.5],
            "pvalue": [0.01, 0.2],
            "significant": [True, False],
            "pvalue_adj": [0.02, 0.4],
            "significant_adj": [True, False],
        },
        index=["GENEA", "GENEB"],
    ).to_csv(tables / "etest.csv")

    pd.DataFrame(
        {
            "method": ["wilcoxon_cell_level_exploratory"] * 3,
            "contrast": ["GENEA_vs_NT"] * 3,
            "names": ["TP53", "MDM2", "GAPDH"],
            "scores": [5.0, 4.0, 1.0],
            "logfoldchanges": [1.2, -0.8, 0.1],
            "pvals": [1e-6, 1e-4, 0.2],
            "pvals_adj": [1e-5, 1e-3, 0.5],
        }
    ).to_csv(tables / "de_GENEA.csv", index=False)

    pd.DataFrame(
        {
            "gene_target": ["GENEA", "GENEB"],
            "n_guides": [3, 2],
            "guides_consistent": [True, False],
            "majority_direction": ["down", "up"],
            "median_target_log2fc": [-1.0, 0.5],
            "directions": ["g1:down|g2:down|g3:down", "g1:up|g2:none"],
        }
    ).to_csv(tables / "gene_guide_consistency.csv", index=False)

    pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3", "c4"],
            "gene_target": ["GENEA", "GENEA", "NT", "NT"],
            "cell_state": ["stress", "proliferative", "epithelial", "epithelial"],
            "phase": ["G1", "S", "G1", "G1"],
        }
    ).to_csv(tables / "cell_annotations.csv", index=False)

    report = {
        "sample_id": "fake",
        "config": {"control": "NT", "replicate_col": None, "min_cells_per_pert": 10},
        "de_note": "Single-sample run: exploratory Wilcoxon.",
        "mixscape_skipped": "skipped for test",
        "gene_target_counts": {"GENEA": 40, "GENEB": 25, "GENEC": 12, "NT": 100},
    }
    (tmp_path / "report.json").write_text(json.dumps(report))
    return tmp_path


def test_summarize_perturbation_five_questions(tmp_path):
    out = _write_fake_results(tmp_path)
    tables = load_result_tables(out)
    s = summarize_perturbation("GENEA", tables, control="NT", min_cells=10)

    assert s["counts"]["n_cells"] == 40
    assert s["counts"]["n_guides"] == 3
    assert s["counts"]["n_replicates"] == 1
    assert s["counts"]["qc_status"] == "guides_consistent"

    assert s["effect_vs_control"]["effect_size"] == 3.0
    assert s["effect_vs_control"]["pvalue_adj"] == 0.02
    assert s["effect_vs_control"]["significant_adj"] is True
    assert s["effect_vs_control"]["confidence_interval"] == "N/A"

    assert isinstance(s["consistency"]["across_guides"], dict)
    assert s["consistency"]["across_guides"]["guides_consistent"] is True
    assert "single-sample" in s["consistency"]["across_samples"].lower()

    genes = s["top_affected"]["genes"]
    assert isinstance(genes, list) and genes[0]["gene"] == "TP53"
    assert "N/A" in s["top_affected"]["pathways"]

    assert any("Wilcoxon" in w or "exploratory" in w for w in s["warnings"])
    assert any("Pathway" in lim for lim in s["limitations"])
    assert s["figures"]["volcano"] == "figures/volcano_GENEA.png"


def test_build_summaries_sorted_and_na_fields(tmp_path):
    out = _write_fake_results(tmp_path)
    tables = load_result_tables(out)
    summaries = build_perturbation_summaries(tables, control="NT")
    names = [s["perturbation"] for s in summaries]
    assert names[0] == "GENEA"
    assert "GENEC" in names

    genec = next(s for s in summaries if s["perturbation"] == "GENEC")
    assert genec["effect_vs_control"]["pvalue_adj"] == "N/A"
    assert genec["top_affected"]["genes"] == "N/A"
    assert genec["counts"]["n_guides"] == "N/A"


def test_write_analysis_report_html_json(tmp_path):
    out = _write_fake_results(tmp_path)
    report = write_analysis_report(out, control="NT")
    assert (out / "report.html").exists()
    html = (out / "report.html").read_text()
    assert "GENEA" in html
    assert "1. Counts" in html or "Counts &amp; QC" in html
    assert "perturbations" in report
    assert "GENEA" in report["perturbations"]
    assert report["perturbation_summary"]["n_perturbations"] == 3
    # CLI-style rebuild should be idempotent
    report2 = write_analysis_report(out)
    assert report2["perturbation_summary"]["n_with_de"] == 1
