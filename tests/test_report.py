"""Tests for per-perturbation report aggregation from tiny fake tables."""

from __future__ import annotations

import json

import pandas as pd

from perturbseq.report import (
    CORE_STATISTICAL_CAVEATS,
    build_analysis_verdict,
    build_perturbation_summaries,
    build_summary_table,
    load_result_tables,
    summarize_perturbation,
    write_analysis_report,
    write_summary_tables,
)


def _write_fake_results(tmp_path):
    tables = tmp_path / "tables"
    figures = tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    (figures / "volcano_GENEA.png").write_bytes(b"fake")
    (figures / "umap.png").write_bytes(b"fake")
    (figures / "edistance.png").write_bytes(b"fake")

    pd.DataFrame(
        {
            "gene_target": ["GENEA", "GENEB", "GENEC"],
            "edistance": [3.0, 1.5, 0.2],
            "n_cells": [40, 25, 12],
            "edistance_ci_low": [2.5, 1.1, 0.05],
            "edistance_ci_high": [3.6, 1.9, 0.4],
            "n_bootstrap": [100, 100, 100],
            "ci_level": [0.95, 0.95, 0.95],
        }
    ).to_csv(tables / "edistance.csv", index=False)

    pd.DataFrame(
        {
            "distance": [3.0, 1.5],
            "pvalue": [0.01, 0.2],
            "significant": [True, False],
            "pvalue_adj": [0.02, 0.4],
            "significant_adj": [True, False],
            "edistance_ci_low": [2.5, 1.1],
            "edistance_ci_high": [3.6, 1.9],
            "n_bootstrap": [100, 100],
            "ci_level": [0.95, 0.95],
        },
        index=["GENEA", "GENEB"],
    ).to_csv(tables / "etest.csv")

    pd.DataFrame(
        {
            "method": ["wilcoxon_cell_level_exploratory"] * 3,
            "evidence_level": ["exploratory"] * 3,
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
            "guide_merge": ["umi", "umi"],
            "weighted_target_log2fc": [-1.1, 0.2],
            "weighted_direction": ["down", "none"],
            "weighted_direction_vote": ["down", "up"],
            "weight_sum": [100.0, 50.0],
            "n_guides_weighted": [3, 2],
            "guide_weights": ["g1:40|g2:30|g3:30", "g1:40|g2:10"],
        }
    ).to_csv(tables / "gene_guide_consistency.csv", index=False)

    pd.DataFrame(
        {
            "guide_id": ["GENEA|design_1", "GENEA|design_2", "GENEB|design_1"],
            "gene_target": ["GENEA", "GENEA", "GENEB"],
            "is_control": [False, False, False],
            "interpretation": [
                "target_knockdown_detected",
                "target_knockdown_detected",
                "inconsistent_guides",
            ],
            "assignment_confidence": ["high", "high", "medium"],
            "target_effect_direction": ["down", "down", "up"],
            "n_cells": [20, 18, 15],
            "target_expr_log2fc_vs_control": [-1.2, -0.9, 0.8],
        }
    ).to_csv(tables / "guide_qc.csv", index=False)

    pd.DataFrame(
        {
            "level": ["gene", "guide"],
            "entity": ["GENEB", "GENEB|design_1"],
            "gene_target": ["GENEB", "GENEB"],
            "warning": ["inconsistent_guides", "low_guide_umi"],
            "detail": ["g1:up|g2:none", "median_guide_umi=2.00 < 5.0"],
        }
    ).to_csv(tables / "qc_warnings.csv", index=False)

    pd.DataFrame(
        {
            "gene_target": ["GENEA", "GENEB"],
            "guide_merge": ["umi", "umi"],
            "weighted_target_log2fc": [-1.1, 0.2],
            "weighted_direction": ["down", "none"],
            "weighted_direction_vote": ["down", "up"],
            "weight_sum": [100.0, 50.0],
            "n_guides_weighted": [3, 2],
            "guide_weights": ["g1:40|g2:30|g3:30", "g1:40|g2:10"],
        }
    ).to_csv(tables / "gene_guide_weighted.csv", index=False)

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
        "config": {
            "control": "NT",
            "replicate_col": None,
            "min_cells_per_pert": 10,
            "etest_power_min_cells": 50,
        },
        "de": {
            "skipped": False,
            "reason": None,
            "detail": "Single-sample run: exploratory Wilcoxon.",
            "evidence_level": "exploratory",
        },
        "experimental_design": {"evidence_level": "exploratory", "replicate_aware": False},
        "mixscape": {"skipped": True, "reason": "user_skip", "detail": "skipped for test"},
        "edistance": {"skipped": False, "reason": None},
        "de_note": "Single-sample run: exploratory Wilcoxon.",
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
    # n_cells=40 < etest_power_min_cells=50 → low_power suppresses reported significance
    assert s["effect_vs_control"]["low_power"] is True
    assert s["effect_vs_control"]["significant_adj_reported"] is False
    ci = s["effect_vs_control"]["confidence_interval"]
    assert isinstance(ci, dict)
    assert ci["low"] == 2.5 and ci["high"] == 3.6
    assert ci["source"] == "edistance_bootstrap_cells"
    assert s["effect_vs_control"]["evidence_level"] == "exploratory"

    assert isinstance(s["consistency"]["across_guides"], dict)
    assert s["consistency"]["across_guides"]["guides_consistent"] is True
    assert "single-sample" in s["consistency"]["across_samples"].lower()

    genes = s["top_affected"]["genes"]
    assert isinstance(genes, list) and genes[0]["gene"] == "TP53"
    assert "N/A" in s["top_affected"]["pathways"]

    assert any("Wilcoxon" in w or "exploratory" in w for w in s["warnings"])
    assert any("Pathway" in lim for lim in s["limitations"])
    assert s["figures"]["volcano"] == "figures/volcano_GENEA.png"
    assert s["figures"]["umap"] == "figures/umap.png"
    assert s["figures"]["edistance"] == "figures/edistance.png"


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


def test_summary_table_csv_and_md(tmp_path):
    out = _write_fake_results(tmp_path)
    tables = load_result_tables(out)
    summaries = build_perturbation_summaries(tables, control="NT")
    df = build_summary_table(summaries)
    assert "GENEA" in set(df["perturbation"])
    assert "edistance" in df.columns
    assert "top_gene" in df.columns
    csv_path, md_path = write_summary_tables(out, summaries)
    assert csv_path.exists() and md_path.exists()
    loaded = pd.read_csv(csv_path)
    assert len(loaded) == len(summaries)
    md = md_path.read_text()
    assert "GENEA" in md
    assert "| perturbation |" in md or "perturbation" in md


def test_html_contains_caveats_checklist_and_download_links(tmp_path):
    out = _write_fake_results(tmp_path)
    write_analysis_report(out, control="NT")
    html = (out / "report.html").read_text()
    assert 'id="statistical-caveats"' in html
    assert "<details" in html and "checklist" in html
    assert any(c[:40] in html for c in CORE_STATISTICAL_CAVEATS)
    assert 'href="summary_table.csv"' in html
    assert 'href="summary_table.md"' in html
    assert (out / "summary_table.csv").exists()
    assert (out / "summary_table.md").exists()
    assert "figures/volcano_GENEA.png" in html
    assert "figures/umap.png" in html
    assert "figures/edistance.png" in html
    assert "QC warnings" in html


def test_html_maps_qc_warnings_and_interpretations(tmp_path):
    out = _write_fake_results(tmp_path)
    report = write_analysis_report(out, control="NT")
    geneb = report["perturbations"]["GENEB"]
    assert geneb["guide_qc"]["status"] == "guides_inconsistent"
    assert any(w["warning"] == "inconsistent_guides" for w in geneb["guide_qc"]["warnings"])
    assert any(i["interpretation"] == "inconsistent_guides" for i in geneb["guide_qc"]["interpretations"])
    assert isinstance(geneb["guide_qc"]["weighted_summary"], dict)
    assert geneb["guide_qc"]["weighted_summary"]["weighted_direction_vote"] == "up"

    html = (out / "report.html").read_text()
    assert "inconsistent_guides" in html
    assert "low_guide_umi" in html
    assert "Guide interpretations" in html
    assert "QC warnings" in html
    assert "Weighted guide summary" in html
    assert "target_knockdown_detected" in html  # GENEA interpretation on card


def test_build_analysis_verdict_from_fake_tables(tmp_path):
    """Verdict helper: significant (non-low-power) + inconsistent guides + Wilcoxon caveat."""
    summaries = [
        {
            "perturbation": "GENEA",
            "counts": {"qc_status": "guides_consistent"},
            "effect_vs_control": {
                "significant_adj_reported": True,
                "low_power": False,
                "de_method": "wilcoxon_cell_level_exploratory",
            },
            "consistency": {"across_guides": {"guides_consistent": True}},
        },
        {
            "perturbation": "GENEB",
            "counts": {"qc_status": "guides_inconsistent"},
            "effect_vs_control": {
                "significant_adj_reported": False,
                "low_power": True,
                "de_method": "wilcoxon_cell_level_exploratory",
            },
            "consistency": {"across_guides": {"guides_consistent": False}},
        },
        {
            "perturbation": "GENEC",
            "counts": {"qc_status": "ok"},
            "effect_vs_control": {
                "significant_adj_reported": False,
                "low_power": True,
                "de_method": "N/A",
            },
            "consistency": {"across_guides": "N/A"},
        },
    ]
    report = {"config": {"replicate_col": None}, "de": {"detail": "Single-sample run: exploratory Wilcoxon."}}
    verdict = build_analysis_verdict(summaries, report)
    assert verdict["n_perturbations"] == 3
    assert verdict["n_significant_reported"] == 1
    assert verdict["n_low_power"] == 2
    assert verdict["n_guides_inconsistent"] == 1
    assert "Wilcoxon" in verdict["summary_sentence"] or "exploratory" in verdict["summary_sentence"]
    assert "1 have E-test significance reported" in verdict["summary_sentence"]
    assert "1 gene(s) show inconsistent" in verdict["summary_sentence"]

    out = _write_fake_results(tmp_path)
    written = write_analysis_report(out, control="NT")
    assert "verdict" in written
    assert written["summary_sentence"] == written["verdict"]["summary_sentence"]
    assert written["verdict"]["n_guides_inconsistent"] >= 1
    html = (out / "report.html").read_text()
    assert 'id="analysis-verdict"' in html
    assert written["summary_sentence"] in html
