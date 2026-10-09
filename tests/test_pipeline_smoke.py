from perturbseq.cli import main


def test_pipeline_smoke_without_pertpy(tmp_path):
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
            "--n-perms",
            "10",
        ]
    )
    assert rc == 0
    assert (out / "sample1.tertiary.h5ad").exists()
    assert (out / "report.json").exists()
    assert (out / "run_manifest.json").exists()
    assert (out / "run.log").exists()
    assert (out / "stages" / "1_input_validation" / "status.json").exists()
    assert (out / "figures" / "qc_histograms.png").exists()
    assert (out / "figures" / "guide_composition.png").exists()
    assert (out / "tables" / "guide_qc.csv").exists()
    assert (out / "tables" / "gene_guide_consistency.csv").exists()
    assert (out / "tables" / "qc_warnings.csv").exists()
    assert (out / "tables" / "composition_audit.csv").exists()
    import json

    report = json.loads((out / "report.json").read_text())
    assert "guide_qc" in report
    assert "n_warnings" in report["guide_qc"]
    assert report["mixscape"]["skipped"] is True
    assert report["mixscape"]["reason"] == "user_skip"
    assert report["edistance"]["skipped"] is True
    assert report["edistance"]["reason"] == "user_skip"
    assert isinstance(report.get("de"), dict) and "skipped" in report["de"]
    assert "matrix_provenance" in report
    assert "X_umap" in report["matrix_provenance"]
    assert report["matrix_provenance"]["pca_source"] == "log1p_hvg"
    assert isinstance(report["matrix_provenance"]["n_hvg"], int)
    assert isinstance(report["matrix_provenance"]["n_pcs"], int)
    assert report["matrix_provenance"]["n_hvg"] > 0
    assert report["matrix_provenance"]["n_pcs"] > 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "success"
    assert "versions" in manifest
    assert "environment" in manifest
    assert "random_seed" in manifest
