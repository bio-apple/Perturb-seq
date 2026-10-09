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
    assert (out / "figures" / "qc_histograms.png").exists()
    assert (out / "figures" / "guide_composition.png").exists()
