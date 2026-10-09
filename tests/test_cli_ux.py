"""CLI UX: config --generate and run --dry-run."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from perturbseq.cli import main
from perturbseq.io import resolve_sample_files
from perturbseq.pipeline import PipelineConfig, dry_run_plan, pipeline_config_from_mapping
from perturbseq.repro import generate_config_template


def test_config_generate_flag_and_action(tmp_path: Path):
    out = tmp_path / "pipeline.template.yaml"
    assert main(["config", "--generate", "-o", str(out)]) == 0
    text = out.read_text()
    assert "sample_id" in text
    assert "input_dir" in text
    assert "control_patterns" in text
    assert text.lstrip().startswith("#")
    loaded = yaml.safe_load(text)
    assert isinstance(loaded, dict)
    assert "sample_id" in loaded

    out2 = tmp_path / "via_action.yaml"
    assert main(["config", "generate", "-o", str(out2)]) == 0
    assert out2.is_file()
    assert "n_mads" in out2.read_text()


def test_generate_config_template_covers_pipeline_fields():
    from dataclasses import fields

    text = generate_config_template()
    for f in fields(PipelineConfig):
        if f.name == "extra":
            continue
        assert f.name in text, f"missing field {f.name} in template"


def test_dry_run_validates_and_does_not_write_results(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
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
            "--dry-run",
        ]
    )
    assert rc == 0
    captured = capsys.readouterr()
    assert "Dry-run OK" in captured.out
    assert "1_input_validation" in captured.out
    assert "3a_perturbation_modeling" in captured.out
    assert "skip mixscape" in captured.out
    assert "skip edistance" in captured.out
    assert "Resolved parameters:" in captured.out
    assert not (out / "report.json").exists()
    assert not (out / "run_manifest.json").exists()


def test_dry_run_plan_helper(tmp_path: Path):
    demo = tmp_path / "demo"
    assert main(["write-demo", "--output-dir", str(demo)]) == 0
    cfg = pipeline_config_from_mapping(
        {
            "input_dir": demo,
            "output_dir": tmp_path / "out",
            "skip_de": True,
        }
    )
    plan = dry_run_plan(cfg)
    assert plan["sample_id"] == "sample1"
    assert "matrix" in plan["input_files"]
    de_stage = next(s for s in plan["stages"] if s["stage"] == "3b_statistical_inference")
    assert "skip de" in de_stage["notes"]


def test_missing_files_error_suggests_commands(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="write-demo") as exc_info:
        resolve_sample_files(tmp_path / "missing", "sample1")
    msg = str(exc_info.value)
    assert "--dry-run" in msg
    assert "--sample-id" in msg
