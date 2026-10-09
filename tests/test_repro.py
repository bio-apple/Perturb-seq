"""Unit tests for checksum / config / stage resume helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from perturbseq.pipeline import pipeline_config_from_mapping
from perturbseq.repro import (
    STAGE_ORDER,
    RunTracker,
    deep_merge,
    load_yaml_config,
    sha256_file,
    sha256_paths,
)


def test_sha256_file_streams_and_is_stable(tmp_path: Path):
    path = tmp_path / "blob.bin"
    path.write_bytes(b"perturb-seq" * 1000)
    a = sha256_file(path, chunk_size=32)
    b = sha256_file(path, chunk_size=4096)
    assert a == b
    assert len(a) == 64
    path.write_bytes(b"changed")
    assert sha256_file(path) != a


def test_sha256_paths_missing_is_empty(tmp_path: Path):
    present = tmp_path / "ok.txt"
    present.write_text("hi")
    missing = tmp_path / "nope.txt"
    out = sha256_paths([present, missing])
    assert out[str(present)]
    assert out[str(missing)] == ""


def test_load_yaml_config_and_pipeline_mapping(tmp_path: Path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(
        "\n".join(
            [
                "sample_id: demo",
                f"input_dir: {tmp_path / 'in'}",
                f"output_dir: {tmp_path / 'out'}",
                "n_mads: 3.5",
                "skip_mixscape: true",
                "control_patterns:",
                '  - "^nt$"',
                '  - "^negctrl"',
                "random_state: 7",
            ]
        )
    )
    data = load_yaml_config(cfg)
    assert data["n_mads"] == 3.5
    merged = deep_merge(data, {"n_perms": 11, "n_mads": None})
    assert merged["n_mads"] == 3.5
    assert merged["n_perms"] == 11
    pc = pipeline_config_from_mapping(merged)
    assert pc.sample_id == "demo"
    assert pc.skip_mixscape is True
    assert pc.random_state == 7
    assert pc.control_patterns == ("^nt$", "^negctrl")


def test_run_tracker_status_resume_and_invalidate(tmp_path: Path):
    out = tmp_path / "results"
    tracker = RunTracker(out, resume=False, seed=0)
    stage = STAGE_ORDER[0]
    inputs = {"a.txt": "abc"}
    params = {"n_mads": 5.0}
    payload = tracker.begin_stage(stage, inputs=inputs, params=params)
    artifact = tracker.stage_dir(stage) / "out.txt"
    artifact.write_text("ok")
    tracker.finish_stage(
        stage,
        payload,
        status="success",
        outputs={"out": str(artifact)},
    )
    status = tracker.load_status(stage)
    assert status is not None
    assert status["status"] == "success"
    assert status["output_checksums"]["out"] == sha256_file(artifact)
    assert (out / "run_manifest.json").is_file()
    assert (out / "run.log").is_file()

    resumed = RunTracker(out, resume=True, seed=0)
    assert resumed.should_skip(stage, params=params, input_checksums=inputs) is True
    # Param change invalidates this and downstream.
    assert resumed.should_skip(stage, params={"n_mads": 9.0}, input_checksums=inputs) is False
    assert resumed.load_status(stage) is None
    assert resumed.load_status(STAGE_ORDER[1]) is None


def test_resume_skips_when_outputs_present(tmp_path: Path):
    out = tmp_path / "results"
    tracker = RunTracker(out, resume=False, seed=1)
    stage = "2_preprocessing_qc"
    params = {"x": 1}
    inputs = {"in": "hash1"}
    payload = tracker.begin_stage(stage, inputs=inputs, params=params)
    ckpt = tracker.checkpoint_path(stage)
    ckpt.write_bytes(b"fake-h5ad")
    tracker.finish_stage(stage, payload, status="success", outputs={"adata": str(ckpt)})

    resumed = RunTracker(out, resume=True, seed=1)
    assert resumed.should_skip(stage, params=params, input_checksums=inputs)
    # Missing output forces rerun.
    ckpt.unlink()
    assert not resumed.should_skip(stage, params=params, input_checksums=inputs)


def test_pipeline_config_requires_paths():
    with pytest.raises(ValueError, match="input_dir"):
        pipeline_config_from_mapping({"sample_id": "x"})


def test_resume_skips_completed_stages_on_demo(tmp_path: Path):
    """Tiny demo: first run completes; second --resume skips unchanged stages."""
    from perturbseq.cli import main

    demo = tmp_path / "demo"
    out = tmp_path / "results"
    assert main(["write-demo", "--output-dir", str(demo)]) == 0
    argv = [
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
    assert main(argv) == 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "success"
    assert "random_seed" in manifest
    assert manifest["inputs"]
    assert all(stage in manifest["stages"] for stage in STAGE_ORDER[:2])

    # Simulate interrupted final stage: drop report stage status + h5ad.
    report_status = out / "stages" / "5_report" / "status.json"
    h5ad = out / "sample1.tertiary.h5ad"
    assert report_status.is_file()
    report_status.unlink()
    if h5ad.is_file():
        h5ad.unlink()

    assert main([*argv, "--resume"]) == 0
    assert h5ad.is_file()
    log = (out / "run.log").read_text()
    assert "Resuming: skip stage 1_input_validation" in log
    assert "Resuming: skip stage 2_preprocessing_qc" in log
    status5 = json.loads((out / "stages" / "5_report" / "status.json").read_text())
    assert status5["status"] == "success"
