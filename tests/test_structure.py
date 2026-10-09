"""Import surface and pipeline stage wiring smoke checks."""

from pathlib import Path

from perturbseq import PIPELINE_STAGES, PipelineConfig
from perturbseq.analysis import run_deseq2_or_wilcoxon, run_mixscape
from perturbseq.annotate import DEFAULT_CONTROL_PATTERNS, annotate_guides, parse_gene_target
from perturbseq.guides import annotate_guides as annotate_guides_new
from perturbseq.perturbation import run_mixscape as run_mixscape_pert
from perturbseq.pipeline import (
    pipeline_config_from_mapping,
    stage_input_validation,
    stage_perturbation_modeling,
    stage_preprocessing_qc,
    stage_report,
    stage_robustness,
    stage_statistical_inference,
)
from perturbseq.preprocess import preprocess_rna as preprocess_rna_alias
from perturbseq.preprocessing import preprocess_rna
from perturbseq.report import init_report, write_html_summary, write_report_json
from perturbseq.robustness import integrate_evidence
from perturbseq.statistics import run_deseq2_or_wilcoxon as run_de


def test_target_modules_importable():
    import perturbseq.guide_qc as guide_qc
    import perturbseq.guides as guides
    import perturbseq.io as io
    import perturbseq.parallel as parallel
    import perturbseq.perturbation as perturbation
    import perturbseq.preprocessing as preprocessing
    import perturbseq.qc as qc
    import perturbseq.report as report
    import perturbseq.robustness as robustness
    import perturbseq.statistics as statistics

    assert guides.annotate_guides is annotate_guides_new
    assert hasattr(guide_qc, "run_guide_qc")
    assert hasattr(io, "validate_sample_inputs")
    assert hasattr(qc, "sample_qc_summary")
    assert preprocess_rna is preprocess_rna_alias
    assert run_mixscape_pert is run_mixscape
    assert run_de is run_deseq2_or_wilcoxon
    assert hasattr(perturbation, "run_edistance")
    assert hasattr(preprocessing, "preprocess_rna")
    assert hasattr(report, "finalize_outputs")
    assert hasattr(robustness, "integrate_evidence")
    assert hasattr(statistics, "check_experimental_design")
    assert hasattr(statistics, "run_de_contrasts")
    assert hasattr(parallel, "parallel_map")
    assert hasattr(parallel, "resolve_n_jobs")
    assert parse_gene_target("NTC_01") == "NT"
    assert DEFAULT_CONTROL_PATTERNS


def test_pipeline_stages_documented():
    assert PIPELINE_STAGES == (
        "1_input_validation",
        "2_preprocessing_qc",
        "3a_perturbation_modeling",
        "3b_statistical_inference",
        "4_robustness",
        "5_report",
    )
    for fn in (
        stage_input_validation,
        stage_preprocessing_qc,
        stage_perturbation_modeling,
        stage_statistical_inference,
        stage_robustness,
        stage_report,
    ):
        assert callable(fn)


def test_report_and_robustness_helpers(tmp_path: Path):
    report = init_report("sample1")
    integrate_evidence(
        report,
        guide_summary={
            "n_genes_consistent_guides": 2,
            "n_genes_inconsistent_guides": 0,
            "n_genes_consistency_unknown": 0,
            "n_warnings": 0,
        },
        mixscape_ok=False,
        n_de_contrasts=1,
        design={"replicate_aware": False},
    )
    assert report["confidence_flags"]["integrated_confidence"] in {"high", "moderate", "low"}
    json_path = write_report_json(report, tmp_path)
    html_path = write_html_summary(report, tmp_path)
    assert json_path.exists()
    assert html_path.exists()
    assert "Perturb-seq" in html_path.read_text()


def test_compat_aliases_share_callables():
    assert annotate_guides is annotate_guides_new
    assert run_mixscape is run_mixscape_pert
    assert run_deseq2_or_wilcoxon is run_de


def test_pipeline_config_from_mapping(tmp_path: Path):
    cfg = pipeline_config_from_mapping(
        {
            "input_dir": tmp_path / "in",
            "output_dir": tmp_path / "out",
            "sample_id": "s1",
            "control_patterns": ["^nt$", "^ntc$"],
        }
    )
    assert isinstance(cfg, PipelineConfig)
    assert cfg.sample_id == "s1"
    assert cfg.control_patterns == ("^nt$", "^ntc$")
