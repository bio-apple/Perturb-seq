import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq.assign import control_cell_mask, qc_metrics
from perturbseq.config import ConfigError, control_spec, validate_config, verify_species_consistency
from perturbseq.distance_tests import write_edistance_results
from perturbseq.io import InputMismatch, load_sample, matrix_only_as_no_guide
from perturbseq.perturb_effect import (
    descriptive_effects,
    run_mixscale,
    run_mixscape,
    tutorial_differential_expression,
)
from perturbseq.plotting import FigureStore, draw_descriptive_heatmap
from perturbseq.preprocess import _drop_mitochondrial_hvg


def _base(organism: str = "human") -> dict:
    return {
        "experiment": {
            "sample_id": "toy",
            "organism": organism,
            "crispr_mode": "CRISPR_KO",
        },
        "labels": {
            "nt_control": "Non-Targeting",
            "multi_guide_policy": "drop",
            "guides_per_cell": 1,
        },
        "qc": {
            "thresholds_confirmed": True,
            "min_umi": 1,
            "max_umi": 10,
            "min_genes": 1,
            "max_mito_frac": 0.2,
            "min_cells_per_target": 1,
            "min_cells_per_guide": 1,
        },
        "input": {"guide_separator": "|", "barcode_policy": "strict"},
        "output": {"output_dir": "results"},
        "preprocess": {},
        "perturb": {},
        "statistics": {"mode": "descriptive"},
        "augur": {},
        "plotting": {"theme": "whitegrid", "formats": ["png"]},
        "runtime": {"seed": 0},
    }


def test_only_human_and_mouse_pass_config_validation():
    human = _base("Human")
    validate_config(human)
    assert human["experiment"]["organism"] == "human"
    assert human["qc"]["mito_gene_prefix"] == "MT-"
    mouse = _base("MOUSE")
    validate_config(mouse)
    assert mouse["experiment"]["organism"] == "mouse"
    assert mouse["qc"]["mito_gene_prefix"] == "mt-"


@pytest.mark.parametrize(
    "organism",
    [None, "", "   ", "organelle", "细胞器", "human,mouse", "rat", "unknown", "other"],
)
def test_organism_other_than_human_or_mouse_fails(organism):
    config = _base()
    if organism is None:
        del config["experiment"]["organism"]
    else:
        config["experiment"]["organism"] = organism
    with pytest.raises(ConfigError, match="只能填 human 或 mouse"):
        validate_config(config)


def test_configured_mito_prefix_must_match_species():
    config = _base("human")
    config["qc"]["mito_gene_prefix"] = "mt-"
    with pytest.raises(ConfigError, match="MT-"):
        validate_config(config)


def test_human_ensg_and_mt_pass_and_ignore_chloroplast():
    config = _base("human")
    validate_config(config)
    ids = ["ENSG00000198888", "ENSG00000198763", "ATCG00010"]
    symbols = ["MT-ND1", "ACTB", "pt-rbcL"]
    verify_species_consistency(config, ids, symbols)
    assert config["experiment"]["organism"] == "human"
    assert config["experiment"]["resolved_organism"] == "human"
    assert config["experiment"]["resolved_mito_prefix"] == "MT-"
    assert config["experiment"]["species_check"] == "passed"
    assert any("叶绿体" in item for item in config["experiment"]["other_organelle_prefixes"])
    var = pd.DataFrame({"gene_id": ids, "gene_symbol": symbols}, index=symbols)
    adata = AnnData(np.array([[10.0, 40.0, 50.0]]), var=var)
    fraction = qc_metrics(adata, config["qc"]["mito_gene_prefix"])["pct_counts_mito"].iloc[0]
    assert fraction == pytest.approx(10 / 100)
    adata.var["highly_variable"] = [True, True, False]
    adata.var["dispersions_norm"] = [5.0, 1.0, 0.0]
    _drop_mitochondrial_hvg(adata, config, 2)
    assert not bool(adata.var.loc["MT-ND1", "highly_variable"])


def test_mouse_ensmusg_and_mt_pass():
    config = _base("mouse")
    validate_config(config)
    ids = ["ENSMUSG00000064341", "ENSMUSG00000000001"]
    symbols = ["mt-Nd1", "Actb"]
    verify_species_consistency(config, ids, symbols)
    assert config["experiment"]["organism"] == "mouse"
    assert config["experiment"]["resolved_mito_prefix"] == "mt-"
    assert config["experiment"]["species_check"] == "passed"
    var = pd.DataFrame({"gene_id": ids, "gene_symbol": symbols}, index=symbols)
    adata = AnnData(np.array([[5.0, 15.0]]), var=var)
    fraction = qc_metrics(adata, "mt-")["pct_counts_mito"].iloc[0]
    assert fraction == pytest.approx(5 / 20)
    adata.var["highly_variable"] = [True, True]
    adata.var["dispersions_norm"] = [4.0, 1.0]
    _drop_mitochondrial_hvg(adata, config, 2)
    assert not bool(adata.var.loc["mt-Nd1", "highly_variable"])


def test_human_rejects_mouse_gene_id_or_symbol():
    config = _base("human")
    validate_config(config)
    with pytest.raises(ConfigError, match="ENSMUSG00000064341"):
        verify_species_consistency(config, ["ENSMUSG00000064341", "ENSG00000198888"], ["ACTB"])
    assert config["experiment"]["organism"] == "human"
    assert config["experiment"]["species_check"] == "failed"
    fresh = _base("human")
    validate_config(fresh)
    with pytest.raises(ConfigError, match="mt-Nd1"):
        verify_species_consistency(fresh, ["ENSG00000198888"], ["mt-Nd1", "MT-ND1"])
    assert fresh["experiment"]["organism"] == "human"
    assert fresh["experiment"]["species_check"] == "failed"


def test_mouse_rejects_human_gene_id_or_symbol():
    config = _base("mouse")
    validate_config(config)
    with pytest.raises(ConfigError, match="ENSG00000198888"):
        verify_species_consistency(config, ["ENSG00000198888"], ["Actb"])
    assert config["experiment"]["organism"] == "mouse"
    assert config["experiment"]["species_check"] == "failed"
    fresh = _base("mouse")
    validate_config(fresh)
    with pytest.raises(ConfigError, match="MT-ND1"):
        verify_species_consistency(fresh, ["ENSMUSG00000064341"], ["MT-ND1"])
    assert fresh["experiment"]["organism"] == "mouse"
    assert fresh["experiment"]["species_check"] == "failed"


def test_default_crispr_mode_is_knockout_and_aliases_keep_direction():
    missing = _base()
    del missing["experiment"]["crispr_mode"]
    validate_config(missing)
    assert missing["experiment"]["crispr_mode"] == "CRISPR_KO"
    assert missing["experiment"]["resolved_expected_direction"] == "not_required"
    for raw, mode, direction in (
        ("KO", "CRISPR_KO", "not_required"),
        ("i", "CRISPRi", "down"),
        ("a", "CRISPRa", "up"),
    ):
        config = _base()
        config["experiment"]["crispr_mode"] = raw
        validate_config(config)
        assert config["experiment"]["crispr_mode"] == mode
        assert config["experiment"]["resolved_expected_direction"] == direction
    other = _base()
    other["experiment"]["crispr_mode"] = "other"
    with pytest.raises(ConfigError, match="crispr_mode_description"):
        validate_config(other)


def test_missing_nt_skips_dependent_steps(tmp_path):
    config = _base()
    del config["labels"]["nt_control"]
    validate_config(config)
    assert config["labels"]["nt_control"] == ""
    adata = AnnData(np.ones((2, 2)))
    effects, guides, status = descriptive_effects(adata, config, branch="main")
    assert effects.empty and guides.empty
    assert list(effects.columns) == [] or "pvalue" in effects.columns
    assert status["reason"].iloc[0] == "no_zero_guide_cells"
    assert status["status"].iloc[0] == "skipped"
    assert "no_nt_control" not in set(status["reason"].astype(str))
    run_mixscape(adata, config)
    run_mixscale(adata, config)
    assert adata.uns["mixscape_skip_reason"] == "no_zero_guide_cells"
    assert adata.uns["mixscale_skip_reason"] == "no_zero_guide_cells"
    write_edistance_results(tmp_path, adata, config)
    skipped = tmp_path / "E-distance" / "figures" / "edistance_etest.skipped.txt"
    assert skipped.read_text().strip() == "no_zero_guide_cells"
    store = FigureStore(tmp_path, config, {"sample": "toy"})
    draw_descriptive_heatmap(store, effects)
    counted = AnnData(np.ones((4, 2)))
    counted.obs["target_gene"] = ["A", "A", "B", "B"]
    counted.obs["perturbation_status"] = "targeting"
    counted.obs["guide_id"] = ["g1", "g1", "g2", "g2"]
    counted.layers["counts"] = np.ones((4, 2))
    _frames, de_status = tutorial_differential_expression(counted, config, branch="main")
    assert set(de_status["reason"]) == {"single_sample_de_not_estimable"}


def test_guide_exploratory_and_nonstrict_barcodes_still_fail():
    exploratory = _base()
    exploratory["statistics"]["mode"] = "guide_exploratory"
    with pytest.raises(ConfigError, match="unsupported_inference_backend"):
        validate_config(exploratory)
    loose = _base()
    loose["input"]["barcode_policy"] = "intersection"
    with pytest.raises(ConfigError, match="strict"):
        validate_config(loose)


def test_designed_nt_does_not_use_zero_guide_cells():
    config = _base()
    validate_config(config)
    spec = control_spec(config)
    assert spec["temporary"] is False
    assert spec["control_source"] == "designed_nt"
    obs = pd.DataFrame(
        {
            "assignment_status": ["single_nt", "no_guide", "single_targeting", "multi_guide"],
            "n_guides": [1, 0, 1, 2],
            "perturbation_status": ["nt", "no_guide", "targeting", "multi_guide"],
            "target_gene": ["Non-Targeting", pd.NA, "Actb", pd.NA],
        }
    )
    mask = control_cell_mask(obs, spec)
    assert mask.tolist() == [True, False, False, False]


def test_zero_guide_cells_are_temporary_control_for_log2fc():
    config = _base()
    del config["labels"]["nt_control"]
    validate_config(config)
    spec = control_spec(config)
    assert spec["control_source"] == "zero_guide_cells"
    assert spec["control_label"] == "zero_guide"
    assert "Non-Targeting" not in spec["control_label"]
    obs = pd.DataFrame(
        {
            "assignment_status": ["no_guide", "no_guide", "single_targeting", "multi_guide"],
            "n_guides": [0, 0, 1, 2],
            "perturbation_status": ["no_guide", "no_guide", "targeting", "multi_guide"],
            "target_gene": [pd.NA, pd.NA, "Actb", pd.NA],
            "guide_id": ["", "", "g1", ""],
            "guide_name": ["", "", "g1", ""],
        }
    )
    mask = control_cell_mask(obs, spec)
    assert mask.tolist() == [True, True, False, False]
    counts = np.array(
        [[8.0, 2.0], [8.0, 2.0], [2.0, 8.0], [1.0, 1.0]],
        dtype=float,
    )
    adata = AnnData(counts.copy())
    adata.obs = obs
    adata.var = pd.DataFrame(
        {"gene_id": ["gA", "gB"], "gene_symbol": ["Actb", "Gapdh"]},
        index=["Actb", "Gapdh"],
    )
    adata.layers["counts"] = counts.copy()
    effects, _guides, status = descriptive_effects(adata, config, branch="main")
    assert "no_nt_control" not in set(status["reason"].astype(str))
    assert set(status.loc[status["target_gene"] == "Actb", "status"]) == {"descriptive"}
    actb = effects[effects["target_gene"] == "Actb"]
    assert not actb.empty
    assert actb["pvalue"].isna().all()
    assert actb["fdr"].isna().all()
    assert actb["ci_lower"].isna().all()
    assert actb["ci_upper"].isna().all()
    assert int(actb["n_nt_cells"].iloc[0]) == 2


def _write_gz(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _toy_sample(folder: Path, assignment_rows: str) -> dict[str, Path]:
    _write_gz(
        folder / "features.tsv.gz",
        "ENSMUSG00000000001\tActb\tGene Expression\n"
        "ENSMUSG00000000002\tGapdh\tGene Expression\n"
        "gA\tGuideA\tCRISPR Guide Capture\n",
    )
    _write_gz(folder / "barcodes.tsv.gz", "A\nB\nC\n")
    _write_gz(
        folder / "matrix.mtx.gz",
        "%%MatrixMarket matrix coordinate integer general\n"
        "3 3 3\n"
        "1 1 5\n"
        "2 2 4\n"
        "3 1 2\n",
    )
    (folder / "reference.csv").write_text(
        "id,name,feature_type,target_gene_name,target_gene_id\n"
        "gA,GuideA,CRISPR Guide Capture,Actb,ENSMUSG00000000001\n"
    )
    (folder / "assignments.csv").write_text(
        "cell_barcode,num_features,feature_call,num_transcripts\n" + assignment_rows
    )
    return {
        "matrix": folder / "matrix.mtx.gz",
        "barcodes": folder / "barcodes.tsv.gz",
        "features": folder / "features.tsv.gz",
        "feature_reference": folder / "reference.csv",
        "assignments": folder / "assignments.csv",
    }


def test_mouse_does_not_run_human_cell_cycle():
    from perturbseq.pipeline import _score_cell_cycle
    from perturbseq.plotting import score_human_cell_cycle_after_log1p

    config = _base("mouse")
    validate_config(config)
    symbols = ["MCM5", "PCNA", "MKI67", "TOP2A"]
    adata = AnnData(np.ones((4, len(symbols))))
    adata.var_names = symbols
    adata.var["gene_symbol"] = symbols
    adata.uns["log1p"] = {"base": None}
    _score_cell_cycle(adata, config)
    assert adata.uns["cell_cycle_status"] == "skipped_no_mouse_cell_cycle_genes"
    assert "小鼠" in str(adata.uns["cell_cycle_note"])
    assert "phase" not in adata.obs.columns
    assert "cell_cycle_s_found" not in adata.uns
    assert score_human_cell_cycle_after_log1p(adata, organism="mouse") == "skipped_not_human"
    assert "phase" not in adata.obs.columns


def test_matrix_only_is_no_guide_and_assignment_only_still_fails(tmp_path):
    kept = _toy_sample(tmp_path / "kept", "A,1,GuideA,2\n")
    adata, annotated, checks = load_sample(
        kept,
        nt_control="",
        feature_type_gex="Gene Expression",
        feature_type_crispr="CRISPR Guide Capture",
        multi_guide_policy="drop",
        guides_per_cell=1,
        min_guide_umi=None,
        mito_prefix="mt-",
    )
    assert checks["n_matrix_only_as_no_guide"] == 2
    assert set(adata.obs_names) == {"A", "B", "C"}
    status = annotated.set_index("cell_barcode")["assignment_status"]
    assert status.loc["A"] == "single_targeting"
    assert status.loc["B"] == "no_guide"
    assert status.loc["C"] == "no_guide"
    assert int(annotated.set_index("cell_barcode").loc["B", "n_guides"]) == 0

    missing = _toy_sample(tmp_path / "missing", "A,1,GuideA,2\nZ,1,GuideA,2\n")
    with pytest.raises(InputMismatch, match="仅分配") as caught:
        load_sample(
            missing,
            nt_control="",
            feature_type_gex="Gene Expression",
            feature_type_crispr="CRISPR Guide Capture",
            multi_guide_policy="drop",
            guides_per_cell=1,
            min_guide_umi=None,
            mito_prefix="mt-",
        )
    assert caught.value.only_assignment == ["Z"]
    assert "B" in caught.value.only_matrix
    assert matrix_only_as_no_guide(caught.value.only_matrix, caught.value.only_assignment) == 0
