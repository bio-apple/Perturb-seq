import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from perturbseq.guides import (
    CONTROL_KEY,
    INTRON_CLASS,
    NT_CLASS,
    PERT_CLASS,
    annotate_assignments,
    assignment_summary,
    guide_frequency,
    guide_target_mapping,
    segment_feature_call,
    unmatched_records,
)
from perturbseq.io import load_analysis_matrix


def test_segment_keeps_pipe_inside_guide_names():
    names = ["NEK4|design_3", "KSR2|design_1", "STRADB|design_1", "Non-Targeting_Control"]
    ordered = sorted(names, key=len, reverse=True)
    assert segment_feature_call("NEK4|design_3", ordered) == ["NEK4|design_3"]
    assert segment_feature_call("STRADB|design_1|KSR2|design_1", ordered) == [
        "STRADB|design_1",
        "KSR2|design_1",
    ]
    assert segment_feature_call("Non-Targeting_Control|NEK4|design_3", ordered) == [
        "Non-Targeting_Control",
        "NEK4|design_3",
    ]


def test_shared_nontargeting_name_is_not_a_perturbation():
    reference = pd.DataFrame(
        {
            "id": ["EditR-NTC-01", "VSGHSM_1|Positive_Control", "VSGHSM_2"],
            "name": ["Non-Targeting_Control", "ACTB|design_1", "NEK4|design_3"],
            "target_gene_id": ["", "ENSG0001", "ENSG0002"],
            "target_gene_name": ["", "ACTB", "NEK4"],
        }
    )
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3", "c4"],
            "num_features": [1, 1, 2, 0],
            "feature_call": [
                "Non-Targeting_Control",
                "ACTB|design_1",
                "Non-Targeting_Control|NEK4|design_3",
                "",
            ],
            "num_transcripts": ["10", "20", "3|4", ""],
        }
    )
    annotated = annotate_assignments(assignments, reference)
    assert annotated.loc[0, "icm_class"] == NT_CLASS
    assert annotated.loc[0, "perturbation"] == CONTROL_KEY
    assert annotated.loc[0, "target_gene_name"] == ""
    assert annotated.loc[1, "icm_class"] == PERT_CLASS
    assert annotated.loc[1, "perturbation"] == "ACTB"
    assert annotated.loc[1, "target_gene_id"] == "ENSG0001"
    assert annotated.loc[2, "icm_class"] == ""
    assert list(annotated.loc[2, "guide_parts"]) == ["Non-Targeting_Control", "NEK4|design_3"]
    _hist, frequency = guide_frequency(annotated)
    actb = frequency.set_index("guide").loc["ACTB|design_1"]
    assert actb["n_cells"] == 1
    assert actb["transcript_sum"] == 20


def test_missing_target_columns_skip_gene_labels():
    reference = pd.DataFrame(
        {
            "id": ["g1"],
            "name": ["NEK4|design_3"],
            "gene_id": ["ENSG"],
            "gene_name": ["NEK4"],
        }
    )
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["c1"],
            "num_features": [1],
            "feature_call": ["NEK4|design_3"],
            "num_transcripts": ["5"],
        }
    )
    annotated = annotate_assignments(assignments, reference)
    assert annotated.loc[0, "icm_class"] == PERT_CLASS
    assert annotated.loc[0, "perturbation"] == ""
    assert annotated.loc[0, "target_gene_name"] == ""


def test_intron_control_is_not_a_perturbation():
    reference = pd.DataFrame(
        {
            "id": ["i1"],
            "name": ["IntronCtrl_ALAS1_1"],
            "target_gene_id": ["ENSG"],
            "target_gene_name": ["ALAS1"],
        }
    )
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["c1"],
            "num_features": [1],
            "feature_call": ["IntronCtrl_ALAS1_1"],
            "num_transcripts": ["8"],
        }
    )
    annotated = annotate_assignments(assignments, reference)
    assert annotated.loc[0, "icm_class"] == INTRON_CLASS
    assert annotated.loc[0, "perturbation"] == ""
    assert annotated.loc[0, "target_gene_name"] == ""


def test_guide_statistics_keep_unmatched_barcodes_and_guides():
    reference = pd.DataFrame(
        {
            "id": ["EditR-NTC-01", "g-actb", "g-missing"],
            "name": ["Non-Targeting_Control", "ACTB|design_1", "MISSING|design_1"],
            "target_gene_id": ["", "ENSG0001", "ENSG0009"],
            "target_gene_name": ["", "ACTB", "MISSING"],
        }
    )
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["keep", "multi", "unknown", "ghost"],
            "num_features": [1, 2, 1, 0],
            "feature_call": [
                "Non-Targeting_Control",
                "ACTB|design_1|Non-Targeting_Control",
                "NOT_A_GUIDE",
                "",
            ],
            "num_transcripts": ["10", "3|4", "1", ""],
        }
    )
    annotated = annotate_assignments(assignments, reference)
    mapping = guide_target_mapping(reference).set_index("guide_id")
    assert mapping.loc["EditR-NTC-01", "target_gene_name"] == ""
    assert mapping.loc["g-actb", "target_gene_name"] == "ACTB"
    assert mapping.loc["g-actb", "guide"] == "ACTB|design_1"

    unmatched = unmatched_records(annotated, reference, pd.Index(["keep", "matrix-only"]))
    kinds = unmatched["record_type"].value_counts().to_dict()
    assert kinds["barcode_not_in_assignment"] == 1
    assert kinds["barcode_not_in_matrix"] == 3
    assert kinds["guide_not_in_reference"] == 1
    assert kinds["guide_not_in_assignment"] == 1
    assert set(unmatched.loc[unmatched["record_type"] == "guide_not_in_reference", "guide"]) == {"NOT_A_GUIDE"}
    assert set(unmatched.loc[unmatched["record_type"] == "guide_not_in_assignment", "guide"]) == {"MISSING|design_1"}

    summary = assignment_summary(annotated, reference, unmatched, pd.Index(["keep", "matrix-only"])).set_index("item")
    assert int(summary.loc["no_guide", "n"]) == 1
    assert int(summary.loc["single_guide", "n"]) == 1
    assert int(summary.loc["multiple_guides", "n"]) == 1
    assert int(summary.loc["ambiguous", "n"]) == 1
    assert int(summary.loc["passing_cells", "n"]) == 4
    assert int(summary.loc["analysis_cells", "n"]) == 1


def test_matrix_only_barcodes_are_dropped(tmp_path: Path):
    sample = "tiny"
    genes = ["G1", "G2"]
    barcodes = ["keep", "drop"]
    matrix = sparse.csr_matrix(np.array([[1, 5], [2, 6]], dtype=np.float32))
    _write_tiny_sample(tmp_path, sample, genes, barcodes, matrix)
    assignments = pd.DataFrame(
        {
            "cell_barcode": ["keep"],
            "num_features": [0],
            "feature_call": [""],
            "num_transcripts": [""],
        }
    )
    adata, kept = load_analysis_matrix(tmp_path, sample, assignments)
    assert list(adata.obs_names) == ["keep"]
    assert list(kept["cell_barcode"]) == ["keep"]
    assert adata.n_vars == 1
    assert float(adata.X[0, 0]) == 1


def _write_tiny_sample(directory: Path, sample: str, genes: list[str], barcodes: list[str], matrix) -> None:
    prefix = directory / f"{sample}.scRNA"
    feature_lines = ["GENE\tG1\tGene Expression", "GUIDE\tNEK4|design_3\tCRISPR Direct Capture"]
    with gzip.open(f"{prefix}.filtered.features.tsv.gz", "wt") as handle:
        handle.write("\n".join(feature_lines) + "\n")
    with gzip.open(f"{prefix}.filtered.barcodes.tsv.gz", "wt") as handle:
        handle.write("\n".join(barcodes) + "\n")
    with gzip.open(f"{prefix}.filtered.matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, sparse.coo_matrix(matrix))
    pd.DataFrame({"id": ["g"], "name": ["NEK4|design_3"]}).to_csv(
        f"{prefix}.feature_barcode_reference.csv", index=False
    )
    pd.DataFrame(
        {"cell_barcode": ["keep"], "num_features": [0], "feature_call": [""], "num_transcripts": [""]}
    ).to_csv(f"{prefix}.positive_cell_guide_assignments.csv", index=False)
