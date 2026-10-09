from perturbseq.annotate import parse_gene_target, split_feature_calls


def test_parse_gene_target_controls_and_guides():
    assert parse_gene_target("NegCtrl_0001") == "NT"
    assert parse_gene_target("NTC_01") == "NT"
    assert parse_gene_target("non-targeting") == "NT"
    assert parse_gene_target("Non-Targeting_Control") == "NT"
    assert parse_gene_target("PDCD10_4") == "PDCD10"
    assert parse_gene_target("ATM/design_3") == "ATM"
    assert parse_gene_target("IFNGR2_1") == "IFNGR2"
    assert parse_gene_target("NEK4|design_3") == "NEK4"


def test_split_feature_calls_illumina_names():
    known = {
        "NEK4|design_3",
        "STRADB|design_1",
        "KSR2|design_1",
        "Non-Targeting_Control",
        "PIM1|design_3",
    }
    assert split_feature_calls("NEK4|design_3", known) == ["NEK4|design_3"]
    assert split_feature_calls("STRADB|design_1|KSR2|design_1", known) == ["STRADB|design_1", "KSR2|design_1"]
    assert split_feature_calls(
        "Non-Targeting_Control|Non-Targeting_Control|PIM1|design_3", known
    ) == ["Non-Targeting_Control", "Non-Targeting_Control", "PIM1|design_3"]
