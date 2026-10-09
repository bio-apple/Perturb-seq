from perturbseq.annotate import annotate_guides, filter_singlets
from perturbseq.demo import write_demo_dragen
from perturbseq.io import load_dragen_sample
from perturbseq.qc import add_qc_metrics, filter_cells


def test_load_demo_and_annotate(tmp_path):
    write_demo_dragen(tmp_path, n_cells=60, n_genes=40, seed=1)
    rna, crispr, assignments, feature_ref = load_dragen_sample(tmp_path, "sample1")
    assert rna.n_obs == 60
    assert rna.n_vars == 40
    assert crispr.n_vars == 6
    assert "CRISPR" in "".join(crispr.var["feature_types"].astype(str))
    rna = annotate_guides(rna, assignments, feature_ref)
    assert {"guide_id", "gene_target", "perturbation", "num_features"} <= set(rna.obs.columns)
    assert "NT" in set(rna.obs["gene_target"])
    assert (rna.obs["gene_target"] == "IFNGR2").any()
    n_multi = int((rna.obs["num_features"] != 1).sum())
    singlets = filter_singlets(rna, singlet_only=True)
    assert singlets.n_obs == rna.n_obs - n_multi
    rna = add_qc_metrics(rna)
    filtered, log = filter_cells(rna, n_mads=8, min_cells=1)
    assert filtered.n_obs > 10
    assert log["n_cells_after_qc"] == filtered.n_obs
