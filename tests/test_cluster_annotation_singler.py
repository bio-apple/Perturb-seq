import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq.plotting import SingleRReference, cluster_annotations, singler_reference_name


def _reference(organism: str) -> SingleRReference:
    rng = np.random.default_rng(0)
    matrix = rng.normal(0, 0.1, size=(30, 6))
    matrix[:15, :3] += 5
    matrix[15:, 3:] += 5
    name = singler_reference_name(organism)
    return SingleRReference(
        name=name,
        organism=organism,
        matrix=matrix,
        labels=("TypeA", "TypeA", "TypeA", "TypeB", "TypeB", "TypeB"),
        features=tuple(f"g{i}" for i in range(1, 31)),
    )


def _adata() -> AnnData:
    expression = np.zeros((2, 30))
    expression[0, :15] = 4
    expression[1, 15:] = 4
    genes = [f"g{i}" for i in range(1, 31)]
    var = pd.DataFrame({"gene_symbol": genes}, index=genes)
    obs = pd.DataFrame({"leiden": ["0", "1"], "cell_type": ["ALL", "ALL"]})
    adata = AnnData(expression, obs=obs, var=var)
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    return adata


def test_mouse_singler_rejects_human_reference_and_returns_cluster_labels():
    assert singler_reference_name("mouse") == "MouseRNAseqData"
    assert singler_reference_name("human") == "HumanPrimaryCellAtlasData"
    adata = _adata()
    human = _reference("human")
    with pytest.raises(ValueError, match="不一致"):
        cluster_annotations(adata, "mouse", reference=human)
    sneaky = human._replace(organism="mouse")
    with pytest.raises(ValueError, match="不一致"):
        cluster_annotations(adata, "mouse", reference=sneaky)

    labels, summary = cluster_annotations(adata, "mouse", reference=_reference("mouse"))
    assert list(labels) == ["TypeA", "TypeB"]
    assert list(summary["singler_label"]) == ["TypeA", "TypeB"]
    assert list(summary["reference"]) == ["MouseRNAseqData", "MouseRNAseqData"]
    assert list(adata.obs["singler_label"]) == ["TypeA", "TypeB"]
    assert list(adata.obs["cell_type"]) == ["ALL", "ALL"]
