"""Parallel helper and DE outer-loop smoke tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from perturbseq.parallel import parallel_map, resolve_n_jobs
from perturbseq.statistics import run_de_contrasts


def test_resolve_n_jobs():
    assert resolve_n_jobs(1) == 1
    assert resolve_n_jobs(0) == 1
    assert resolve_n_jobs(4) == 4
    assert resolve_n_jobs(-1) >= 1


def test_parallel_map_n_jobs_gt_1():
    def square(x: int) -> int:
        return x * x

    items = list(range(8))
    sequential = parallel_map(square, items, n_jobs=1)
    parallel = parallel_map(square, items, n_jobs=2)
    assert sequential == [x * x for x in items]
    assert parallel == sequential


def test_run_de_contrasts_n_jobs():
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 40
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(float)
    # Make group A slightly higher on first genes for a detectable contrast.
    X[:20, :5] += 5.0
    obs = pd.DataFrame(
        {
            "de_group": ["A"] * 20 + ["B"] * 20 + ["NT"] * 20,
        }
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    adata = AnnData(X=X, obs=obs, var=var)

    tables_1, errors_1 = run_de_contrasts(
        adata, ["A", "B"], reference="NT", groupby="de_group", min_cells=5, n_jobs=1
    )
    tables_2, errors_2 = run_de_contrasts(
        adata, ["A", "B"], reference="NT", groupby="de_group", min_cells=5, n_jobs=2
    )
    assert not errors_1 and not errors_2
    assert [g for g, _ in tables_1] == [g for g, _ in tables_2] == ["A", "B"]
    for (_, t1), (_, t2) in zip(tables_1, tables_2, strict=True):
        assert not t1.empty and not t2.empty
        assert list(t1.columns) == list(t2.columns)
