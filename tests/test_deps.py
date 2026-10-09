"""pertpy optional-dependency messaging."""

from __future__ import annotations

import sys
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from perturbseq._deps import (
    PERTPY_INSTALL,
    PERTPY_MISSING_MSG,
    PERTPY_SKIP_FLAGS,
    is_pertpy_import_error,
    pertpy_available,
    require_pertpy,
    warn_if_pertpy_missing,
)


def test_require_pertpy_message_when_missing():
    """Force ImportError for pertpy and check user-facing guidance."""
    real_import = __import__

    def _block_pertpy(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "pertpy" or (isinstance(name, str) and name.startswith("pertpy.")):
            raise ImportError("No module named 'pertpy'", name="pertpy")
        return real_import(name, globals, locals, fromlist, level)

    saved = sys.modules.pop("pertpy", None)
    try:
        with patch("builtins.__import__", side_effect=_block_pertpy):
            with pytest.raises(ImportError, match="pertpy is required") as excinfo:
                require_pertpy()
        msg = str(excinfo.value)
        assert "pip install" in msg
        assert PERTPY_SKIP_FLAGS in msg
        assert "pertpy" in msg.lower()
        assert ".[de]" in msg or "pertpy[de]" in msg or PERTPY_INSTALL in msg
    finally:
        if saved is not None:
            sys.modules["pertpy"] = saved


def test_is_pertpy_import_error_detects_name_and_message():
    assert is_pertpy_import_error(ImportError("No module named 'pertpy'", name="pertpy"))
    wrapped = ImportError(PERTPY_MISSING_MSG)
    wrapped.__cause__ = ImportError("No module named 'pertpy'", name="pertpy")
    assert is_pertpy_import_error(wrapped)
    assert not is_pertpy_import_error(ImportError("No module named 'numpy'", name="numpy"))
    assert not is_pertpy_import_error(ValueError("something else"))


def test_warn_if_pertpy_missing_respects_skip_flags():
    with patch("perturbseq._deps.pertpy_available", return_value=False):
        assert warn_if_pertpy_missing(need_mixscape=False, need_distance=False, need_deseq2=False) is None
        msg = warn_if_pertpy_missing(need_mixscape=True, need_distance=True, need_deseq2=False)
        assert msg is not None
        assert "Mixscape" in msg
        assert "E-distance" in msg
        assert PERTPY_SKIP_FLAGS in msg


def test_run_mixscape_surfaces_install_guidance():
    from perturbseq.perturbation import run_mixscape

    adata = AnnData(
        X=np.zeros((20, 5)),
        obs=pd.DataFrame(
            {
                "perturbation": ["NT"] * 10 + ["GENE"] * 10,
                "gene_target": ["NT"] * 10 + ["GENE"] * 10,
            },
            index=[f"c{i}" for i in range(20)],
        ),
    )
    with patch("perturbseq.perturbation.require_pertpy", side_effect=ImportError(PERTPY_MISSING_MSG)):
        with pytest.raises(ImportError, match="pip install") as excinfo:
            run_mixscape(adata)
        assert PERTPY_SKIP_FLAGS in str(excinfo.value)


def test_pertpy_available_returns_bool():
    assert isinstance(pertpy_available(), bool)
