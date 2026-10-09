"""Compat alias for ``preprocessing`` — thin re-exports only.

Owns: nothing (no logic). Prefer ``perturbseq.preprocessing``.
Does NOT own: normalize/PCA/UMAP/Leiden (see ``preprocessing``).
"""

from perturbseq.preprocessing import preprocess_rna

__all__ = ["preprocess_rna"]
