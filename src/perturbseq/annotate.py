"""Backward-compatible re-exports; prefer ``perturbseq.guides``."""

from perturbseq.guides import (
    DEFAULT_CONTROL_PATTERNS,
    annotate_guides,
    build_feature_target_map,
    filter_singlets,
    parse_gene_target,
    split_feature_calls,
)

__all__ = [
    "DEFAULT_CONTROL_PATTERNS",
    "annotate_guides",
    "build_feature_target_map",
    "filter_singlets",
    "parse_gene_target",
    "split_feature_calls",
]
