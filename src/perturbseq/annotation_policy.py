"""Sample-type policy for whether full cell-cycle / state annotation runs.

Clustering / UMAP always run upstream regardless of this decision.
Kept free of scanpy/pertpy imports so unit tests stay lightweight.
"""

from __future__ import annotations

from typing import Any

SAMPLE_TYPES = ("cell_line", "primary", "mixed", "unknown")


def resolve_cell_annotation_decision(
    sample_type: str = "cell_line",
    *,
    skip_cell_annotation: bool | None = None,
) -> dict[str, Any]:
    """Decide whether to run cell-cycle / state annotation.

    Precedence: explicit ``skip_cell_annotation`` override (True=skip, False=run)
    always wins over ``sample_type``.

    Default policy (when override is None)::

        cell_line → skip  (homogeneous lines; Leiden/UMAP still available)
        primary   → run
        mixed     → run
        unknown   → skip  (Perturb-seq demos are usually cell lines; note recorded)

    Returns a dict with keys: ``skip``, ``reason``, ``detail``, ``sample_type``,
    ``source`` (``override`` | ``sample_type``).
    """
    st = (sample_type or "unknown").strip().lower()
    if st not in SAMPLE_TYPES:
        raise ValueError(f"sample_type must be one of {SAMPLE_TYPES}, got {sample_type!r}")

    if skip_cell_annotation is True:
        return {
            "skip": True,
            "reason": "user_skip",
            "detail": "Cell annotation skipped (--skip-cell-annotation / skip_cell_annotation=true).",
            "sample_type": st,
            "source": "override",
        }
    if skip_cell_annotation is False:
        return {
            "skip": False,
            "reason": "user_run",
            "detail": "Cell annotation forced (--run-cell-annotation / skip_cell_annotation=false).",
            "sample_type": st,
            "source": "override",
        }

    if st == "cell_line":
        return {
            "skip": True,
            "reason": "sample_type=cell_line",
            "detail": (
                "Homogeneous cell-line Perturb-seq: skip full cell-cycle/state annotation; "
                "Leiden/UMAP clustering still ran. Use --run-cell-annotation to force."
            ),
            "sample_type": st,
            "source": "sample_type",
        }
    if st in ("primary", "mixed"):
        return {
            "skip": False,
            "reason": f"sample_type={st}",
            "detail": f"Sample type {st!r}: run cell-cycle / state annotation.",
            "sample_type": st,
            "source": "sample_type",
        }
    # unknown — skip with note (least surprising for typical Perturb-seq cell-line demos)
    return {
        "skip": True,
        "reason": "sample_type=unknown",
        "detail": (
            "sample_type=unknown: skipping cell annotation (Perturb-seq default assumes "
            "cell-line-like demos). Set --sample-type primary|mixed or --run-cell-annotation "
            "if tissue / mixed populations need labels."
        ),
        "sample_type": st,
        "source": "sample_type",
    }
