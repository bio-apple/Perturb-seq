"""Optional heavy dependencies (pertpy / PyDESeq2) with explicit failure messages."""

from __future__ import annotations

from types import ModuleType

PERTPY_INSTALL = 'pip install -e ".[de]"  # or: pip install "pertpy[de]>=1.3"'
PERTPY_SKIP_FLAGS = "--skip-mixscape --skip-distance"
PERTPY_MISSING_MSG = (
    "pertpy is required for Mixscape, E-distance, and PyDESeq2. "
    f"Install with: {PERTPY_INSTALL}. "
    f"To continue without those steps, pass {PERTPY_SKIP_FLAGS} "
    "(Wilcoxon DE without replicates still works without pertpy)."
)


def is_pertpy_import_error(exc: BaseException) -> bool:
    """True when *exc* (or its cause chain) is a missing-pertpy ImportError."""
    seen: set[int] = set()
    cur: BaseException | None = exc
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        if isinstance(cur, ImportError):
            name = getattr(cur, "name", None)
            if name == "pertpy" or (isinstance(name, str) and name.startswith("pertpy.")):
                return True
            text = str(cur).lower()
            if "pertpy" in text:
                return True
        cur = cur.__cause__ or cur.__context__
    return False


def pertpy_available() -> bool:
    try:
        import pertpy  # noqa: F401
    except ImportError:
        return False
    return True


def require_pertpy() -> ModuleType:
    """Import pertpy or raise ImportError with install / skip-flag guidance."""
    try:
        import pertpy as pt
    except ImportError as exc:
        raise ImportError(PERTPY_MISSING_MSG) from exc
    return pt


def warn_if_pertpy_missing(*, need_mixscape: bool, need_distance: bool, need_deseq2: bool) -> str | None:
    """Return a user-facing warning string when pertpy is needed but absent."""
    if not (need_mixscape or need_distance or need_deseq2):
        return None
    if pertpy_available():
        return None
    parts = []
    if need_mixscape:
        parts.append("Mixscape")
    if need_distance:
        parts.append("E-distance")
    if need_deseq2:
        parts.append("PyDESeq2 (replicate-aware DE)")
    return (
        f"pertpy is not installed; requested step(s) will fail or degrade: {', '.join(parts)}. "
        f"{PERTPY_MISSING_MSG}"
    )


__all__ = [
    "PERTPY_INSTALL",
    "PERTPY_MISSING_MSG",
    "PERTPY_SKIP_FLAGS",
    "is_pertpy_import_error",
    "pertpy_available",
    "require_pertpy",
    "warn_if_pertpy_missing",
]
