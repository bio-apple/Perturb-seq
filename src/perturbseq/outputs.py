"""Lean output policy: keep checkpoints and PNGs only."""

from __future__ import annotations

from pathlib import Path


def is_kept_output(path: Path, root: Path) -> bool:
    """True for checkpoints/*.h5ad and any PNG under the output root."""
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    if not path.is_file():
        return False
    suffix = path.suffix.lower()
    if suffix == ".png":
        return True
    if suffix == ".h5ad" and relative.parts and relative.parts[0] == "checkpoints":
        return True
    return False


def cleanup_lean_outputs(folder: Path) -> int:
    """Delete pipeline products that are neither checkpoints nor PNGs.

    Does not touch paths outside ``folder``. Returns the number of files removed.
    Empty directories left behind (except ``checkpoints/``) are pruned.
    """
    root = Path(folder)
    if not root.is_dir():
        return 0
    removed = 0
    for path in sorted(root.rglob("*"), reverse=True):
        if not path.is_file():
            continue
        if is_kept_output(path, root):
            continue
        path.unlink(missing_ok=True)
        removed += 1
    _prune_empty_dirs(root, keep_named={"checkpoints"})
    return removed


def _prune_empty_dirs(root: Path, *, keep_named: set[str]) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        if not path.is_dir():
            continue
        if path.name in keep_named and path.parent.resolve() == root.resolve():
            continue
        try:
            next(path.iterdir())
        except StopIteration:
            path.rmdir()
        except OSError:
            pass
