import os
import tempfile
from pathlib import Path


def _cache_dir(name: str) -> str:
    candidates = [
        Path.home() / ".cache" / "perturbseq" / name,
        Path(tempfile.gettempdir()) / "perturbseq" / name,
    ]
    for path in candidates:
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        return str(path)
    raise RuntimeError(f"无法创建缓存目录 {name}")


for _key, _name in (("MPLCONFIGDIR", "mpl"), ("NUMBA_CACHE_DIR", "numba")):
    os.environ.setdefault(_key, _cache_dir(_name))

from perturbseq.pipeline import main  # noqa: E402

if __name__ == "__main__":
    main()
