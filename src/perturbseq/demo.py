from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from perturbseq.io import SAMPLE_FILES


def write_demo_dragen(output_dir: Path, n_cells: int = 120, n_genes: int = 80, seed: int = 0) -> Path:
    """Write a tiny DRAGEN-like MEX + CRISPR assignment bundle for local testing."""
    rng = np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample = "sample1"
    genes = [f"GENE{i:03d}" for i in range(n_genes - 8)] + [f"MT-GENE{i}" for i in range(8)]
    guides = [
        ("NegCtrl_0001", "NT"),
        ("IFNGR2_1", "IFNGR2"),
        ("IFNGR2_2", "IFNGR2"),
        ("STAT1_1", "STAT1"),
        ("JAK2_1", "JAK2"),
        ("PDCD10_1", "PDCD10"),
    ]
    barcodes = [f"CELL{i:04d}" for i in range(n_cells)]
    n_features = n_genes + len(guides)
    matrix = rng.negative_binomial(4, 0.4, size=(n_features, n_cells)).astype(np.float32)
    assignments = []
    for i, barcode in enumerate(barcodes):
        roll = rng.random()
        if roll < 0.08:
            g1, g2 = rng.choice(len(guides), size=2, replace=False)
            assignments.append((barcode, 2, f"{guides[g1][0]}|{guides[g2][0]}", f"{rng.integers(5, 40)}|{rng.integers(5, 40)}"))
            continue
        if roll < 0.16:
            assignments.append((barcode, 0, "", ""))
            continue
        gi = int(i % len(guides))
        guide_id, target = guides[gi]
        umi = int(rng.integers(8, 60))
        assignments.append((barcode, 1, guide_id, str(umi)))
        if target != "NT":
            module = slice(0, 12)
            matrix[module, i] = matrix[module, i] + rng.integers(4, 12, size=12)
        matrix[n_genes + gi, i] = umi
    mtx_path = output_dir / SAMPLE_FILES["matrix"].format(sample=sample)
    tmp = output_dir / "_matrix.mtx"
    mmwrite(tmp, sparse.csr_matrix(matrix))
    import gzip
    import shutil

    with tmp.open("rb") as src, gzip.open(mtx_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp.unlink()

    pd.Series(barcodes).to_csv(
        output_dir / SAMPLE_FILES["barcodes"].format(sample=sample),
        index=False,
        header=False,
        compression="gzip",
    )
    feature_rows = [[f"ENSG{i:03d}", genes[i], "Gene Expression"] for i in range(n_genes)]
    feature_rows += [[guide_id, guide_id, "CRISPR Direct Capture"] for guide_id, _ in guides]
    pd.DataFrame(feature_rows).to_csv(
        output_dir / SAMPLE_FILES["features"].format(sample=sample),
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
    ref = pd.DataFrame(
        {
            "id": [g[0] for g in guides],
            "name": [g[0] for g in guides],
            "read": "R2",
            "target": "5p",
            "position": -20,
            "sequence": ["A" * 20] * len(guides),
            "feature_type": "CRISPR Direct Capture",
        }
    )
    ref.to_csv(output_dir / SAMPLE_FILES["feature_ref"].format(sample=sample), index=False)
    assign_df = pd.DataFrame(assignments, columns=["cell_barcode", "num_features", "feature_call", "num_transcripts"])
    assign_df.to_csv(output_dir / SAMPLE_FILES["assignments"].format(sample=sample), index=False)
    return output_dir
