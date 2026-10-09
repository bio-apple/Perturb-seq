from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from perturbseq.io import SAMPLE_FILES


def write_demo_dragen(output_dir: Path, n_cells: int = 120, n_genes: int = 80, seed: int = 0) -> Path:
    """Write a tiny DRAGEN-like MEX + CRISPR assignment bundle for local testing.

    Plants target gene symbols in the expression matrix so guide QC can score
    target-gene LFC. IFNGR2 has two guides with conflicting target effects
    (strong knockdown vs opposing/weak) so ``gene_guide_consistency`` flags
    multi-guide inconsistency — a hands-on pitfall for the tutorial.
    """
    rng = np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample = "sample1"
    guides = [
        ("NegCtrl_0001", "NT"),
        ("IFNGR2_1", "IFNGR2"),
        ("IFNGR2_2", "IFNGR2"),
        ("STAT1_1", "STAT1"),
        ("JAK2_1", "JAK2"),
        ("PDCD10_1", "PDCD10"),
    ]
    # Real target symbols first so guide_qc can resolve target_in_matrix / LFC.
    target_symbols = ["IFNGR2", "STAT1", "JAK2", "PDCD10"]
    n_mt = 8
    n_filler = max(0, n_genes - len(target_symbols) - n_mt)
    genes = target_symbols + [f"GENE{i:03d}" for i in range(n_filler)] + [f"MT-GENE{i}" for i in range(n_mt)]
    if len(genes) != n_genes:
        raise ValueError(f"n_genes={n_genes} too small for demo layout (need ≥{len(target_symbols) + n_mt})")
    barcodes = [f"CELL{i:04d}" for i in range(n_cells)]
    n_features = n_genes + len(guides)
    matrix = rng.negative_binomial(4, 0.4, size=(n_features, n_cells)).astype(np.float32)
    # NT-like baseline for planted targets (overwritten per guide below).
    for ti in range(len(target_symbols)):
        matrix[ti, :] = rng.integers(10, 16, size=n_cells).astype(np.float32)
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
            # Generic module shift on filler genes (leave planted targets for LFC logic).
            start = len(target_symbols)
            end = min(start + 12, n_genes - n_mt)
            width = end - start
            if width > 0:
                matrix[start:end, i] = matrix[start:end, i] + rng.integers(4, 12, size=width)
            if target in target_symbols:
                ti = target_symbols.index(target)
                if guide_id == "IFNGR2_1":
                    # Strong knockdown → target_effect_direction "down"
                    matrix[ti, i] = float(rng.integers(1, 3))
                elif guide_id == "IFNGR2_2":
                    # Conflicting / weak guide → "up" vs NT baseline (inconsistency)
                    matrix[ti, i] = float(rng.integers(22, 30))
                else:
                    # Single-guide targets: mild consistent knockdown
                    matrix[ti, i] = float(rng.integers(1, 4))
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
