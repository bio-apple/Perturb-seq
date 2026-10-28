ScTypeDB_full.tsv — official ScType marker database (Ianevski et al.),
converted from ScTypeDB_full.xlsx for offline use (tsv; no Excel dependency).
Source: https://github.com/IanevskiAleksandr/sc-type
License: GPL-3.0 (upstream).

Pipeline usage (perturbseq.sctype; coexist with SingleR):
  config: sctype.enabled (default true), sctype.tissue (required to run),
          sctype.db_path (default: this directory's ScTypeDB_full.tsv),
          sctype.low_score_fraction (default 0.25 → Unknown clusters)
  obs column: sctype_classification
  figure: UMAP/umap_sctype.png (written in preprocess)
  If tissue is unset/null, ScType is skipped (no error exit).

Tissues in ScTypeDB_full.tsv:
  Adrenal, Brain, Eye, Heart, Hippocampus, Immune system, Intestine,
  Kidney, Liver, Lung, Muscle, Pancreas, Placenta, Spleen, Stomach, Thymus

Required columns: tissueType, cellName, geneSymbolmore1, geneSymbolmore2
Custom DB: same columns; formats tsv/csv/xlsx via sctype.db_path.
