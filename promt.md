从 DRAGEN scRNA CRISPR / Perturb-seq 二级分析结果出发

用 scanpy + pertpy 做三级分析。流程对齐 sc-best-practices · perturbation modeling、pertpy（Heumos et al. 2026）以及 Illumina Connected Multiomics Perturb-seq walkthrough。





# Data quality control
- Cells with < 200 expressed genes, classified as blank controls or containing a large fraction of mitochondrial genes (over 10%) were filtered.
- Genes expressed in less than three cells were filtered.
- reported that at least 30 cells are required to capture each perturbation phenotype

