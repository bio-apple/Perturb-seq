从 DRAGEN scRNA CRISPR / Perturb-seq 二级分析结果出发:https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode

输入:

sample1.scRNA.filtered.matrix.mtx.gz
sample1.scRNA.filtered.barcodes.tsv.gz
sample1.scRNA.filtered.features.tsv.gz
sample1.scRNA.feature_barcode_reference.csv
sample1.scRNA.positive_cell_guide_assignments.csv

用 scanpy + pertpy+Seurat 做三级分析，流程对齐 
sc-best-practices perturbation modeling:https://www.sc-best-practices.org/conditions/perturbation-modeling/
Illumina Connected Multiomics Perturb-seq walkthrough:https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq

三级分析内容必须包括：

| 分析模块 | 推荐输出 | 优先级 |
| --- | --- | --- |
| **细胞 QC** | QC violin、scatter、过滤前后细胞数统计 | 必做 |
| **gRNA QC** | guide 细胞数柱状图、每细胞 guide 数分布、assignment 汇总 | 必做 |
| **UMAP / clustering** | cluster UMAP、marker 表、细胞状态注释 | 必做 |
| **Perturbation 映射** | 按 guide、target gene、NTC 着色的 UMAP | 必做 |
| **细胞组成变化** | 各扰动在 cluster 中的比例及计数 | 推荐 |
| **差异表达** | 每个扰动对照 NTC 的 DEG 表、火山图、热图 | 有合适对照时 |
| **通路分析** | GO / GSEA 富集结果、富集图 | 有可靠 DEG 或排序结果时 |


只针对单样本

Data quality control
Cells with <200 expressed genes, classified as blank controls or containing a large fraction of mitochondrial genes (over 10%) were filtered
Genes expressed in less than three cells were filtered
reported that at least 30cells are required to capture each perturbation phenotype