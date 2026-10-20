# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤、描述性 CPM log2FC、Mixscape / Mixscale 和扰动距离。一次运行一个生物学样本。

文档索引见 [docs/README.md](docs/README.md)。

## 环境搭建

已安装 conda。在仓库根目录执行：

```bash
conda create -y -n perturbseq-tertiary python=3.12 pip
conda activate perturbseq-tertiary
pip install -r requirements.txt
pip install -e . --no-deps
conda install -y -c conda-forge r-base=4.4.3
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install(version = "3.20", ask = FALSE, update = FALSE); BiocManager::install(c("SingleR", "celldex"), version = "3.20", ask = FALSE, update = FALSE)'
```

簇注释用 R 4.4.3 + Bioconductor 3.20 的 SingleR / celldex；`singler.enabled: false` 可跳过。celldex 首次下载参考需联网。

## 最小运行示例

```bash
conda activate perturbseq-tertiary
python -m perturbseq \
  --sample-id SAMPLE \
  --organism human \
  --input-dir /path/to/sample \
  --output-dir /path/to/output \
  --guide-separator '|' \
  --nt-control NT_NAME
```

`NT_NAME` 换成参考表里的非靶向 guide 名。去掉 `--nt-control` 时使用 `zero_guide` 临时对照。

## 输出概览

- `qc/`：mRNA / gRNA 统计
- `UMAP/`：表达 UMAP、簇注释、细胞周期、Mixscape signature
- `DGE/Mixscape/`：Mixscape 门槛、摘要、class / 2x2 / fraction
- `DGE/mixscale/`：Mixscale ridge / scatter、门槛
- `DGE/descriptive/`：描述性 CPM log2FC 与相关效应表
- `E-distance/`：对当前对照的能量距离与 guide 伪批量距离

图只写 PNG（无 PDF、无同名图旁 txt）。细节见 [docs/](docs/README.md)。


