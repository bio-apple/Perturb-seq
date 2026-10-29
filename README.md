# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤，并跑 Mixscape 与 Mixscale。一次运行一个生物学样本。

步骤：`io` → `assign` → `preprocess` → **Mixscape ‖ Mixscale**（进程池并发）。跑完后 lean cleanup 只保留 PNG。

## 环境搭建

前置要求：conda、Python >= 3.11（推荐 3.12）。R 4.4.3 + Bioconductor 3.20 仅 SingleR 簇注释需要（可选）。

```bash
conda create -y -n perturbseq-tertiary python=3.12 pip
conda activate perturbseq-tertiary
pip install -r requirements.txt
pip install -e . --no-deps

# SingleR 簇注释（可选；singler.enabled: false 可跳过）
conda install -y -c conda-forge r-base=4.4.3
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install(version = "3.20", ask = FALSE, update = FALSE); BiocManager::install(c("SingleR", "celldex"), version = "3.20", ask = FALSE, update = FALSE)'
```

## 最小运行示例

```bash
conda activate perturbseq-tertiary
python -m perturbseq \
  --sample-id SAMPLE \
  --organism human \
  --input-dir /path/to/sample \
  --output-dir /path/to/output \
  --guide-separator '|' \
  --nt-control 'NEG_CONTROL/Non-Targeting/NegCtrl'
```

完整参数与配置见 [docs/analysis.md](docs/analysis.md)，产物清单见 [docs/outputs.md](docs/outputs.md)。
