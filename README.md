# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤，并按 CRISPR 模式跑 Mixscape 或 Mixscale。一次运行一个生物学样本。

步骤：`io` → `assign` → `preprocess` → 分析步。`--steps all` 时：**KO → Mixscape**；**CRISPRi / CRISPRa → Mixscale**；**`other` → 仅 `io,assign,preprocess`**（不做 Mixscape / Mixscale；显式写这两步会拒绝）。KO/i/a 若命令行同时指定两者，仍可并发。图分目录 `Mixscape/`、`mixscale/` 等。成功跑完后磁盘**只保留 PNG**（不写 `checkpoints/`，不支持 `--resume`；须从 `io` 连续跑）。细则见 [docs/](docs/README.md)、[docs/analysis.md](docs/analysis.md)。

线程默认：`runtime.n_threads = max(1, CPU数//3)`，用于 OpenMP / BLAS / `scanpy.n_jobs` 与并行步进程数上限。

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

`NT_NAME` 换成参考表里的非靶向 guide 名。去掉 `--nt-control`，或指定了但数据中没有可用 NT 细胞时，使用 `zero_guide` 临时对照（结果里会标明 `temporary` / `zero_guide_cells`）。

| 常用参数 | 说明 |
| --- | --- |
| `--organism` | 必填 `human` / `mouse` → `experiment.organism` |
| `--crispr-mode` | `KO` / `i` / `a` / `other`（默认 `KO`） |
| `--steps` | 默认 `all`：KO→mixscape；i/a→mixscale；other→仅至 preprocess（须从 `io` 起；KO/i/a 显式同写两者时可并行） |
| 基因数 / mito | 默认 200–6000、mito≥20% 剔除；可用 `--min-genes` / `--max-genes` / `--max-mito-frac` 覆盖 |
| QC 指标 | `nCount` / `nFeature` / `percent.mt`（`nCount` 只作图，不按 UMI 绝对阈值过滤细胞） |

## 输出概览

成功结束后只保留各步 PNG：

| 目录 | 内容 |
| --- | --- |
| `qc/` | mRNA / gRNA / QC 统计图 |
| `UMAP/` | 表达 UMAP、簇注释、细胞周期、Mixscape signature |
| `Mixscape/` | Mixscape class / 2x2 等 |
| `mixscale/` | Mixscale ridge / scatter |

图只写 PNG（`plotting.dpi` 默认 300）。`plotting.save_plot_data` 默认 `false`。配置与规则见 [docs/analysis.md](docs/analysis.md)，产物清单见 [docs/outputs.md](docs/outputs.md)。
