# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤，并跑 Mixscape 与 Mixscale。一次运行一个生物学样本。

步骤：`io` → `assign` → `preprocess` → **Mixscape ‖ Mixscale**（进程池并发；图目录 `Mixscape/`、`mixscale/`）。跑完后 lean cleanup：**主要保留 PNG**；`--output-dir` 与 `--input-dir` 相同时额外保护 `{sample}.scRNA.*` 与已解析的 DRAGEN 输入。不写 `checkpoints/`，不支持 `--resume`；须从 `io` 连续跑。细则见 [docs/](docs/README.md)、[docs/analysis.md](docs/analysis.md)、[docs/outputs.md](docs/outputs.md)。

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

簇注释：**SingleR**（R 4.4.3 + Bioconductor 3.20 / celldex；`singler.enabled: false` 可跳过；首次下载参考需联网）与纯 Python **ScType** 并存。ScType 仅 CLI（不在 `config.yaml`），参数见下表。

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

| 常用参数 | 说明 |
| --- | --- |
| `--organism` | 必填 `human` / `mouse` → `experiment.organism` |
| `--crispr-mode` | `KO` / `i` / `a`（默认 `KO`）。写入预期方向元数据与图标签（KO / KD / OE）：KO→`not_required`，i→`down`，a→`up`。步骤固定跑 Mixscape 与 Mixscale |
| `--nt-control` | 默认 `NEG_CONTROL/Non-Targeting/NegCtrl`：在参考表 `id` / `name` / `target_gene_name` / `target_gene_id` 上做忽略大小写的 `re.search`；可用 `/` 或 `,` 自定义正则。未匹配 → 自动 `zero_guide`；也可填 `zero_guide`。匹配到但 QC 后无可用 NT 细胞 → 回退临时 `zero_guide` |
| `--sctype-tissue` | 启用 ScType；未传则跳过。choices：`Adrenal`、`Brain`、`Eye`、`Heart`、`Hippocampus`、`Immune system`、`Intestine`、`Kidney`、`Liver`、`Lung`、`Muscle`、`Pancreas`、`Placenta`、`Spleen`、`Stomach`、`Thymus`（来自 `data/sctype/ScTypeDB_full.tsv`，GPL-3.0；[上游](https://github.com/IanevskiAleksandr/sc-type)）。产物：`obs['sctype_classification']`、`UMAP/umap_sctype.png` |
| `--sctype-db-path` | 可选自定义 marker DB（tsv/csv/xlsx）；须含列 `tissueType`、`cellName`、`geneSymbolmore1`、`geneSymbolmore2` |
| `--sctype-low-score-fraction` | 可选；低分簇 → `Unknown`（默认 0.25） |
| 基因数 / mito | 默认 200–6000、mito≥20% 剔除；可用 `--min-genes` / `--max-genes` / `--max-mito-frac` 覆盖 |
| QC 指标 | `nCount` / `nFeature` / `percent.mt`（`nCount` 只作图，不按 UMI 绝对阈值过滤细胞） |

## 输出概览

成功结束后磁盘以各步 PNG 为主（中间 TSV/JSON/YAML 会被 lean cleanup 删除）：

| 目录 | 内容 |
| --- | --- |
| `qc/` | assign：过滤前 / 描述性 QC 统计图 |
| `UMAP/` | QC-pass：HVG/PCA/marker、表达 UMAP、SingleR / ScType、细胞周期；Mixscape signature |
| `Mixscape/` | Mixscape class / 2x2 等 |
| `mixscale/` | Mixscale ridge / scatter 等 |

图只写 PNG（`plotting.dpi` 默认 300）。`plotting.save_plot_data` 默认 `false`。嵌入用 **QC-pass** 细胞；Mixscape / Mixscale 用更严的 `analysis_eligible` ∪ 对照。配置与规则见 [docs/analysis.md](docs/analysis.md)，产物清单见 [docs/outputs.md](docs/outputs.md)。
