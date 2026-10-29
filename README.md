# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤，并跑 Mixscape 与 Mixscale。一次运行一个生物学样本。

步骤：`io` → `assign` → `preprocess` → **Mixscape ‖ Mixscale**（进程池并发）。跑完后 lean cleanup 只保留 PNG。不写 `checkpoints/`，不支持 `--resume`。细则见 [docs/](docs/README.md)、[docs/analysis.md](docs/analysis.md)、[docs/outputs.md](docs/outputs.md)。

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

## 命令行参数

| 参数 | 说明 |
| --- | --- |
| `--sample-id` | 样本名，用于拼接 DRAGEN 输入文件名 |
| `--organism` | 必填 `human` / `mouse` |
| `--input-dir` / `--output-dir` | DRAGEN 输入目录 / 结果输出目录 |
| `--guide-separator` | 多条 guide 之间的分隔符（如 `|`） |
| `--crispr-mode` | `KO`（默认）/ `i` / `a`。图标签对应 KO / KD / OE |
| `--nt-control` | 默认 `NEG_CONTROL/Non-Targeting/NegCtrl`，在参考表四列上做忽略大小写的 `re.search`；可用 `/` 或 `,` 自定义正则。未匹配 → 自动 `zero_guide`；也可显式填 `zero_guide`。匹配到但 QC 后无可用 NT 细胞 → 回退临时 `zero_guide` |
| `--max-mito-frac` | 覆盖默认 `0.15`（15%） |
| `--min-genes` / `--max-genes` | 覆盖默认 200 / 6000 |
| `--config` | 软件默认 YAML；缺省 `config/config.yaml` |
| `--sctype-tissue` | 启用 ScType；choices：`Adrenal`、`Brain`、`Eye`、`Heart`、`Hippocampus`、`Immune system`、`Intestine`、`Kidney`、`Liver`、`Lung`、`Muscle`、`Pancreas`、`Placenta`、`Spleen`、`Stomach`、`Thymus`（来自 `data/sctype/ScTypeDB_full.tsv`，GPL-3.0；[上游](https://github.com/IanevskiAleksandr/sc-type)） |
| `--sctype-db-path` | marker DB 路径（tsv/csv/xlsx）；默认 `data/sctype/ScTypeDB_full.tsv`；须含列 `tissueType`、`cellName`、`geneSymbolmore1`、`geneSymbolmore2` |
| `--sctype-low-score-fraction` | 低分簇 → `Unknown`（默认 0.25） |

QC 过滤只用基因数与线粒体比例，不做文库 UMI 绝对阈值过滤。线程默认 `max(1, CPU//3)`。

## 输出概览

| 目录 | 内容 |
| --- | --- |
| `qc/` | 过滤前 / 描述性 QC 统计图 |
| `UMAP/` | HVG/PCA/marker、表达 UMAP、SingleR / ScType、细胞周期；Mixscape signature |
| `Mixscape/` | Mixscape class / 2x2 等 |
| `mixscale/` | Mixscale ridge / scatter 等 |

图只写 PNG（`plotting.dpi` 默认 300）。嵌入用 QC-pass 细胞；Mixscape / Mixscale 用更严的 `analysis_eligible` ∪ 对照。产物清单见 [docs/outputs.md](docs/outputs.md)。
