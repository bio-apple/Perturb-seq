# Perturb-seq 三级分析流程

[English](README.md) | 简体中文

从 DRAGEN scRNA CRISPR / Perturb-seq **二级分析**结果出发，用 **scanpy + pertpy** 做三级分析。流程对齐 [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/)、[pertpy](https://pertpy.scverse.org/en/latest/index.html)（[Heumos et al. 2026](https://www.nature.com/articles/s41592-025-02909-7)）以及 [Illumina Connected Multiomics Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq)。

> 英文主 README：[README.md](README.md)（安装包元数据指向该文件）。本页覆盖安装 / 运行 / 输出；统计局限与 flag 对照以 [STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md)、[CLI_YAML.md](docs/CLI_YAML.md) 为准（中英文 README 不重复展开）。

## 输入（DRAGEN 二级分析）

每个样本（示例 `sample1`）需要：

| 文件 | 作用 |
| --- | --- |
| `{sample}.scRNA.filtered.matrix.mtx.gz` | 过滤后细胞 × 特征计数（基因 + CRISPR） |
| `{sample}.scRNA.filtered.barcodes.tsv.gz` | 细胞 barcode |
| `{sample}.scRNA.filtered.features.tsv.gz` | 基因/gRNA ID、名称、feature type |
| `{sample}.scRNA.feature_barcode_reference.csv` | gRNA 参考（`id,name,...,feature_type`） |
| `{sample}.scRNA.positive_cell_guide_assignments.csv` | DRAGEN GMM guide calling：`cell_barcode,num_features,feature_call,num_transcripts` |

字段定义见 [DRAGEN v4.5 CRISPR mode](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode)。`feature_call` / `num_transcripts` 在多 guide 时用 `|` 分隔。

把上述文件放到一个目录，例如 `data/raw/`。

## 流程（先看数据，再定阈值）

文档入口：**[docs/tutorial.ipynb](docs/tutorial.ipynb)**（交互 notebook，合成 demo 可离线）· **[docs/tutorial.md](docs/tutorial.md)**（CLI 对照）· **[docs/CLI_YAML.md](docs/CLI_YAML.md)**（CLI ↔ YAML）· **[docs/STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md)**（统计局限）· **[docs/ANALYSIS_DEPENDENCIES.md](docs/ANALYSIS_DEPENDENCIES.md)**（步骤依赖 / 矩阵层）· **[docs/SC_BEST_PRACTICES_MAP.md](docs/SC_BEST_PRACTICES_MAP.md)**（流水线步骤 ↔ sc-best-practices 章节）。

```
DRAGEN MEX + guide assignment
  → 拆分 Gene Expression / CRISPR Direct Capture
  → 注释 guide_id、gene_target、NT vs perturbed
  → QC（MAD；推荐 --perturbation-aware-qc）→ 默认 singlet（num_features == 1）
  → normalize / log1p（保留 layers['counts']）+ HVG / PCA / Leiden / UMAP
  → 细胞注释（周期 + 状态）→ Mixscape（KO/NP；基因组规模可 auto-skip）
  → E-distance / E-test（含 bootstrap CI / low_power；pre/post Mixscape）
  → DE：优先 pseudobulk + PyDESeq2（真重复→inferential；否则 exploratory）；Wilcoxon 仅 <2 个 pseudobulk 单元时回退
  → 扰动空间（--perturbation-space：pca_silhouette | kmeans | lr_classifier）
  → Guide QC + 报告图（guide-consistency / target-validation PNG）
```

对照默认映射 `NegCtrl*` / `NT`/`NTC` / `non-targeting` → `NT`（`--control-patterns` 可改）。统计边界与文献只维护一份：[STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md)。

### 常见坑（短表）

完整讨论：[STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md) · [tutorial.md](docs/tutorial.md#常见坑--common-pitfalls)。

| 坑 | 简述 |
| --- | --- |
| 同基因多 guide 不一致 | 勿合并成单一「基因 KO」；看 `gene_guide_consistency.csv` |
| 多 guide on-target 失败 | `potential_low_efficiency` ≠「无表型」 |
| Mixscape 大量 NP | 转录组分类 ≠ 已证实无效编辑 |
| 无真重复谈 FDR | 仅 exploratory；`evidence_level=inferential` 需 ≥2 真重复 |
| UMAP/Leiden 当效应 | 不是扰动证据；优先 E-distance / 有重复时的 DE |
| 过滤改变命中 | 看 `composition_audit.csv`；推荐 `--perturbation-aware-qc` |
| 默认信任 DRAGEN assignment | 三级流程默认不重 call guide |

## 安装

**默认推荐**（含 `pertpy[de]`：Mixscape / E-distance / E-test / PyDESeq2）：

```bash
# pip（约束见 constraints.txt，与 pyproject / environment.yml 对齐）
pip install -e ".[de,dev]" -c constraints.txt

# 或 conda（conda-forge 装核心科学栈，再 editable 装 .[de,dev]）
conda env create -f environment.yml
conda activate perturbseq-tertiary
```

**轻量 / 仅核心**（scanpy / anndata / Leiden 等，**不含** pertpy）：

```bash
pip install -e ".[dev]" -c constraints.txt
# 生产跑数但不装开发工具：pip install -e . -c constraints.txt
```

安装后可用 `python -m perturbseq …` 或入口脚本 `perturbseq` / `perturbseq-tertiary`。

本地 lint hook（可选）：`pre-commit install`（配置见 `.pre-commit-config.yaml`，规则对齐 `pyproject.toml` 的 ruff/mypy）。

### 已测试版本 / 兼容性

本流程曾用 **pertpy==1.4.0** 与 **scanpy==1.12.4** 做端到端验证（conda 环境 `perturbseq-tertiary`，Python 3.12.15；版本来自 `importlib.metadata`）。下表中的支持范围为 `pyproject.toml` / `constraints.txt` / `environment.yml` 声明——排查安装问题时优先对照这些文件，而不是已漂移的本地环境。

| 包 | 已测试（E2E） | 支持范围 |
| --- | --- | --- |
| Python | 3.12.15 | `>=3.10` |
| scanpy | 1.12.4 | `>=1.10,<2` |
| pertpy | 1.4.0 | `>=1.3,<2`（`pertpy` / `de` / `all` extras） |
| anndata | 0.13.4 | `>=0.10,<0.14` |

### 无 pertpy 时的降级路径

| 能力 | 有 `.[de]`（`pertpy[de]`） | 仅核心（无 pertpy） |
| --- | --- | --- |
| 输入校验、guide 注释、QC / MAD 过滤、singlet | ✅ | ✅ |
| normalize / HVG / PCA / Leiden / UMAP、细胞注释 | ✅ | ✅ |
| Guide QC（assignment / 靶效应 / 毒性一致性） | ✅ | ✅ |
| 探索性 DE（无真重复：pseudobulk PyDESeq2 / Wilcoxon 回退） | ✅（PyDESeq2 需 pertpy） | ✅ 仅 Wilcoxon；可加 `--skip-de` |
| Mixscape（KO/NP 分类、`X_pert`） | ✅ | ❌ → `mixscape.reason=missing_pertpy`；请加 `--skip-mixscape` |
| E-distance / E-test | ✅ | ❌ → `edistance.reason=missing_pertpy`；请加 `--skip-distance` |
| Pseudobulk + PyDESeq2（有/无 `replicate_col`） | ✅ | ❌ → 无法走 PyDESeq2；无 pertpy 时 Wilcoxon 回退 |
| 扰动空间聚类（`--perturbation-space`） | ✅ `pca_silhouette`（默认）/ `kmeans` / `lr_classifier` | 成功时 `skipped: false` + `method`；失败 → `reason=failed` |

缺 pertpy 且未加 `--skip-mixscape --skip-distance` 时，流程会发出明确警告，并在 `report.json` 写入显式状态块（如 `"mixscape": {"skipped": true, "reason": "missing_pertpy", ...}`），而不是静默省略字段或裸 `null`。无 pertpy 的推荐命令：

```bash
python -m perturbseq run \
  --input-dir data/demo --output-dir results/demo --sample-id sample1 \
  --skip-mixscape --skip-distance
```

## 运行

零基础逐步教程：交互 **[docs/tutorial.ipynb](docs/tutorial.ipynb)**（推荐）或 CLI **[docs/tutorial.md](docs/tutorial.md)**。

```bash
jupyter notebook docs/tutorial.ipynb
# 或: jupyter lab docs/tutorial.ipynb
```

先用合成 DRAGEN 文件检查环境：

```bash
python -m perturbseq write-demo --output-dir data/demo
python -m perturbseq run --input-dir data/demo --output-dir results/demo --sample-id sample1
```

YAML（默认 `configs/default.yaml`；CLI 覆盖同名键；完整对照：**[CLI_YAML.md](docs/CLI_YAML.md)**）+ `--resume`：

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1 \
  --resume
```

真实数据：`--input-dir data/raw --output-dir results/sample1 --sample-id sample1 --control NT`。

常用参数（细节见 [CLI_YAML.md](docs/CLI_YAML.md)）：`--config`、`--resume`、`--dry-run`、`--random-state`、`--n-jobs`（默认 `1`；`-1` = 全 CPU；并行 DE/E-test 外层循环）、`--keep-multiplets`、`--perturbation-aware-qc`、`--replicate-col`（≥2 真重复 → `evidence_level=inferential`）、`--de-prefer-pseudobulk`（默认 true，无真重复也优先 exploratory PyDESeq2）、`--pseudo-replicates`、`--n-bootstrap`、`--etest-power-min-cells`、`--skip-mixscape` / `--skip-distance` / `--skip-de` / `--skip-cell-annotation`、`--mixscape-mode` / `--mixscape-targets` / `--mixscape-top-n`、`--perturbation-type`、`--perturbation-space`、`--guide-merge`、`--on-target-*`、`--guide-reassign`、`--de-covariates`、`--report-plot-top-n`。

独立子命令：`annotate`、`guide-qc`、`report`（从已有 `results/` 重建 HTML/JSON / 报告图，无需重跑全流程）。

## 输出

```
results/sample1/
  sample1.tertiary.h5ad
  report.json / report.html     # verdict + 每扰动卡片；DE 带 evidence_level
  run_manifest.json / run.log
  stages/<stage>/status.json    # inputs / params_hash（供 --resume）
  figures/                      # QC、UMAP、E-distance、火山图、guide_consistency_*、target_validation_*
  tables/
    composition_audit.csv, qc_filter_by_guide.csv
    edistance.csv, etest.csv, distances.csv, de_*.csv
    guide_qc.csv, gene_guide_consistency.csv, qc_warnings.csv
    gene_guide_weighted.csv     # guide_merge ≠ none 时
    perturbation_clusters.csv
    pre_mixscape/ / post_mixscape/
```

阶段：`1_input_validation` → `2_preprocessing_qc` → `3a_perturbation_modeling` → `3b_statistical_inference` → `4_robustness` → `5_report`。章节对照：[SC_BEST_PRACTICES_MAP.md](docs/SC_BEST_PRACTICES_MAP.md)。矩阵层 / 可选步骤 `skipped`+`reason`：[ANALYSIS_DEPENDENCIES.md](docs/ANALYSIS_DEPENDENCIES.md)。模块对照：[structure.md](docs/structure.md)。

关键 `obs`：`guide_id`、`gene_target`、`perturbation`、`num_features`、`guide_umi`、`mixscape_class*`。Guide QC 三层（assignment / 靶效应 / 毒性）、on-target、`guide_merge`、`guide_reassign` 的解读见 [STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md)。

## 方法选型

| 步骤 | 选择 | 不选的 |
| --- | --- | --- |
| Guide 归属 | DRAGEN GMM；可选 `--guide-reassign` | 默认不重 call；CatchR / Cell Ranger FB 未安装 |
| 无效扰动细胞 | Mixscape（Papalexi 2021；KO 假设更强） | 把所有 targeting 细胞都当 KO |
| 效应大小 | E-distance / E-test（Peidli 2024；可 `X_pert`） | 只看 UMAP 是否分开 |
| 差异表达 | 优先 pseudobulk + PyDESeq2；Wilcoxon 仅回退 | 静默把细胞级 Wilcoxon 当人口 FDR |
| 未见扰动预测 | 不作为默认步骤 | scGen / 基础模型作主路径 |

## 项目结构

```
src/perturbseq/   # io, guides, qc, preprocessing, perturbation, statistics, report, report_plots, …
configs/default.yaml
docs/             # tutorial.ipynb / tutorial.md / CLI_YAML / STATISTICAL_CAVEATS / …
tests/
```

Lint：

```bash
pip install -e ".[dev]" -c constraints.txt
ruff check src tests
python -m mypy src/perturbseq/qc.py src/perturbseq/perturbation.py src/perturbseq/statistics.py src/perturbseq/_deps.py
```
