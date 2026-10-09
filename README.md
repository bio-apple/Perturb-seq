# Perturb-seq 三级分析流程

从 DRAGEN scRNA CRISPR / Perturb-seq **二级分析**结果出发，用 **scanpy + pertpy** 做三级分析。流程对齐 [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/)、[pertpy](https://pertpy.scverse.org/en/latest/index.html)（[Heumos et al. 2026](https://www.nature.com/articles/s41592-025-02909-7)）以及 [Illumina Connected Multiomics Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq)。

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

步骤依赖、描述性 vs 推断性、矩阵层见 **[docs/ANALYSIS_DEPENDENCIES.md](docs/ANALYSIS_DEPENDENCIES.md)**。

```
DRAGEN MEX + guide assignment
  → 拆分 Gene Expression / CRISPR Direct Capture
  → 注释 guide_id、gene_target、NT vs perturbed
  → QC（n_counts / n_genes / %MT，MAD 自适应）                         [filtering]
  → 默认保留 num_features == 1 的单 guide 细胞                         [filtering]
  → normalize / log1p（保留 layers['counts']）+ HVG / PCA / Leiden     [descriptive]
  → UMAP                                                               [visualization only]
  → 细胞注释：细胞周期 + 状态打分（非组织细胞类型）                      [descriptive]
  → Mixscape：signature → KO/NP；可按 KO+NT 子集                       [filtering/classification]
  → E-distance / E-test（report 标注 pca_source；pre/post Mixscape）   [inferential]
  → 差异表达：有重复则 pseudobulk + PyDESeq2；否则探索性 Wilcoxon      [inferential/exploratory]
  → 按扰动平均 PCA 轮廓聚类
  → Guide QC：按 guide 评估 assignment / 靶效应 / 毒性与一致性
```

**统计边界**：细胞不是独立生物学重复。单样本、无 `replicate` 时，Wilcoxon 只能当探索，不能当组间结论。Mixscape 把“拿到 targeting guide 但转录组仍像对照”的细胞标成 NP；NP 与真正无效应无法仅靠转录组区分。全局 Leiden/UMAP **不是**扰动效应证据。QC / singlet / Mixscape KO 过滤会改变样本组成，见 `tables/composition_audit.csv`。同基因多 guide 的效应方向不一致时，不要把该基因的细胞简单合并当作单一效力结论。

对照识别默认把 `NegCtrl*`、`NT`/`NTC`、`non-targeting` 等映射为 `NT`。guide 名如 `PDCD10_4`、`ATM/design_3` 会解析成靶基因。可用 `--control-patterns` 改。

## 安装

```bash
conda env create -f environment.yml
conda activate perturbseq-tertiary
```

核心依赖：`scanpy`、`pertpy>=1.3`（Mixscape / Distance / PyDESeq2）。无 pertpy 时仍可跑 QC 与可视化，但需加 `--skip-mixscape --skip-distance`。

## 运行

先用合成 DRAGEN 文件检查环境：

```bash
python -m perturbseq write-demo --output-dir data/demo
python -m perturbseq run --input-dir data/demo --output-dir results/demo --sample-id sample1
```

用 YAML 统一参数（默认读 `configs/default.yaml`；CLI 覆盖同名键）：

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1 \
  --skip-mixscape --skip-distance
```

中断后从上一成功阶段恢复（输入 checksum 与阶段参数未变则跳过）：

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1 \
  --resume
```

真实二级分析结果：

```bash
python -m perturbseq run \
  --input-dir data/raw \
  --output-dir results/sample1 \
  --sample-id sample1 \
  --control NT
```

常用参数：

- `--config`：YAML 参数文件（避免 notebook/脚本默认值分叉）
- `--resume`：跳过已成功且 inputs/params 未变的阶段；变更会失效下游
- `--random-state`：随机种子（写入 manifest）
- `--keep-multiplets`：保留 0 或 ≥2 条 guide 的细胞
- `--replicate-col`：`obs` 中的生物学重复列；有重复才走 PyDESeq2
- `--skip-mixscape` / `--skip-cell-annotation` / `--skip-distance` / `--skip-de`
- `--perturbation-type KO`：CRISPRi 可改为 `KD`

对已有 h5ad 只补注释：

```bash
python -m perturbseq annotate \
  --h5ad results/A549_2/A549_2.tertiary.h5ad \
  --output-dir results/A549_2 \
  --inplace
```

对已有 tertiary h5ad 只跑 guide QC（无需重载 MTX）：

```bash
python -m perturbseq guide-qc \
  --h5ad results/A549_2/A549_2.tertiary.h5ad \
  --output-dir results/A549_2 \
  --control NT
```

从已有 `results/` 重建解读报告（无需重跑全流程）：

```bash
python -m perturbseq report --output-dir results/A549_2
```

## 输出

```
results/sample1/
  sample1.tertiary.h5ad
  report.json         含 perturbations：每个扰动的结构化摘要（机器可读）
  report.html         按「用户会问的问题」组织的可读报告
  run_manifest.json   可审计总账：输入 sha256、环境/版本、种子、参数、各阶段状态
  run.log             完整运行日志
  stages/<stage>/     status.json（inputs/outputs/checksums/params/success|failed）
                      + adata.h5ad 检查点（供 --resume）
  provenance_manifest.json  轻量 provenance（与 run_manifest 互补）
  figures/            QC、guide 组成、UMAP（可视化标注）、E-distance、火山图、guide_qc
  tables/
    composition_audit.csv   各过滤阶段细胞组成
    edistance.csv           主路径（Mixscape 成功时为 post-Mixscape）
    pre_mixscape/ / post_mixscape/   Mixscape 前后 E-distance（若跑 Mixscape）
    de_*.csv、etest.csv、perturbation_clusters.csv、
    guide_qc.csv、gene_guide_consistency.csv、qc_warnings.csv
```

阶段名与编排一致：`1_input_validation` → `2_preprocessing_qc` → `3a_perturbation_modeling` → `3b_statistical_inference` → `4_robustness` → `5_report`。同输入 + 同 config 重跑时，用 `run_manifest.json` / `stages/*/status.json` 对照每步 checksum 与参数，而不仅是再出一个 `.h5ad`。

`report.json` 另含 `matrix_provenance`、`mixscape_status`、`composition_notes`（见 [docs/ANALYSIS_DEPENDENCIES.md](docs/ANALYSIS_DEPENDENCIES.md)）。`report.html` / `report.json` → `perturbations` 对每个扰动回答：①细胞/guide/重复数与 QC；②相对 NT 的效应量、CI、校正 p；③跨 guide/样本/细胞状态一致性；④顶层受影响基因（通路富集可选，缺依赖则 N/A）；⑤局限与需人工复核的警告。CSV / `.h5ad` / JSON 仍是下游接口。

`obs` 中的关键列：`guide_id`、`gene_target`、`perturbation`、`num_features`、`guide_umi`、`mixscape_class`、`mixscape_class_global`（若 Mixscape 成功）。`layers['counts']` 为原始 UMI；`X` 为 log-norm；UMAP/Leiden 仅描述结构。

Guide QC 三层证据（写在 `tables/guide_qc.csv` 与 `report.json` → `guide_qc`）：

1. **assignment confidence**：细胞数、guide UMI、`num_features` / singlet 比例、跨 `sample_id` 检出率
2. **transcriptional effect**：靶基因相对 NT 的 log2FC（靶基因须在表达矩阵中）
3. **cytotoxicity / state**：`n_counts` / `n_genes` / `%MT` 相对对照，以及 `phase` / `cell_state`（若有）

`gene_guide_consistency.csv` 检查同基因多 guide 效应方向是否一致；`qc_warnings.csv` 标记细胞过少、UMI/检出偏低、效应离群、同基因 guide 不一致。`interpretation` 列用于区分「guide 不足 / 技术差」与「guide 合格但无靶效应」。

## 方法选型（为什么是这些，而不是预测模型）

| 步骤 | 选择 | 不选的 |
| --- | --- | --- |
| Guide 归属 | 沿用 DRAGEN GMM（Illumina 对该试剂盒的推荐） | 默认不再重跑 assignment；若要对比方法见 crispat |
| 无效扰动细胞 | Mixscape（Papalexi 2021；pertpy 实现） | 把所有 targeting 细胞都当 KO |
| 效应大小 | E-distance / E-test（Peidli 2024） | 只看 UMAP 是否分开 |
| 差异表达 | 有重复：pseudobulk + PyDESeq2（Squair 2021） | 把细胞当独立样本出组间 p 值 |
| 未见扰动预测 | 不作为默认步骤 | scGen/基础模型在独立基准上常不优于简单基线 |

PerturBase 适合查公开 Perturb-seq 数据集与可视化对照，不是本仓库的计算内核。

## 项目结构

对齐工程阶段图（Input → QC → 扰动建模 ∥ 统计推断 → 稳健性 → 报告）：

```
src/perturbseq/
  io.py              # 1. DRAGEN MEX / guide 参考 / metadata 校验
  guides.py          # 2. guide 解析、NT 识别；re-export guide_qc
  guide_qc.py        #    guide 级 QC / 一致性（实现保留，由 guides 汇总）
  qc.py              # 2. 细胞 & 样本 QC
  preprocessing.py   # 2. 标准化 / PCA / UMAP（preprocess.py 兼容别名）
  perturbation.py    # 3a. Mixscape、E-distance、guide consistency
  statistics.py      # 3b. pseudobulk DE、FDR、实验设计检查
  robustness.py      # 4. 敏感性 / 效应一致性 / confidence flags
  report.py          # 5. H5AD / JSON / HTML / provenance
  pipeline.py        # 薄编排器（阶段 1–5；3a∥3b 概念并行）
  annotate.py        # → guides 兼容 re-export
  analysis.py        # → perturbation + statistics 兼容 re-export
  cell_annotation.py # 细胞周期 / 状态注释（可选步骤）
  repro.py           # sha256 / YAML / RunTracker / resume / run_manifest
  cli.py             # python -m perturbseq …（--config / --resume）
configs/default.yaml # 集中默认参数
tests/
docs/structure.md    # 阶段与模块对照
```
