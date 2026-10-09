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

```
DRAGEN MEX + guide assignment
  → 拆分 Gene Expression / CRISPR Direct Capture
  → 注释 guide_id、gene_target、NT vs perturbed
  → QC（n_counts / n_genes / %MT，MAD 自适应，不套固定 10%）
  → 默认保留 num_features == 1 的单 guide 细胞（与 ICM 一致；多 guide 按实验设计决定）
  → normalize / log1p / HVG / PCA / neighbors / UMAP / Leiden
  → Mixscape：perturbation signature → KO/NP 分类 → LDA；并尝试 Mixscale 连续打分
  → E-distance / E-test（相对 NT）
  → 差异表达：有生物学重复则 pseudobulk + PyDESeq2；单样本仅探索性 Wilcoxon
  → 按扰动平均 PCA 轮廓聚类（相似通路/复合物）
```

**统计边界**：细胞不是独立生物学重复。单样本、无 `replicate` 时，Wilcoxon 只能当探索，不能当组间结论。Mixscape 把“拿到 targeting guide 但转录组仍像对照”的细胞标成 NP；NP 与真正无效应无法仅靠转录组区分。

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

真实二级分析结果：

```bash
python -m perturbseq run \
  --input-dir data/raw \
  --output-dir results/sample1 \
  --sample-id sample1 \
  --control NT
```

常用参数：

- `--keep-multiplets`：保留 0 或 ≥2 条 guide 的细胞
- `--replicate-col`：`obs` 中的生物学重复列；有重复才走 PyDESeq2
- `--skip-mixscape` / `--skip-distance` / `--skip-de`
- `--perturbation-type KO`：CRISPRi 可改为 `KD`

## 输出

```
results/sample1/
  sample1.tertiary.h5ad
  report.json
  figures/          QC、guide 组成、UMAP、E-distance、火山图
  tables/           edistance.csv、etest.csv、de_*.csv、perturbation_clusters.csv
```

`obs` 中的关键列：`guide_id`、`gene_target`、`perturbation`、`num_features`、`guide_umi`、`mixscape_class`、`mixscape_class_global`（若 Mixscape 成功）。

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

```
src/perturbseq/     加载、注释、QC、预处理、pertpy 分析、CLI
tests/              解析、IO、无 pertpy 的冒烟测试
```
