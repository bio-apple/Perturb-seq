# 分析规则

写出哪些文件、图何时跳过，见 [outputs.md](outputs.md)。

## 目录

1. [命令行与配置](#命令行与配置)
2. [物种](#物种)
3. [CRISPR 模式](#crispr-模式)
4. [对照](#对照)
5. [QC](#qc)
6. [分配](#分配)
7. [Barcode](#barcode)
8. [统计](#统计)
9. [Mixscape 与 Mixscale](#mixscape-与-mixscale)
10. [预处理与簇注释](#预处理与簇注释)
11. [靶基因表达](#靶基因表达)

---

## 命令行与配置

`config/config.yaml` 只放软件默认参数。样本名、物种、目录、guide 分隔符、非靶向对照和 CRISPR 模式由命令行传入。不传 `--config` 时读这份默认文件；命令行覆盖同名字段。

出现代码不读取的配置键会失败，错误为 `未识别的配置键` + 键路径。`labels.flag_guides_above` 一旦写出，会以「已停用」拒绝启动。

| 参数 | 说明 |
| --- | --- |
| `--crispr-mode` | 默认 `CRISPR_KO`，也可 `i` / `a` / `other` |
| `--steps` | 默认 `all`：`io` → `assign` → `preprocess` → `perturb` → `augur` → `enrich` → `report` |
| `--resume` | 从上一 checkpoint 继续 |
| `--validate-only` | 只检查参考表和 barcode |
| `--batch-key` / `--cell-type-key` | 空则不做 Harmony；`cell_type` 保持 `ALL` |

输入文件由 `--sample-id` 与 `--input-dir` 拼出：

```text
{sample_id}.scRNA.filtered.matrix.mtx.gz
{sample_id}.scRNA.filtered.barcodes.tsv.gz
{sample_id}.scRNA.filtered.features.tsv.gz
{sample_id}.scRNA.feature_barcode_reference.csv
{sample_id}.scRNA.positive_cell_guide_assignments.csv
```

最小命令见仓库根目录 README。

---

## 物种

`--organism` 必填，只能是 `human` 或 `mouse`（大小写不敏感，加载后小写）。

| 物种 | 基因号 | 线粒体前缀 |
| --- | --- | --- |
| human | `ENSG` | `MT-` |
| mouse | `ENSMUSG` | `mt-` |

读入后做基因号 / 符号一致性核对，不根据基因改写物种。矛盾或不一致时停止，并写出数量与例子。叶绿体等前缀记在 `other_organelle_prefixes`，不计入 percent.mt。

结果写入 `resolved_config.yaml` 与 `run_info.json`：`resolved_organism`、`resolved_mito_prefix`、`species_check`、`organism_evidence`。验证失败仍写 `resolved_config.yaml`，然后停止。

---

## CRISPR 模式

| 输入 | 解析为 | 预期方向 |
| --- | --- | --- |
| 缺省 / `KO` | `CRISPR_KO` | `not_required` |
| `i` | `CRISPRi` | `down` |
| `a` | `CRISPRa` | `up` |
| `other` | `other` | 必须同时给 `--crispr-description` 与 `--expected-direction`（`down` / `up` / `bidirectional` / `unknown`） |

模式不由矩阵或分配文件改写，也不由 guide 名猜测。

---

## 对照

### 设计好的 NT

`--nt-control` 写参考表里的非靶向 guide 名。写了且数据里确有这些细胞时，对照就是这组细胞；zero-guide 不并入。

### 临时对照 `zero_guide`

去掉 `--nt-control` 时：

- 对照来源：`control_source=zero_guide_cells`，标签 `zero_guide`
- 细胞：分配表里 0 条 guide，以及仅矩阵、记为 `no_guide` 的细胞（须通过 QC，且不是 `low_library_no_guide`）
- 之后仍做描述性 CPM log2FC、Mixscape、Mixscale、对对照的 E-distance

`multi_guide` 与 `single_targeting` 不是对照。输入目录另有一份 `*.positive_cell_guide_assignments.csv` 时停止。

### 临时对照警告与 QC 对比

`run` 一开始写入：

- 运行日志
- `tables/temporary_control_warning.txt`
- `run_info.json` 的 `temporary_control_warning`

警告要点：当前是临时 zero-guide；这些细胞可能在 UMI、基因数、线粒体、环境 RNA、细胞状态上与靶向细胞不同；后续效应相对这组对照。

`report` 另画 `qc/figures/temporary_control_qc.png`（UMI / 基因数 / percent.mt），表在 `qc/tables/temporary_control_qc.tsv`。对比只作描述，不因分布差异中止。设计好的 NT 不写警告，并删掉已有警告文件与对比图。

`low_library_no_guide` 不进入临时对照，也不进入设计 NT。真正检出 NT guide 的 `single_nt` 仍是设计对照。

---

## QC

### 四个绝对阈值

命令行须同时给出 `--min-umi`、`--max-umi`、`--min-genes`、`--max-mito-frac`，或全部留空。部分给出会失败。

全部留空时，在分配步骤用全部细胞计算，写入 `tables/qc_thresholds.json`，`source=median_mad`：

- UMI / 基因数：`log1p` → median ± `qc.n_mad`（默认 3）× 1.4826 × MAD → `expm1`
- 线粒体比例：原分数上同一公式，限制在 \[0, 1\]
- MAD 为 0 时改用标准差

命令行给齐四个数时 `source=command_line`。

### 核糖体比例

按基因符号前缀：人类 `RPL`/`RPS`，小鼠 `Rpl`/`Rps`（大小写不敏感）。这是核糖体蛋白，不是线粒体；`MRPL`/`Mrpl` 不计；不用 rRNA 基因集。比例 0–1，图上 `percent.ribo` 为 0–100。

未给 `--max-ribo-frac` 时，上限为 median + `qc.n_mad` × 1.4826 × MAD（`ribo_frac_source=median_mad`）。单独给出（0–1）时用该值。超过则 QC 排除，原因 `high_ribo`，主过滤阶段仍是 `qc_failed`。

### gRNA 捕获效率与文库大小

- 效率 = CRISPR guide UMI 合计 / Gene Expression `total_counts`（guide UMI 为 0 则效率为 0）
- Spearman（QC 前全部细胞）→ `qc/tables/grna_capture_library_correlation.tsv`（此 p 不进效应表）
- 图：`qc/figures/grna_capture_vs_library_size.png`

避免文库偏浅、未捕获 guide 的细胞被当作阴性对照：

1. 在已检出 ≥1 条 guide 的细胞上，对 `log1p(total_counts)` 取 median − n_mad × 1.4826 × MAD，再 `expm1`
2. 无 guide 且 `total_counts` 低于该下限 → `qc_pass=false`，阶段 `low_library_no_guide`
3. 仍记在账本，分配状态仍是 `no_guide`，但不进对照

无细胞有 guide，或有 guide 细胞文库方差为 0 时跳过该排除（`no_cells_with_guide` / `library_size_variance_zero`）。

---

## 分配

`labels.guides_per_cell` 默认 1：下游只保留恰好 1 条 guide 的细胞。

| 状态 | 含义 |
| --- | --- |
| `no_guide` | 0 条 |
| `single_targeting` | 1 条靶向 |
| `single_nt` | 1 条设计 NT |
| `multi_guide` | 多于 1 条 |

`labels.multi_guide_policy` 只能是 `drop`。`--guide-separator` 必须非空；guide 名内部含分隔符时按参考库最长匹配切分。

过滤顺序：

```text
low_library_no_guide → no_guide → multi_guide → low_guide_umi → qc_failed → low_cells_target → retained
```

规则跳过时，`low_library_no_guide` 新增排除数为 0。

`guides_per_cell` 图：`qc/figures/guides_per_cell.png`。纵轴 `Cells`；分母含 `no_guide` 与 `multi_guide`；柱上标计数与一位小数百分比。不画 guide 数 > 5 的红线。

---

## Barcode

`input.barcode_policy` 必须是 `strict`。

| 情况 | 行为 |
| --- | --- |
| 仅 assignment-only | `io` / `--validate-only` 失败，写 `barcode_differences.tsv`、`input_validation.json` |
| 仅 matrix-only，且 assignment-only = 0 | 记为 `no_guide` 并保留；合格者可进 `zero_guide` |
| 两侧都有差异 | 失败 |
| 集合一致 | 差异表只有表头 |

---

## 统计

`statistics.mode` 只能是 `descriptive`。`guide_exploratory` → `unsupported_inference_backend`。`allow_random_split_replicates: true` 被拒绝。

### 描述性 log2FC

按靶基因加总原始计数 → CPM：

```text
log2((CPM_target + pseudocount_cpm) / (CPM_control + pseudocount_cpm))
```

两侧计数都为 0 → `log2fc_descriptive` 缺失，`effect_status=not_detected_both`。

| 表 | 分支 |
| --- | --- |
| `effect_main.tsv` | 全部可分析细胞 |
| `effect_sensitivity.tsv` | Mixscape `perturbed` + 同一批对照 |

`pvalue` / `fdr` / `ci_lower` / `ci_upper` 为 NA，`inference_status=not_run_single_sample_descriptive`。细胞水平 Wilcoxon、Mixscale、E-test 的 p 分别留在 `DGE/plot_data/`、`obs['mixscale_score']`、`E-distance/tables/edistance_etest.tsv`。

### Bootstrap 区间（描述性不确定性）

配置：`bootstrap_n`（默认 200）、`bootstrap_ci`（默认 0.95）、种子 `runtime.seed`。

| 侧 | 单位 |
| --- | --- |
| 靶基因 ≥2 条 guide | 按 guide 有放回 |
| 靶基因仅 1 条 guide | 按细胞有放回 |
| 设计 NT ≥2 条 guide | 按 guide |
| 设计 NT 1 条 / `zero_guide` | 按细胞 |

每次靶基因与对照同时重抽。新列：`log2fc_bootstrap_ci_*`、`bootstrap_unit`、`bootstrap_n`、`bootstrap_ci_level`、`uncertainty_meaning`（固定 `descriptive_estimate_uncertainty_not_significance_test`）、`bootstrap_status`。

区间表示样本内描述性不确定性，不能当作显著性检验或生物学重复置信区间。

### 单样本 DE 与 Enrichr

edgeR / PyDESeq2：`single_sample_de_not_estimable`，四张表只有表头。

Enrichr：`|log2fc_descriptive| ≥ plotting.volcano_lfc`（默认 1.0）的固定名单，需联网。adjusted P → `enrichment_fdr`。无入选基因 → `no_effect_ranked_genes`；断网 → `enrichr_failed`。

无细胞类型注释时 `cell_type=ALL`，`cell_type_scope=unstratified`。

---

## Mixscape 与 Mixscale

Mixscape 是敏感性分支。标签：`nt` / `perturbed` / `non_perturbed` / `ambiguous` / `not_evaluated`。未过 DE 门槛保持 `skipped_no_de_genes`。`mixscape_fraction` 分母为该靶全部靶向细胞。

两套门槛分开读：

| | Mixscape | Mixscale |
| --- | --- | --- |
| 对齐 | Seurat Mixscape vignette | Mixscale 教程 |
| `logfc_threshold` | 0.25 | 0.2 |
| `min_de_genes` | 5 | 5 |
| `max_de_genes` | 100 | 100 |
| `pval_cutoff` | 0.05 | 0.05 |
| 其他 | `n_dims=40`，`n_control_neighbors=20`，`iter_num=10` | 不读 Mixscape 的 0.25 / `n_dims` |

仅通过对应门槛的靶基因进入 2×2 / 山脊图 / 散点图。2×2 还需 NP、NT、KO 与扰动分数齐全。Mixscale 把软件填的假分数 `1` 置为缺失；散点图再丢掉有限分数全为 1 的靶基因。

散点图：`nbin=10`；横轴为该靶基因内分数分位数箱；纵轴为该基因 log 归一化均值；NT 不分面。

---

## 预处理与簇注释

1. 从 `layers['counts']` 归一化到 `normalize_target_sum`（10000）→ `log1p`
2. HVG：Seurat flavor，`n_hvg=2000`；去掉本物种线粒体前缀后补非线粒体基因
3. PCA 在标准化 HVG 上，`n_pcs=30`；`X` 保持 log 归一化
4. Leiden：igraph，分辨率来自配置
5. `--batch-key` 空 → 不跑 Harmony；有键且水平 < 2 → `not_run_fewer_than_two_batch_levels`

### 细胞周期

在已有 `uns['log1p']` 后，用 [scib `get_cell_cycle_genes`](https://scib.readthedocs.io/en/latest/api/scib.preprocessing.score_cell_cycle.html)（Tirosh / Regev）对 human / mouse 打分，写入 `phase` / `S_score` / `G2M_score`。不改 `X`，不做回归。

| `cell_cycle_status` | 含义 |
| --- | --- |
| `scored_after_log1p` | 成功 |
| `skipped_missing_cell_cycle_genes` | 矩阵缺基因 |
| `skipped_scib_unavailable` | 未装 scib |

### 簇注释

对 Leiden 簇均值表达跑 R SingleR + celldex（小鼠 `MouseRNAseqData`，人类 `HumanPrimaryCellAtlasData`，`label.main`）。标签 → `obs['singler_label']`；`cell_type` 不改。`singler.enabled: false` → 跳过，原因 `singler_disabled`。

---

## 靶基因表达

`target_expression_effect.tsv`：靶基因自身 RNA 的描述性 CPM log2FC。

| `target_expression_status` | 原因 |
| --- | --- |
| `not_evaluable` | `gene_not_in_rna` / `nt_expression_below_threshold` / `insufficient_cells` |

图与表的路径见 [outputs.md](outputs.md)。
