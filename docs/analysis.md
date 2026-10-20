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
10. [相对对照的标准化](#相对对照的标准化)
11. [预处理与簇注释](#预处理与簇注释)
12. [靶基因表达](#靶基因表达)

---

## 命令行与配置

`config/config.yaml` 只放软件默认参数。样本名、物种（`--organism` → `experiment.organism`）、目录、guide 分隔符、非靶向对照和 CRISPR 模式（`--crispr-mode`：`KO` / `i` / `a` / `other`）由命令行传入。不传 `--config` 时读这份默认文件；命令行覆盖同名字段。

出现代码不读取的配置键会失败，错误为 `未识别的配置键` + 键路径。`labels.flag_guides_above` 一旦写出，会以「已停用」拒绝启动。

| 参数 | 说明 |
| --- | --- |
| `--organism` | 必填 `human` / `mouse`，写入 `experiment.organism` |
| `--crispr-mode` | `KO` / `i` / `a` / `other`，默认 `KO`；内部映射见 [CRISPR 模式](#crispr-模式) |
| `--steps` | 默认 `all`：`io` → `assign` → `preprocess` → `mixscape` → `mixscale` → `perturb` → `report` |
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

`--organism` 必填，写入配置键 `experiment.organism`；只能是 `human` 或 `mouse`（大小写不敏感，加载后小写）。

| 物种 | 基因号 | 线粒体前缀 |
| --- | --- | --- |
| human | `ENSG` | `MT-` |
| mouse | `ENSMUSG` | `mt-` |

读入后做基因号 / 符号一致性核对，不根据基因改写物种。矛盾或不一致时停止，并写出数量与例子。叶绿体等前缀记在 `other_organelle_prefixes`，不计入 percent.mt。

结果写入 `resolved_config.yaml` 与 `run_info.json`：`resolved_organism`、`resolved_mito_prefix`、`species_check`、`organism_evidence`。验证失败仍写 `resolved_config.yaml`，然后停止。

---

## CRISPR 模式

对外与命令行统一写 `KO` / `i` / `a` / `other`（`--crispr-mode`，配置键 `experiment.crispr_mode`）。解析后内部枚举为 `CRISPR_KO` / `CRISPRi` / `CRISPRa` / `other`（见下表）；完整名别名仍可接受，但错误信息与默认参数对外只写简写。

| CLI / 配置输入 | 内部枚举 | 预期方向 |
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

### UMI 绝对阈值

命令行须同时给出 `--min-umi`、`--max-umi`，或全部留空。只给其中一个会失败。默认配置不设 UMI 绝对阈值。

- 成对给出时写入 `tables/qc_thresholds.json`，`source=command_line`，`umi_filter=absolute_cutoffs`
- 全部留空时不做 UMI 绝对阈值过滤：`source=skipped`，`umi_filter=disabled`

不再用 median ± MAD 推导 UMI 上下限。

### 表达基因数

默认剔除表达基因数 ≤ 200 或 ≥ 6000 的细胞：`qc.min_genes=200`，`qc.max_genes=6000`。保留区间是开区间 `(min_genes, max_genes)`。可用 `--min-genes` / `--max-genes` 覆盖。不并入 UMI 成套绝对阈值。超过上限的原因是 `high_genes`，低于等于下限的原因是 `low_genes`。靶基因 / guide 细胞数门槛默认 `min_cells_per_target=20`、`min_cells_per_guide=10`。

### 线粒体比例

默认剔除线粒体 reads 比例 ≥ 20% 的细胞：`qc.max_mito_frac=0.2`。可用 `--max-mito-frac` 覆盖。失败原因是 `high_mito`。

> **方法学警示**：20% 是软件默认固定阈值，不是按本样本分布自适应的。肿瘤、高代谢或线粒体相关扰动研究应先看 `qc_violin` / `qc_scatter` 的线粒体分布，再决定是否用 `--max-mito-frac` 调整；过度过滤会丢掉真实生物学变异。

### 核糖体比例

按基因符号前缀：人类 `RPL`/`RPS`，小鼠 `Rpl`/`Rps`（大小写不敏感）。这是核糖体蛋白，不是线粒体；`MRPL`/`Mrpl` 不计；不用 rRNA 基因集。比例 0–1，图上 `percent.ribo` 为 0–100。

须显式给出 `--max-ribo-frac` / `qc.max_ribo_frac`（0–1）才做 ribo 过滤（`ribo_frac_source=command_line`）。未给出时不做 ribo 过滤（`ribo_frac_source=disabled`）。超过则 QC 排除，原因 `high_ribo`，主过滤阶段仍是 `qc_failed`。不再用 MAD 推导上限。

### gRNA 捕获效率与文库大小

- 效率 = CRISPR guide UMI 合计 / Gene Expression `total_counts`（guide UMI 为 0 则效率为 0）
- Spearman（QC 前全部细胞）→ `qc/tables/grna_capture_library_correlation.tsv`（此 p 不进效应表）

原 `low_library_no_guide`（按有 guide 细胞的 MAD 文库下限排除浅库无 guide 细胞）已移除，不再做该排除；过滤账本仍保留该阶段名，新增排除数为 0。

---

## 分配

`labels.targets_per_cell` 默认 1：下游只保留恰好 1 个唯一靶基因的细胞。判定按靶基因数，不按 guide 条数；同一靶基因的多条 guide（双 guide 文库）仍是 `single_targeting`。

| 状态 | 含义 |
| --- | --- |
| `no_guide` | 0 条 guide |
| `single_targeting` | 恰好 1 个靶向基因（可有多条 guide） |
| `single_nt` | 恰好 1 个设计 NT 标签（可有多条 NT guide） |
| `multi_guide` | 唯一靶基因数 > 1（不同靶，或靶向与 NT 混打） |

`labels.multi_guide_policy` 只能是 `drop`。`--guide-separator` 必须非空；guide 名内部含分隔符时按参考库最长匹配切分。

过滤顺序：

```text
low_library_no_guide → no_guide → multi_guide → low_guide_umi → qc_failed → low_cells_target → retained
```

`low_library_no_guide` 已禁用（不再用 MAD），该阶段新增排除数恒为 0。

`guides_per_cell` 图：`qc/figures/guides_per_cell.png` 仍按 guide 条数分布作图（描述捕获），与资格判定用的 `n_targets` 不同。纵轴 `Cells`；分母含 `no_guide` 与 `multi_guide`；柱上标计数与一位小数百分比。

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
| `DGE/descriptive/tables/effect_main.tsv` | 全部可分析细胞 |
| `DGE/descriptive/tables/effect_sensitivity.tsv` | Mixscape `perturbed` + 同一批对照 |

`pvalue` / `fdr` / `ci_lower` / `ci_upper` 为 NA，`inference_status=not_run_single_sample_descriptive`。细胞水平 Wilcoxon、Mixscale、E-test 的 p 分别留在 `DGE/Mixscape/plot_data/`（2x2 火山）、`obs['mixscale_score']`、`E-distance/tables/edistance_etest.tsv`。

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

### 单样本 DE

edgeR / PyDESeq2：`single_sample_de_not_estimable`，四张表只有表头。

无细胞类型注释时 `cell_type=ALL`，`cell_type_scope=unstratified`。

---

## Mixscape 与 Mixscale

流水线里是两个独立步骤（实现分别为 `perturbseq.mixscape` / `perturbseq.mixscale`），checkpoint 为 `mixscape.h5ad` → `mixscale.h5ad`。可单独重跑，例如 `--steps mixscape` 或 `--steps mixscale --resume`（Mixscale 需要上一拍的 `X_pert`）。之后的 `perturb` 只做描述性效应、伪批量距离与 E-test，不再内嵌这两步。

Mixscape 是敏感性分支。标签：`nt` / `perturbed` / `non_perturbed` / `ambiguous` / `not_evaluated`。未过 DE 门槛保持 `skipped_no_de_genes`。`mixscape_fraction` 分母为该靶全部靶向细胞。

两套门槛分开读：

| | Mixscape | Mixscale |
| --- | --- | --- |
| 对齐 | Seurat Mixscape vignette | Mixscale 教程 |
| `logfc_threshold` | **0.25**（`mixscape.logfc_threshold`；勿与 Mixscale 的 0.2 混用） | **0.2**（`mixscale.logfc_threshold`） |
| `min_de_genes` | 5 | 5 |
| `max_de_genes` | 100 | 100 |
| `pval_cutoff` | 0.05 | 0.05 |
| 其他 | `n_dims=30`（与 `preprocess.n_pcs` 一致），`n_control_neighbors=20`，`iter_num=10` | 不读 Mixscape 的 0.25 / `n_dims` |

两套门槛独立：改 Mixscape 不会改 Mixscale，反之亦然。`perturb.mixscape_logfc` 仅在 `mixscape.logfc_threshold` 缺省时回退，须与 Mixscape 一致。

仅通过对应门槛的靶基因进入 2×2 / 山脊图 / 散点图。2×2 还需 NP、NT、KO 与扰动分数齐全。Mixscale 把软件填的假分数 `1` 置为缺失；散点图再丢掉有限分数全为 1 的靶基因。

> **口径**：Mixscape 的 KO/NP 是样本内分类（编辑成功与否的启发式），Mixscale 分数是连续扰动强度。二者都不是生物学重复水平上的推断；描述性 CPM log2FC 与 E-distance 见下文。

散点图：`nbin=10`；横轴为该靶基因内分数分位数箱；纵轴为该基因 log 归一化均值；NT 不分面。

---

## 相对对照的标准化

评估单基因扰动时，不能把「全细胞 HVG z-score」和「相对对照的扰动签名」混为一谈。本流水线按路径分开处理；**不再额外写一层「相对 NT 均值/标准差的全局基因 Z-score」**，以免和 Mixscape 的 `X_pert` 冲突。

| 步骤 | 输入矩阵 | 相对对照？ | 口径 |
| --- | --- | --- | --- |
| 预处理 `X` | 原始 counts → library-size（默认 1e4）→ `log1p` | 否 | 全基因组；`layers['counts']` 保留原始计数 |
| 表达 PCA / 邻居 / Leiden / UMAP | HVG 切片上 `sc.pp.scale`（全分析细胞均值/标准差）→ PCA | 否 | 描述性表达空间；可被细胞周期、文库深度主导；**不用于**扰动分类或描述性 log2FC |
| Mixscape DE 门槛 | `X`（log1p）上 Wilcoxon vs 对照 | 是（对比组） | 对照 = 设计 NT 或临时 `zero_guide` |
| Mixscape `X_pert` | pertpy `perturbation_signature`：每个细胞减去 k 个最近对照邻域的平均表达（线性空间取均值再 `log1p`） | 是 | 对齐 Seurat `CalcPerturbSig`；签名在**未 scale 的 log1p** 上算 |
| Mixscape 分类 GMM | `X_pert` 上 DE 基因，`scale=True` | 是（签名空间内再基因 scale） | 对齐 Seurat `RunMixscape` 使用 `scale.data`；不是相对池化 NT 的第二套全局 Z |
| Mixscape signature UMAP | `X_pert` → PCA（`zero_center=True`，不基因 scale）→ UMAP | 是 | 对齐 vignette `ScaleData(do.scale=F, do.center=T)` + `RunPCA`/`RunUMAP`；坐标在 `X_umap_mixscape` |
| Mixscale | `X_pert`，`scale=True`；投影相对对照分布 z-score，对照细胞分数为 0 | 是 | 连续扰动强度，不是编辑成功与否 |
| 描述性效应 | 原始 counts 加总 → CPM log2FC vs 对照 | 是（比值） | `log2((CPM_t+α)/(CPM_c+α))`；**不是**基因级 Z-score |
| E-distance | 预处理 `X_pca`（表达空间） | 对比组是对照，嵌入不是签名 | 对齐 pertpy Distance Tests；**不用** `X_pert`；置换单位是细胞 → 样本内距离检验，不是生物学重复 FDR |
| 热图着色 | 展示时对所选细胞做行 z-score | 仅作图 | 不改 `X` / `X_pert` |

无设计 NT 时，上表「对照」一律为临时 `zero_guide`；图与 `temporary_control_warning` 会标明。不要把表达 UMAP/Leiden 当作「已相对对照标准化后的扰动结构」。

---

## 预处理与簇注释

PCA、邻居图、Leiden 聚类和 UMAP 都建立在 top `preprocess.n_hvg`（默认 **2000**）高变基因上。

1. 从 `layers['counts']` 对**全部基因**做文库大小归一化（`normalize_target_sum`，默认 10000）→ `log1p`。完整矩阵留给后续效应表；`layers['counts']` 仍是原始计数。
2. 高变基因：Seurat flavor，取 top 2000；去掉本物种线粒体前缀后用后面的非线粒体基因补满。
3. 只在这 2000 个 HVG 上做 z-score（相对**全部分析细胞**，不是相对 NT/`zero_guide`），再算 PCA（`n_pcs` 默认 30，与 `mixscape.n_dims` 一致）；`X` 保持全基因 log 归一化稀疏矩阵。
4. 邻居图、UMAP、Leiden（igraph，`leiden_resolution` 默认 **0.5**）都用这份 HVG-PCA（若跑了 Harmony 则用 `X_pca_harmony`）。这是描述性表达嵌入，不是扰动签名空间。
5. `--batch-key` 空 → 不跑 Harmony；有键且水平 < 2 → `not_run_fewer_than_two_batch_levels`。

参数记在 `adata.uns['preprocess']`。

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
