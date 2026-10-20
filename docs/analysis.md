# 分析规则

写出哪些文件、图何时跳过，见 [outputs.md](outputs.md)。

## 目录

1. [流水线步骤](#流水线步骤)
2. [命令行与配置](#命令行与配置)
3. [物种](#物种)
4. [CRISPR 模式](#crispr-模式)
5. [对照](#对照)
6. [QC](#qc)
7. [分配](#分配)
8. [Barcode](#barcode)
9. [统计](#统计)
10. [Mixscape 与 Mixscale](#mixscape-与-mixscale)
11. [相对对照的标准化](#相对对照的标准化)
12. [预处理与簇注释](#预处理与簇注释)
13. [扰动距离](#扰动距离)
14. [靶基因表达](#靶基因表达)

---

## 流水线步骤

`--steps` 默认 `all`，顺序固定：

```text
io → assign → preprocess → mixscape → mixscale → perturb → report
```

| 步骤 | 作用 | checkpoint |
| --- | --- | --- |
| `io` | 读入矩阵 / 参考 / 分配，barcode 校验 | `checkpoints/io.h5ad` |
| `assign` | QC、资格过滤、对照标签 | `checkpoints/assign.h5ad` |
| `preprocess` | 归一化、HVG、PCA（默认 30）、邻居 / Leiden / UMAP；可选 Harmony、细胞周期 | `checkpoints/preprocess.h5ad` |
| `mixscape` | DE 门槛、`X_pert`、KO/NP 分类、signature UMAP（`n_dims` 默认 30） | `checkpoints/mixscape.h5ad` |
| `mixscale` | 独立 DE 门槛与连续扰动分数（依赖 `X_pert`） | `checkpoints/mixscale.h5ad` |
| `perturb` | 描述性 CPM log2FC、伪批量、E-test（不再内嵌 Mixscape/Mixscale） | `checkpoints/perturb.h5ad` |
| `report` | 汇总表与 PNG | `adata_processed.h5ad` 等 |

可子集重跑，例如 `--steps mixscape`，或 `--steps mixscale --resume`（从上一 checkpoint 载入）。`preprocess.n_pcs` 与 `mixscape.n_dims` 必须同为 30（或同改、且一致），否则配置校验失败。

---

## 命令行与配置

`config/config.yaml` 只放软件默认参数。样本名、物种、目录、guide 分隔符、非靶向对照和 CRISPR 模式由命令行传入。不传 `--config` 时读默认文件；命令行覆盖同名字段。

出现代码不读取的配置键会失败（`未识别的配置键` + 路径）。`labels.flag_guides_above` 一旦写出 → 已停用拒绝；`labels.guides_per_cell` → 须改为 `labels.targets_per_cell`。

| 参数 | 说明 |
| --- | --- |
| `--organism` | 必填 `human` / `mouse` → `experiment.organism` |
| `--crispr-mode` | `KO` / `i` / `a` / `other`，默认 `KO` |
| `--steps` | 见上；逗号分隔子集 |
| `--resume` | 从上一 checkpoint 继续 |
| `--validate-only` | 只检查参考表和 barcode |
| `--batch-key` / `--cell-type-key` | 空则不做 Harmony；`cell_type` 保持 `ALL` |
| `--min-umi` / `--max-umi` | 成对给出才做 UMI 绝对阈值；都不填则跳过 |
| `--min-genes` / `--max-genes` | 覆盖默认 200 / 6000 |
| `--max-mito-frac` | 覆盖默认 `0.2` |
| `--max-ribo-frac` | 显式给出才做 ribo 过滤 |

输入文件由 `--sample-id` 与 `--input-dir` 拼出：

```text
{sample_id}.scRNA.filtered.matrix.mtx.gz
{sample_id}.scRNA.filtered.barcodes.tsv.gz
{sample_id}.scRNA.filtered.features.tsv.gz
{sample_id}.scRNA.feature_barcode_reference.csv
{sample_id}.scRNA.positive_cell_guide_assignments.csv
```

图：只写 PNG，dpi=`plotting.dpi`（默认 300）。`plotting.formats` 键可存在（默认 `[png]`），实现不读该键切换格式。

---

## 物种

`--organism` 必填，写入 `experiment.organism`；只能是 `human` 或 `mouse`（大小写不敏感，加载后小写）。

| 物种 | 基因号 | 线粒体前缀 |
| --- | --- | --- |
| human | `ENSG` | `MT-` |
| mouse | `ENSMUSG` | `mt-` |

读入后核对基因号 / 符号与声明物种一致，不根据基因改写物种。矛盾时停止并写出数量与例子。叶绿体等前缀记在 `other_organelle_prefixes`，不计入 percent.mt。

结果写入 `resolved_config.yaml` 与 `run_info.json`：`resolved_organism`、`resolved_mito_prefix`、`species_check`、`organism_evidence`。验证失败仍写 `resolved_config.yaml`，然后停止。

---

## CRISPR 模式

对外写 `KO` / `i` / `a` / `other`（`--crispr-mode` → `experiment.crispr_mode`）。解析后内部为 `CRISPR_KO` / `CRISPRi` / `CRISPRa` / `other`；完整名别名可接受，但错误信息与默认参数对外只用简写。

| CLI / 配置 | 内部枚举 | 预期方向 |
| --- | --- | --- |
| 缺省 / `KO` | `CRISPR_KO` | `not_required` |
| `i` | `CRISPRi` | `down` |
| `a` | `CRISPRa` | `up` |
| `other` | `other` | 须同时给 `--crispr-description` 与 `--expected-direction`（`down` / `up` / `bidirectional` / `unknown`） |

模式不由矩阵或分配文件改写，也不由 guide 名猜测。

---

## 对照

### 设计 NT

`--nt-control` 写参考表里的非靶向 guide 名。写了且数据里确有这些细胞时，对照就是这组；zero-guide 不并入。

### 临时对照 `zero_guide`

去掉 `--nt-control` 时：

- `control_source=zero_guide_cells`，标签 `zero_guide`
- 细胞：分配表里 0 条 guide，以及仅矩阵、记为 `no_guide` 的细胞（须通过 QC）
- 之后仍做描述性 CPM log2FC、Mixscape、Mixscale、对对照的 E-distance

`multi_guide` 与 `single_targeting` 不是对照。输入目录另有一份 `*.positive_cell_guide_assignments.csv` 时停止。

### 警告与 QC 对比

`run` 一开始写入：运行日志、`tables/temporary_control_warning.txt`、`run_info.json` 的 `temporary_control_warning`。

> **警示**：临时 zero-guide 可能在 UMI、基因数、线粒体、环境 RNA、细胞状态上与靶向细胞不同；后续效应都相对这组对照，不宜直接当作设计 NT 的生物学结论。

`report` 另画 `qc/figures/temporary_control_qc.png`（UMI / 基因数 / percent.mt），表在 `qc/tables/temporary_control_qc.tsv`。对比只作描述，不因分布差异中止。设计 NT 不写警告，并删掉已有警告文件与对比图。

`low_library_no_guide` 不进入临时对照，也不进入设计 NT（该规则已禁用，见 QC）。真正检出 NT guide 的 `single_nt` 仍是设计对照。

---

## QC

### UMI

须同时给出 `--min-umi` 与 `--max-umi`（或配置成对），或全部留空。只给其中一个会失败。默认配置不设绝对阈值。

- 成对给出 → `tables/qc_thresholds.json`：`source=command_line`，`umi_filter=absolute_cutoffs`
- 全部留空 → 不做 UMI 绝对过滤：`source=skipped`，`umi_filter=disabled`

无 MAD / 分布自适应 UMI 上下限。

### 表达基因数

默认剔除 `n_genes ≤ 200` 或 `≥ 6000`：`qc.min_genes` / `qc.max_genes`。保留开区间 `(min_genes, max_genes)`。`--min-genes` / `--max-genes` 可覆盖。原因码：`high_genes` / `low_genes`。靶基因 / guide 细胞数门槛默认 `min_cells_per_target=20`、`min_cells_per_guide=10`。

### 线粒体

默认剔除 mito 比例 ≥ 20%：`qc.max_mito_frac=0.2`。`--max-mito-frac` 可覆盖。原因码 `high_mito`。

> **警示**：20% 是软件固定默认，不是按本样本自适应。肿瘤、高代谢或线粒体相关扰动应先看 `qc_violin` / `qc_scatter`，再决定是否调低/调高；过度过滤会丢掉真实生物学变异。

### 核糖体

按基因符号前缀：人类 `RPL`/`RPS`，小鼠 `Rpl`/`Rps`（大小写不敏感）。核糖体蛋白，不是线粒体；`MRPL`/`Mrpl` 不计。比例 0–1，图上 `percent.ribo` 为 0–100。

须显式 `--max-ribo-frac` / `qc.max_ribo_frac` 才过滤（`ribo_frac_source=command_line`）；未给出则 `ribo_frac_source=disabled`。超过 → `high_ribo`，主过滤阶段仍为 `qc_failed`。无 MAD 上限。

### gRNA 捕获

效率 = CRISPR guide UMI 合计 / Gene Expression `total_counts`（guide UMI 为 0 则效率为 0）。Spearman（QC 前全部细胞）→ `qc/tables/grna_capture_library_correlation.tsv`（此 p 不进效应表）。不再画捕获效率散点图。

### `low_library_no_guide`

已禁用（原按有 guide 细胞 MAD 文库下限排除浅库无 guide 的规则已移除）。过滤账本仍保留该阶段名，新增排除数恒为 0。

---

## 分配

`labels.targets_per_cell` 默认 1：下游只保留恰好 1 个唯一靶基因的细胞。按靶基因数判定，不按 guide 条数；同一靶的多条 guide（双 guide 文库）仍是 `single_targeting`。

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

`qc/figures/guides_per_cell.png` 按 **guide 条数** 分布（描述捕获），与资格用的 `n_targets` 不同。纵轴 `Cells`；分母含 `no_guide` 与 `multi_guide`；柱上标计数与一位小数百分比。

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

> **警示**：单样本描述性模式。效应表无组间推断 p/FDR；bootstrap 区间是样本内不确定性，不能当作生物学重复置信区间或显著性检验。

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

### Bootstrap

`bootstrap_n`（默认 200）、`bootstrap_ci`（默认 0.95）、种子 `runtime.seed`。

| 侧 | 单位 |
| --- | --- |
| 靶基因 ≥2 条 guide | 按 guide 有放回 |
| 靶基因仅 1 条 guide | 按细胞有放回 |
| 设计 NT ≥2 条 guide | 按 guide |
| 设计 NT 1 条 / `zero_guide` | 按细胞 |

每次靶基因与对照同时重抽。列：`log2fc_bootstrap_ci_*`、`bootstrap_unit`、`bootstrap_n`、`bootstrap_ci_level`、`uncertainty_meaning`（固定 `descriptive_estimate_uncertainty_not_significance_test`）、`bootstrap_status`。

### 单样本 DE 占位

edgeR / PyDESeq2：`single_sample_de_not_estimable`，四张表只有表头。无细胞类型注释时 `cell_type=ALL`，`cell_type_scope=unstratified`。

---

## Mixscape 与 Mixscale

两个独立步骤（`perturbseq.mixscape` / `perturbseq.mixscale`），checkpoint `mixscape.h5ad` → `mixscale.h5ad`。`perturb` 只做描述性效应、伪批量距离与 E-test。

Mixscape 是敏感性分支。标签：`nt` / `perturbed` / `non_perturbed` / `ambiguous` / `not_evaluated`。未过 DE 门槛 → `skipped_no_de_genes`。`mixscape_fraction` 分母为该靶全部靶向细胞。

| | Mixscape | Mixscale |
| --- | --- | --- |
| 对齐 | Seurat Mixscape vignette | Mixscale 教程 |
| `logfc_threshold` | **0.25**（`mixscape.logfc_threshold`） | **0.2**（`mixscale.logfc_threshold`） |
| `min_de_genes` / `max_de_genes` | 5 / 100 | 5 / 100 |
| `pval_cutoff` | 0.05 | 0.05 |
| 其他 | `n_dims=30`（= `preprocess.n_pcs`），`n_control_neighbors=20`，`iter_num=10` | 不读 Mixscape 的 0.25 / `n_dims` |

两套门槛独立。`perturb.mixscape_logfc` 仅在 `mixscape.logfc_threshold` 缺省时回退，须与 Mixscape 一致。

仅通过对应门槛的靶基因进入 2×2 / 山脊图 / 散点图。2×2 还需 NP、NT、KO 与扰动分数齐全。Mixscale 把软件填的假分数 `1` 置为缺失；散点图再丢掉有限分数全为 1 的靶基因。散点：`nbin=10`；横轴为该靶内分数分位数箱；纵轴为该基因 log 归一化均值；NT 不分面。

> **口径**：Mixscape 的 KO/NP 是样本内分类（编辑成功启发式），Mixscale 分数是连续扰动强度。二者都不是生物学重复水平上的推断。

---

## 相对对照的标准化

评估单基因扰动时，勿把「全细胞 HVG z-score」与「相对对照的扰动签名」混为一谈。本流水线按路径分开；**不再另写一层「相对 NT 均值/标准差的全局基因 Z-score」**，以免和 Mixscape 的 `X_pert` 冲突。

| 步骤 | 输入 | 相对对照？ | 口径 |
| --- | --- | --- | --- |
| 预处理 `X` | counts → library-size（默认 1e4）→ `log1p` | 否 | 全基因组；`layers['counts']` 保留原始计数 |
| 表达 PCA / 邻居 / Leiden / UMAP | HVG 上 `sc.pp.scale`（全分析细胞）→ PCA | 否 | 描述性表达空间；**不用于**扰动分类或描述性 log2FC |
| Mixscape DE 门槛 | `X`（log1p）上 Wilcoxon vs 对照 | 是（对比组） | 对照 = 设计 NT 或临时 `zero_guide` |
| Mixscape `X_pert` | pertpy `perturbation_signature`：减 k 个最近对照邻域均值（线性空间取均值再 `log1p`） | 是 | 对齐 Seurat `CalcPerturbSig`；在未 scale 的 log1p 上算 |
| Mixscape 分类 GMM | `X_pert` 上 DE 基因，`scale=True` | 是（签名内再基因 scale） | 对齐 `RunMixscape` 的 `scale.data` |
| Mixscape signature UMAP | `X_pert` → PCA（`zero_center=True`，不基因 scale）→ UMAP | 是 | 坐标在 `X_umap_mixscape` |
| Mixscale | `X_pert`，`scale=True`；对照细胞分数为 0 | 是 | 连续扰动强度 |
| 描述性效应 | 原始 counts 加总 → CPM log2FC vs 对照 | 是（比值） | **不是**基因级 Z-score |
| E-distance | 预处理 `X_pca` | 对比组是对照，嵌入不是签名 | 不用 `X_pert`；细胞置换 → 样本内距离检验 |
| 热图着色 | 展示时行 z-score | 仅作图 | 不改 `X` / `X_pert` |

无设计 NT 时，上表「对照」一律为临时 `zero_guide`。不要把表达 UMAP/Leiden 当作「已相对对照标准化后的扰动结构」。

---

## 预处理与簇注释

PCA、邻居、Leiden、UMAP 建立在 top `preprocess.n_hvg`（默认 **2000**）HVG 上。

1. 从 `layers['counts']` 对**全部基因**做文库大小归一化（默认 10000）→ `log1p`；`layers['counts']` 仍是原始计数。
2. HVG：Seurat flavor，top 2000；去掉本物种线粒体前缀后用非线粒体基因补满。
3. 只在这 2000 HVG 上 z-score（相对**全部分析细胞**），再 PCA（`n_pcs` 默认 30，与 `mixscape.n_dims` 一致）；`X` 保持全基因 log 归一化稀疏矩阵。
4. 邻居、UMAP、Leiden（igraph，`leiden_resolution` 默认 **0.5**）用这份 HVG-PCA（有 Harmony 则用 `X_pca_harmony`）。
5. `--batch-key` 空 → 不跑 Harmony；有键且水平 < 2 → `not_run_fewer_than_two_batch_levels`。

参数记在 `adata.uns['preprocess']`。

### 细胞周期

在已有 `uns['log1p']` 后，用 scib `get_cell_cycle_genes`（Tirosh / Regev）对 human / mouse 打分 → `phase` / `S_score` / `G2M_score`。不改 `X`，不做回归。

| `cell_cycle_status` | 含义 |
| --- | --- |
| `scored_after_log1p` | 成功 |
| `skipped_missing_cell_cycle_genes` | 矩阵缺基因 |
| `skipped_scib_unavailable` | 未装 scib |

### 簇注释

对 Leiden 簇均值表达跑 R SingleR + celldex（小鼠 `MouseRNAseqData`，人类 `HumanPrimaryCellAtlasData`，`label.main`）。标签 → `obs['singler_label']`；`cell_type` 不改。`singler.enabled: false` → 跳过整张 `umap_cluster_annotation`（原因 `singler_disabled`）。

---

## 扰动距离

| 产物 | 条件 | 口径 |
| --- | --- | --- |
| `E-distance/` E-test | 有对照与 `X_pca` | 表达 PCA 上相对对照；**细胞置换**，样本内距离检验，不是生物学重复 FDR |
| guide 伪批量距离热图 / MDS | 配置里显式写出 `perturb.min_guide_replicates`；默认 yaml **不设**此键 → `min_guide_replicates_not_set` 跳过 | guide 不是生物学重复；标题会标明非 E-test |

---

## 靶基因表达

`target_expression_effect.tsv`：靶基因自身 RNA 的描述性 CPM log2FC。

| `target_expression_status` | 原因 |
| --- | --- |
| `not_evaluable` | `gene_not_in_rna` / `nt_expression_below_threshold` / `insufficient_cells` |

路径见 [outputs.md](outputs.md)。
