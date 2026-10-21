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
10. [差异表达（WilcoxonTest）](#差异表达wilcoxontest)
11. [Mixscape 与 Mixscale](#mixscape-与-mixscale)
12. [相对对照的标准化](#相对对照的标准化)
13. [预处理与簇注释](#预处理与簇注释)
14. [扰动距离](#扰动距离)
15. [靶基因表达](#靶基因表达)

---

## 流水线步骤

`--steps` 默认 `all`。依赖关系：

```text
io → assign → preprocess → { dge ‖ mixscape ‖ mixscale } → perturb
```

`dge` / `mixscape` / `mixscale` 只依赖 `preprocess`，彼此无边；同一次 `run` 里只要这三者中有 **≥2 个** 待跑，就会用进程池**并发**执行（各进程读临时 preprocess 副本，写各自 PNG 目录；进程间用临时目录交换，**不**写出 `checkpoints/`）。仅指定其中一步时仍串行。`perturb` 在三者之后：以 preprocess 为底，合并本轮 Mixscape 标签与 Wilcoxon `uns`（sensitivity 需要 Mixscape；Mixscale 不并入 perturb）。

不支持 `--resume` / 不落盘 h5ad：`--steps` 须从 `io` 起连续跑（通常 `--steps all`）。

### 数据分析流程图

```mermaid
flowchart TD
  IN["输入<br/>matrix / barcodes / features<br/>feature_barcode_reference<br/>positive_cell_guide_assignments"] --> IO

  IO["io<br/>读入 + barcode 校验"] --> ASG

  ASG["assign<br/>QC · 资格过滤 · 对照解析"] --> CTRL{"对照"}
  CTRL -->|"--nt-control 且数据中有可用 NT 细胞"| NT["设计 NT"]
  CTRL -->|"未指定，或指定了但无可用 NT 细胞"| ZG["临时 zero_guide"]
  NT --> PRE
  ZG --> PRE

  PRE["preprocess<br/>归一化 · HVG · PCA/neighbors<br/>Leiden/UMAP · 可选 Harmony/细胞周期"] --> PAR

  subgraph PAR["并行（≥2 步时 ProcessPool；受 runtime.n_threads 上限）"]
    DGE["dge<br/>pertpy WilcoxonTest"]
    MX["mixscape<br/>Seurat Mixscape vignette"]
    MS["mixscale<br/>satijalab Mixscale<br/>自算 X_pert_mixscale"]
  end

  DGE --> OUT_W["内存 wilcoxon uns"]
  MX --> OUT_MX["Mixscape/<br/>layer X_pert"]
  MS --> OUT_MS["mixscale/"]

  PRE --> PT
  MX -.->|"合并 mixscape 标签"| PT
  DGE -.->|"合并 wilcoxon uns"| PT
  PT["perturb<br/>描述性 CPM · E-test"] --> OUT_DGE["DGE/"]
  PT --> OUT_ED["E-distance/"]
  PT --> OUT["**/*.png only"]
```

| 步骤 | 作用 | 输出目录 |
| --- | --- | --- |
| `io` | 读入矩阵 / 参考 / 分配（账本在内存） | — |
| `assign` | QC、资格过滤、对照标签；写出 `qc/` 图 | `qc/` |
| `preprocess` | 归一化、HVG、PCA（30）、UMAP 等；写出 HVG/PCA/簇注释/细胞周期图 | `qc/`、`UMAP/` |
| `dge` | pertpy `WilcoxonTest`（与另两者独立；无独立 PNG） | — |
| `mixscape` | Seurat Mixscape：分类 + `X_pert`；写出 Mixscape / 相关 UMAP | **`Mixscape/`**、`UMAP/` |
| `mixscale` | Mixscale：自算签名 + 分数；写出 ridge / scatter / DE 计数 / DE 热图 | **`mixscale/`** |
| `perturb` | 描述性效应、E-test；写出描述性图与 E-test 图 | `DGE/`、`E-distance/` |

图在**对应分析步骤结束时**写出。跑完做 lean cleanup：磁盘只留 `**/*.png`（细节见 [outputs.md](outputs.md)）。

### 并行触发与限制

| 情形 | 行为 |
| --- | --- |
| `--steps all`（默认） | preprocess 之后三步并发，再跑 perturb |
| `--steps` 含 ≥2 个 `dge`/`mixscape`/`mixscale` 且从 `io` 连续跑到该处 | 一并并发 |
| 只指定三者之一（仍须从 `io` 起） | 单进程串行 |
| 共享输入 | 各 worker 读**临时** preprocess.h5ad（非输出目录） |
| 输出隔离 | `mixscape`→`Mixscape/`；`mixscale`→`mixscale/`；日志带 `[step]` 前缀，可能交错 |
| 线程 | 默认 `runtime.n_threads = max(1, CPU//3)`；限制 OpenMP/BLAS/`scanpy.n_jobs` 与 ProcessPool `max_workers` |
| 资源 | 并发时峰值内存约 ×3；matplotlib 用 `Agg` |

`preprocess.n_pcs` 与 `mixscape.n_dims` 须同为 30。

---

## 命令行与配置

`config/config.yaml` 只放软件默认参数。样本名、物种、目录、guide 分隔符、非靶向对照和 CRISPR 模式由命令行传入。不传 `--config` 时读默认文件；命令行覆盖同名字段。

出现代码不读取的配置键会失败（`未识别的配置键` + 路径）。`labels.flag_guides_above` 写出即拒绝；`labels.guides_per_cell` 须改为 `labels.targets_per_cell`。

| 参数 | 说明 |
| --- | --- |
| `--organism` | 必填 `human` / `mouse` → `experiment.organism` |
| `--crispr-mode` | `KO` / `i` / `a` / `other`，默认 `KO` |
| `--steps` | 见上；须从 `io` 起；含 ≥2 个 `dge`/`mixscape`/`mixscale` 时并发 |
| `--validate-only` | 只检查参考表和 barcode |
| `--batch-key` / `--cell-type-key` | 空则不做 Harmony；`cell_type` 保持 `ALL` |
| `--min-genes` / `--max-genes` | 覆盖默认 200 / 6000 |
| `--max-mito-frac` | 覆盖默认 `0.2` |

输入文件由 `--sample-id` 与 `--input-dir` 拼出：

```text
{sample_id}.scRNA.filtered.matrix.mtx.gz
{sample_id}.scRNA.filtered.barcodes.tsv.gz
{sample_id}.scRNA.filtered.features.tsv.gz
{sample_id}.scRNA.feature_barcode_reference.csv
{sample_id}.scRNA.positive_cell_guide_assignments.csv
```

图：只写 PNG，dpi=`plotting.dpi`（默认 300）。`plotting.formats` 键可存在（默认 `[png]`），实现固定写 PNG。

---

## 物种

`--organism` 必填，写入 `experiment.organism`；只能是 `human` 或 `mouse`（大小写不敏感，加载后小写）。

| 物种 | 基因号 | 线粒体前缀 |
| --- | --- | --- |
| human | `ENSG` | `MT-` |
| mouse | `ENSMUSG` | `mt-` |

读入后核对基因号 / 符号与声明物种一致，不根据基因改写物种。矛盾时停止并写出数量与例子。叶绿体等前缀记在 `other_organelle_prefixes`，不计入 percent.mt。

物种核对结果写入配置对象；跑程中可能短暂写出 `resolved_config.yaml` / `run_info.json`，成功结束后 lean cleanup 会删掉。验证失败时 `resolved_config.yaml` 仍可能留在输出目录（未完成 cleanup）。

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

单一来源：`control_spec(config)`。`assign` 用账本调用 `resolve_control_spec` 写入 `config['control']`、`adata.uns`；后续 dge / Mixscape / Mixscale / perturb 只读该解析结果（同一次运行内存传递）。

### 设计 NT

`--nt-control` 写参考表里的非靶向 guide 名。写了**且** QC 通过后账本里确有对应 NT 细胞时，对照就是这组；zero-guide 不并入。

### 临时对照 `zero_guide`

以下任一情况使用临时对照：

1. 未指定 `--nt-control` / `labels.nt_control` 为空  
2. 指定了设计 NT，但数据中没有可用（QC 通过）的 NT 细胞 → **自动回退**到 zero_guide，并在 note / 日志中写明 `fallback_from`

临时对照时：

- `control_source=zero_guide_cells`，标签 `zero_guide`，`temporary_control=true`
- 细胞：分配表里 0 条 guide，以及仅矩阵、记为 `no_guide` 的细胞（须通过 QC）
- Wilcoxon / 描述性 CPM log2FC / Mixscape / Mixscale / E-distance 均相对该对照继续跑（无对照细胞时各步显式 skip，原因如 `no_zero_guide_cells`）

`multi_guide` 与 `single_targeting` 不是对照。输入目录另有一份 `*.positive_cell_guide_assignments.csv` 时停止。

### 警告与 QC 对比

`run` 开始时按配置写出临时警告；`assign` 解析（含数据级回退后）再写一次：日志与 `run_info.json` 的 `temporary_control_warning`（跑程中也可能短暂写出 `tables/temporary_control_warning.txt`，结束 cleanup 后不保留）。

> **警示**：临时 zero-guide 可能在 UMI、基因数、线粒体、环境 RNA、细胞状态上与靶向细胞不同；后续效应都相对这组对照，不宜直接当作设计 NT 的生物学结论。

`assign` 步写出 `qc/temporary_control_qc.png`（UMI / 基因数 / percent.mt）。对比只作描述，不因分布差异中止。设计 NT 不写警告，并去掉已有警告文件与对比图。对照掩码在 pandas CoW 下会拷贝为可写数组，避免 `output array is read-only` 导致对比图被跳过。

`low_library_no_guide` 阶段在过滤账本中保留，新增排除数为 0；该状态不进入临时对照，也不进入设计 NT。真正检出 NT guide 的 `single_nt` 仍是设计对照。

---

## QC

QC **指标**（均计算并画在 `qc/`）：文库大小（`nCount` / `total_counts`）、表达基因数（`nFeature` / `n_genes`）、线粒体比例（`percent.mt` / `mito_frac`）。

QC **过滤**只用基因数与线粒体比例；**不做**文库 UMI / nCount 绝对阈值过滤（无 `min_umi` / `max_umi` / `min_counts` / `max_counts`）。

### 表达基因数

默认剔除 `n_genes ≤ 200` 或 `≥ 6000`：`qc.min_genes` / `qc.max_genes`。保留开区间 `(min_genes, max_genes)`。`--min-genes` / `--max-genes` 可覆盖。原因码：`high_genes` / `low_genes`。靶基因 / guide 细胞数门槛默认 `min_cells_per_target=20`、`min_cells_per_guide=10`。

### 线粒体

默认剔除 mito 比例 ≥ 20%：`qc.max_mito_frac=0.2`。`--max-mito-frac` 可覆盖。原因码 `high_mito`。

> **警示**：20% 是软件固定默认，不是按本样本自适应。肿瘤、高代谢或线粒体相关扰动应先看 `qc_metrics_panel`（含 percent.mt violin/直方图）/ `qc_scatter_panel`，再决定是否调低/调高；过度过滤会丢掉真实生物学变异。

### 文库大小（nCount）

`total_counts` 写入账本并出现在 `qc_metrics_panel` / 散点图；不按绝对上下限剔除细胞。guide 级 `qc.min_guide_umi`（单条 guide 的 UMI 是否过线）与文库大小过滤无关。

### GEX marker 图

`preprocess` 后写出 `marker_detection_heatmap` 与 `dotplot_markers`（白底 Library-D / Seurat 点图风格）。基因列表优先读 `plotting.markers`；未配置时按 `experiment.organism` 使用默认 housekeeping（人类：`ACTB,GAPDH,B2M,MALAT1,NEAT1,UBC,PPIA,EEF1A1`；小鼠对应首字母小写形式）。仅使用矩阵中存在的基因；全缺失则 skip。行标签优先 `singler_label`（否则 `cell_type`）拼成 `leiden: annotation`；无注释时退回裸 Leiden id。点图：点大小 = Percent Expressed，颜色 = Average Expression（右侧双图例）。

### gRNA 捕获

效率 = CRISPR guide UMI 合计 / Gene Expression `total_counts`（guide UMI 为 0 则效率为 0）。Spearman（QC 前全部细胞）记入 `adata.uns['qc_thresholds'].grna_capture_library_correlation`（此 p 不进效应表）。

### `low_library_no_guide`

过滤账本保留该阶段名；当前实现不应用该排除规则，新增排除数恒为 0。

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

`qc/guides_per_cell.png` 按 **guide 条数** 分布（描述捕获），与资格用的 `n_targets` 不同。纵轴 `Cells`；分母含 `no_guide` 与 `multi_guide`；柱上标计数与一位小数百分比。

---

## Barcode

`input.barcode_policy` 必须是 `strict`。

| 情况 | 行为 |
| --- | --- |
| 仅 assignment-only | `io` / `--validate-only` 失败，写 `barcode_differences.tsv`、`input_validation.json` |
| 仅 matrix-only，且 assignment-only = 0 | 记为 `no_guide` 并保留；合格者可进 `zero_guide` |
| 两侧都有差异 | 失败 |
| 集合一致 | 差异表只有表头 |

成功 `run` 结束后这些校验文件随 lean cleanup 删除；`--validate-only` 保留。

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

| 表（`perturb` 步内存 `adata.uns`；图在 `DGE/`） | 分支 |
| --- | --- |
| `effect_main` | 全部可分析细胞 |
| `effect_sensitivity` | Mixscape `perturbed` + 同一批对照 |

`pvalue` / `fdr` / `ci_lower` / `ci_upper` 为 NA，`inference_status=not_run_single_sample_descriptive`。主 DE 见下一节 WilcoxonTest；Mixscape 2×2、Mixscale 分数、E-test 分别在 Mixscape 图 / `obs['mixscale_score']` / `uns['edistance_etest']`。

### Bootstrap

`bootstrap_n`（默认 200）、`bootstrap_ci`（默认 0.95）、种子 `runtime.seed`。

| 侧 | 单位 |
| --- | --- |
| 靶基因 ≥2 条 guide | 按 guide 有放回 |
| 靶基因仅 1 条 guide | 按细胞有放回 |
| 设计 NT ≥2 条 guide | 按 guide |
| 设计 NT 1 条 / `zero_guide` | 按细胞 |

每次靶基因与对照同时重抽。列：`log2fc_bootstrap_ci_*`、`bootstrap_unit`、`bootstrap_n`、`bootstrap_ci_level`、`uncertainty_meaning`（固定 `descriptive_estimate_uncertainty_not_significance_test`）、`bootstrap_status`。

无细胞类型注释时 `cell_type=ALL`，`cell_type_scope=unstratified`。

---

## 差异表达（WilcoxonTest）

步骤 `dge`（`perturbseq.differential_expression`）对齐 pertpy [Differential gene expression](https://pertpy.scverse.org/en/stable/tutorials/notebooks/differential_gene_expression.html) 中 **Simple statistical tests → `WilcoxonTest`**：对表达矩阵做两组比较（Mann-Whitney U），**不经过** Mixscape / Mixscale。

| 项 | 口径 |
| --- | --- |
| API | `pertpy.tools.WilcoxonTest.compare_groups` |
| 矩阵 | preprocess 后的 `adata.X`（log1p） |
| 对比 | 每个 `target_gene` vs 当前对照（`control_spec`：设计 NT 或临时 `zero_guide`） |
| `log_fc` | pertpy：`log2(mean_log1p_g1) - log2(mean_log1p_g0)`，**无** pseudocount；零均值 → ±inf/NaN，写出前改为 NA（p 值不动）。**不是** counts 空间的 log2FC |
| 与 Mixscape | **独立**：不读 `mixscape_status` / `X_pert`；Mixscape 自有 DE 门槛只服务分类 |
| 输出 | 内存 `wilcoxon_de_status` / meta（全基因表不落盘）；无独立 PNG |

> **警示**：单生物学样本下细胞不是独立生物学重复；`p_value` / `adj_p_value` 是细胞级秩和检验，不能当作多样本组间 FDR。

随 `--steps all`（从 `io` 起）一并运行。

---

## Mixscape 与 Mixscale

两个**互不依赖**的步骤（`perturbseq.mixscape` / `perturbseq.mixscale`），与 `dge` 一样只依赖 `preprocess`；默认全量跑时三者并发（见上文「并行触发与限制」）。结果目录分别为 `Mixscape/` 与 `mixscale/`。Mixscale **自行**算 `layers['X_pert_mixscale']`（对齐教程 CalcPerturbSig），不读 Mixscape 的 `X_pert` 或分类标签。

| | Mixscape | Mixscale | Wilcoxon (`dge`) |
| --- | --- | --- | --- |
| 参考 | [Mixscape vignette](https://satijalab.org/seurat/articles/mixscape_vignette) | [Mixscale](https://satijalab.github.io/Mixscale/) | [pertpy WilcoxonTest](https://pertpy.scverse.org/en/stable/tutorials/notebooks/differential_gene_expression.html) |
| 输出目录 | `Mixscape/` | `mixscale/` | （无独立 PNG） |
| 签名层 | `X_pert` | `X_pert_mixscale` | — |
| `logfc_threshold`（\|logFC\|） | **0.1** | **0.1** | — |
| `min_de_genes` | 5 | 5 | — |
| `pval_cutoff` | 0.05 | 0.05 | — |
| 主要产物 | KO/NP 分类、signature UMAP、2×2、DE 计数 | `mixscale_score`、DE 计数、DE 热图 | 靶 vs 对照 DE 表 |

Mixscape 标签：`nt` / `perturbed` / `non_perturbed` / `ambiguous` / `not_evaluated`；未过 DE 门槛 → `skipped_no_de_genes`。仅通过对应门槛的靶进入 2×2 / 山脊 / 散点。Mixscape 2×2（`Mixscape/{靶}_mixscape_2x2.png`）为后验 / 扰动分数 / 细胞级火山 / Top 基因热图；热图列块为 **NP | NT | KO**（CRISPRi/a 时 KO 块显示 KD/OE）。Mixscale 假分数 `1` 置缺失。

Mixscale 过 DE 门后，用**本步** `uns['mixscale_de_gate']` / `uns['mixscale_de_genes']`（非 Mixscape）写出与 Mixscape 对称的计数图与热图：`mixscale/target_de_counts.png`、`target_in_de_gene_counts.png`、`{靶}_mixscale_heatmap.png`；「靶基因 symbol 精确出现在自身 DE 列表」的热图再复制到 `mixscale/target_in_de/`（对标 `Mixscape/target_in_de/`）。

**Mixscale 分组（现状）**：官方 [Mixscale_DoHeatmap](https://satijalab.github.io/Mixscale/) 只有 NT vs 全部靶向细胞并按 `mixscale_score` 排序，不做 NP/KO 分类。本流水线在**同跑且已合并 Mixscape 标签**时，热图 / ridge / scatter 优先按 **NT | NP | KO**（默认；CRISPRi→KD、CRISPRa→OE、other→PT）分块或分面，块内仍按 `mixscale_score` 排序；无 Mixscape 分类时回退为文档等价的 NT | perturbed。并行步结束后若两者皆有，会重画 `mixscale/` 图以带上 NT/NP/KO。

**无 ridge / scatter / DE 热图**：`passed=True` 的靶为 0 时，`mixscale_score_ridge` / `mixscale_scatter` / DE 热图记 `no_targets_passed_de_gate` 并跳过 PNG。临时 `zero_guide`、弱效应或细胞数偏少时可能全部靶 `n_de_genes=0`。可调 `mixscale.logfc_threshold` / `min_de_genes`。

> **口径**：计算上三者独立；作图在合并后可对齐 Mixscape 标签。KO/NP 与 Mixscale 分数都不是生物学重复推断；Wilcoxon p 为细胞级。

---

## 相对对照的标准化

评估单基因扰动时，勿把「全细胞 HVG z-score」与「相对对照的扰动签名」混为一谈。本流水线按路径分开：表达嵌入用全分析细胞 scale；扰动分类与 Mixscale 用相对对照的签名层。

| 步骤 | 输入 | 相对对照？ | 口径 |
| --- | --- | --- | --- |
| 预处理 `X` | counts → library-size（默认 1e4）→ `log1p` | 否 | 全基因组；`layers['counts']` 保留原始计数 |
| 表达 PCA / 邻居 / Leiden / UMAP | HVG 上 `sc.pp.scale`（全分析细胞）→ PCA | 否 | 描述性表达空间；**不用于**扰动分类或描述性 log2FC |
| Mixscape DE 门槛 | `X`（log1p）上 Wilcoxon vs 对照 | 是（对比组） | 对照 = 设计 NT 或临时 `zero_guide` |
| Mixscape `X_pert` | pertpy `perturbation_signature`：减 k 个最近对照邻域均值（线性空间取均值再 `log1p`） | 是 | 对齐 Seurat `CalcPerturbSig`；在未 scale 的 log1p 上算 |
| Mixscape 分类 GMM | `X_pert` 上 DE 基因，`scale=True` | 是（签名内再基因 scale） | 对齐 `RunMixscape` 的 `scale.data` |
| Mixscape signature UMAP | `X_pert` → PCA（`zero_center=True`，不基因 scale）→ UMAP | 是 | 坐标在 `X_umap_mixscape` |
| Mixscale | 自算 `X_pert_mixscale`，`scale=True`；对照细胞分数为 0 | 是 | 与 Mixscape 步骤独立；连续扰动强度 |
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

对 Leiden 簇均值表达跑 R SingleR + celldex（小鼠 `MouseRNAseqData`，人类 `HumanPrimaryCellAtlasData`，`label.main`）。标签 → `obs['singler_label']`；`cell_type` 不改。`singler.enabled: false` → 跳过组合图 `umap_annotation_cell_cycle_panel`（原因 `singler_disabled`）。

---

## 扰动距离

| 产物 | 条件 | 口径 |
| --- | --- | --- |
| `E-distance/` E-test | 有对照与 `X_pca` | 表达 PCA 上相对对照；**细胞置换**，样本内距离检验，不是生物学重复 FDR |

实现用 pertpy 的 energy-distance 统计量，但**不**物化 N×N pairwise 矩阵（按组流式置换）。`edistance.n_perms` 默认 1000；超大对照（常见于临时 zero-guide）按 `edistance.max_contrast_cells`（默认 4000）随机下采样。

图为 `E-distance/edistance_etest.png`。

---

## 靶基因表达

`uns['target_expression_effect']`（`perturb` 步内存）：靶基因自身 RNA 的描述性 CPM log2FC；图在 `DGE/`。

| `target_expression_status` | 原因 |
| --- | --- |
| `not_evaluable` | `gene_not_in_rna` / `nt_expression_below_threshold` / `insufficient_cells` |

路径与图见 [outputs.md](outputs.md)。
