# 分析规则

图路径与输出规则见 [outputs.md](outputs.md)。

## 目录

1. [流水线步骤](#流水线步骤)
2. [命令行与配置](#命令行与配置)
3. [物种](#物种)
4. [CRISPR 模式](#crispr-模式)
5. [对照](#对照)
6. [QC](#qc)
7. [分配](#分配)
8. [Barcode](#barcode)
9. [Mixscape 与 Mixscale](#mixscape-与-mixscale)
10. [预处理与簇注释](#预处理与簇注释)

---

## 流水线步骤

```text
io → assign → preprocess → mixscape ‖ mixscale
```

Mixscape 与 Mixscale 以进程池并发，各读临时 preprocess 副本，互不依赖。不支持 `--resume`，须从 `io` 起连续跑。两者都跑完且已有 Mixscape 标签时，会重画 `mixscale/` 图按 NT|NP|KO 分面。

| 步骤 | 作用 | 输出目录 |
| --- | --- | --- |
| `io` | 读入矩阵 / 参考 / 分配 | — |
| `assign` | QC、资格过滤、对照解析；写出 `qc/` 图 | `qc/` |
| `preprocess` | 归一化、HVG、PCA、UMAP / Leiden / SingleR / ScType | `UMAP/` |
| `mixscape` | `analysis_eligible` ∪ 对照上：签名 → pertpy 分类 | `Mixscape/`、`UMAP/` |
| `mixscale` | 同范围独立签名 → pertpy 分数 | `mixscale/` |

跑完做 lean cleanup：只留 `**/*.png`（细节见 [outputs.md](outputs.md)）。

---

## 命令行与配置

`config/config.yaml` 放软件默认参数；样本、路径、物种等由命令行传入，覆盖同名字段。不识别的配置键会报错。

| 参数 | 说明 |
| --- | --- |
| `--sample-id` | 样本名 |
| `--organism` | 必填 `human` / `mouse` |
| `--input-dir` / `--output-dir` | 输入 / 输出目录 |
| `--guide-separator` | guide 分隔符（如 `|`） |
| `--crispr-mode` | `KO`（默认）/ `i` / `a` |
| `--min-genes` / `--max-genes` | 覆盖默认 200 / 6000 |
| `--max-mito-frac` | 覆盖默认 `0.15` |
| `--nt-control` | 默认 `NEG_CONTROL/Non-Targeting/NegCtrl`；未匹配 → 自动 `zero_guide`；可显式填 `zero_guide` |
| `--sctype-tissue` | 启用 ScType；未传则跳过 |
| `--sctype-db-path` | 自定义 marker DB；默认 `data/sctype/ScTypeDB_full.tsv` |
| `--sctype-low-score-fraction` | 低分簇 → `Unknown`（默认 0.25） |
| `--config` | 默认 YAML；缺省 `config/config.yaml` |

常用 YAML 键（完整见 `config/config.yaml`）：

| 键 | 作用 |
| --- | --- |
| `qc.min_genes` / `max_genes` / `max_mito_frac` | QC 过滤 |
| `qc.min_cells_per_target` | 靶向细胞数 < 该值 → `low_cells_target`（默认 20） |
| `qc.min_cells_per_guide` | 仅统计用，不排除细胞（默认 10） |
| `preprocess.n_hvg` / `n_pcs` / `leiden_resolution` | HVG / PCA / Leiden |
| `mixscape.*` / `mixscale.*` | pertpy 参数（`logfc_threshold` 0.1；互不读取） |
| `singler.enabled` | `false` 跳过 SingleR |
| `runtime.n_threads` | 默认 `max(1, CPU//3)` |

---

## 物种

`--organism` 必填 `human` 或 `mouse`。

| 物种 | 基因号 | 线粒体前缀 |
| --- | --- | --- |
| human | `ENSG` | `MT-` |
| mouse | `ENSMUSG` | `mt-` |

读入后核对基因号与声明物种一致，矛盾时停止。叶绿体等前缀不计入 percent.mt。

---

## CRISPR 模式

| CLI | 内部 | 预期方向 | 图标签 |
| --- | --- | --- | --- |
| `KO`（默认） | `CRISPR_KO` | `not_required` | KO |
| `i` | `CRISPRi` | `down` | KD |
| `a` | `CRISPRa` | `up` | OE |

---

## 对照

单一来源：`assign` 解析后写入 `config['control']`，Mixscape / Mixscale 只读该结果。

### 设计 NT

`--nt-control` 默认在参考表 `id` / `name` / `target_gene_name` / `target_gene_id` 上做 `re.search`（忽略大小写）。匹配到的 guide 且 QC 后有对应细胞时，对照为这组设计 NT。

### 临时对照 `zero_guide`

以下情况使用临时对照：(1) 显式 `--nt-control zero_guide`；(2) 关键词未匹配 → 自动；(3) 匹配到但无可用 NT 细胞 → 自动回退。

临时对照细胞 = 0 条 guide 且 QC 通过的细胞。Mixscape / Mixscale 相对该对照继续跑；无对照细胞时显式 skip。

> **警示**：临时 zero-guide 可能在 UMI、基因数、线粒体等指标上与靶向细胞不同，后续效应都相对这组对照，不宜直接当作设计 NT 的生物学结论。`assign` 步会写出 `qc/temporary_control_qc.png` 对比。

---

## QC

QC 指标：`nCount`、`nFeature`、`percent.mt`（均画在 `qc/`）。过滤只用基因数与线粒体比例，不做文库 UMI 绝对阈值过滤。

- **基因数**：剔除 `n_genes <= 200` 或 `>= 6000`。原因码 `high_genes` / `low_genes`
- **线粒体**：剔除 `mito_frac >= 0.15`。原因码 `high_mito`。15% 是固定默认，非自适应——过度过滤会丢掉真实生物学变异
- **靶细胞数**：靶向细胞数 < `min_cells_per_target`（默认 20）→ `low_cells_target`，不进 `analysis_eligible`

---

## 分配

`labels.targets_per_cell` 默认 1：只保留恰好 1 个唯一靶基因的细胞。同一靶的多条 guide 仍算单靶。

| 状态 | 含义 |
| --- | --- |
| `no_guide` | 0 条 guide |
| `single_targeting` | 1 个靶向基因 |
| `single_nt` | 1 个设计 NT 标签 |
| `multi_guide` | 唯一靶基因数 > 1 |

过滤顺序：`no_guide → multi_guide → qc_failed → low_cells_target → retained`

---

## Barcode

`input.barcode_policy` 为 `strict`：assignment-only barcode 导致失败；matrix-only（无 assignment）记为 `no_guide` 并保留；两侧都有差异时失败。

---

## Mixscape 与 Mixscale

两个互不依赖的步骤，均只依赖 `preprocess`，以进程池并发。细胞范围：`analysis_eligible` ∪ 对照；表达嵌入仍用全部 QC-pass 细胞。

输入对齐 pertpy 1.x。差异基因由各自在未缩放 `log1p` 的 `.X` 上完成。Mixscale 自行算签名层 `X_pert_mixscale`，不读 Mixscape 的结果。

| | Mixscape | Mixscale |
| --- | --- | --- |
| 输出目录 | `Mixscape/` | `mixscale/` |
| 签名层 | `X_pert` | `X_pert_mixscale` |
| `logfc_threshold` | 0.1 | 0.1 |
| `min_de_genes` | 5 | 5 |
| `pval_cutoff` | 0.05 | 0.05 |
| 主要产物 | KO/NP 分类、signature UMAP、2x2 | `mixscale_score`、DE 热图 |

Mixscape 标签：`nt` / `perturbed` / `non_perturbed` / `ambiguous` / `not_evaluated`。2x2 图含后验 / 火山 / 热图，热图对齐 Seurat MixscapeHeatmap（`max.genes=20`）。Mixscale 热图对齐 Mixscale_DoHeatmap（列块 NT | 靶，按 score 排序）。

> KO/NP 与 Mixscale 分数都不是生物学重复推断。

---

## 预处理与簇注释

1. 从 `layers['counts']` 对全部基因做文库归一化（默认 1e4）→ `log1p`
2. HVG：Seurat flavor，top 2000，去掉线粒体基因后补满
3. HVG 上 z-score → PCA（30 维）；`X` 保持全基因 log 归一化
4. 邻居、UMAP、Leiden（`leiden_resolution` 默认 0.5）用 HVG-PCA

细胞范围 = QC-pass 全部细胞，非仅 `analysis_eligible`。

**细胞周期**：scib Tirosh/Regev 基因集打分 → `phase` / `S_score` / `G2M_score`，不改 `X`。

**SingleR**：R + celldex，标签 → `obs['singler_label']`。`singler.enabled: false` 跳过。

**ScType**：CLI `--sctype-tissue` 启用，标签 → `obs['sctype_classification']`。
