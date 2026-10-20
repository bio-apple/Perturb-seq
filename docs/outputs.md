# 输出目录

路径相对于 `--output-dir`。图只写 PNG，dpi 为 `plotting.dpi`（默认 300）。无 PDF、无同名图旁 `.txt`、无 HTML。成功写出 PNG 时会删除同名旧 `.pdf` / `.txt`。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [跳过文件写在哪](#跳过文件写在哪)
3. [校验与分配](#校验与分配)
4. [预处理、效应、距离](#预处理效应距离)
5. [report 图](#report-图)
6. [旧布局兼容](#旧布局兼容)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | mRNA / gRNA 统计（`figures`、`tables`） |
| `UMAP/` | 表达 UMAP（全细胞 HVG scale）、簇注释、细胞周期、Mixscape signature（`X_pert`，相对对照） |
| `DGE/Mixscape/` | Mixscape 门槛、摘要、class / fraction / 2x2 / target_in_de（`figures`、`tables`、`plot_data`） |
| `DGE/mixscale/` | Mixscale ridge / scatter、门槛表（`figures`、`tables`、`plot_data`） |
| `DGE/descriptive/` | 描述性 CPM log2FC、guide 一致性、伪批量图与效应表（非 Mixscape / Mixscale 推断） |
| `E-distance/` | 对当前对照的 E-test；guide 伪批量距离热图 / MDS 与对应表 |
| 输出根目录 | checkpoint、`adata_processed.h5ad`、`resolved_config.yaml`、`run_info.json`、摘要表 |

分桶由 `FigureStore` 按图名与 section 决定。大小写约定：`DGE/Mixscape/`（首字母大写）与 `DGE/mixscale/`（全小写）；读取时也接受旧的 `DGE/mixscape/`。

---

## 跳过文件写在哪

| 图名类型 | 跳过说明路径 |
| --- | --- |
| QC / HVG / PCA | `qc/tables/<图名>.skipped.txt` |
| Mixscape 集合 | `DGE/Mixscape/plot_data/<图名>.skipped.txt` |
| Mixscale 集合 | `DGE/mixscale/plot_data/<图名>.skipped.txt` |
| 描述性 DE（含 `volcano`、`target_effect_forest`） | `DGE/descriptive/plot_data/<图名>.skipped.txt` |
| 以 `umap_` 开头 | `UMAP/tables/<图名>.skipped.txt` |
| `perturbation_distance_heatmap`、`perturbation_embedding` | `E-distance/tables/<图名>.skipped.txt` |
| E-test 缺 `X_pca` | `E-distance/figures/edistance_etest.skipped.txt`（`missing_X_pca`） |
| 其余 | 输出根下 `plot_data/` |

---

## 校验与分配

| 路径 | 何时有 | 何时缺 |
| --- | --- | --- |
| `resolved_config.yaml`、`run_info.json` | `run` 或 `--validate-only`；含物种核对、`temporary_control_warning` | 配置加载失败无文件；物种失败仍写 `resolved_config.yaml` 后停止 |
| `tables/temporary_control_warning.txt` | zero-guide 临时对照进入 `run` | 设计 NT；`--validate-only` |
| `barcode_differences.tsv`、`input_validation.json` | `io` / `--validate-only` | 集合一致时差异表仅表头 |
| `tables/guide_reference_conflicts.tsv` | 参考表冲突 | 参考表通过 |
| `tables/cell_calls.tsv`、`io_checks.json`、`multi_guide_composition.tsv`、`checkpoints/io.h5ad` | `io` 通过 | assignment-only 时提前退出 |
| `cell_status.tsv`、`filter_stage_counts.tsv`、`tables/qc_thresholds.json`、`checkpoints/assign.h5ad` | `assign` 完成（固定阈值：mito≥20%、表达基因数 200–6000；UMI / ribo 仅在显式给出时过滤；`low_library_no_guide` 已禁用） | 无 `io` checkpoint |
| `qc/tables/qc_metrics.tsv` | UMI、基因数、mito、ribo、guide UMI、捕获效率 | 同上 |
| `qc/tables/grna_capture_library_correlation.tsv` | Spearman；含低文库排除是否执行 | 缺列 / 细胞 < 2 / 方差 0 |
| `qc/tables/assignment_summary.tsv` | 四类分配计数 | 同上 |
| `qc/tables/cells_per_guide.tsv`、`cells_per_target.tsv` | guide / 靶基因细胞数 | 同上 |

过滤顺序：`low_library_no_guide` → `no_guide` → `multi_guide` → `low_guide_umi` → `qc_failed` → `low_cells_target` → `retained`。

`low_library_no_guide` 阶段仍在账本，但规则已禁用（原 MAD 路径已移除），新增排除数为 0。

---

## 预处理、效应、距离

| 路径 | 何时有 | 何时缺 |
| --- | --- | --- |
| `checkpoints/preprocess.h5ad` | 全基因归一化；PCA（默认 `n_pcs=30`）/ neighbors / Leiden / UMAP 在 top 2000 HVG 上；可选 Harmony；scib 细胞周期；`uns['preprocess']` 记参数 | 无分析细胞 |
| `checkpoints/mixscape.h5ad` | Mixscape 分类、`X_pert`、signature UMAP（`n_dims` 默认 30，与 `preprocess.n_pcs` 一致） | `mixscape` 未跑 |
| `checkpoints/mixscale.h5ad` | Mixscale 分数（依赖 `X_pert`） | `mixscale` 未跑；缺 `X_pert` 时跳过并写 checkpoint |
| `checkpoints/perturb.h5ad` | 描述性效应、伪批量与 E-distance 之后 | `perturb` 未跑 |
| `DGE/Mixscape/tables/mixscape_de_gate.tsv` | Mixscape 门槛（默认 logfc 0.25） | `no_zero_guide_cells`；`mixscape.enabled=false` |
| `DGE/Mixscape/tables/mixscape_summary.tsv` | 每靶 gene 的 KO/NP/ambiguous 比例 | Mixscape 未跑 |
| `DGE/mixscale/tables/mixscale_de_gate.tsv` | Mixscale 自己的 DE 门槛 | Mixscale 未跑 |
| `DGE/descriptive/tables/effect_main.tsv`、`effect_sensitivity.tsv` | 描述性 CPM log2FC + bootstrap 区间；`pvalue`/`fdr`/`ci_*` 为 NA | 无对照 / 整支跳过时仅表头 |
| `DGE/descriptive/tables/guide_effects.tsv`、`target_expression_effect.tsv` | guide 效应；靶基因自身 RNA | 无可计算靶基因 |
| `DGE/descriptive/tables/de_edger_*.tsv`、`de_pydeseq2_*.tsv` | 仅表头 | `single_sample_de_not_estimable` |
| `comparison_status.tsv` | effect / edgeR / PyDESeq2 状态 | `perturb` 未到 |
| `pseudobulk/guide_metadata.tsv`、`guide_counts.h5ad` | 按 `guide_id` 加总 | 无 guide |
| `E-distance/tables/edistance_etest.tsv`、图 | 对当前对照的 E-test（`X_pca`，细胞置换） | `no_zero_guide_cells` / `missing_X_pca` |
| `E-distance/tables/perturbation_distance.tsv`、`perturbation_embedding.tsv` | 需已配置 `perturb.min_guide_replicates` | `min_guide_replicates_not_set` / `insufficient_cells` |
| `qc_summary.tsv`、`perturbation_summary.tsv`、`plot_manifest.tsv`、`adata_processed.h5ad` | `report` 完成；`n_significant_de_genes` 为 NA | 无 `cells_per_target` 时摘要空 |

效应表 bootstrap 列与重抽样规则见 [analysis.md § 统计](analysis.md#统计)。

---

## report 图

图目录只写 PNG；不写 PDF，也不写同名图旁 `.txt`。方法说明在 `plot_manifest.tsv` 的 `method` 列。跳过原因仍写在 `*.skipped.txt`（tables / plot_data / 少数 figures 旁的 skip 标记），与分析表 TSV、`temporary_control_warning.txt` 等非画图产物无关。

### QC

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `qc/figures/guides_per_cell.png` | 每细胞 guide 条数分布（描述用；资格看 `n_targets`） | 无 `n_guides` |
| `qc/figures/assignment_status_bar.png` | 四类分配 | 无 `assignment_status` |
| `qc/figures/cells_per_target_bar.png` | 靶基因细胞数 | 表空 |
| `qc/figures/target_cell_counts.png` | 函数可画；`report` 当前不调用 | 未调用 |
| `qc/figures/qc_scatter.png` | nCount vs mito / nFeature（全 annotated 账本优先；过滤前） | 缺 QC 列 |
| `qc/figures/qc_violin.png` | nFeature、nCount、percent.mt；有则加 percent.ribo（同上） | 缺 QC 列 |
| `qc/figures/temporary_control_qc.png` | 仅 zero-guide：对照 vs targeting 的 UMI / 基因数 / mt | 设计 NT；缺列或空组 |
| `qc/figures/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | HVG / PCA；`pca_scatter` 为 PC1–PC2（有 Leiden 则着色），图例在图下方多列，显示轴取 0.5–99.5% 百分位窗口（不改 `X_pca`） | 无预处理对象 |

### UMAP

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `umap_mixscape_class` | 表达 UMAP（`mixscape_class`） | 无 `X_umap` 或列 |
| `umap_cell_type.png` | 真实细胞类型 | 仅 `ALL` → `unstratified_cell_type` |
| `umap_cluster_annotation.png` | Leiden + SingleR（左 Leiden，右注释；不再单独出 `umap_leiden`） | `singler_disabled`；无 leiden；无 Rscript |
| `umap_nt_perturbed_cell_cycle.png` | Perturbed / 对照（设计 NT 或 zero-guide cells）/ Phase（scib human/mouse） | 缺基因 / 无 scib / 无 phase |
| `umap_mixscape.png` | signature 空间 | `signature_umap_not_run`；非有限坐标 |
| `umap_mixscape_class_global.png` | 已分类细胞 | 无有限坐标或已分类细胞 |

### DGE / E-distance / 摘要

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `DGE/Mixscape/figures/mixscape_class_global_bar.png` | KO / NP / 对照计数 | Mixscape 未分出类 |
| `DGE/Mixscape/figures/{靶基因}_mixscape_2x2.png` | 火山 + 热图（门槛 0.25） | 缺 NP/NT/KO 或分数 |
| `DGE/Mixscape/figures/mixscape_fraction_bar.png` | 每靶比例 | 摘要空 |
| `DGE/Mixscape/figures/target_de_counts.png` | DE 计数三柱 | 无门槛表 |
| `DGE/Mixscape/figures/target_in_de_gene_counts.png` | 自身在 DE 名单的靶基因 | `no_targets_in_own_de` |
| `DGE/Mixscape/figures/target_in_de/{靶基因}_mixscape_2x2.png` | 靶基因在自身 DE 名单里的复制件 | 无命中 |
| `DGE/mixscale/figures/mixscale_score_ridge.png` | 山脊图（门槛 0.2） | `mixscale_not_run` / `de_gate_not_run` |
| `DGE/mixscale/figures/mixscale_scatter.png` | 分数箱 vs 靶基因表达 | `no_scored_de_targets` 等 |
| `DGE/descriptive/figures/de_gene_heatmap.png` | 描述性 log2FC（行 = gene_id） | 描述性行 < 2 |
| `DGE/descriptive/figures/target_expression_effect.png` | 仅 `branch=main` | 无效应列 |
| `DGE/descriptive/figures/guide_consistency.png`、`pseudobulk_*.png` | guide 一致性 / 伪批量 | 表空；guide < 3 |
| `DGE/descriptive/plot_data/volcano.skipped.txt` | `inference_not_available` | 单样本不画单独火山图 |
| `DGE/descriptive/plot_data/target_effect_forest.skipped.txt` | `inference_not_available` | 同上 |
| `E-distance/figures/edistance_etest.png` | 表达 `X_pca` 上相对对照的 E-test；着色为样本内校正显著；非生物学重复 | `no_zero_guide_cells` / `missing_X_pca` |
| `E-distance/…/perturbation_distance_heatmap.png` 等 | guide 伪批量两两距离（标题标明非 E-test） | 无 `min_guide_replicates` |
| `figures/summary/perturbation_summary.png` | 摘要条形图；descriptive 模式不含显著计数 / FDR 阈值列 | 表空或数值全 NA |

`report` 不调用 `save_target_log2fc_bar`，因此不会重画 `{靶基因}_log2fc.png`。

---

## 旧布局兼容

新跑会写入上表路径，并尽量删除同名旧产物。`report` / E-test 读表时仍按顺序尝试旧路径：

| 新路径 | 仍可读的旧路径 |
| --- | --- |
| `DGE/Mixscape/tables/mixscape_*.tsv` | `DGE/tables/`、`DGE/mixscape/tables/`、根 `tables/` |
| `DGE/descriptive/tables/effect_*.tsv` 等 | 输出根目录同名文件 |
| `E-distance/tables/perturbation_*.tsv` | 输出根目录同名文件 |

仅重跑 `report` 即可把图迁到新桶；Mixscape / Mixscale 表对应 `--steps mixscape` / `mixscale`；描述性效应表需再跑 `perturb`（或全流程）才会落到 `DGE/descriptive/tables/`。
