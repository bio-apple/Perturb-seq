# 输出目录

路径相对于 `--output-dir`。图只写 PNG，dpi 为 `plotting.dpi`（默认 300）。无 PDF、无 HTML。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [跳过文件写在哪](#跳过文件写在哪)
3. [校验与分配](#校验与分配)
4. [预处理、效应、距离、富集](#预处理效应距离富集)
5. [report 图](#report-图)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | mRNA / gRNA 统计（`figures`、`tables`） |
| `DGE/` | Mixscape、Mixscale、描述性 log2FC、Augur、Enrichr（`figures/<section>`、`tables`、`plot_data`） |
| `UMAP/` | 表达 UMAP、簇注释、细胞周期、Mixscape signature |
| `E-distance/` | 对当前对照的 E-test；guide 伪批量距离热图 / MDS |
| 输出根目录 | 效应表、checkpoint、`adata_processed.h5ad`、`resolved_config.yaml`、`run_info.json` |

分桶由 `FigureStore` 按图名决定。

---

## 跳过文件写在哪

| 图名类型 | 跳过说明路径 |
| --- | --- |
| QC / HVG / PCA | `qc/tables/<图名>.skipped.txt` |
| DGE 集合（含 `volcano`） | `DGE/plot_data/<图名>.skipped.txt` |
| 以 `umap_` 开头 | `UMAP/tables/<图名>.skipped.txt` |
| `perturbation_distance_heatmap`、`perturbation_embedding` | `E-distance/tables/<图名>.skipped.txt` |
| E-test 缺 `X_pca` | `E-distance/figures/edistance_etest.skipped.txt`（`missing_X_pca`） |
| `target_effect_forest` | `plot_data/target_effect_forest.skipped.txt` |
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
| `cell_status.tsv`、`filter_stage_counts.tsv`、`tables/qc_thresholds.json`、`checkpoints/assign.h5ad` | `assign` 完成（含 UMI/mito、表达基因数 100–6000、ribo、`low_library_no_guide`） | 无 `io` checkpoint |
| `qc/tables/qc_metrics.tsv` | UMI、基因数、mito、ribo、guide UMI、捕获效率 | 同上 |
| `qc/tables/grna_capture_library_correlation.tsv` | Spearman；含低文库排除是否执行 | 缺列 / 细胞 < 2 / 方差 0 |
| `qc/tables/assignment_summary.tsv` | 四类分配计数 | 同上 |
| `qc/tables/cells_per_guide.tsv`、`cells_per_target.tsv` | guide / 靶基因细胞数 | 同上 |

过滤顺序：`low_library_no_guide` → `no_guide` → `multi_guide` → `low_guide_umi` → `qc_failed` → `low_cells_target` → `retained`。

`low_library_no_guide` 仍在账本，状态 `no_guide`，`qc_pass=false`，不进对照。

---

## 预处理、效应、距离、富集

| 路径 | 何时有 | 何时缺 |
| --- | --- | --- |
| `checkpoints/preprocess.h5ad` | 归一化、HVG、PCA、可选 Harmony；scib 细胞周期 | 无分析细胞 |
| `DGE/tables/mixscape_de_gate.tsv` | Mixscape 门槛（默认 logfc 0.25） | `no_zero_guide_cells`；`mixscape.enabled=false` |
| `DGE/tables/mixscape_summary.tsv` | 每靶 gene 的 KO/NP/ambiguous 比例 | Mixscape 未跑 |
| `effect_main.tsv`、`effect_sensitivity.tsv` | 描述性 CPM log2FC + bootstrap 区间；`pvalue`/`fdr`/`ci_*` 为 NA | 无对照 / 整支跳过时仅表头 |
| `guide_effects.tsv`、`target_expression_effect.tsv` | guide 效应；靶基因自身 RNA | 无可计算靶基因 |
| `DGE/tables/de_edger_*.tsv`、`de_pydeseq2_*.tsv` | 仅表头 | `single_sample_de_not_estimable` |
| `comparison_status.tsv` | effect / edgeR / PyDESeq2 状态 | `perturb` 未到 |
| `pseudobulk/guide_metadata.tsv`、`guide_counts.h5ad` | 按 `guide_id` 加总 | 无 guide |
| `E-distance/tables/edistance_etest.tsv`、图 | 对当前对照的 E-test（`X_pca`，细胞置换） | `no_zero_guide_cells` / `missing_X_pca` |
| `perturbation_distance.tsv`、`perturbation_embedding.tsv` | 需已配置 `perturb.min_guide_replicates` | `min_guide_replicates_not_set` / `insufficient_cells` |
| `augur_scores.tsv`、`checkpoints/augur.h5ad` | Augur | 细胞不足为失败行 |
| `enrichment_results.tsv`、`tables/enrichment_status.txt` | Enrichr（需联网） | `disabled` / 无基因 / `enrichr_failed` 等 |
| `qc_summary.tsv`、`perturbation_summary.tsv`、`plot_manifest.tsv`、`adata_processed.h5ad` | `report` 完成；`n_significant_de_genes` 为 NA | 无 `cells_per_target` 时摘要空 |

效应表 bootstrap 列与重抽样规则见 [analysis.md § 统计](analysis.md#统计)。

---

## report 图

### QC

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `qc/figures/guides_per_cell.png` | 每细胞 guide 数；含 no_guide / multi | 无 `n_guides` |
| `qc/figures/assignment_status_bar.png` | 四类分配 | 无 `assignment_status` |
| `qc/figures/cells_per_guide_bar.png` | guide 细胞数（虚线 = 门槛） | 表空 |
| `qc/figures/cells_per_target_bar.png` | 靶基因细胞数 | 表空 |
| `qc/figures/target_cell_counts.png` | 函数可画；`report` 当前不调用 | 未调用 |
| `qc/figures/multi_guide_composition_bar.png` | multi-guide 组成 | 无 multi |
| `qc/figures/qc_scatter.png` | nCount vs mito / nFeature | 缺 QC 列 |
| `qc/figures/qc_violin.png` | nFeature、nCount、percent.mt；有则加 percent.ribo | 缺 QC 列 |
| `qc/figures/grna_capture_vs_library_size.png` | 捕获效率 vs 文库大小 | `missing_fields` / `no_finite_points` |
| `qc/figures/temporary_control_qc.png` | 仅 zero-guide：对照 vs targeting 的 UMI / 基因数 / mt | 设计 NT；缺列或空组 |
| `qc/figures/cell_filter_waterfall.png` | 过滤瀑布 | 阶段表空 |
| `qc/figures/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | HVG / PCA | 无预处理对象 |

### UMAP

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `umap_assignment` / `umap_target` / `umap_status` / `umap_mixscape_class` / `umap_leiden` | 表达 UMAP | 无 `X_umap` 或列 |
| `umap_cell_type.png` | 真实细胞类型 | 仅 `ALL` → `unstratified_cell_type` |
| `umap_cluster_annotation.png` | Leiden + SingleR（需 R / celldex，首次下载需联网） | `singler_disabled`；无 leiden；无 Rscript |
| `umap_nt_perturbed_cell_cycle.png` | Perturbed / Status / Phase（scib human/mouse） | 缺基因 / 无 scib / 无 phase |
| `umap_mixscape.png` | signature 空间 | `signature_umap_not_run`；非有限坐标 |
| `umap_mixscape_class_global.png` | 已分类细胞 | 无有限坐标或已分类细胞 |

### DGE / E-distance / 摘要

| 路径 | 内容 | 会缺 |
| --- | --- | --- |
| `mixscape/mixscape_class_global_bar.png` | KO / NP / 对照计数 | Mixscape 未分出类 |
| `mixscape/{靶基因}_mixscape_2x2.png` | 火山 + 热图（门槛 0.25） | 缺 NP/NT/KO 或分数 |
| `mixscape/mixscape_fraction_bar.png` | 每靶比例 | 摘要空 |
| `mixscale/mixscale_score_ridge.png` | 山脊图（门槛 0.2） | `mixscale_not_run` / `de_gate_not_run` |
| `mixscale/mixscale_scatter.png` | 分数箱 vs 靶基因表达 | `no_scored_de_targets` 等 |
| `target_in_de/{靶基因}_mixscape_2x2.png` | 靶基因在自身 DE 名单里的复制件 | 无命中 |
| `plot/target_de_counts.png` | DE 计数三柱 | 无门槛表 |
| `plot/target_in_de_gene_counts.png` | 自身在 DE 名单的靶基因 | `no_targets_in_own_de` |
| `de/de_gene_heatmap.png` | 描述性 log2FC（行 = gene_id） | 描述性行 < 2 |
| `de/target_expression_effect.png` | 仅 `branch=main` | 无效应列 |
| `de/guide_consistency.png`、`pseudobulk_*.png` | guide 一致性 / 伪批量 | 表空；guide < 3 |
| `DGE/plot_data/volcano.skipped.txt` | `inference_not_available` | 单样本不画单独火山图 |
| `plot_data/target_effect_forest.skipped.txt` | `inference_not_available` | 同上 |
| `E-distance/…/perturbation_distance_heatmap.png` 等 | guide 伪批量两两距离（标题标明非 E-test） | 无 `min_guide_replicates` |
| `summary/augur_score_bar.png` | Augur AUC | 表空或细胞不足 |
| `summary/enrichment_dotplot.png` | Enrichr | `no_candidate_genes` |
| `summary/perturbation_summary.png` | 摘要条形图 | 表空或数值全 NA |

`report` 不调用 `save_target_log2fc_bar`，因此不会重画 `{靶基因}_log2fc.png`。
