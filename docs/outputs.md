# 输出目录

路径相对于命令行 `--output-dir`。图只保存 PNG，dpi 用 `plotting.dpi`（默认 300）。没有 PDF，没有 HTML 报告。

`FigureStore.skip` 按图名写原因，不看调用时的 section：mRNA/gRNA 统计和高变基因、PCA 写 `qc/tables/<图名>.skipped.txt`，`_DGE_FIGURES` 里的图写 `DGE/plot_data/<图名>.skipped.txt`，文件名以 `umap_` 开头的图写 `UMAP/tables/<图名>.skipped.txt`，`perturbation_distance_heatmap` 和 `perturbation_embedding` 写 `E-distance/tables/<图名>.skipped.txt`。其余图名写到输出根下的 `plot_data/`。`volcano` 在 DGE 集合里，所以单独火山图的跳过说明是 `DGE/plot_data/volcano.skipped.txt`。`target_effect_forest` 不在分桶名单里，跳过说明在 `plot_data/target_effect_forest.skipped.txt`。E-test 缺 `X_pca` 时，`distance_tests` 直接写 `E-distance/figures/edistance_etest.skipped.txt`，内容 `missing_X_pca`。

下面是当前代码的输出，不绑定某一次样本运行。`config/config.yaml` 只含分析软件默认参数。样本、物种、目录和对照来自命令行。命令行、分配、barcode 和各步骤的判定见 [analysis.md](analysis.md)。

mRNA 和 gRNA 的统计图、表在 `qc/`（`qc/figures` 与 `qc/tables`），包括 `guides_per_cell`、分配柱、每个 guide / 靶基因的细胞数、过滤瀑布、高变基因和 PCA。差异基因门槛、Mixscape 2×2、Mixscale 山脊图和散点图、靶基因 DE 计数、靶基因出现在自己 DE 名单里的计数、描述性 log2FC 热图和靶基因效应图、guide 伪批量图、Augur 和 Enrichr 在 `DGE/`（`DGE/figures/<section>`、`DGE/tables`、`DGE/plot_data`）。文件名以 `umap_` 开头的图在 `UMAP/`（`UMAP/figures` 与 `UMAP/tables`），包括簇注释、细胞周期 UMAP 和 Mixscape signature UMAP。每个通过 Mixscape 差异基因门槛的靶基因对当前对照（设计好的 NT，或 `zero_guide`）的能量距离在 `E-distance/`（`E-distance/figures` 与 `E-distance/tables`）。guide 伪批量的两两距离热图和 MDS 也进这个目录，标题写明那是 guide 伪批量距离。描述性 CPM 效应表、原始输入、checkpoint、h5ad 在输出根目录，不进 `DGE/`。分桶由 `FigureStore` 按图名决定。

## 代码会写的图和表

### 校验和分配

| 路径 | 何时有 | 何时缺 |
| --- | --- | --- |
| `resolved_config.yaml`、`run_info.json` | 进入 `run` 或 `--validate-only` 就会写配置；`run_info.json` 在每个步骤之后更新。`biological_sample_count` 为 1。`experiment.organism` 只能是输入的 `human` 或 `mouse`。`resolved_organism` 是这个输入值，不是从基因号推断的。`resolved_mito_prefix` 是采用的线粒体前缀（人类 `MT-`，小鼠 `mt-`）。`species_check` 为 `passed` 或 `failed`。`organism_evidence` 写下基因号和符号的数量与例子。`other_organelle_prefixes` 记下叶绿体等非线粒体细胞器。进入 `run` 后 `temporary_control_warning` 有值：zero-guide 临时对照时是完整警告，设计好的 NT 对照时是空字符串 | 配置加载失败时没有。物种一致性验证失败时仍写 `resolved_config.yaml`（`species_check` 为 `failed`），然后停止 |
| `tables/temporary_control_warning.txt` | `run` 一开始，且当前对照是临时 zero-guide cells。正文与运行日志、`run_info.json` 的 `temporary_control_warning` 相同：这些细胞可能在 UMI、基因数、线粒体比例、环境 RNA 和细胞状态上与靶向细胞不同，后续 log2FC、Mixscape、Mixscale、E-distance 都以这组临时对照为参照。分布对比只作描述 | 设计好的 NT 对照不写该文件。`--validate-only` 不进入 `run`，因此也不写 |
| `barcode_differences.tsv`、`input_validation.json` | `--validate-only` 或 `io`。`input.barcode_policy` 必须是 `strict`。其它取值在加载配置时会被拒绝，也不会改成两边求交。规则不看文件名或样本名。assignment-only 使 `status=failed` 并停止。仅有 matrix-only 时 `status=passed`，这些细胞记为 `no_guide`，`n_matrix_only_as_no_guide` 是保留数量；没有设计 NT 时，它们和分配表里 0 条 guide 的细胞一起作为 `zero_guide` 临时对照。两侧都有差异时仍然失败 | 集合一致时差异表只有表头 |
| `tables/guide_reference_conflicts.tsv` | guide 参考表里同一 id 多个靶基因，或非 NT guide 的 `target_gene_name` 为空 | 参考表通过时不写 |
| `tables/cell_calls.tsv`、`tables/io_checks.json`、`tables/multi_guide_composition.tsv`、`checkpoints/io.h5ad` | `io` 通过之后。matrix-only 细胞在 `cell_calls` 里是 `no_guide` | 存在 assignment-only 时 `io` 在这些文件之前退出 |
| `cell_status.tsv`、`filter_stage_counts.tsv`、`tables/qc_thresholds.json`、`checkpoints/assign.h5ad` | `assign` 完成。`tables/qc_thresholds.json` 记录四个 QC 绝对阈值。命令行同时给出 `--min-umi`、`--max-umi`、`--min-genes`、`--max-mito-frac` 时用这四个数；只给其中一部分会失败；全部留空时 `source` 为 `median_mad`：UMI 和基因数先 `log1p`，界是 median ± `qc.n_mad`（默认 3）× 1.4826 × MAD，再 `expm1`；线粒体比例在原分数上用同一公式，并限制在 0 到 1。MAD 为 0 时改用标准差。命令行给齐四个数时 `source` 为 `command_line`。同一文件还写 `max_ribo_frac` 和 `ribo_frac_source`（未给 `--max-ribo-frac` 时为 `median_mad`，单独给出时为 `command_line`），以及 `low_library_no_guide`（是否执行、原因、`min_total_counts`、来源 `median_mad` 或 `skipped`） | 没有 `io` 的 checkpoint 时这一步不会开始 |
| `qc/tables/qc_metrics.tsv` | `total_counts`、`n_genes_by_counts`、`pct_counts_mito`、`pct_counts_ribo`、`guide_umi_total`、`grna_capture_efficiency`、`n_guides` | 同上 |
| `qc/tables/grna_capture_library_correlation.tsv` | QC 通过前全部细胞上，gRNA 捕获效率与文库大小的 Spearman 系数、细胞数和 p 值。同时记下 `low_library_no_guide` 是否执行、原因和下限。这个 p 不写入 `effect_main.tsv` | 缺少捕获效率列时 `status=skipped`，`reason=missing_capture_columns`。细胞少于 2 或方差为 0 时也跳过，原因分别是 `fewer_than_2_cells`、`zero_variance` |
| `qc/tables/assignment_summary.tsv` | 四类分配状态的细胞数：`no_guide`、`single_targeting`、`single_nt`、`multi_guide` | 同上 |
| `qc/tables/cells_per_guide.tsv`、`qc/tables/cells_per_target.tsv` | guide 和靶基因的细胞数 | 同上 |

过滤阶段按这个顺序记账：`low_library_no_guide`、`no_guide`、`multi_guide`、`low_guide_umi`、`qc_failed`、`low_cells_target`，剩下的是 `retained`。`low_library_no_guide` 仍记在分配账本里，分配状态保持 `no_guide`，但 `qc_pass` 为 false，不进入设计 NT，也不进入 `zero_guide`。规则跳过时该阶段新增排除数为 0，原因在 `tables/qc_thresholds.json`。

### 预处理、效应、距离、富集

| 路径 | 何时有 | 何时缺 |
| --- | --- | --- |
| `checkpoints/preprocess.h5ad` | `preprocess` 完成。归一化到配置的 `normalize_target_sum` 再 `log1p`。高变基因去掉该物种的线粒体前缀（人类 `MT-`，小鼠 `mt-`）。叶绿体等其他细胞器不从高变基因里按线粒体排除，也不计入 percent.mt。PCA 在标准化后的高变基因上，`X` 保持 log 归一化。没有外部注释时 `cell_type` 为 `ALL`，`cell_type_scope` 为 `unstratified`。效应表的 `cell_type` 列也是 `ALL`。`batch_key` 为空时不调用 Harmony。配置了 `batch_key` 且有效水平少于 2 时，`uns['batch_correction']` 为 `not_run_fewer_than_two_batch_levels`。已有 `uns['log1p']` 时用 scib 的 human/mouse 细胞周期基因打分，写入 `phase`、`S_score`、`G2M_score`，不改 `X`，也不做回归。成功时 `uns['cell_cycle_status']` 为 `scored_after_log1p`；矩阵缺基因时为 `skipped_missing_cell_cycle_genes`；未装 scib 时为 `skipped_scib_unavailable` | 没有 `assign` 的分析细胞时不会开始 |
| `DGE/tables/mixscape_de_gate.tsv` | Mixscape 门槛，`logfc_threshold` 来自 `mixscape`（默认 0.25）。少于 `min_de_genes` 时 `passed` 为 false，`reason` 为 `skipped_no_de_genes`。对照是设计好的 NT，或没有设计 NT 时的 `zero_guide` | 临时对照细胞数为 0 时不跑，原因 `no_zero_guide_cells`。`mixscape.enabled` 为 false，或 Mixscape 在写表前失败时也没有这张表 |
| `DGE/tables/mixscape_summary.tsv` | 每个靶基因的 `perturbed`、`non_perturbed`、`ambiguous` 计数和比例。分母是该靶基因的全部靶向细胞，组为空时比例写成 0.0 | Mixscape 没跑 |
| `effect_main.tsv`、`effect_sensitivity.tsv` | 输出根目录。对照是设计好的 NT，或 `zero_guide`。主效应按靶基因把原始计数加总，再换成 CPM，每个基因一行：`log2((CPM_target+pseudocount_cpm)/(CPM_control+pseudocount_cpm))`，`pseudocount_cpm` 来自配置。靶基因和对照的原始计数都是 0 时，`log2fc_descriptive` 为缺失，`effect_status` 为 `not_detected_both`。`pvalue`、`fdr`、`ci_lower`、`ci_upper` 为 NA，`inference_status` 为 `not_run_single_sample_descriptive`。另有样本内 bootstrap 百分位区间 `log2fc_bootstrap_ci_lower`、`log2fc_bootstrap_ci_upper`、`bootstrap_unit`（该靶基因的 `guide` 或 `cell`）、`bootstrap_n`、`bootstrap_ci_level`、`uncertainty_meaning`（固定 `descriptive_estimate_uncertainty_not_significance_test`）、`bootstrap_status`。至少 2 条不同 guide 时按 guide 有放回重抽样，抽到的条数等于观察到的 guide 数，合并后再算 CPM log2FC；只有 1 条 guide 时按该靶基因的细胞有放回重抽样。对照每次一起重抽：设计好的 NT 至少 2 条 guide 时按 guide，只有 1 条时按细胞；`zero_guide` 临时对照按细胞。次数和水平来自 `statistics.bootstrap_n`（默认 200）和 `statistics.bootstrap_ci`（默认 0.95），种子是 `runtime.seed`。点估计缺失或任一侧抽样单位少于 2 个时区间为缺失，`bootstrap_status` 写原因，缺失不写成 0。这些区间表示单个样本内部描述性估计的不确定性，不能当作显著性检验，也不能当作生物学重复的置信区间。细胞水平 Wilcoxon 的 `adj_p_value`、Mixscale 分数、E-test 的 p 分别留在 `DGE/plot_data/`、`obs['mixscale_score']` 和 `E-distance/tables/edistance_etest.tsv`，不写入这两张表。`branch=main` 用全部可分析细胞；`branch=sensitivity` 只用 Mixscape `perturbed` 加对照。这两个分支名写在 `perturb` 步骤里 | 临时对照细胞数为 0 时该对比跳过，原因 `no_zero_guide_cells`。对照细胞少于 `statistics.min_nt_cells`、文库大小为 0，或靶基因细胞少于 `qc.min_cells_per_target` 时，该靶基因只有 skipped 行。整支对比跳过时效应表只有表头，不填区间 |
| `guide_effects.tsv`、`target_expression_effect.tsv` | 输出根目录。对照是设计好的 NT 或 `zero_guide`。靶基因自身 RNA 在基因缺失、对照表达过低、或对照 / 该靶基因细胞数为 0 时，`target_expression_status` 为 `not_evaluable`，原因分别是 `gene_not_in_rna`、`nt_expression_below_threshold`、`insufficient_cells` | 没有对照细胞或没有可计算的靶基因时对应表为空 |
| `DGE/tables/de_edger_{main,sensitivity}.tsv`、`de_pydeseq2_{main,sensitivity}.tsv` | pertpy 差异表达教程的伪批量 edgeR 和 PyDESeq2。单样本时四张表只有表头。描述性 CPM log2FC 不写进这些列 | `comparison_status.tsv` 里两条模型的 `status` 为 `skipped`，`reason` 为 `single_sample_de_not_estimable` |
| `comparison_status.tsv` | 输出根目录。汇总 effect、edgeR 和 PyDESeq2 | `perturb` 没跑到 |
| `pseudobulk/guide_metadata.tsv`、`pseudobulk/guide_counts.h5ad` | 按真实 `guide_id` 加总计数 | 没有带 guide id 的分析细胞时矩阵为空 |
| `E-distance/tables/edistance_etest.tsv`、`E-distance/figures/edistance_etest.png` | pertpy `DistanceTest("edistance")`。对照是设计好的 NT 或 `zero_guide`。分组是 `mixscape_de_gate.tsv` 里 `passed` 为 true 的靶基因，加上该对照。嵌入是已有的 `X_pca`。置换单位是这一个样本里的细胞。`n_perms` 1000，`padj_threshold` 0.0015，校正为 holm-sidak。`inference` 为 `within_sample_distance_test`。这些 p 值留在这个目录 | 临时对照细胞数为 0 时只写 `E-distance/figures/edistance_etest.skipped.txt`，内容 `no_zero_guide_cells`。没有 `X_pca` 时同一文件的内容是 `missing_X_pca`。没有门槛表时改为检验每一个有细胞的 `target_gene` |
| `perturbation_distance.tsv`、`perturbation_embedding.tsv` | 仅当配置里已经有 `perturb.min_guide_replicates`，而且函数返回非空表。代码不替缺失的键填默认值 | 键不存在时写 `plot_data/perturbation_distance_heatmap.skipped.txt`，内容 `min_guide_replicates_not_set`。键存在但靶基因不够时，同一文件的内容是 `insufficient_cells` |
| `augur_scores.tsv`、`checkpoints/augur.h5ad` | `augur` 完成。没有注释时按 `ALL` 跑 | 细胞数不够时表里是失败行 |
| `enrichment_results.tsv`、`tables/enrichment_status.txt` | `enrich`。基因名单是 `effect_status=descriptive` 且 `\|log2fc_descriptive\| ≥ plotting.volcano_lfc`（默认 1.0）。`gseapy` 访问 Enrichr 服务器，需要联网。Enrichr 的 adjusted P 写在 `enrichment_fdr` | `enrichment.enabled` 为 false、没装 `gseapy`、物种不是 human/mouse、没有效应表、没有靶基因达到 `min_genes`，或 Enrichr 连不上。后两种原因分别是 `no_effect_ranked_genes` 和 `enrichr_failed` |
| `qc_summary.tsv`、`perturbation_summary.tsv`、`plot_manifest.tsv`、`adata_processed.h5ad` | `report` 完成。摘要里的 `n_significant_de_genes` 写成 NA | 没有 `qc/tables/cells_per_target.tsv` 时摘要为空 |

### report 会画的图

| 路径 | 内容 | 会缺的情况 |
| --- | --- | --- |
| `qc/figures/guides_per_cell.png` | 纵轴 `Cells`。分母是账本里的全部细胞，含 `no_guide` 和 `multi_guide`。柱上是整数和一位小数百分比。这张图不画 guide 数大于 5 的红线。作图表在 `qc/tables/guides_per_cell.tsv` | 没有 `cell_status.tsv` 的 `n_guides` |
| `qc/figures/assignment_status_bar.png` | 四类分配状态的细胞数，柱上是整数 | 没有 `assignment_status` |
| `qc/figures/cells_per_guide_bar.png` | guide 细胞数，虚线是 `min_cells_per_guide`。画出的 guide 不超过 40 个时，柱端标整数 | 表为空 |
| `qc/figures/cells_per_target_bar.png` | 每个靶基因的细胞数，虚线是 `min_cells_per_target`。不超过 40 个时标基因名和整数 | 表为空 |
| `qc/figures/target_cell_counts.png` | `draw_target_cell_counts` 会写这张图、同名 txt 和 `qc/tables/target_cell_counts.tsv`。细胞数为 0 的靶基因不画柱，并标 y=20 和 y=40。`report` 当前不调用这个函数 | 未调用时保持缺图 |
| `qc/figures/multi_guide_composition_bar.png` | multi-guide 的组成，柱上是整数细胞数 | 没有 multi-guide 细胞 |
| `qc/figures/qc_scatter.png` | nCount 对 percent.mt、nCount 对 nFeature | 缺少 QC 列 |
| `qc/figures/qc_violin.png` | `nFeature_RNA`、`nCount_RNA`、`percent.mt`（0–100）。账本里有核糖体比例时再加 `percent.ribo`（0–100）。文件名是 `qc_violin.png` | 缺少 QC 列 |
| `qc/figures/grna_capture_vs_library_size.png` | 标题是 gRNA 捕获效率与文库大小的关系。横轴是文库大小 `total_counts`，纵轴是 gRNA 捕获效率。没检测到 guide 的细胞单独标出。执行了低文库排除时画 `low_library_no_guide` 下限。作图表在 `qc/tables/grna_capture_vs_library_size.tsv`。Spearman 结果在 `qc/tables/grna_capture_library_correlation.tsv` | 缺少 `total_counts` 或 `grna_capture_efficiency` 时 skip，原因 `missing_fields`。没有有限点时原因 `no_finite_points` |
| `qc/figures/temporary_control_qc.png` | 只在 zero-guide 临时对照时由 `report` 绘制。同一张图比较通过 QC 的 zero-guide cells（临时对照）和 `perturbation_status=targeting` 的可分析细胞：UMI（`total_counts`）、基因数（`n_genes_by_counts`）、`percent.mt`（0–100）。作图表 `qc/tables/temporary_control_qc.tsv` 的列是 `group`、`metric`、`n_cells`、`median`、`mad`、`q1`、`q3`。组名是 `zero-guide cells（临时对照）` 和 `targeting`。对比只作描述，分布有差异时分析继续 | 设计好的 NT 对照不生成这张图和这张表。缺少 QC 列时 skip，原因 `missing_qc_columns`。缺少分组所需列时原因 `missing_fields`。没有通过 QC 的 zero-guide 细胞时原因含 `no_zero_guide_cells`，没有可分析靶向细胞时含 `no_targeting_cells`。跳过说明在 `qc/tables/temporary_control_qc.skipped.txt` |
| `qc/figures/cell_filter_waterfall.png` | 各过滤阶段剩下的细胞，柱上是整数 | 阶段表为空 |
| `qc/figures/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | 高变基因和 PCA。作图表在 `qc/tables/` | 没有预处理后的对象 |
| `UMAP/figures/umap_assignment.png`、`umap_target.png`、`umap_status.png`、`umap_mixscape_class.png`、`umap_leiden.png` | 表达 UMAP。作图数据在 `UMAP/tables/` | 没有 `X_umap` 或对应列 |
| `UMAP/figures/umap_cell_type.png` | 有真实细胞类型注释时才画 | 标签只有 `ALL` 时 skip，原因 `unstratified_cell_type` |
| `UMAP/figures/umap_cluster_annotation.png` | Leiden（igraph）和 SingleR 簇标签并排。先对簇内 log 归一化表达取平均，再由 `Rscript` 调用 SingleR 与 celldex。小鼠参考 `MouseRNAseqData`，人类参考 `HumanPrimaryCellAtlasData`，标签字段 `label.main`。参考必须和 `experiment.organism` 是同一个物种。图上写 SingleR、参考名和 Leiden 分辨率。说明 txt 与 PNG 放在一起，表在 `UMAP/tables/`。标签在 `obs['singler_label']`，`cell_type` 保持原值。第一次运行会从 ExperimentHub 下载参考，需要联网 | `singler.enabled` 为 false 时 skip，原因 `singler_disabled`，不调用 R。没有 leiden 或 `X_umap` 时原因 `missing_fields`。找不到 `Rscript` 时原因是 `没有找到 Rscript，无法运行 SingleR。` |
| `UMAP/figures/umap_nt_perturbed_cell_cycle.png` | 同一套表达 UMAP 的三个面板：Perturbed、Perturbation Status、Cell Cycle Phase。对照 / targeting 来自 `perturbation_status`。人类与小鼠在 `X` 已 `log1p` 时用 scib 细胞周期基因打分（`phase` / `S_score` / `G2M_score`）。PNG 旁有 txt，没有单独的 tsv | 缺基因（`skipped_missing_cell_cycle_genes`）、未装 scib、没有 `log1p` 或没有 `X_umap` / `phase` 时不画 |
| `UMAP/figures/umap_mixscape.png` | signature 空间，按 `mixscape_class` 着色 | 没有 `X_umap_mixscape` 时 `signature_umap_not_run`。未进入签名的细胞坐标是缺失值，任一缺失就 skip，原因 `non_finite_coordinates` |
| `UMAP/figures/umap_mixscape_class_global.png` | signature 空间里已分类的细胞：去掉 `not_evaluated` 和 `skipped_no_de_genes`，也去掉缺失坐标 | 没有有限坐标或没有已分类细胞 |
| `DGE/figures/mixscape/mixscape_class_global_bar.png` | 已分类细胞的计数：KO、NP、对照，柱上是整数。`CRISPR_KO` 时扰动类显示为 KO | Mixscape 没跑，或没有分出类 |
| `DGE/figures/mixscape/{靶基因}_mixscape_2x2.png` | 通过 Mixscape 门槛，并且 NP、NT、KO 和扰动分数都有的靶基因。火山图是细胞水平 Wilcoxon，竖线是 Mixscape `logfc_threshold`（默认 0.25）。`plotting.volcano_lfc`（默认 1.0）不画到这张 2×2 上。右下角热图最多 20 个基因、每组最多 300 个细胞，各组不抽成同样多，列序 NP、NT、KO，组间空白列宽是 `max(6, round(0.08 × 最大组细胞数))`。KO 的扰动分数颜色是 `#EE6A50`。作图数据在 `DGE/plot_data/`，其中 `adj_p_value` 属于这张火山图 | 缺 NP、NT 或 KO，或没有分数。记录在 `DGE/plot_data/mixscape_target_panels.tsv`，原因 `missing_ko_np_nt_or_score` |
| `DGE/figures/mixscape/mixscape_fraction_bar.png` | 每个靶基因的 perturbed、non_perturbed、ambiguous 比例 | 摘要为空 |
| `DGE/figures/mixscale/mixscale_score_ridge.png` | 对照，以及通过 Mixscale 门槛并且有有限分数的靶基因。门槛是 `min_de_genes=5`、`logfc_threshold=0.2`、`max_de_genes=100`、`pval_cutoff=0.05`、`fine_mode` false | 没有分数时 `mixscale_not_run`；没有门槛时 `de_gate_not_run` |
| `DGE/figures/mixscale/mixscale_scatter.png` | 通过 Mixscale 门槛且有真实分数的靶基因，每个一个分面。`nbin=10`。横轴是该靶基因内的扰动分数分位数箱，纵轴是该靶基因自己的 log 归一化表达均值。每个分面一种颜色，无图例，NT 不占分面。纵轴可在精确匹配失败后使用唯一的不区分大小写 `gene_symbol`；DE 名单是否含靶基因仍区分大小写。有限分数全部是 1 的靶基因不画。`DGE/plot_data/mixscale_scatter.tsv` 只有箱均值 | `mixscale_not_run`、`de_gate_not_run`，或 `no_scored_de_targets` |
| `DGE/figures/target_in_de/{靶基因}_mixscape_2x2.png` | 从 `DGE/figures/mixscape/` 复制 symbol 与自己的 Mixscape DE 名单完全一致（区分大小写）的靶基因。同时写 `DGE/plot_data/target_in_de.tsv`。表里的 `logfc` 来自这条 Wilcoxon 名单 | 没有任何靶基因命中时目录里没有新的复制件 |
| `DGE/figures/plot/target_de_counts.png` | 三根蓝柱：靶基因总数、有差异基因的靶基因（不少于 `min_de_genes`）、差异基因名单里含对应靶基因 | 没有 Mixscape 门槛表 |
| `DGE/figures/plot/target_in_de_gene_counts.png` | 差异基因名单包含自身 symbol 的靶基因，柱高是未截断的 `n_de_genes` | 没有这样的靶基因时 skip，原因 `no_targets_in_own_de` |
| `DGE/figures/de/de_gene_heatmap.png` | `effect_main.tsv` 里 `effect_status=descriptive` 的 log2FC，按行做 z-score。行键是 `gene_id`。共用 symbol 的不同基因各占一行；该 symbol 在图里出现多次时，标签是 `symbol (gene_id)`。标题写明相对当前对照的描述性 log2FC | 没有描述性行，或基因少于 2 个 |
| `DGE/figures/de/target_expression_effect.png` | 只保留 `branch=main`。纵轴用表里的 `pseudocount_cpm` 和当前对照（设计好的 NT 或 `zero_guide`）。数据来自 `target_expression_effect.tsv` | 表里没有效应列 |
| `DGE/figures/de/guide_consistency.png`、`pseudobulk_library_size.png`、`pseudobulk_pca.png` | guide 间的靶基因 CPM log2FC，以及 guide 伪批量。PCA 把每条 guide 缩放到 10000 再 `log1p`，标题是 `guide pseudobulk, not replicates`。作图表在 `DGE/plot_data/` | 对应表为空。伪批量可用 guide 少于 3 时 skip，原因 `fewer_than_3_guides` |
| `DGE/plot_data/volcano.skipped.txt` | `report` 直接 skip，原因 `inference_not_available`。2×2 里的火山图随 Mixscape 面板走 | 单样本描述性模式下不画单独火山图 |
| `plot_data/target_effect_forest.skipped.txt` | 同样 skip，原因 `inference_not_available` | 单样本描述性模式下不画 |
| `E-distance/figures/perturbation_distance_heatmap.png`、`perturbation_embedding.png` | 仅在距离表已经写成时由 `report` 绘制。热图标题是 `pairwise energy distance, not the E-test`，MDS 标题是 `MDS of pairwise distances, not the E-test`。这是 guide 伪批量之间的两两能量距离。作图表在 `E-distance/tables/` | 当前配置没有 `perturb.min_guide_replicates`，`perturb` 把原因写在 `plot_data/perturbation_distance_heatmap.skipped.txt`。`report` 读不到距离表时，再在 `E-distance/tables/` 写一份 skip，原因 `missing_fields` |
| `DGE/figures/summary/augur_score_bar.png` | Augur AUC。作图表在 `DGE/plot_data/` | 分数表为空或细胞不足 |
| `DGE/figures/summary/enrichment_dotplot.png` | 上面那份固定基因名单的 Enrichr。横轴是 adjusted P | 没有入选基因时 skip，原因 `no_candidate_genes` |
| `DGE/figures/summary/perturbation_summary.png` | 摘要条形图。显著性计数在摘要表里是 NA，因此这张图不拿它当柱高 | 摘要表为空，或数值列全是 NA |

`report` 不调用 `save_target_log2fc_bar`。函数仍在 `plotting.py` 里，下一次 `report` 不会因为它重画 `{靶基因}_log2fc.png`。
