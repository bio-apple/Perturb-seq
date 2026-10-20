# Perturb-seq 单样本三级分析

从一份 DRAGEN scRNA CRISPR 二级结果出发，完成 guide 分配、细胞过滤、描述性 CPM log2FC、Mixscape / Mixscale 敏感性分析、扰动距离和 Augur。一次运行只接受一个生物学样本，`experiment.biological_sample_count` 在加载配置时写成 1。图只写 PNG，dpi 用 `plotting.dpi`（默认 300）。没有 PDF，也没有 HTML 报告。

各目录的输出契约见 [docs/outputs.md](docs/outputs.md)。

## 运行

环境是 `/opt/miniconda3/envs/perturbseq-tertiary/`（Python 3.12）。入口是 `python -m perturbseq`（`perturbseq.__main__` 调用 `perturbseq.cli.main`）。`pyproject.toml` 里的 console script `perturbseq-tertiary` 指向同一个函数。必须显式传入 `--config`。程序没有默认配置。

```bash
/opt/miniconda3/envs/perturbseq-tertiary/bin/python -m perturbseq --config config/config.yaml --steps all
/opt/miniconda3/envs/perturbseq-tertiary/bin/python -m perturbseq --config config/config.yaml --steps io,assign,preprocess
/opt/miniconda3/envs/perturbseq-tertiary/bin/python -m perturbseq --config config/config.yaml --steps perturb,augur --resume
/opt/miniconda3/envs/perturbseq-tertiary/bin/python -m perturbseq --config config/config.yaml --validate-only
```

`--steps` 默认 `all`。步骤名是 `io`、`assign`、`preprocess`、`perturb`、`augur`、`enrich`、`report`。`--resume` 从上一步的 checkpoint 接着跑。`--validate-only` 只检查参考表和 barcode，不进入分配之后的步骤。

`experiment.sample_id` 必填。`experiment.organism` 必填，只能是 `human` 或 `mouse`。大小写不敏感，加载后写成小写。`organelle`、`细胞器`、`other`、`human,mouse`、空字符串、缺省或其他值都会失败，错误说明只能填 `human` 或 `mouse`。细胞器不是物种。人类线粒体前缀是 `MT-`，小鼠是 `mt-`；配置里写了 `mito_gene_prefix` 就必须与物种一致，没写就按物种填上。

读入基因号和基因符号后做一致性核对，不根据基因改写物种。`human` 对 `ENSG` 和 `MT-`，`mouse` 对 `ENSMUSG` 和 `mt-`。基因号与符号互相矛盾，或与输入物种不一致，就停止，并写出各自的数量和例子。叶绿体等其他细胞器前缀记在 `other_organelle_prefixes`，不计入 percent.mt。核对结果在 `resolved_config.yaml` 和 `run_info.json`：`resolved_organism` 是输入物种，`resolved_mito_prefix` 是采用的线粒体前缀，`species_check` 是 `passed` 或 `failed`，`organism_evidence` 是依据。验证失败时仍会写出 `resolved_config.yaml`，然后停止。

`experiment.crispr_mode` 可以不写，缺省是 `CRISPR_KO`。也可以写 `KO`、`i`、`a`、`other`，分别对应 `CRISPR_KO`、`CRISPRi`、`CRISPRa`、`other`。`other` 仍要填写 `crispr_mode_description` 和 `expected_target_direction`。`labels.nt_control` 可以不写或留空。没有非靶向对照时，QC、分配、预处理、聚类和 UMAP 照常做；描述性 log2FC、Mixscape、Mixscale、对 NT 的 E-distance 跳过，原因是 `no_nt_control`。输入目录里如果还有另一份 `*.positive_cell_guide_assignments.csv`，流程会停止。

`qc.thresholds_confirmed` 必须为 `true`，并且 `min_umi`、`max_umi`、`min_genes`、`max_mito_frac` 已经写成数值。`config/config.yaml` 里的四个数是该样本事先记下的阈值，运行时直接使用这些数。

`statistics.mode` 默认 `descriptive`，配置里也只接受这个值。写成 `guide_exploratory` 时，加载配置失败，错误文本是 `unsupported_inference_backend`。`allow_random_split_replicates: true` 同样被拒绝。

默认测试数据是 `config/config.yaml`：`experiment.sample_id` 为 `A549_2`，`input.input_dir` 为 `/Users/yfan/Desktop/test_data`，物种 `human`，NT 名 `Non-Targeting_Control`。按样本名寻找这五个文件：

- `A549_2.scRNA.filtered.matrix.mtx.gz`
- `A549_2.scRNA.filtered.barcodes.tsv.gz`
- `A549_2.scRNA.filtered.features.tsv.gz`
- `A549_2.scRNA.feature_barcode_reference.csv`
- `A549_2.scRNA.positive_cell_guide_assignments.csv`

仓库里没有 `config/hefei_third.yaml`。合肥样本不是默认数据。要跑别的样本，另写配置并用 `--config` 指向它。

## 分配

`labels.guides_per_cell` 默认是 1。进入下游的细胞必须恰好有 1 条 guide。

分配状态互斥，只有四类：`no_guide`、`single_targeting`、`single_nt`、`multi_guide`。0 条在 `no_guide` 阶段排除。多于 1 条保持 `multi_guide`，在 `multi_guide` 阶段排除。`labels.multi_guide_policy` 只能是 `drop`。配置里如果出现 `labels.flag_guides_above`，加载配置时会拒绝启动。

`input.guide_separator` 必须是已核实的非空字符串。A549_2 是 `|`。guide 名内部也含这个字符时，按参考库最长匹配切分。

过滤阶段按这个顺序记账：`no_guide`、`multi_guide`、`low_guide_umi`、`qc_failed`、`low_cells_target`，剩下的是 `retained`。

`guides_per_cell` 图写到 `qc/figures/guides_per_cell.png`，表在 `qc/tables/guides_per_cell.tsv`。纵轴是 `Cells`。分母是分配账本里的全部细胞，包含 `no_guide` 和 `multi_guide`。柱上同时标整数计数和一位小数百分比。这张图不画 guide 数大于 5 的红线。

## Barcode

`input.barcode_policy` 必须是 `strict`。其它取值在加载配置时会被拒绝。矩阵 barcode 集合和分配 barcode 集合只要有差异，`io` 和 `--validate-only` 就停止，并写出 `barcode_differences.tsv` 和 `input_validation.json`。集合一致时差异表只有表头。

A549_2 的 strict 校验记录是仅矩阵 11397 个 barcode、仅分配 0 个，`status` 为 `failed`。因此用 `config/config.yaml` 做 `strict` 时，全程跑不完。工作区里现在没有 `results/`。

## 统计

`statistics.mode` 是 `descriptive`。guide 只用于描述 guide 之间的效应范围。主效应按靶基因把原始计数加总，再换成 CPM。每个基因一行：

`log2((CPM_target + pseudocount_cpm) / (CPM_NT + pseudocount_cpm))`

`statistics.pseudocount_cpm` 从配置读取，`config/config.yaml` 里是 `1.0`。靶基因和 NT 的原始计数都是 0 时，该基因的 `log2fc_descriptive` 为缺失，`effect_status` 为 `not_detected_both`。

`effect_main.tsv` 和 `effect_sensitivity.tsv` 写在输出目录根下。`branch=main` 用全部可分析细胞。`branch=sensitivity` 只保留 Mixscape 标成 `perturbed` 的靶向细胞，再加上同一批 NT。这两个分支名写死在 `perturb` 步骤里。两张表的 `pvalue`、`fdr`、`ci_lower`、`ci_upper` 为 NA，`inference_status` 是 `not_run_single_sample_descriptive`。细胞水平 Wilcoxon 的 `adj_p_value`、Mixscale 分数、E-test 的 p 分别留在 `DGE/plot_data/`、`obs['mixscale_score']` 和 `E-distance/tables/edistance_etest.tsv`。`report` 把 `qc_summary.tsv` 和 `perturbation_summary.tsv` 的 `n_significant_de_genes` 写成 NA。

单样本时 edgeR 和 PyDESeq2 都走 pertpy 差异表达教程的伪批量模型，状态是 `skipped`，原因是 `single_sample_de_not_estimable`。`DGE/tables/de_edger_main.tsv`、`de_edger_sensitivity.tsv`、`de_pydeseq2_main.tsv`、`de_pydeseq2_sensitivity.tsv` 只有表头。描述性 CPM log2FC 只在效应表的 `log2fc_descriptive`。

Enrichr 用 `|log2fc_descriptive| ≥ plotting.volcano_lfc` 的固定基因名单（默认 1.0）。返回的 adjusted P 记在 `enrichment_fdr`，是这个固定名单上的基因集统计。没有任何靶基因达到 `enrichment.min_genes` 时，`tables/enrichment_status.txt` 写 `no_effect_ranked_genes`。

没有细胞类型注释时，`obs['cell_type']` 写成 `ALL`，`cell_type_scope` 为 `unstratified`。效应表的 `cell_type` 列也是 `ALL`。Mixscape 在 `cell_type_key` 为空时不按细胞类型拆分，日志会写 `Mixscape scope: ALL cells, not stratified`。

## Mixscape 与 Mixscale

Mixscape 是敏感性分支。分类标签是 `nt`、`perturbed`、`non_perturbed`、`ambiguous`、`not_evaluated`。没过 Mixscape 差异基因门槛的靶基因保持 `skipped_no_de_genes`，不改成 NP。`mixscape_fraction` 的分母是该靶基因的全部靶向细胞；组为空时比例写成 0.0。全局类柱状图和 `umap_mixscape_class_global` 只保留已分类细胞，去掉 `not_evaluated` 和 `skipped_no_de_genes`。`X_umap_mixscape` 对未进入签名的细胞是缺失值；`draw_umap` 看到任一缺失坐标就跳过 `umap_mixscape`，原因 `non_finite_coordinates`。

两套门槛分开读，互不借用：

- Mixscape 对齐 Seurat Mixscape vignette。`mixscape.n_dims` 默认 40，矩阵比 40 小才截断；`n_control_neighbors` 默认 20，作为 `perturbation_signature` 的 `n_neighbors`，NT 细胞更少时才减小。`logfc_threshold` 0.25，`min_de_genes` 5，`pval_cutoff` 0.05，检验是 Wilcoxon，`iter_num` 10，`fine_mode` 固定为 false。热图 `mixscape_heatmap_max_genes` 20、`mixscape_heatmap_max_cells` 300：NP、NT、KO 各组各自最多抽这么多细胞，人数不相等时保持原样。热图在 2×2 右下角，列序 NP、NT、KO，组间空白列宽是 `max(6, round(0.08 × 最大组细胞数))`。扰动分数图里 KO 的颜色是 `#EE6A50`。2×2 火山图的竖线是 Mixscape 的 `logfc_threshold` 0.25。`plotting.volcano_lfc`（默认 1.0）用于 Enrichr 名单和单独火山图的跳过逻辑，不画到这张 2×2 上。只有一个生物学样本，调用里不传 replicate 拆分。配置里的 `mixscape.ambiguous_enabled`、`ambiguous_lower`、`ambiguous_upper`，以及 `statistics.main_branch`、`statistics.sensitivity_branch`，当前代码不读取。后验分类沿用软件自己的 0.5，并记在 `uns['mixscape_parameters']`。
- Mixscale 对齐 Mixscale 教程。`min_de_genes` 5，`logfc_threshold` 0.2，`max_de_genes` 100，`fine_mode` false，`pval_cutoff` 0.05。不读取 Mixscape 的 0.25，也不把 `CalcPerturbSig` 的 `n_dims` 套到 Mixscale。

只有通过对应差异基因门槛的靶基因才进入 Mixscape 2×2 和 Mixscale 山脊图 / 散点图。2×2 还要求该靶基因同时有 NP、NT 和 KO 细胞以及扰动分数。软件把未评分的 perturbation 填成 1 时，`run_mixscale` 在靶基因不在已评分集合、或投影到的差异基因少于 `min_de_genes` 时把分数写成缺失。散点图再丢掉有限分数全部是 1 的靶基因。

Mixscale 散点图 `nbin=10`。横轴是该靶基因内部扰动分数的分位数箱，纵轴是这个靶基因自己在 `adata.X` 里的 log 归一化表达均值。每个分面一种颜色，无图例。NT 不占分面。通过 Mixscale 门槛且有真实分数的靶基因全部画出。纵轴先按 `gene_symbol` 或 var 名精确匹配；两边都没有时，恰好一个不区分大小写的 `gene_symbol` 可以用。差异基因名单里是否包含靶基因仍然区分大小写。

## 预处理与簇注释

预处理从 `layers['counts']` 重新归一化到 `preprocess.normalize_target_sum`（10000），再 `log1p`。高变基因用 Seurat flavor，`n_hvg` 2000。线粒体基因从高变集合里拿掉，用后面的非线粒体基因补上。人类用 `MT-`，小鼠用 `mt-`。配置里的前缀和物种不一致时，加载配置失败。叶绿体等其他细胞器基因不计入 percent.mt，也不按线粒体从高变基因里排除。PCA 做在标准化后的高变基因上，`n_pcs` 30，`X` 保持 log 归一化。

`preprocess.batch_key` 为空时不调用 Harmony。A549_2 的 `batch_key` 是空。只有配置了 `batch_key`、而且有效水平少于 2 时，`uns['batch_correction']` 才写成 `not_run_fewer_than_two_batch_levels`。邻居数 `n_neighbors` 15，Leiden 用 igraph，`leiden_resolution` 1.0。

细胞周期只在人类、并且 `uns['log1p']` 已经存在时打分。打分不改 `X`，也不做回归。小鼠的状态是 `skipped_not_human`。

`umap_cluster_annotation` 是 Leiden 和状态 marker 的相对 z-score。状态是 S、G2M、Epithelial、Stress、Interferon、Hypoxia。最高分 `z ≥ 0.5` 且与第二名差距 `≥ 0.25` 才赋值，否则 `unassigned`。标签写在图和 `UMAP/tables/umap_cluster_annotation.tsv` 里。`obs['cell_type']` 仍是 `ALL` 或外部注释，簇状态不写进这一列。

## CRISPR 模式

`experiment.crispr_mode` 缺省时是 `CRISPR_KO`，`resolved_expected_direction` 为 `not_required`。`KO` 与 `CRISPR_KO` 相同。`i` 与 `CRISPRi` 解析为 `down`，`a` 与 `CRISPRa` 解析为 `up`。`other` 仍是独立模式，必须同时有描述和方向（`down`、`up`、`bidirectional` 或 `unknown`），不会当成 KO。`config/config.yaml` 已经把 A549_2 写成 `CRISPR_KO`，文件头说明这是按 Edit-R NTC 和 `gene|design_N` 命名记下的假设，运行时保持这个值。矩阵和分配文件本身没有被代码用来确认这是敲除、CRISPRi 还是 CRISPRa。

`target_expression_effect.tsv` 记录靶基因自身 RNA 的描述性 CPM log2FC。基因不在矩阵里（`gene_not_in_rna`）、NT 检出比例或平均 CPM 低于 `target_expression` 的门槛（`nt_expression_below_threshold`）、或 NT / 该靶基因细胞数为 0（`insufficient_cells`）时，`target_expression_status` 是 `not_evaluable`。

## 图与输出

`plotting.formats` 只有 `png`。分桶由 `FigureStore` 按图名决定：

- `qc/figures` 与 `qc/tables`：mRNA QC 和 gRNA 统计，包括 `guides_per_cell`、分配柱、每个 guide / 靶基因的细胞数、过滤瀑布、高变基因和 PCA。
- `DGE/figures/<section>`、`DGE/tables`、`DGE/plot_data`：差异基因门槛、Mixscape 2×2、Mixscale 山脊图和散点图、靶基因 DE 计数、靶基因出现在自己 DE 名单里的计数、描述性 log2FC 热图和靶基因效应图、guide 伪批量图、Augur 和 Enrichr。
- `UMAP/figures` 与 `UMAP/tables`：文件名以 `umap_` 开头的图，包括簇注释、细胞周期 UMAP 和 Mixscape signature UMAP。
- `E-distance/figures` 与 `E-distance/tables`：每个通过 Mixscape 差异基因门槛的靶基因对 `labels.nt_control` 的 pertpy `DistanceTest("edistance")`，用对象里已有的 `X_pca`，置换单位是细胞，校正是 holm-sidak。guide 伪批量的两两距离热图和 MDS 也进这个目录；标题写明那是 guide 伪批量距离。当前配置没有 `perturb.min_guide_replicates`，热图因此跳过，`plot_data/perturbation_distance_heatmap.skipped.txt` 的内容是 `min_guide_replicates_not_set`。E-test 的 p 不写入效应表。

描述性效应表在输出根目录。主效应图 `target_expression_effect` 只画 `branch=main`，纵轴是 `log2((CPM_target+1)/(CPM_NT+1))`。guide 伪批量 PCA 把每条 guide 缩放到 10000 再 `log1p`，标题是 `guide pseudobulk, not replicates`；guide 少于 3 条时跳过，原因 `fewer_than_3_guides`。描述性热图的行是 `gene_id`。共用一个 symbol 的不同基因保持各行；这个 symbol 在图中出现多次时，标签写成 `symbol (gene_id)`。

已有的 54 张 `{靶基因}_mixscape_2x2.png` 来自旧运行。ACTB 那张图的脚注写着 `lfc=1.0`。当前代码传给 2×2 的是 Mixscape `logfc_threshold` 0.25，签名宽度默认 `n_dims` 40。旧的 `mixscale_scatter.txt` 把 Eif6 和 Plk5 记在 `omitted_missing_gene`，文件里也还没有当前代码会写的 `y_axis_symbol_rule` 行。工作区的 `results/` 里现在没有这些 PNG。

## 测试

```bash
/opt/miniconda3/envs/perturbseq-tertiary/bin/python -m pytest
```

`tests/test_a549_assignments.py` 在 `/Users/yfan/Desktop/test_data` 存在时核对 A549_2 的 guide 参考表和分配计数。

