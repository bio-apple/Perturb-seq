# 分析规则

这里是命令行、配置和各步骤的判定。写出哪些文件、图何时跳过，见 [outputs.md](outputs.md)。

## 命令行与配置

`config/config.yaml` 只放分析软件的默认参数，例如 Mixscape、Mixscale、预处理、Augur 和作图。样本名、物种、输入目录、输出目录、guide 分隔符、非靶向对照和 CRISPR 模式在命令行输入。不传 `--config` 时读取这份默认参数文件。命令行会覆盖文件里的同名字段。

加载配置时，映射里出现分析不读取的键会失败。错误文本是 `未识别的配置键` 加上键路径，多个键用顿号连接。`labels.flag_guides_above` 仍在允许列表里；一旦写上，另外以「已停用，不能作为剔除阈值」拒绝启动。

最小命令见仓库根目录的 README。`--crispr-mode` 缺省是 `CRISPR_KO`，也可以写 `i`、`a` 或 `other`。`--steps` 默认 `all`。步骤名是 `io`、`assign`、`preprocess`、`perturb`、`augur`、`enrich`、`report`。`--resume` 从上一步的 checkpoint 接着跑。`--validate-only` 只检查参考表和 barcode，不进入分配之后的步骤。

程序没有内置样本。五个 DRAGEN 文件由 `--sample-id` 和 `--input-dir` 拼出来：

- `{sample_id}.scRNA.filtered.matrix.mtx.gz`
- `{sample_id}.scRNA.filtered.barcodes.tsv.gz`
- `{sample_id}.scRNA.filtered.features.tsv.gz`
- `{sample_id}.scRNA.feature_barcode_reference.csv`
- `{sample_id}.scRNA.positive_cell_guide_assignments.csv`

`--batch-key` 和 `--cell-type-key` 不填时分别为空。空的批次列不做 Harmony。空的细胞类型列使 `cell_type` 保持 `ALL`。

## 物种

`--organism` 必填，只能是 `human` 或 `mouse`。大小写不敏感，加载后写成小写。`organelle`、`细胞器`、`other`、`human,mouse`、空字符串或其他值都会失败，错误说明只能填 `human` 或 `mouse`。细胞器不是物种。人类线粒体前缀是 `MT-`，小鼠是 `mt-`，由物种决定。

读入基因号和基因符号后做一致性核对，不根据基因改写物种。`human` 对 `ENSG` 和 `MT-`，`mouse` 对 `ENSMUSG` 和 `mt-`。基因号与符号互相矛盾，或与输入物种不一致，就停止，并写出各自的数量和例子。叶绿体等其他细胞器前缀记在 `other_organelle_prefixes`，不计入 percent.mt。核对结果在 `resolved_config.yaml` 和 `run_info.json`：`resolved_organism` 是输入物种，`resolved_mito_prefix` 是采用的线粒体前缀，`species_check` 是 `passed` 或 `failed`，`organism_evidence` 是依据。验证失败时仍会写出 `resolved_config.yaml`，然后停止。

## CRISPR 模式

`--crispr-mode` 可以不写，缺省是 `CRISPR_KO`。也可以写 `KO`、`i`、`a`、`other`，分别对应 `CRISPR_KO`、`CRISPRi`、`CRISPRa`、`other`。`CRISPRi` 的预期方向是 `down`，`CRISPRa` 是 `up`，`CRISPR_KO` 是 `not_required`。`other` 必须同时有 `--crispr-description` 和 `--expected-direction`（`down`、`up`、`bidirectional` 或 `unknown`）。矩阵和分配文件不会被用来改写这个模式。

`resolved_expected_direction` 由上面的 CRISPR 模式决定。代码不根据 guide 名字猜测敲除、干扰或激活。

## 对照

`--nt-control` 写该样本参考表里的非靶向 guide 名。可以不写。写了并且数据里真有这些细胞时，对照就是这组设计好的非靶向细胞，zero-guide 细胞不并进去。没写时，分配表里 0 条 guide 的细胞，以及仅在矩阵里、记为 `no_guide` 的细胞，是临时对照：`control_source=zero_guide_cells`，标签是 `zero_guide`，不会写成某个 Non-Targeting 名字。有了这组临时对照之后，描述性 CPM log2FC、Mixscape、Mixscale、对对照的 E-distance 都会做。`multi_guide` 和 `single_targeting` 不是对照。输入目录里如果还有另一份 `*.positive_cell_guide_assignments.csv`，流程会停止。

临时对照在 `run` 一开始写入运行日志、`tables/temporary_control_warning.txt` 和 `run_info.json` 的 `temporary_control_warning`。全文是：当前对照是临时 zero-guide cells，这不是设计好的非靶向对照。zero-guide 细胞可能在 UMI、基因数、线粒体比例、环境 RNA 和细胞状态上与靶向细胞不同。后续 log2FC、Mixscape、Mixscale、E-distance 都是相对这组临时对照。UMI、基因数和线粒体比例的对比只作描述，分布差异不会中止分析。`report` 用通过 QC 的 zero-guide 细胞（`assignment_status=no_guide`、`n_guides=0`、`qc_pass`）和 `perturbation_status=targeting` 的可分析细胞画 `qc/figures/temporary_control_qc.png`：同一张图比较 UMI（`total_counts`）、基因数（`n_genes_by_counts`）和 `percent.mt`（0–100）。`qc/tables/temporary_control_qc.tsv` 给出每组细胞数、中位数、MAD 和四分位数，组名是 `zero-guide cells（临时对照）` 和 `targeting`。缺少 QC 列或某一组没有细胞时，跳过原因写在 `qc/tables/temporary_control_qc.skipped.txt`。设计好的 NT 不写这份警告，并删掉已有的警告文件和 `qc/figures/temporary_control_qc.png`。

`low_library_no_guide` 的细胞不进入这组临时对照，也不进入设计好的 NT。真正检出 NT guide 的 `single_nt` 细胞仍是设计对照。已有 guide 的细胞不会因为捕获效率低被改成 NT。

## QC 阈值

四个 QC 绝对阈值不在默认参数文件里。命令行同时给出 `--min-umi`、`--max-umi`、`--min-genes`、`--max-mito-frac` 时使用这四个数。只给其中一部分会失败。全部留空时，在分配步骤用全部细胞计算，并写进 `tables/qc_thresholds.json`，`source` 为 `median_mad`。UMI 和基因数先做 `log1p`，界是 median ± `qc.n_mad`（默认 3）× 1.4826 × MAD，再 `expm1` 回到原尺度；线粒体比例在原分数上用同一公式，并限制在 0 到 1。MAD 为 0 时改用标准差。命令行给齐四个数时 `source` 为 `command_line`。

核糖体蛋白基因比例是另一项上限，不是第五个必须和上面四个一起给出的绝对阈值。按基因符号前缀计，人类 `RPL`/`RPS`，小鼠 `Rpl`/`Rps`，大小写不敏感。这是核糖体蛋白基因，不是线粒体，也不使用 rRNA 基因集。`MRPL`/`Mrpl` 不以这些前缀开头，不计入。比例是 0–1 的分数，图上的 `percent.ribo` 显示 0–100。未单独给 `--max-ribo-frac` 时，上限是全部细胞的 median + `qc.n_mad` × 1.4826 × MAD，`tables/qc_thresholds.json` 的 `ribo_frac_source` 为 `median_mad`。单独给出 `--max-ribo-frac`（0 到 1）时用这个数，`ribo_frac_source` 为 `command_line`。超过上限的细胞在 QC 阶段排除，原因记为 `high_ribo`，主过滤阶段仍是 `qc_failed`。

gRNA 捕获效率是该细胞 CRISPR guide UMI 合计除以 Gene Expression 的 `total_counts`。guide UMI 合计为 0 时效率为 0。文库大小仍保留为 `total_counts`。相关在 QC 通过前的全部细胞上计算 Spearman，系数、细胞数和 p 值写在 `qc/tables/grna_capture_library_correlation.tsv`。这个 p 不写入 `effect_main.tsv`。

这条检查用来避免文库偏浅、没捕获到 guide 的细胞被当作阴性对照。下限来自已经检测到至少 1 条 guide 的细胞：对它们的 `log1p(total_counts)` 取 median − `qc.n_mad` × 1.4826 × MAD，再 `expm1`，与 UMI 阈值同一套 MAD。没检测到 guide、且 `total_counts` 低于该下限的细胞，`qc_pass` 为 false，主过滤阶段是 `low_library_no_guide`。它们留在分配账本里，分配状态仍是 `no_guide`，但 `control_cell_mask` 不包括它们。没有细胞检测到 guide，或这些细胞的文库大小方差为 0 时，跳过这条排除，不让流程崩溃，也不把全部 `no_guide` 改成设计好的 NT。`tables/qc_thresholds.json` 的 `low_library_no_guide.reason` 分别是 `no_cells_with_guide` 和 `library_size_variance_zero`。图是 `qc/figures/grna_capture_vs_library_size.png`，标题写明捕获效率与文库大小的关系，并标出没检测到 guide 的细胞。

## 分配

`labels.guides_per_cell` 默认是 1。进入下游的细胞必须恰好有 1 条 guide。

分配状态互斥，只有四类：`no_guide`、`single_targeting`、`single_nt`、`multi_guide`。0 条在 `no_guide` 阶段排除。多于 1 条保持 `multi_guide`，在 `multi_guide` 阶段排除。`labels.multi_guide_policy` 只能是 `drop`。配置里如果出现 `labels.flag_guides_above`，加载配置时会拒绝启动。

`--guide-separator` 必须是已核实的非空字符串。guide 名内部也含这个字符时，按参考库最长匹配切分。

过滤阶段按这个顺序记账：`low_library_no_guide`、`no_guide`、`multi_guide`、`low_guide_umi`、`qc_failed`、`low_cells_target`，剩下的是 `retained`。规则跳过时 `low_library_no_guide` 阶段新增排除数为 0。

`guides_per_cell` 图写到 `qc/figures/guides_per_cell.png`，表在 `qc/tables/guides_per_cell.tsv`。纵轴是 `Cells`。分母是分配账本里的全部细胞，包含 `no_guide` 和 `multi_guide`。柱上同时标整数计数和一位小数百分比。这张图不画 guide 数大于 5 的红线。

## Barcode

`input.barcode_policy` 必须是 `strict`。其它取值在加载配置时会被拒绝，也不会改成两边求交。规则不看文件名或样本名：

- 分配表有、矩阵没有（assignment-only）：`io` 和 `--validate-only` 失败，并写出 `barcode_differences.tsv` 和 `input_validation.json`。
- 矩阵有、分配没有（matrix-only），且 assignment-only 为 0：这些细胞记为 `no_guide` 并保留，流程不中止。没有设计 NT 时，其中通过 QC、且不是 `low_library_no_guide` 的细胞，和分配表里同样合格的 0 条 guide 细胞一起作为 `zero_guide` 临时对照。
- 两侧都有差异：仍然失败。
- 集合一致时差异表只有表头。

## 统计

`statistics.mode` 默认 `descriptive`，配置里也只接受这个值。写成 `guide_exploratory` 时，加载配置失败，错误文本是 `unsupported_inference_backend`。`allow_random_split_replicates: true` 同样被拒绝。

`statistics.mode` 是 `descriptive`。guide 只用于描述 guide 之间的效应范围。主效应按靶基因把原始计数加总，再换成 CPM。每个基因一行：

`log2((CPM_target + pseudocount_cpm) / (CPM_control + pseudocount_cpm))`

对照是设计好的 NT，或没有设计 NT 时的 `zero_guide` 细胞。`statistics.pseudocount_cpm` 从配置读取。靶基因和对照的原始计数都是 0 时，该基因的 `log2fc_descriptive` 为缺失，`effect_status` 为 `not_detected_both`。

`effect_main.tsv` 和 `effect_sensitivity.tsv` 写在输出目录根下。`branch=main` 用全部可分析细胞。`branch=sensitivity` 只保留 Mixscape 标成 `perturbed` 的靶向细胞，再加上同一批对照细胞。这两个分支名写在 `perturb` 步骤里。细胞水平 Wilcoxon 的 `adj_p_value`、Mixscale 分数、E-test 的 p 分别留在 `DGE/plot_data/`、`obs['mixscale_score']` 和 `E-distance/tables/edistance_etest.tsv`。`report` 把 `qc_summary.tsv` 和 `perturbation_summary.tsv` 的 `n_significant_de_genes` 写成 NA。

有描述性 log2FC 的靶基因同时写出样本内 bootstrap 百分位区间。`statistics.bootstrap_n` 是正整数，默认 200。`statistics.bootstrap_ci` 是 (0, 1) 内的数，默认 0.95。随机种子用 `runtime.seed`。靶基因在可分析细胞里有至少 2 条不同 guide 时，重抽样单位是 guide：有放回抽取，每次抽到的 guide 数等于观察到的 guide 数，把抽到的 guide 原始计数合并后再算同一条 CPM log2FC。只有 1 条 guide 时，重抽样单位改为该靶基因的细胞，有放回抽取，次数等于细胞数。对照按自己的 guide 数决定单位，并且在每一次抽样里和靶基因一起重抽。设计好的 NT 有至少 2 条 guide 时按 guide 重抽样，只有 1 条 guide 时按细胞重抽样。`zero_guide` 临时对照没有 guide，按细胞重抽样。

新列是 `log2fc_bootstrap_ci_lower`、`log2fc_bootstrap_ci_upper`、`bootstrap_unit`、`bootstrap_n`、`bootstrap_ci_level`、`uncertainty_meaning`、`bootstrap_status`。`bootstrap_unit` 记该靶基因的单位，取值 `guide` 或 `cell`。`uncertainty_meaning` 固定为 `descriptive_estimate_uncertainty_not_significance_test`。点估计缺失时 `bootstrap_status` 为 `point_estimate_missing`。靶基因、对照的抽样单位都少于 2 个时为 `fewer_than_2_target_and_control_sampling_units`；只有靶基因一侧少于 2 个时为 `fewer_than_2_target_sampling_units`；只有对照一侧少于 2 个时为 `fewer_than_2_control_sampling_units`。百分位算不出来时为 `bootstrap_interval_undefined`，算出时为 `computed`。这些情况下区间两列保持缺失，不会写成 0。`pvalue`、`fdr`、`ci_lower`、`ci_upper` 保持为 NA，`inference_status` 保持 `not_run_single_sample_descriptive`。敏感性分支整支跳过时，`effect_sensitivity.tsv` 只有表头。

这些区间是单个样本内部重抽样得到的描述性不确定性。它们不能当作显著性检验，也不能当作生物学重复的置信区间。

单样本时 edgeR 和 PyDESeq2 都走 pertpy 差异表达教程的伪批量模型，状态是 `skipped`，原因是 `single_sample_de_not_estimable`。`DGE/tables/de_edger_main.tsv`、`de_edger_sensitivity.tsv`、`de_pydeseq2_main.tsv`、`de_pydeseq2_sensitivity.tsv` 只有表头。描述性 CPM log2FC 只在效应表的 `log2fc_descriptive`。

Enrichr 用 `|log2fc_descriptive| ≥ plotting.volcano_lfc` 的固定基因名单（默认 1.0）。`gseapy` 把名单提交到 Enrichr 服务器，这一步需要联网。返回的 adjusted P 记在 `enrichment_fdr`，是这个固定名单上的基因集统计。没有任何靶基因达到 `enrichment.min_genes` 时，`tables/enrichment_status.txt` 写 `no_effect_ranked_genes`。连不上服务器时原因是 `enrichr_failed`。

没有细胞类型注释时，`obs['cell_type']` 写成 `ALL`，`cell_type_scope` 为 `unstratified`。效应表的 `cell_type` 列也是 `ALL`。Mixscape 在 `cell_type_key` 为空时不按细胞类型拆分，日志会写 `Mixscape scope: ALL cells, not stratified`。

## Mixscape 与 Mixscale

Mixscape 是敏感性分支。分类标签是 `nt`、`perturbed`、`non_perturbed`、`ambiguous`、`not_evaluated`。没过 Mixscape 差异基因门槛的靶基因保持 `skipped_no_de_genes`，不改成 NP。`mixscape_fraction` 的分母是该靶基因的全部靶向细胞；组为空时比例写成 0.0。全局类柱状图和 `umap_mixscape_class_global` 只保留已分类细胞，去掉 `not_evaluated` 和 `skipped_no_de_genes`。`X_umap_mixscape` 对未进入签名的细胞是缺失值；`draw_umap` 看到任一缺失坐标就跳过 `umap_mixscape`，原因 `non_finite_coordinates`。

两套门槛分开读，互不借用：

- Mixscape 对齐 Seurat Mixscape vignette。`mixscape.n_dims` 默认 40，矩阵比 40 小才截断；`n_control_neighbors` 默认 20，作为 `perturbation_signature` 的 `n_neighbors`，NT 细胞更少时才减小。`logfc_threshold` 0.25，`min_de_genes` 5，`pval_cutoff` 0.05，检验是 Wilcoxon，`iter_num` 10，`fine_mode` 固定为 false。热图 `mixscape_heatmap_max_genes` 20、`mixscape_heatmap_max_cells` 300：NP、NT、KO 各组各自最多抽这么多细胞，人数不相等时保持原样。热图在 2×2 右下角，列序 NP、NT、KO，组间空白列宽是 `max(6, round(0.08 × 最大组细胞数))`。扰动分数图里 KO 的颜色是 `#EE6A50`。2×2 火山图的竖线是 Mixscape 的 `logfc_threshold` 0.25。`plotting.volcano_lfc`（默认 1.0）用于 Enrichr 名单和单独火山图的跳过逻辑，不画到这张 2×2 上。只有一个生物学样本，调用里不传 replicate 拆分。配置里出现未识别的键会报错并停止。`mixscape.ambiguous_enabled`、`mixscape.ambiguous_lower`、`mixscape.ambiguous_upper`、`statistics.main_branch`、`statistics.sensitivity_branch` 已从默认配置删除。后验分类沿用软件自己的 0.5，并记在 `uns['mixscape_parameters']`。
- Mixscale 对齐 Mixscale 教程。`min_de_genes` 5，`logfc_threshold` 0.2，`max_de_genes` 100，`fine_mode` false，`pval_cutoff` 0.05。不读取 Mixscape 的 0.25，也不把 `CalcPerturbSig` 的 `n_dims` 套到 Mixscale。

只有通过对应差异基因门槛的靶基因才进入 Mixscape 2×2 和 Mixscale 山脊图 / 散点图。2×2 还要求该靶基因同时有 NP、NT 和 KO 细胞以及扰动分数。软件把未评分的 perturbation 填成 1 时，`run_mixscale` 在靶基因不在已评分集合、或投影到的差异基因少于 `min_de_genes` 时把分数写成缺失。散点图再丢掉有限分数全部是 1 的靶基因。

Mixscale 散点图 `nbin=10`。横轴是该靶基因内部扰动分数的分位数箱，纵轴是这个靶基因自己在 `adata.X` 里的 log 归一化表达均值。每个分面一种颜色，无图例。NT 不占分面。通过 Mixscale 门槛且有真实分数的靶基因全部画出。纵轴先按 `gene_symbol` 或 var 名精确匹配；两边都没有时，恰好一个不区分大小写的 `gene_symbol` 可以用。差异基因名单里是否包含靶基因仍然区分大小写。

## 预处理与簇注释

预处理从 `layers['counts']` 重新归一化到 `preprocess.normalize_target_sum`（10000），再 `log1p`。高变基因用 Seurat flavor，`n_hvg` 2000。线粒体基因从高变集合里拿掉，用后面的非线粒体基因补上。线粒体前缀由 `--organism` 决定，人类是 `MT-`，小鼠是 `mt-`。叶绿体等其他细胞器基因不计入 percent.mt，也不按线粒体从高变基因里排除。PCA 做在标准化后的高变基因上，`n_pcs` 30，`X` 保持 log 归一化。

`--batch-key` 为空时不调用 Harmony。只有传入了批次列、而且有效水平少于 2 时，`uns['batch_correction']` 才写成 `not_run_fewer_than_two_batch_levels`。邻居数和 Leiden 分辨率从 `config/config.yaml` 读取，Leiden 用 igraph。

细胞周期在预处理后、已有 `uns['log1p']` 时打分。基因集来自 [scib `get_cell_cycle_genes` / `score_cell_cycle`](https://scib.readthedocs.io/en/latest/api/scib.preprocessing.score_cell_cycle.html)（Tirosh / Regev lab 列表）：人类为大写符号，小鼠为首字母大写符号。用 `scanpy.tl.score_genes_cell_cycle` 写入 `phase`、`S_score`、`G2M_score`，`use_raw=False`，不改 `X`，也不做回归。`human` 与 `mouse` 都打分；缺基因时 `uns['cell_cycle_status']` 为 `skipped_missing_cell_cycle_genes`，未装 scib 时为 `skipped_scib_unavailable`，成功时为 `scored_after_log1p`。

`umap_cluster_annotation` 用 SingleR 给每个 Leiden 簇一个标签。先对簇内 log 归一化表达取平均，再由 `Rscript` 跑 SingleR。小鼠参考是 celldex 的 `MouseRNAseqData`，人类参考是 `HumanPrimaryCellAtlasData`，标签字段是 `label.main`。参考必须和 `experiment.organism` 是同一个物种。标签写入 `obs['singler_label']`、图和 `UMAP/tables/umap_cluster_annotation.tsv`。`obs['cell_type']` 仍是 `ALL` 或外部注释，不改成 SingleR 标签。`singler.enabled` 为 `false` 时跳过这一步，原因 `singler_disabled`。

## 靶基因表达

`target_expression_effect.tsv` 记录靶基因自身 RNA 的描述性 CPM log2FC。基因不在矩阵里（`gene_not_in_rna`）、NT 检出比例或平均 CPM 低于 `target_expression` 的门槛（`nt_expression_below_threshold`）、或 NT / 该靶基因细胞数为 0（`insufficient_cells`）时，`target_expression_status` 是 `not_evaluable`。

图的分桶、跳过原因和各表的列见 [outputs.md](outputs.md)。
