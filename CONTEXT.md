# Perturb-seq 三级分析

针对一个 Sample，在 DRAGEN 二级分析之后描述基因表达结构与 guide 构成，并用 Mixscape 与 E-distance 描述扰动效应。基因层面的差异和 E-test 都是这一个文库内部的细胞级证据。

## 样本与输入

**Sample**:
一个 DRAGEN 输出前缀及其五份二级输出，对应一个文库。没有生物学重复轴。
_Avoid_: 批次, lane, study

**标准测试 Sample**:
前缀为 `A549_2` 的那一个 Sample。参考表同时有 `target_gene_id` 和 `target_gene_name`。本分析用它核对结果。
_Avoid_: sample1

**Guide call**:
DRAGEN 用 GMM 判定一个 passing cell 表达哪些 CRISPR feature，记在 positive cell guide assignments 里。三级分析沿用这一判定。条数看 `num_features`。
_Avoid_: 重新归属, guide assignment, 按 `|` 的个数数 guide

**Guide**:
参考表里的一条 CRISPR feature。`feature_call` 先对 `id`，对不上再对 `name`。标准测试 Sample 里，单 guide 细胞的 `feature_call` 全部等于 `name`，而 `name` 内部可以含 `|`。多 guide 细胞的 `feature_call` 是若干个 `name` 用 `|` 接起来。
_Avoid_: 按每个 `|` 切开

**人类**:
本分析的物种。线粒体基因前缀为 `MT-`，核糖体基因前缀为 `RPS` 或 `RPL`。`target_gene_name` 保持参考表里的大小写。
_Avoid_: 小鼠

**敲除**:
本分析解释 KO / NP 时所假定的实验类型：CRISPR 敲除。
_Avoid_: CRISPRi, CRISPRa, 敲低

## 细胞集

**Passing cell**:
出现在 guide assignment 里的细胞。Guide 频率在全部 passing cell 上计算，包括零条 guide 和多条 guide。只在表达矩阵里、没有 assignment 行的 barcode 不是 Passing cell。
_Avoid_: 单 guide 细胞, filtered 矩阵里的全部 barcode

**分析细胞**:
Passing cell 里能在表达矩阵中找到的细胞。结果里的细胞表装的是这群细胞。Cluster、Mixscape 和扰动标签只标在其中的单 guide 细胞上。
_Avoid_: filtered 矩阵里的全部 barcode

**零 guide 细胞**:
Guide call 里一条 guide 都没有的 passing cell（`num_features` 为 0）。只进入 guide 频率。
_Avoid_: 非靶向对照

**多 guide 细胞**:
Guide call 里叫到两条及以上 guide 的 passing cell。只进入 guide 频率。
_Avoid_: 扰动

**单 guide 细胞**:
Guide call 中恰好叫到一条 guide 的细胞（`num_features` 为 1）。基因表达 QC、Cluster、Marker、组成、分类、Mixscape、E-distance 和探索性差异表达都在这群细胞上。
_Avoid_: 单细胞, singlet

## 分类

**非靶向对照**:
`feature_call` 命中非靶向名字的单 guide 细胞。命中不区分大小写，子串为 `NegCtrl`、`Non-Targeting`、`NonTargeting`、`NTC`。Mixscape、E-distance 和探索性差异表达的对照都是这一类。一个 Sample 里没有这一类时，三者都没有对照可比。
_Avoid_: NP, 内含子对照, 零 guide 细胞

**内含子对照**:
`feature_call` 命中 `IntronCtrl` 的单 guide 细胞。留在基因表达的 Cluster 里，不进入对照，也不进入扰动。
_Avoid_: 非靶向对照, 扰动

**受扰动细胞**:
单 guide 细胞里，既不是非靶向对照也不是内含子对照的细胞。Mixscape 再把它们分成 KO 和 NP。
_Avoid_: KO

**KO**:
Mixscape 对受扰动细胞的标签，表示扰动签名落在已扰动分量。按敲除来读。
_Avoid_: 受扰动细胞, 敲除在实验上已经成功

**NP**:
Mixscape 对受扰动细胞的标签，表示扰动签名落在未扰动分量。按敲除来读时，它是「转录组上看不到敲除」。
_Avoid_: 非靶向对照

## 表达结构

**基因表达**:
Sample 矩阵里属于基因表达的那些 feature。QC 分布、聚类、marker 和细胞组成都在这一部分上计算。
_Avoid_: 全矩阵, 把 guide 算进聚类

**CRISPR direct capture**:
同一矩阵里属于 guide 的那些 feature。只用于 guide 频率，以及把 guide 标回细胞。
_Avoid_: 把 guide 放进聚类或 marker

**Guide 频率**:
全部 Passing cell 上，`num_features` 的分布，以及每个 Guide 的细胞数和转录本总和。同一个 name 在一个细胞里只计一次。
_Avoid_: 按 `|` 的个数计数

**表达空间**:
单 guide 细胞上，log 标准化之后的高变基因所做的 PCA。高变基因最多 2000 个，基因总数更少就全用。Cluster 和 E-distance 都在这一个 PCA 上。UMAP 只是这张邻域图的展示。
_Avoid_: X_pert, 在 UMAP 坐标上聚类

**Cluster**:
表达空间的 PCA 邻域图上划出的一群单 guide 细胞。
_Avoid_: 细胞类型, UMAP 分群

**Marker**:
把某个 Cluster 和其他单 guide 细胞分开的基因。比较用全部基因表达，去掉全零基因。p 值只给这张表排序。
_Avoid_: 探索性差异表达, 通路基因

**组成**:
单 guide 细胞在各 Cluster 中的占比。
_Avoid_: 扰动组成

**扰动组成**:
某个扰动或非靶向对照内部，各 Cluster 所占的细胞比例。内含子对照不单独成组。
_Avoid_: 组成

## 扰动效应

**靶基因**:
feature barcode 参考表里的 `target_gene_name`，同一行的 `target_gene_id` 是它的稳定标识。两列都在，且这一行两格都非空，靶基因才成立。同一个 `target_gene_name` 下有多个 id 时，仍是一个靶基因。
_Avoid_: gene_name, gene_id, 从 guide 名解析基因

**扰动**:
一个靶基因。打到该基因的单 guide 细胞都属于这一个扰动。参考表 `id` 带 `Positive_Control` 的 guide 也算，只要 `target_gene_name` 是这个基因。参考表缺 `target_gene_id` 或 `target_gene_name` 任一列时，没有扰动；描述层仍可做。Guide 身份只出现在 guide 频率里。
_Avoid_: guide, 受扰动细胞, 阳性对照类

**E-distance**:
一个扰动的全部单 guide 细胞（含 NP）与非靶向对照，在表达空间里的能量距离。为零表示两组分不开，越大表示效应越强。
_Avoid_: log2FC, p 值, X_pert

**E-test**:
在一个 Sample 内、对进入该扰动 E-distance 的那些细胞置换标签，用来参照 E-distance 的检验。与探索性差异表达同一证据层，不是跨文库可重复的显著性。
_Avoid_: 生物学重复显著性

**E-test 校正**:
一个 Sample 内、对实际跑过的全部 E-test 做的 BH 校正。与探索性差异表达同一证据层。
_Avoid_: 生物学重复 FDR, 探索性校正

**探索性差异表达**:
一个扰动的 KO 细胞与非靶向对照之间、在全部基因表达上的细胞级 Wilcoxon 比较，去掉全零基因。NP 不在其中。不构成生物学重复上的推断。
_Avoid_: pseudobulk, 生物学重复显著性

**探索性校正**:
同一扰动内部、跨基因的 BH 校正，只给探索性差异表达排序。
_Avoid_: 生物学重复 FDR, E-test 校正
