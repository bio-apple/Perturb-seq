# 输出目录

路径相对于 `--output-dir`。图只写 PNG（`plotting.dpi` 默认 300）。

## Lean 输出

成功跑完后删除非 PNG 产物，只保留 `**/*.png`。

**保护例外**：已解析的 DRAGEN 输入文件始终保留；`--input-dir` 与 `--output-dir` 相同时额外保护 `{sample_id}.scRNA.*`。

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | 过滤前 / 描述性 QC 图 |
| `UMAP/` | HVG/PCA/marker、UMAP、注释；Mixscape signature |
| `Mixscape/` | Mixscape 图 |
| `mixscale/` | Mixscale ridge / scatter / DE 热图 |

## 各步写出的图

### assign → `qc/`

| 路径 | 内容 |
| --- | --- |
| `qc/guides_per_cell.png` | 每细胞 guide 条数 |
| `qc/cells_per_target_bar.png` | 按 `target_gene_name` 统计细胞数 |
| `qc/cells_per_guide_name_bar.png` | 按 `name` 统计细胞数 |
| `qc/qc_metrics_panel.png` | 2x3 QC：violin + 直方图（nFeature / nCount / percent.mt） |
| `qc/qc_scatter_panel.png` | 1x2 散点：log10(nCount) vs percent.mt / log10(nFeature) |
| `qc/temporary_control_qc.png` | 仅临时 zero-guide 对照时 |

### preprocess → `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `UMAP/variable_features.png`、`pca_elbow.png` | HVG / PCA |
| `UMAP/marker_detection_heatmap.png`、`dotplot_markers.png` | marker 热图 / 点图 |
| `UMAP/umap_cell_type.png` | 细胞类型（有注释时） |
| `UMAP/umap_sctype.png` | ScType 注释（需 `--sctype-tissue`） |
| `UMAP/umap_annotation_cell_cycle_panel.png` | Leiden + SingleR + Cell Cycle 组合 |

### mixscape → `Mixscape/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `UMAP/umap_mixscape_class.png` 等 | Mixscape UMAP |
| `Mixscape/mixscape_class_global_bar.png` | 分类计数 |
| `Mixscape/{靶基因}_mixscape_2x2.png` | 火山 + 热图 |
| `Mixscape/target_de_counts.png` | 全部靶 vs 进入 GMM 的靶 |

### mixscale → `mixscale/`

| 路径 | 内容 |
| --- | --- |
| `mixscale/mixscale_score_ridge.png`、`mixscale_scatter.png` | 分数山脊 / 散点 |
| `mixscale/target_de_counts.png`、`target_in_de_gene_counts.png` | DE 计数 |
| `mixscale/{靶基因}_mixscale_heatmap.png` | DE 热图（`target_in_de/` 为靶基因在自身 DE 列表的子集） |
