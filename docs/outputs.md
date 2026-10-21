# 输出目录

路径相对于 `--output-dir`。图只写 PNG，dpi 为 `plotting.dpi`（默认 300）。

**精简输出**：成功跑完 `run` 后，磁盘上只保留

| 保留 | 说明 |
| --- | --- |
| `**/*.png` | 各步分析目录下的 PNG（直接写在目录根，如 `qc/*.png`、`DGE/*.png`） |

跑程中可能短暂写出 TSV / JSON / YAML 等中间文件；结束时的 lean cleanup 只留下 PNG。`--validate-only` 不做 cleanup（校验文件留在输出目录）。不写出 `checkpoints/*.h5ad`，不支持 `--resume`；需从 `io` 起连续运行（通常 `--steps all`）。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [各步写出的图](#各步写出的图)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | assign / preprocess QC 图 |
| `UMAP/` | 表达 UMAP、簇注释、细胞周期；Mixscape signature UMAP |
| `Mixscape/` | Mixscape 图 |
| `mixscale/` | Mixscale ridge / scatter / DE 计数 / DE 热图 |
| `DGE/` | 描述性效应图 |
| `E-distance/` | E-test 散点 |

`dge` / `mixscape` / `mixscale` 互不依赖；Wilcoxon（`dge`）无独立 PNG。并行步之间用**临时目录**交换 AnnData，跑完即删，不落盘到输出目录。

---

## 各步写出的图

跳过原因只打日志 / 记在内存。`plotting.save_plot_data` 默认 `false`（不写图旁 TSV）。

### assign → `qc/`

| 路径 | 内容 |
| --- | --- |
| `qc/guides_per_cell.png` | 每细胞 guide 条数 |
| `qc/cells_per_target_bar.png` | 靶基因细胞数 |
| `qc/qc_metrics_panel.png` | 2×3 QC：上行 Seurat 风格 violin（nFeature / nCount / percent.mt）+ jitter；下行 Library-D 直方图（同三指标）。轴标签英文；Identity=`sample_id` |
| `qc/qc_scatter_panel.png` | 1×2 QC 散点：nCount vs percent.mt；nCount vs nFeature |
| `qc/temporary_control_qc.png` | 仅临时 zero-guide 对照时；组标签 `zero-guide` / `targeting`（纯英文） |

### preprocess → `qc/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `qc/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | HVG / PCA |
| `qc/marker_detection_heatmap.png`、`dotplot_markers.png` | marker；白底；行=注释簇名；dotplot 为 Seurat 式 Percent/Average 双图例 |
| `UMAP/umap_cell_type.png` | 真实细胞类型（有注释时） |
| `UMAP/umap_annotation_cell_cycle_panel.png` | 2 行组合：上行 Leiden + SingleR；下行 Perturbed / 对照 / Cell Cycle |

### mixscape → `Mixscape/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `UMAP/umap_mixscape_class.png` / `umap_mixscape.png` / `umap_mixscape_class_global.png` | Mixscape UMAP |
| `Mixscape/mixscape_class_global_bar.png`、`mixscape_fraction_bar.png` | 分类计数 / 比例 |
| `Mixscape/{靶基因}_mixscape_2x2.png`、`target_in_de/...` | 火山 + 热图 |
| `Mixscape/target_de_counts.png`、`target_in_de_gene_counts.png` | DE 计数 |

### mixscale → `mixscale/`

| 路径 | 内容 |
| --- | --- |
| `mixscale/mixscale_score_ridge.png`、`mixscale_scatter.png` | 分数山脊 / 分箱散点；有 Mixscape 标签时按 NT/NP/KO（或 KD/OE）分行或分线 |
| `mixscale/target_de_counts.png`、`target_in_de_gene_counts.png` | Mixscale DE 门计数（对标 Mixscape 同名图） |
| `mixscale/{靶基因}_mixscale_heatmap.png`、`target_in_de/...` | Mixscale DE 热图：有 Mixscape 时列块 **NT \| NP \| KO**；否则 NT \| perturbed；`target_in_de/` 为「靶基因在自身 DE 列表」子集 |

### perturb → `DGE/` + `E-distance/`

| 路径 | 内容 |
| --- | --- |
| `DGE/de_gene_heatmap.png` | 描述性 log2FC |
| `DGE/guide_consistency.png`、`target_expression_effect.png` | guide / 靶基因自身表达 |
| `E-distance/edistance_etest.png` | 细胞级 E-test |
