# 输出目录

路径相对于 `--output-dir`。图只写 PNG，dpi 为 `plotting.dpi`（默认 300）。

**精简输出**：成功跑完 `run` 后，磁盘上只保留

| 保留 | 说明 |
| --- | --- |
| `checkpoints/*.h5ad` | `--resume` 所需；表 / 门槛 / 效应等在对应 checkpoint 的 `uns` |
| `**/*.png` | 各步分析目录下的 PNG（直接写在目录根，如 `qc/*.png`、`DGE/*.png`） |

跑程中可能短暂写出 TSV / JSON / YAML 等中间文件；结束时的 lean cleanup 只留下上述两类。`--validate-only` 不做 cleanup（校验文件留在输出目录）。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [Checkpoint 里有什么](#checkpoint-里有什么)
3. [各步写出的图](#各步写出的图)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `checkpoints/` | 逐步 `*.h5ad`（resume + 表格数据） |
| `qc/` | assign / preprocess QC 图 |
| `UMAP/` | 表达 UMAP、簇注释、细胞周期；Mixscape signature UMAP |
| `Mixscape/` | Mixscape 图 |
| `mixscale/` | Mixscale ridge / scatter |
| `DGE/` | 描述性效应图 |
| `E-distance/` | E-test 散点 |

`dge` / `mixscape` / `mixscale` 互不依赖；Wilcoxon 结果只在 `checkpoints/dge.h5ad`（无独立 PNG 桶）。

---

## Checkpoint 里有什么

| 文件 | 关键 `uns` / 用途 |
| --- | --- |
| `checkpoints/io.h5ad` | `cell_calls`（完整注释账本，供 `--steps assign --resume`） |
| `checkpoints/assign.h5ad` | `qc_thresholds`、`cells_per_target`、对照字段 |
| `checkpoints/preprocess.h5ad` | 归一化 / PCA / UMAP / Leiden 等 |
| `checkpoints/dge.h5ad` | `wilcoxon_de_results`、`wilcoxon_de_status` |
| `checkpoints/mixscape.h5ad` | `mixscape_de_gate`、`mixscape_summary`、分类与 `X_pert` |
| `checkpoints/mixscale.h5ad` | `mixscale_de_gate`、分数与 `X_pert_mixscale` |
| `checkpoints/perturb.h5ad` | `effect_*`、`guide_effects`、`target_expression_effect`、`edistance_etest` 等 |

需要完整对象时读最新相关 checkpoint（通常 `perturb.h5ad`）。`assign --resume` 依赖 `io.h5ad` 的 `uns['cell_calls']`；缺失时请重跑 `io`。

---

## 各步写出的图

跳过原因只打日志 / 记在内存。`plotting.save_plot_data` 默认 `false`（不写图旁 TSV）。

刷新某类图：对该步 `--steps <名> --resume`。

### assign → `qc/`

| 路径 | 内容 |
| --- | --- |
| `qc/guides_per_cell.png` | 每细胞 guide 条数 |
| `qc/cells_per_target_bar.png` | 靶基因细胞数 |
| `qc/qc_metrics_panel.png` | 2×3 QC：上行 Seurat 风格 violin（nFeature / nCount / percent.mt）+ jitter；下行 Library-D 直方图（同三指标）。轴标签英文；Identity=`sample_id` |
| `qc/qc_scatter.png` / `scatter_ncount_nfeature_gex.png` | QC 散点 |
| `qc/temporary_control_qc.png` | 仅临时 zero-guide 对照时；组标签 `zero-guide` / `targeting`（纯英文） |

### preprocess → `qc/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `qc/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | HVG / PCA |
| `qc/marker_detection_heatmap.png`、`dotplot_markers.png` | marker |
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
| `mixscale/mixscale_score_ridge.png`、`mixscale_scatter.png` | 分数山脊 / 箱线 |

### perturb → `DGE/` + `E-distance/`

| 路径 | 内容 |
| --- | --- |
| `DGE/de_gene_heatmap.png` | 描述性 log2FC |
| `DGE/target_expression_effect.png`、`guide_consistency.png` | 靶基因 / guide |
| `E-distance/edistance_etest.png` | E-test 散点 |
