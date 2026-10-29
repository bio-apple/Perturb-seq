# 输出目录

路径相对于 `--output-dir`。图只写 PNG，dpi 为 `plotting.dpi`（默认 300）。

**Lean 输出**：成功跑完 `run` 后，`cleanup_lean_outputs` 删除输出目录内**非 PNG** 的流水线产物（含跑程中短暂写出的 TSV / JSON / YAML 等）。默认只保留：

| 保留 | 说明 |
| --- | --- |
| `**/*.png` | 各步分析目录下的 PNG（直接写在目录根，如 `qc/*.png`、`Mixscape/*.png`、`UMAP/*.png`） |

**保护例外**（不会被 cleanup 删除）：

| 条件 | 保护内容 |
| --- | --- |
| 始终 | 本样本已解析的 DRAGEN 输入文件路径（`resolve_input_paths`） |
| `--output-dir` 与 `--input-dir` 为同一目录 | 额外保护 `{sample_id}.scRNA.*`（避免同目录读写时删掉原始 scRNA 输入） |

非 scRNA 的中间文件仍可能被删。不写出 `checkpoints/*.h5ad`，不支持 `--resume`；需从 `io` 起连续运行。`--validate-only` 不做 cleanup（校验文件留在输出目录）。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [各步写出的图](#各步写出的图)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | assign 步：过滤前 / 描述性 QC 图 |
| `UMAP/` | QC-pass 嵌入相关图（HVG/PCA/marker、UMAP、注释）；Mixscape signature UMAP |
| `Mixscape/` | Mixscape 图 |
| `mixscale/` | Mixscale ridge / scatter / DE 计数 / DE 热图 |

`mixscape` / `mixscale` 互不依赖；并行步之间用**临时目录**交换 AnnData，跑完即删，不落盘到输出目录。

---

## 各步写出的图

跳过原因只打日志 / 记在内存。`plotting.save_plot_data` 默认 `false`（不写图旁 TSV）。

### assign → `qc/`

| 路径 | 内容 |
| --- | --- |
| `qc/guides_per_cell.png` | 每细胞 guide 条数 |
| `qc/cells_per_target_bar.png` | 按参考表 **`target_gene_name`**（账本 `target_gene`）统计细胞数（含 0；不做 analysis_eligible 过滤）。柱下标注靶名。虚线=`qc.min_cells_per_target`；点线/实线/点划线=**min / median / max**。图例另报 total cells、过门槛靶数与其 sum n_cells（`n_cells ≥ min_cells_per_target`） |
| `qc/cells_per_guide_name_bar.png` | 按参考表 **`name`**（账本 `guide_name`）统计细胞数（含 0；不做 analysis_eligible 过滤）。柱下标注 guide 名。虚线=`qc.min_cells_per_guide`（仅参考，不排除细胞） |
| `qc/qc_metrics_panel.png` | 2×3 QC：上行 Seurat 风格 violin（nFeature / nCount / percent.mt）+ jitter；下行 Library-D 直方图（同三指标）。轴标签英文；Identity=`sample_id` |
| `qc/qc_scatter_panel.png` | 1×2 QC 散点：log10(nCount) vs percent.mt；log10(nCount) vs log10(nFeature) |
| `qc/temporary_control_qc.png` | 仅临时 zero-guide 对照时；组标签 `zero-guide` / `targeting`（纯英文） |

### preprocess → `UMAP/`（QC-pass 细胞）

| 路径 | 内容 |
| --- | --- |
| `UMAP/variable_features.png`、`pca_elbow.png` | HVG / PCA |
| `UMAP/marker_detection_heatmap.png`、`dotplot_markers.png` | marker；热图/点图白底；热图紫→黑→黄（检测比例）；dotplot 浅灰→蓝，点大小=Percent Expressed |
| `UMAP/umap_cell_type.png` | 真实细胞类型（有注释时） |
| `UMAP/umap_sctype.png` | ScType 注释（`obs['sctype_classification']`；需命令行 `--sctype-tissue`） |
| `UMAP/umap_annotation_cell_cycle_panel.png` | 2 行组合（QC-pass 嵌入全细胞，非仅 analysis_eligible）：上行 Leiden + SingleR；下行 Perturbed / 对照 / Cell Cycle；`singler.enabled: false` 时跳过 |

### mixscape → `Mixscape/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `UMAP/umap_mixscape_class.png` / `umap_mixscape.png` / `umap_mixscape_class_global.png` | Mixscape UMAP |
| `Mixscape/mixscape_class_global_bar.png` | 分类计数 |
| `Mixscape/{靶基因}_mixscape_2x2.png` | 火山（adj-p Up/Down）+ 热图（对齐 Seurat MixscapeHeatmap：RNA、NT-high、raw p、块内后验排序，最多 Top20） |
| `Mixscape/target_de_counts.png` | 全部靶 vs 进入 Mixscape GMM 的靶。pertpy 不返回基因名，因此没有 `target_in_de_gene_counts` |

### mixscale → `mixscale/`

| 路径 | 内容 |
| --- | --- |
| `mixscale/mixscale_score_ridge.png`、`mixscale_scatter.png` | 分数山脊 / 分箱散点；有 Mixscape 标签时按 NT/NP/KO（或 KD/OE）分行或分线 |
| `mixscale/target_de_counts.png`、`target_in_de_gene_counts.png` | 来自 pertpy `mixscale_de_genes` 的计数（不是 Mixscape 的基因表） |
| `mixscale/{靶基因}_mixscale_heatmap.png`、`target_in_de/...` | 对齐 Mixscale_DoHeatmap：列块 **NT \| 靶**（按 score）；行 = 靶基因 + pertpy 按原始 p 排序的 DE 基因；`target_in_de/` 为「靶基因在自身 DE 列表」子集 |
