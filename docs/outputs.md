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

非 scRNA 的中间文件仍可能被删。不写出 `checkpoints/*.h5ad`，不支持 `--resume`；需从 `io` 起连续运行（通常 `--steps all`）。`--validate-only` 不做 cleanup（校验文件留在输出目录）。无独立 `DGE/` / `E-distance/` / `wilcoxon/` 输出目录（旧布局若残留，作图时会尝试清掉）。

判定规则见 [analysis.md](analysis.md)。

## 目录

1. [目录结构](#目录结构)
2. [各步写出的图](#各步写出的图)

---

## 目录结构

| 目录 | 内容 |
| --- | --- |
| `qc/` | assign / preprocess QC 图 |
| `UMAP/` | 表达 UMAP、SingleR / ScType、细胞周期；Mixscape signature UMAP |
| `Mixscape/` | Mixscape 图（`--steps` 含 `mixscape` 时） |
| `mixscale/` | Mixscale ridge / scatter / DE 计数 / DE 热图（含 `mixscale` 时） |

`mixscape` / `mixscale` 互不依赖；并行步之间用**临时目录**交换 AnnData，跑完即删，不落盘到输出目录。

---

## 各步写出的图

跳过原因只打日志 / 记在内存。`plotting.save_plot_data` 默认 `false`（不写图旁 TSV）。

### assign → `qc/`

| 路径 | 内容 |
| --- | --- |
| `qc/guides_per_cell.png` | 每细胞 guide 条数 |
| `qc/cells_per_target_bar.png` | 各靶基因细胞数（参考/分配全部靶，含 0；不做 analysis_eligible / 门槛过滤）。虚线=`qc.min_cells_per_target`；点线/实线/点划线=全靶 **min / median / max**（`n_cells`）。图例另报 **total cells**（全靶含 0 的 `n_cells` 总和，与柱高之和一致）、**Targets with n_cells ≥ threshold** 个数，以及这些靶的 **sum n_cells**（过门槛后剩余总细胞数；与 assign 过滤一致：`n_cells < min_cells_per_target` 记 low_cells_target，故保留为 ≥） |
| `qc/qc_metrics_panel.png` | 2×3 QC：上行 Seurat 风格 violin（nFeature / nCount / percent.mt）+ jitter；下行 Library-D 直方图（同三指标）。轴标签英文；Identity=`sample_id` |
| `qc/qc_scatter_panel.png` | 1×2 QC 散点：nCount vs percent.mt；nCount vs nFeature |
| `qc/temporary_control_qc.png` | 仅临时 zero-guide 对照时；组标签 `zero-guide` / `targeting`（纯英文） |

### preprocess → `qc/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `qc/variable_features.png`、`pca_elbow.png`、`pca_scatter.png` | HVG / PCA |
| `qc/marker_detection_heatmap.png`、`dotplot_markers.png` | marker；白底；行=注释簇名；dotplot 为 Seurat 式 Percent/Average 双图例 |
| `UMAP/umap_cell_type.png` | 真实细胞类型（有注释时） |
| `UMAP/umap_sctype.png` | ScType 注释（`obs['sctype_classification']`；需 `sctype.enabled` 且设置 `tissue`） |
| `UMAP/umap_annotation_cell_cycle_panel.png` | 2 行组合（QC-pass 嵌入全细胞，非仅 analysis_eligible）：上行 Leiden + SingleR；下行 Perturbed / 对照 / Cell Cycle；`singler.enabled: false` 时跳过 |

### mixscape → `Mixscape/` + `UMAP/`

| 路径 | 内容 |
| --- | --- |
| `UMAP/umap_mixscape_class.png` / `umap_mixscape.png` / `umap_mixscape_class_global.png` | Mixscape UMAP |
| `Mixscape/mixscape_class_global_bar.png` | 分类计数 |
| `Mixscape/{靶基因}_mixscape_2x2.png`、`target_in_de/...` | 火山 + 热图 |
| `Mixscape/target_de_counts.png`、`target_in_de_gene_counts.png` | DE 计数 |

### mixscale → `mixscale/`

| 路径 | 内容 |
| --- | --- |
| `mixscale/mixscale_score_ridge.png`、`mixscale_scatter.png` | 分数山脊 / 分箱散点；有 Mixscape 标签时按 NT/NP/KO（或 KD/OE）分行或分线 |
| `mixscale/target_de_counts.png`、`target_in_de_gene_counts.png` | Mixscale DE 门计数（对标 Mixscape 同名图） |
| `mixscale/{靶基因}_mixscale_heatmap.png`、`target_in_de/...` | Mixscale DE 热图：有 Mixscape 时列块 **NT \| NP \| KO**；否则 NT \| perturbed；`target_in_de/` 为「靶基因在自身 DE 列表」子集 |
