# 模块与流水线阶段

工程结构按职责拆分，逻辑以迁移/薄适配为主，不是整库重写。

| 阶段 | 职责 | 编排 | 领域模块 |
| --- | --- | --- | --- |
| 1 Input & validation | DRAGEN MEX、guide 参考、assignments | `stages/input_validation.py` | `io.py` |
| 2 Preprocessing & QC | guide 归属、细胞/样本 QC、归一化降维 | `stages/preprocessing_qc.py` | `guides.py`, `qc.py`, `preprocessing.py`, `composition.py` |
| 3a Perturbation modeling | Mixscape、E-distance、perturbation-space、guide consistency | `stages/perturbation_modeling.py` | `perturbation.py` (+ `guide_qc.py`) |
| 3b Statistical inference | pseudobulk DE、FDR、设计检查 | `stages/statistical_inference.py` | `statistics.py` |
| 4 Robustness | 敏感性、效应一致性、confidence flags | `stages/robustness.py` | `robustness.py` |
| 5 Report | H5AD / CSV / JSON / HTML / provenance / report PNGs | `stages/report.py` | `report.py`, `report_plots.py` |

分层：`pipeline.py` 只做编排 + checkpoint/`--resume`；`stages/` 调用领域模块，不把算法塞进编排器。

步骤依赖（描述性 vs 推断性、矩阵层、pre/post Mixscape）：见 [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md)。与 sc-best-practices 章节对照：见 [SC_BEST_PRACTICES_MAP.md](SC_BEST_PRACTICES_MAP.md)。生物模块覆盖清单（细胞 QC / gRNA QC / UMAP / 扰动映射 / 组成 / DE / 通路）：见同文档 [§ Tertiary module coverage](SC_BEST_PRACTICES_MAP.md#tertiary-module-coverage-checklist)。与 Seurat Mixscape vignette 对照（默认 Python；Seurat 为可选 R）：见 [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md)。

`pipeline.run_pipeline` 按 1→2→(3a then 3b)→4→5 编排。3a 与 3b 在概念上并行；进程内因 DE 可用 E-distance 排序而先跑建模。

兼容别名（薄 re-export，无独立逻辑）：`annotate.py` → `guides`；`preprocess.py` → `preprocessing`；`analysis.py` → `perturbation` + `statistics`。

`--resume`：每阶段 `status.json` 记录 `input_checksums` 与 `params_hash`（及完整 `params`）；变更则经 `repro.RunTracker` 自该阶段起失效，不复用旧 checkpoint。

并发保留模块（未合并改写实现，由编排器/阶段调用）：

- `guide_qc.py` — guide 级 QC / 一致性（经 `guides` re-export）
- `cell_annotation.py` — 细胞周期 / 状态注释
- `composition.py` — 样本组成审计（可按需挂入 stage 2）
- `repro.py` — YAML / checksum / resume / `RunTracker`
