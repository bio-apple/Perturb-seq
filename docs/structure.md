# 模块与流水线阶段

工程结构按职责拆分，逻辑以迁移/薄适配为主，不是整库重写。

| 阶段 | 职责 | 模块 |
| --- | --- | --- |
| 1 Input & validation | DRAGEN MEX、guide 参考、assignments | `io.py` |
| 2 Preprocessing & QC | guide 归属、细胞/样本 QC、归一化降维 | `guides.py`, `qc.py`, `preprocessing.py`, `composition.py` |
| 3a Perturbation modeling | Mixscape、E-distance、guide consistency | `perturbation.py` (+ `guide_qc.py`) |
| 3b Statistical inference | pseudobulk DE、FDR、设计检查 | `statistics.py` |
| 4 Robustness | 敏感性、效应一致性、confidence flags | `robustness.py` |
| 5 Report | H5AD / CSV / JSON / HTML / provenance | `report.py` |

步骤依赖（描述性 vs 推断性、矩阵层、pre/post Mixscape）：见 [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md)。

`pipeline.py` 按 1→2→(3a then 3b)→4→5 编排。3a 与 3b 在概念上并行；进程内因 DE 可用 E-distance 排序而先跑建模。

兼容别名：`annotate.py` → `guides`；`preprocess.py` → `preprocessing`；`analysis.py` → `perturbation` + `statistics`。

并发保留模块（未合并改写实现，由编排器调用）：

- `guide_qc.py` — guide 级 QC / 一致性（经 `guides` re-export）
- `cell_annotation.py` — 细胞周期 / 状态注释
- `composition.py` — 样本组成审计（可按需挂入 stage 2）
- `repro.py` — YAML / checksum / resume / `RunTracker`
