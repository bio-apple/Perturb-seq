# 文档索引

| 文档 | 内容 |
| --- | --- |
| [analysis.md](analysis.md) | 流水线步骤、CLI / 配置、物种、对照、QC、分配、统计、Mixscape / Mixscale、预处理 |
| [outputs.md](outputs.md) | 精简输出（checkpoint + PNG）、图路径、checkpoint `uns` |

安装与最小命令见仓库根目录 [README.md](../README.md)。

默认步骤：`io` → `assign` → `preprocess` → `{dge ‖ mixscape ‖ mixscale}` → `perturb`（后三者并发，见 [analysis.md](analysis.md)）。

数据分析流程图（Mermaid，含输入 / 对照 / checkpoint / 输出分桶）：见 [analysis.md § 流水线步骤](analysis.md#流水线步骤)。
