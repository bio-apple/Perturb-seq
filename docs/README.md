# 文档索引

| 文档 | 内容 |
| --- | --- |
| [analysis.md](analysis.md) | 流水线步骤、CLI / 配置、物种、对照、QC、分配、Mixscape / Mixscale、预处理与 SingleR / ScType |
| [outputs.md](outputs.md) | lean 输出（PNG + 输入保护规则）、图路径 |

安装与最小命令见仓库根目录 [README.md](../README.md)。ScType 由 CLI 启用（`--sctype-tissue`）；离线库、tissue 列表与自定义 DB 见根目录 [README.md](../README.md)。

步骤：`io` → `assign` → `preprocess` → `mixscape` ‖ `mixscale`。不写 checkpoints，不支持 `--resume`，见 [analysis.md](analysis.md)。

数据分析流程图（Mermaid）：见 [analysis.md § 流水线步骤](analysis.md#流水线步骤)。
