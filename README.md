# TFEB 综合实验设计

本项目协作编写一份英文分子生物学综合实验设计报告：克隆人 TFEB，并研究其在 HeLa 细胞中的营养依赖性亚细胞定位。最终交付为一份英文 PDF，所有实验结果均为预期结果，不代表已完成实验。

作者为蒋思源、宋子研；英文署名暂用 **Siyuan Jiang / Ziyan Song**，学号与实际作者贡献待确认。

## 当前状态

截至 2026-10-23，B-v1 已固定为人 TFEB `NM_007162.3` / `NP_009093.1`、pEGFP-N1 空载体 `U55762.1` 与 C 端 EGFP 融合。序列设计经本地独立路径验收：最终质粒 6134 bp，融合蛋白 722 aa，理论未修饰平均质量约 80.425 kDa。阶段结果以 [B 执行进度](docs/B执行进度.md) 与 [STATUS](docs/STATUS.md) 的实际证据为准。

整篇报告仍是工作稿：独立 A 实验细节、身份信息、真实个人贡献及课程提交要求尚待完善。计算和图表均是设计预测，没有湿实验数据。详见 [课程要求与待办](docs/10-23%20课程要求与待办.md) 和 [B 给 A 的交接](docs/B交接给A.md)。

## B 部分执行入口

蒋思源负责的 B 部分按 [详细执行计划与产出验收](docs/10-23%20B部分执行计划.md) 的 P0–P9 推进。新会话直接使用 [执行指令](docs/10-23%20B部分新会话指令.md)。计划尚未合入 main 时，从 `codex/b-execution-plan` 读取；不要只读初始化版本就开始设计。

用户已明确：**本项目不使用 UF，少量序列分析和报告工作在本地完成；不再处理 MacBook 同步。** 完成 B 所需的 A 内容可以一并补齐，不需要等 A 先完成。当前实现分支为 `codex/b-design`；完整阶段清单保留于 `docs/B执行进度.md`，不以单个 ready 文件或成功编译代替最终验收。

## 文件入口

| 路径 | 用途 |
|---|---|
| `报告/main.tex` | 英文报告主文件 |
| `报告/sections/` | 按章节拆分的正文，便于两人分别编辑 |
| `报告/references.bib` | 共用参考文献库 |
| `资料/原始提纲/TFEB_Experimental_Design_Outline.docx` | 保留不变的原始提纲 |
| `docs/10-23 课程要求与待办.md` | 课程要求来源、缺项及待确认事项 |
| `docs/10-23 模板来源.md` | ElegantPaper 固定版本、许可及编译说明 |
| `设计/construct_spec.json`、`construct.gb`、`construct.fasta` | 唯一设计参数及可导入序列工具的最终构建 |
| `设计/primers.tsv`、`expected_fragments.tsv`、`sanger_coverage.tsv` | 完整引物、鉴定预测及测序覆盖 |
| `资料/B来源索引.tsv`、`资料/序列/` | 来源台账、带版本原始记录与派生参考序列 |
| `scripts/B-sequence/`、`scripts/B-methods/`、`scripts/B-figures/` | 可复算序列、方法与图表；运行前阅读对应说明 |

课程课件不上传公开仓库；本仓库只保留与完成作业有关的要求摘要。

## 两人协作

公开仓库：[siyuanj/tfeb-experimental-design](https://github.com/siyuanj/tfeb-experimental-design)。宋子研尚未注册 GitHub，协作者权限尚未添加；注册并提供账号后再添加。

1. 每次开始工作先读取项目规则与 `docs/STATUS.md`，检查自己的未提交改动，再拉取远端最新内容。
2. `main` 保留稳定版本；从最新 `main` 创建功能分支，通过 Pull Request 合并。
3. 按下表编辑各自负责的 `sections/` 文件；共同文件如 `main.tex` 和 `references.bib` 修改前先协调。
4. 提交只包含本次相关改动。PR 说明改了什么、核查了什么、还有哪些未定项；合并前检查 PDF 与引用。

### 已确定的计划分工

| 作者 | 计划负责范围 | 主要章节文件 |
|---|---|---|
| 蒋思源（B） | CDS、载体与融合方向；引物与 PCR；酶切连接、DH5α 筛选、Sanger；质粒和克隆流程图；对应预期结果、方法及文献 | `methods-cloning.tex`、`results-cloning.tex` |
| 宋子研（A） | 摘要与关键词；TFEB、自噬/溶酶体、mTORC1 背景和假设；HeLa 培养转染、饥饿/Torin1、显微/DAPI、N/C 定量统计、Western blot；对应结果、讨论、局限及文献 | `abstract.tex`、`background.tex`、`methods-cell.tex`、`results-cell.tex`、`discussion.tex` |

文件位于 `报告/sections/`。`design.tex` 汇总两部分方法，`expected-results.tex` 汇总两部分预期结果；`declarations.tex` 与 `references.bib` 共同维护。

B 可以先独立完成克隆设计；固定载体、融合方向、完整序列及预期蛋白大小后，将这些信息交给 A 用于细胞实验与 Western blot 设计。上述是 **planned contributions**，最终 Author Contributions 按真实完成的工作填写。

如 B 设计需要修改摘要、背景、细胞检测接口或讨论中的相关内容，可按用户授权直接补齐，并记录实际改动；宋子研已有修改仍需保留。A 中独立未完成项不阻塞 B，但必须在最终交付时明确列出。

## 编译

使用 ElegantPaper **v0.12**，英文模式，参考文献后端显式设为 `bibend=bibtex`。在 `报告/` 下可使用 Tectonic：

```bash
tectonic -X compile main.tex
```

或使用 XeLaTeX 与 BibTeX：

```bash
xelatex main.tex
bibtex main
xelatex main.tex
xelatex main.tex
```

上述是项目约定的构建方式，实际已验证的环境与结果以 STATUS 为准。模板官方默认 Biber；本项目选用 BibTeX 是为了兼容 Tectonic 与 XeLaTeX 两条构建路径，详见 [模板来源](docs/10-23%20模板来源.md)。

## 本地复算与版本绑定

复用旧课程已安装的本地质粒环境，解释器为 `/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python`，不修改该环境。依赖实测版本记录于 `设计/environment.json`；其他机器应按记录建立等效环境，不能直接照用此机器路径。

序列复算顺序及长度定义见 [序列说明](scripts/B-sequence/README.md)。科学输入与派生文件通过 `docs/协作/sequence-ready.json`、`methods-ready.json`、`report-ready.json` 逐阶段绑定 SHA-256。重新生成序列或方法可能改变文件哈希；先撤回受影响的 ready，再按依赖重验，不能混用旧图、旧 PDF 和新序列。

主线最终结构核查命令为 `python scripts/B-audit/final_audit.py`，只读取科学文件并写出独立核查 JSON；须在三个阶段均就绪后运行。它不代替科学审阅、逐页目视或 Git/PR 交付检查。
