# TFEB 综合实验设计

本项目协作编写一份英文分子生物学综合实验设计报告：克隆人 TFEB，并研究其在 HeLa 细胞中的营养依赖性亚细胞定位。最终交付为一份英文 PDF，所有实验结果均为预期结果，不代表已完成实验。

作者为蒋思源、宋子研；英文署名暂用 **Siyuan Jiang / Ziyan Song**，学号与实际作者贡献待确认。

用户明确本次作业不要求很强的创新性：当前保持 TFEB 构建、HeLa 表达及正常培养/饥饿/Torin1 定位比较。时间过程、营养恢复与突变体等作为后续完善项。已有研究依据与阅读入口见 [文献调研](文献调研/README.md)。

## 当前状态

截至 2026-10-23，本轮仅初始化协作项目、迁移已有提纲并建立 LaTeX 文档结构。当前是工作草稿，尚不能作为最终作业提交；序列、引物、载体、实验条件和身份信息仍需补齐。详见 [课程要求与待办](docs/10-23%20课程要求与待办.md)，本次核查与构建结果见 [STATUS](docs/STATUS.md)。

## 文件入口

| 路径 | 用途 |
|---|---|
| `报告/main.tex` | 英文报告主文件 |
| `报告/sections/` | 按章节拆分的正文，便于两人分别编辑 |
| `报告/references.bib` | 共用参考文献库 |
| `资料/原始提纲/TFEB_Experimental_Design_Outline.docx` | 保留不变的原始提纲 |
| `docs/10-23 课程要求与待办.md` | 课程要求来源、缺项及待确认事项 |
| `docs/10-23 模板来源.md` | ElegantPaper 固定版本、许可及编译说明 |
| `文献调研/` | 文献总报告、检索记录和分类阅读笔记；扩展建议不增加当前必做范围 |

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
