# B 课程要求与产物映射

日期：2026-10-23。来源为 `docs/10-23 课程要求与待办.md` 对课件第2、3、5–7页的摘录；原课件不公开。

| 要求 | 支撑产物 | 正文/验收 |
|---|---|---|
| 基因信息、选择依据、研究基础 | `资料/序列/`、`资料/B来源索引.tsv`、`docs/B设计决策.md` | methods-cloning；background；P1–P3/P8 |
| 正确5′→3′引物、保护碱基、Tm | `设计/primers.tsv`、`primer_checks.json`、`construct_spec.json` | 引物表、初期退火模型和优化边界；P4 |
| 合理PCR程序与体系 | `设计/methods_parameters.json`、方法正文 | 方法反应表；P6 |
| 酶切、载体/插入处理、定向读框 | `construct_metrics.json`、`junctions.tsv`、完整构建/蛋白 | 接头和片段表；P4–P6 |
| 抗生素、克隆筛选、完整鉴定 | `expected_fragments.tsv`、`sanger_primers.tsv`、`sanger_coverage.tsv` | KanR与Neo/G418分开；PCR/胶/Sanger接受规则；P6 |
| 质粒图与预期结果、成功标准 | `报告/images/B/`、可重跑图脚本 | 环图/流程/预测胶，不伪造照片；P7 |
| 后续表达与定位鉴定接口 | `docs/B交接给A.md` | TFEB–EGFP、空载体、理论蛋白质量及局限；P8 |
| 一份英文报告、作者贡献和AI声明、参考文献 | `报告/main.tex`及sections、references.bib、main.pdf | 独立A和身份/提交待办保留；P8–P9 |
| 可复核准确性 | `scripts/B-sequence/`、`设计/verification_summary.json`、`docs/B验收记录.md` | 独立序列路径+错误注入+PDF逐页验收；P5/P9 |

本表是映射，不是阶段通过证据。实际通过状态仅见 `docs/B执行进度.md` 与验收记录。
