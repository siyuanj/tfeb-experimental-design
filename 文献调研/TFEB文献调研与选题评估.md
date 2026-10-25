# TFEB 文献调研与课程选题评估

> 2026-10-25 范围更新：用户已要求继续研究延伸；当前新增 C 端位点 × 标签 × 刺激比较，详见 [机制延伸记录](../docs/10-25%20机制延伸与流程图重绘.md)。以下保留 10 月 23 日调研及当时范围决定，不代表当前仅限基础复现。

调研日期：2026-10-23（America/New_York）。用途：为综合实验设计确定问题边界、文献依据与后续阅读重点。本文为中文研究笔记，英文提交稿仍由 `报告/main.tex` 及其章节文件维护。

## 1. 结论：核心实验已有直接先例，适合定位为验证性课程设计

**人 TFEB 的克隆、荧光融合表达，以及在 HeLa 细胞中比较饥饿和 mTOR 抑制后的核/胞质定位，已经被直接研究。** 2012 年的两篇关键论文已覆盖这套组合；2018 年研究进一步解析营养恢复后的核输出；2025 年已有 HeLa 单细胞动态成像和微流控研究。因此，现有三组定位方案不能写成首次发现。[R03](https://doi.org/10.1038/emboj.2012.32)、[R04](https://doi.org/10.1126/scisignal.2002790)、[R10](https://doi.org/10.1038/s41467-018-05862-6)、[R16](https://doi.org/10.1038/s42003-025-07870-x)

这并不妨碍它成为课程大作业。讲义要求开放性实验设计、完整方法与预期结果，不要求实际完成实验。**用户决定（2026-10-23）：本次作业不需要很强的创新性；当前保留原基础方案，动态与机制扩展作为之后希望完善时的下一步。** 因此，文献调研用于支持现有方案、核准细节和避免错误解释，不以追求创新为由扩大必做范围。这里记录用户的范围决定，不代替教师对最终报告的评价。

本次最影响后续设计的发现是：

- 最贴近当前问题的经典质粒是 **TFEB–EGFP，EGFP 在 C 端**；原提纲拟选的 **EGFP–TFEB，EGFP 在 N 端**也有文献先例，但二者不能混用名称和终止密码子逻辑。[R04](https://doi.org/10.1126/scisignal.2002790)、[R12](https://doi.org/10.1038/s41598-021-00613-y)、[Addgene #38119](https://www.addgene.org/38119/)
- 核定位、靶基因转录和自噬通量是不同读出；仅看到 TFEB 入核，不能声称自噬降解能力已经提高。[R02](https://doi.org/10.1126/science.1204592)、[R15](https://doi.org/10.1016/j.isci.2024.110432)
- 瞬时过表达可能改变定位。2026 年研究在 HeLa 中观察到表达量与核定位正相关；“胞质滞留系统被过量蛋白占满”是作者提出的解释，而不是本项目实测事实。[R18](https://doi.org/10.1038/s41467-026-71001-1)
- S211A、S142A 和再喂养实验均已有研究，可以提高课程设计的解释力，不能仅因加入它们就宣称创新。[R08](https://doi.org/10.1080/15548627.2016.1271514)、[R09](https://doi.org/10.1038/s41467-018-04849-7)、[R10](https://doi.org/10.1038/s41467-018-05862-6)

## 2. 调研范围与证据等级

本次整理 **19 篇原始研究、1 份方法指南和 3 项官方质粒记录**，覆盖 2009–2026 年。重点是营养信号、核质转运、荧光构建、突变体、表达量及功能读出；没有试图罗列全部肿瘤、神经退行性疾病或药物筛选论文。

检索入口为 PubMed、PMC、Europe PMC 和原期刊页面，通过题名、作者、DOI 追溯核验；质粒信息以 Addgene 官方记录为准。优先使用原始研究，方法指南单独列出。检索过程和本地证据见 [检索记录](检索记录.md)。

本文是**定向叙述性调研**：没有预注册系统综述、穷尽性去重统计或统一质量评分。表中“正文”表示读到与问题相关的结果、图注或方法，不表示逐字通读全文、补充材料和原始数据。网页可访问、论文已发表与结论可推广到所有细胞条件，也不是同一件事。

## 3. 大家已经研究到哪一步

### 3.1 建立功能与营养定位机制：2009–2017

| 编号与研究 | 实际研究内容及主要证据 | 对本项目的意义 | 本次阅读层级 |
|---|---|---|---|
| **R01 Sardiello 2009**，Science | 识别协调溶酶体生成和功能的基因网络；TFEB 调节相关基因，过表达提高部分底物清除 | 解释为什么选 TFEB；不是当前三组处理的直接方法来源 | 元数据、摘要 |
| **R02 Settembre 2011**，Science | 结合 TFEB 增益/缺失、LC3 相关检测及降解读出，连接自噬与溶酶体生成；研究饥饿、S142 与 ERK2 | 功能背景；提醒功能结论需要定位以外的证据 | 正文及关键图注 |
| **R03 Settembre 2012**，EMBO J | HeLa 稳定 TFEB–GFP，高内涵成像和核/胞质强度比；研究饥饿、Torin、磷酸化位点和 Rag | 与本项目高度重合；直接支持模型和读出 | 正文、关键图注及指定方法 |
| **R04 Roczniak-Ferguson 2012**，Sci Signal | HeLa/HeLa M、pEGFP-N1 中 TFEB–GFP；比较 EBSS、rapamycin、Torin1；研究 S211–14-3-3、核定位信号及内源 TFEB | 最应先读的实验模板；当前基础假设已有直接证据 | 摘要、正文、Fig. 3/5 及克隆/分级方法 |
| **R05 Martina 2012**，Autophagy | 主要用 ARPE-19 的 TFEB–Flag，结合催化抑制剂和 RPTOR/RICTOR 干预；支持 mTORC1 调节核转运 | 独立机制证据；药物效应存在模型差异 | 摘要、关键结果及图注 |
| **R06 Martina & Puertollano 2013**，J Cell Biol | Rag 参与氨基酸依赖的 TFEB 溶酶体募集；分析 N 端区域与突变 | N 端不是无功能的接头区；选标签方向须有依据 | 关键结果、方法相关段落 |
| **R07 Medina 2015**，Nat Cell Biol | HeLa TFEB–GFP 与内源 TFEB；calcineurin/MCOLN1 干预支持溶酶体 Ca²⁺–去磷酸化途径 | 饥饿响应包含磷酸酶调控；不能仅用一种药物概括全部机制 | 主图图注、定向结果与方法 |
| **R08 Vega-Rubin-de-Celis 2017**，Autophagy | HeLa S211A 仍可分布于核和胞质，并继续响应 Torin1；进一步研究 S122 | 不把 S211A 预设为“永远全核定位、完全不受药物影响” | 正文相关结果及 Fig. 1/4 |

对应原始来源：[R01](https://doi.org/10.1126/science.1174447)、[R02](https://doi.org/10.1126/science.1204592)、[R03](https://doi.org/10.1038/emboj.2012.32)、[R04](https://doi.org/10.1126/scisignal.2002790)、[R05](https://doi.org/10.4161/auto.19653)、[R06](https://doi.org/10.1083/jcb.201209135)、[R07](https://doi.org/10.1038/ncb3114)、[R08](https://doi.org/10.1080/15548627.2016.1271514)。

### 3.2 从“进入细胞核”发展到转运动力学和具体调控机制：2018–2024

| 编号与研究 | 实际研究内容及主要证据 | 对本项目的意义 | 本次阅读层级 |
|---|---|---|---|
| **R09 Li 2018**，Nat Commun | MCF7/HT29、内源 TFEB 和稳定 TFEB–GFP；S142/S138、GSK3β 与核输出信号联系氨基酸及葡萄糖供应 | “饥饿”必须说清缺什么；营养恢复及突变比较已有先例 | 正文、关键图注和相关方法 |
| **R10 Napolitano 2018**，Nat Commun | HeLa/HEK293T，活细胞成像和突变；营养恢复促进 TFEB 经 CRM1/XPO1 输出，S142/S138 与 S211 的作用不同 | 终点 N/C 比是输入、输出和胞质滞留的合成结果 | 正文、Fig. 1/3/4 及相关方法 |
| **R11 Napolitano 2020**，Nature | 细胞遗传干预及疾病模型显示，TFEB 的 mTORC1 磷酸化具有 Rag/FLCN 依赖的底物特异性 | p-S6K 可作通路参考，不能替代直接测 TFEB；血清与氨基酸饥饿不能等同 | 关键结果、图注；未全面审读动物部分 |
| **R12 Park 2021**，Sci Rep | 研究 TFEB 剪接变体；方法中存在全长 TFEB 插入 pEGFP-C1 的 N 端 EGFP 构建 | N 端融合并非无先例；不能据此认定其在当前 HeLa 条件优于 C 端 | 元数据、指定方法及结果片段 |
| **R13 Cui 2023**，Nature | 冷冻电镜与细胞突变解析 mTORC1–TFEB–Rag–Ragulator 复合物，说明底物如何被呈递 | 为 Rag 依赖机制提供结构依据；不需要在课程中复现结构研究 | 摘要、结构结果和讨论相关段 |
| **R14 Malik 2023**，Science | AMPK/FNIP1 遗传与生化证据，支持能量压力通过 FNIP1 影响 TFEB 的定位和转录程序 | 区分能量缺乏和氨基酸缺乏；拓展背景即可 | 关键正文、主图图注 |
| **R15 Zhao 2024**，iScience | C2C12/U2OS，定向招募磷酸酶和内源标签；去磷酸化促入核，但所测转录输出有限 | 入核与完整功能激活不是同义词 | 摘要、关键结果、方法 |

对应原始来源：[R09](https://doi.org/10.1038/s41467-018-04849-7)、[R10](https://doi.org/10.1038/s41467-018-05862-6)、[R11](https://doi.org/10.1038/s41586-020-2444-0)、[R12](https://doi.org/10.1038/s41598-021-00613-y)、[R13](https://doi.org/10.1038/s41586-022-05652-7)、[R14](https://doi.org/10.1126/science.abj5559)、[R15](https://doi.org/10.1016/j.isci.2024.110432)。

### 3.3 最近研究：动态、表达量和仍在更新的机制

| 编号与研究 | 本次核实的发现 | 与课程的关系 | 本次阅读层级 |
|---|---|---|---|
| **R16 Ruolo 2025**，Commun Biol，2025-03-15 | 稳定 TFEB–GFP/H2B-mCherry HeLa，微流控、单细胞成像和模型；研究营养变化、Torin1/AZD8055 与核输出 | 基础定位及再喂养都已有更深入工作。可借鉴时间过程，不必引入微流控或数学建模 | 关键结果、图注、方法与讨论 |
| **R17 Negoita 2026**，Autophagy，2026-02-23 | AMPK 激活伴随 TFEB N/C 端位点去磷酸化；对既往某些“AMPK 直接磷酸化 C 端位点”的解释提出相反证据 | 具体激酶箭头仍需谨慎；不影响经典定位现象已被研究这一结论 | 元数据、全文关键结果与方法 |
| **R18 Xiong 2026**，Nat Commun，2026-04-01 | IKK–β-TrCP2 调节 TFEB 降解；前期 HeLa 体系观察瞬时表达量与核定位相关，较低表达的稳定体系无同样关系 | 直接提醒表达量控制、不能只挑最亮细胞；区分总量与分布 | 摘要、结果首节和相关方法 |
| **R19 Esposito 2026**，Nature，2026-09-30 | 报告营养非依赖的 v-ATPase/TAX1BP1–TBK1/ULK1–FNIP1 通路，联系 TFEB 与溶酶体应激 | 领域仍有新问题，但不能把“不同应激”当作与饥饿完全同一机制；本课程无须复制 | 期刊/PubMed 摘要、主图题名、出版信息；未读完整方法 |

对应原始来源：[R16](https://doi.org/10.1038/s42003-025-07870-x)、[R17](https://doi.org/10.1080/15548627.2026.2629720)、[R18](https://doi.org/10.1038/s41467-026-71001-1)、[R19](https://doi.org/10.1038/s41586-026-11093-3)。R19 的日期还经 [PubMed PMID 42816627](https://pubmed.ncbi.nlm.nih.gov/42816627/) 核对。本次纳入的最近发表论文为这篇 2026-09-30 论文，不表示它是整个领域截至检索日的最后一篇论文。

## 4. 与我们现有方案逐项比较

| 当前构想 | 先前直接证据 | 本项目可以合理主张什么 |
|---|---|---|
| 克隆人 TFEB，构建荧光表达载体 | R04 的 pEGFP-N1-TFEB；官方 #38119 | 完成明确序列、接头和筛选方案；不是首次构建 TFEB 荧光质粒 |
| HeLa 中表达融合蛋白 | R03、R04、R10、R16、R18 | 建立一个有文献依据的验证模型 |
| 正常培养、饥饿、Torin1 三组 | R03、R04；R16 的动态扩展 | 验证指定构建对处理的定位响应 |
| 核标记界定核区并定量核/胞质荧光；本方案拟用 DAPI | R03、R04、R10、R16 支持定位定量，核区标记方式各异 | 比单张代表图更完整；不能把 H2B-mCherry 或活细胞方法都写成 DAPI 固定染色 |
| anti-GFP 全细胞 Western blot | 既有 TFEB 融合表达研究 | 支持蛋白表达、表观大小及裂解产物判断；不能单独证明核转位 |
| 新增 S211A 或营养恢复 | R08、R09、R10 | 增强机制比较或可逆性验证；不能仅凭这些改变主张首创 |

以上比较的核心来源是 [R04](https://doi.org/10.1126/scisignal.2002790)、[R10](https://doi.org/10.1038/s41467-018-05862-6) 与 [R16](https://doi.org/10.1038/s42003-025-07870-x)。

### 一个最直接的实例

R04 的 Figure 3A 就比较了 EBSS 饥饿、rapamycin 和 Torin1；其方法将 TFEB cDNA 亚克隆至 pEGFP-N1。Figure 3B 包含核/胞质比的时间过程，Figure 5 还检验内源 TFEB。因此，现有方案的相似性已经可以根据具体图号判断，不只是“这领域有人研究过”。文中的剂量和时间是**论文已使用条件**，本调研没有把它们直接宣布为本项目的最优条件。[R04 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3437338/)

## 5. 当前基础方案需讲清的细节，以及扩展时的解释边界

### 5.1 标签方向和序列版本

官方 [Addgene #38119](https://www.addgene.org/38119/) 为人 TFEB、pEGFP-N1、C 端 EGFP，页面注明删除 TFEB 终止密码子，历史参考号为 NM_007162.2。它适合作为已有构建的基准；本次未下载或比对全序列，不能直接照抄酶切位点和引物。

另两项记录分别是 [#99957](https://www.addgene.org/99957/)（人 S211A，C 端 Myc/His，无 GFP）和 [#79014](https://www.addgene.org/79014/)（小鼠 Tfeb、N 端 Flag、S142A/S211A）。它们不能当作同一种人 TFEB–EGFP 构建。

**我的判断：**后续优先以有直接 HeLa 定位证据的 C 端荧光构建为参考，比较后再决定。不能只因担心 C 端功能就断定 N 端标签更好，TFEB 的 N 端本身也参与调控；现有资料没有完成两种标签方向的系统性优劣比较。[R06](https://doi.org/10.1083/jcb.201209135)、[R12](https://doi.org/10.1038/s41598-021-00613-y)

### 5.2 “饥饿”和“mTOR 抑制”不是充分的方法描述

应写清究竟去掉氨基酸、血清、葡萄糖还是多种成分，并明确介质、时间与恢复方式。Torin1 作用于 mTOR 催化活性，不能由单药结果独立得出“效应只由 mTORC1 造成”。Rapamycin 在不同文献模型中的效应也不一致，不能与 Torin1 简单互换。[R03](https://doi.org/10.1038/emboj.2012.32)、[R05](https://doi.org/10.4161/auto.19653)、[R09](https://doi.org/10.1038/s41467-018-04849-7)、[R11](https://doi.org/10.1038/s41586-020-2444-0)

### 5.3 定位、转录和降解功能分层判断

| 拟观察结果 | 可直接支持的结论 | 若要作进一步主张，需要增加什么 |
|---|---|---|
| 融合蛋白核/胞质比变化 | 该构建在该条件下重新分布 | 内源 TFEB 定位核验，检验是否与内源行为一致 |
| 全细胞 anti-GFP WB | 融合蛋白表达及表观大小 | 核/质分级和纯度标记才可用分级 WB 核验分布 |
| CLEAR 报告或靶基因 RT-qPCR | 对应转录输出发生变化 | TFEB 遗传干预等支持其特异贡献 |
| 恰当对照下的 LC3 周转/货物降解 | 自噬相关过程或降解功能的证据 | 需按具体检测排除阻断、积累等替代解释 |

TFEB 可在溶酶体受损时作为补偿反应入核，所以“进核”不是“降解已经增强”的充分条件。功能证据层级可先看 [R02](https://doi.org/10.1126/science.1204592) 和 [R15](https://doi.org/10.1016/j.isci.2024.110432)；[M01 自噬检测指南](https://doi.org/10.1080/15548627.2020.1797280)仅作为多种互补检测的补充方法依据，本次没有逐项审阅整份指南。

### 5.4 表达量和成像纳入标准

**当前方案的质量建议**：记录表达强度、预先定义成像与纳入标准，避免只选择最亮细胞，并检查结论是否被极高表达细胞主导。内源 TFEB 染色可作为之后的独立核验，不列为当前新增必做实验。本次不设置没有依据的数值阈值，也没有确认课程拥有稳定细胞株。[R18](https://doi.org/10.1038/s41467-026-71001-1)

### 5.5 后续若增加突变体，需要回答什么问题

S211A 可用于检验该位点参与胞质滞留的作用；S142A/S138A 与核输出问题更直接相关。不同突变不能统一称为“持续激活型”。尤其 S211A 可继续响应 Torin1，并保留相当核输出能力。若选一个突变，应先写清比较假设，再设计预期结果。[R08](https://doi.org/10.1080/15548627.2016.1271514)、[R09](https://doi.org/10.1038/s41467-018-04849-7)、[R10](https://doi.org/10.1038/s41467-018-05862-6)

### 5.6 统计与结果解释

现行 `methods-cell.tex` 已要求背景校正、处理接近背景的胞质信号，并明确同一次实验中的细胞不算独立生物学重复。这些原则无需重新“修复”。后续需落实独立重复、每重复采样、图像分割及主要比较；若增加构建×处理或时间因素，应随设计调整统计方法，而不是机械沿用单因素 ANOVA。这里是统计设计建议，尚未进行样本量或效能计算。

## 6. 当前范围与后续完善

**用户已明确：当前不要求很强的创新性，采用基础方案；下列 B/C 是以后希望完善时的备选操作。** 后续扩展都基于已有文献，不属于本次必做项。

| 方案 | 核心问题 | 优点与新增负担 | 我的判断 |
|---|---|---|---|
| A：基础验证 | 自建 TFEB 荧光构建能否对饥饿/Torin1 产生预期定位响应？ | 克隆、表达、成像连贯；生物学问题较基础 | **当前保留的作业范围** |
| B：动态与可逆性 | 饥饿后 TFEB 的核富集能否在补回营养后逆转，变化过程如何？ | 增加时间点与恢复组；固定细胞分时取样可比较群体平均响应 | 后续若想完善，优先考虑 |
| C：一个位点的机制比较 | S211 位点改变是否影响基础分布及处理响应？ | 增加突变构建和测序、组数与统计复杂度；与分子克隆工作紧密衔接 | 之后再考虑 WT 与一个突变体的比较 |

方案 B 的直接依据为 [R10](https://doi.org/10.1038/s41467-018-05862-6) 和 [R16](https://doi.org/10.1038/s42003-025-07870-x)，C 为 [R04](https://doi.org/10.1126/scisignal.2002790) 和 [R08](https://doi.org/10.1080/15548627.2016.1271514)。

若未来采用 B，候选题目可写为（不替换当前题目）：

> **Construction of a fluorescent TFEB reporter and analysis of its reversible nutrient-dependent nucleocytoplasmic redistribution in HeLa cells**

固定细胞分时取样只能比较群体平均的再分布和营养恢复响应；若要追踪同一细胞的可逆轨迹，或估计核输入/输出速率，需要相应活细胞或专门转运实验。当前作业无需为此增加平台或实验组。上述题目和扩展思路留作未来参考，本轮不继续追求“未被报道”的组合。

## 7. 阅读优先级，以及两位作者的阅读衔接

当前基础方案先精读下表前两篇；其余四篇作为解释方法局限或后续完善时的参考，不要求为当前作业增加相应实验。

| 顺序 | 文献及重点 | 阅读后应能回答 |
|---|---|---|
| 1 | **R04，2012**：Fig. 3/5、构建与成像方法 | 我们三组定位实验与前人有多接近？哪一端加 GFP？ |
| 2 | **R03，2012**：Fig. 2/3/5 | 如何把药物、磷酸化和营养信号连成证据链？ |
| 3 | **R10，2018**：Fig. 1/3/4 | 营养恢复如何影响核输出？终点比值为什么不能证明单一转运机制？ |
| 4 | **R16，2025**：Fig. 1/2 与实验设计 | 时间过程多提供了什么信息？哪些复杂技术可以不纳入课程？ |
| 5 | **R18，2026**：结果第一节与 Supplementary Fig. 1 对应说明 | 如何避免表达量造成的定位混淆？ |
| 6 | **R08，2017**：Fig. 1 及相关结果 | 若选择 S211A，哪些预期合理，哪些是过度简化？ |

按项目已确定的计划分工，蒋思源的 B 部分可优先提取 R04/#38119/R08 的构建依据；宋子研的 A 部分可优先阅读 R03/R10/R16/R18，并用 R01/R02 补背景。两人需要共同固定融合方向和构建信息，避免克隆部分写 EGFP–TFEB，细胞与 WB 部分却沿用 TFEB–EGFP 的假设。上述是阅读建议，不代表作者已完成这些贡献。

## 8. 现有文件状态与下一步

本次回读确认：现行 `报告/references.bib` 已正确使用 Settembre 2012 的 EMBO J 条目和 Napolitano 2018；`methods-cloning.tex` 已分别解释两种融合方向的终止密码子逻辑；`methods-cell.tex` 已区分独立实验与细胞数；讨论也已区分定位与通量。原 Word 的错引和构建解释属于历史问题，**本轮没有再次改写现行正文或原始 Word**。

当前范围已按用户决定保持基础方案，无需再次选择 A/B/C。下一步先精读与三组基础比较最相关的 R04/R03，再固定带版本的 TFEB 序列、载体与融合方向，分别补克隆设计和细胞处理/采样参数。R10/R16 与突变体研究可帮助理解背景，时间过程、再喂养、突变体及额外功能实验继续放在后续完善清单。

未完成的核验包括：完整候选质粒/转录本比对、引物和读框、全部补充材料与原始图像、全部文章的更正/撤稿状态、课程资源和截止日期。本次没有开展实验、序列分析或科研计算；没有形成任何本人的实验结果。

## 9. 完整参考清单

以下 R01–R19 为本次纳入的原始研究；M01 为方法指南，不能计为独立 TFEB 机制发现。题名按出版记录列出，作者采用首作者等的简写，便于后续按正式引文格式补全。

1. **R01.** Sardiello M, et al. **A gene network regulating lysosomal biogenesis and function.** *Science*. 2009;325:473–477. [DOI 10.1126/science.1174447](https://doi.org/10.1126/science.1174447)；[PMID 19556463](https://pubmed.ncbi.nlm.nih.gov/19556463/)。
2. **R02.** Settembre C, et al. **TFEB links autophagy to lysosomal biogenesis.** *Science*. 2011;332:1429–1433. [DOI 10.1126/science.1204592](https://doi.org/10.1126/science.1204592)；[PMID 21617040](https://pubmed.ncbi.nlm.nih.gov/21617040/)。
3. **R03.** Settembre C, et al. **A lysosome-to-nucleus signalling mechanism senses and regulates the lysosome via mTOR and TFEB.** *EMBO J*. 2012;31:1095–1108. [DOI 10.1038/emboj.2012.32](https://doi.org/10.1038/emboj.2012.32)；[PMID 22343943](https://pubmed.ncbi.nlm.nih.gov/22343943/)。
4. **R04.** Roczniak-Ferguson A, et al. **The transcription factor TFEB links mTORC1 signaling to transcriptional control of lysosome homeostasis.** *Sci Signal*. 2012;5:ra42. [DOI 10.1126/scisignal.2002790](https://doi.org/10.1126/scisignal.2002790)；[PMID 22692423](https://pubmed.ncbi.nlm.nih.gov/22692423/)。
5. **R05.** Martina JA, et al. **MTORC1 functions as a transcriptional regulator of autophagy by preventing nuclear transport of TFEB.** *Autophagy*. 2012;8:903–914. [DOI 10.4161/auto.19653](https://doi.org/10.4161/auto.19653)；[PMID 22576015](https://pubmed.ncbi.nlm.nih.gov/22576015/)。
6. **R06.** Martina JA, Puertollano R. **Rag GTPases mediate amino acid-dependent recruitment of TFEB and MITF to lysosomes.** *J Cell Biol*. 2013;200:475–491. [DOI 10.1083/jcb.201209135](https://doi.org/10.1083/jcb.201209135)；[PMC3575543](https://pmc.ncbi.nlm.nih.gov/articles/PMC3575543/)。
7. **R07.** Medina DL, et al. **Lysosomal calcium signalling regulates autophagy through calcineurin and TFEB.** *Nat Cell Biol*. 2015;17:288–299. [DOI 10.1038/ncb3114](https://doi.org/10.1038/ncb3114)；[PMID 25720963](https://pubmed.ncbi.nlm.nih.gov/25720963/)。
8. **R08.** Vega-Rubin-de-Celis S, et al. **Multistep regulation of TFEB by MTORC1.** *Autophagy*. 2017;13:464–472. [DOI 10.1080/15548627.2016.1271514](https://doi.org/10.1080/15548627.2016.1271514)；[PMID 28055300](https://pubmed.ncbi.nlm.nih.gov/28055300/)。
9. **R09.** Li L, et al. **A TFEB nuclear export signal integrates amino acid supply and glucose availability.** *Nat Commun*. 2018;9:2685. [DOI 10.1038/s41467-018-04849-7](https://doi.org/10.1038/s41467-018-04849-7)；[PMID 29992949](https://pubmed.ncbi.nlm.nih.gov/29992949/)。
10. **R10.** Napolitano G, et al. **mTOR-dependent phosphorylation controls TFEB nuclear export.** *Nat Commun*. 2018;9:3312. [DOI 10.1038/s41467-018-05862-6](https://doi.org/10.1038/s41467-018-05862-6)；[PMID 30120233](https://pubmed.ncbi.nlm.nih.gov/30120233/)。
11. **R11.** Napolitano G, et al. **A substrate-specific mTORC1 pathway underlies Birt–Hogg–Dubé syndrome.** *Nature*. 2020;585:597–602. [DOI 10.1038/s41586-020-2444-0](https://doi.org/10.1038/s41586-020-2444-0)；[PMID 32612235](https://pubmed.ncbi.nlm.nih.gov/32612235/)。
12. **R12.** Park JY, et al. **A splicing variant of TFEB negatively regulates the TFEB-autophagy pathway.** *Sci Rep*. 2021;11:21119. [DOI 10.1038/s41598-021-00613-y](https://doi.org/10.1038/s41598-021-00613-y)；[PMID 34702966](https://pubmed.ncbi.nlm.nih.gov/34702966/)。
13. **R13.** Cui Z, et al. **Structure of the lysosomal mTORC1–TFEB–Rag–Ragulator megacomplex.** *Nature*. 2023;614:572–579. [DOI 10.1038/s41586-022-05652-7](https://doi.org/10.1038/s41586-022-05652-7)；[PMID 36697823](https://pubmed.ncbi.nlm.nih.gov/36697823/)。
14. **R14.** Malik N, et al. **Induction of lysosomal and mitochondrial biogenesis by AMPK phosphorylation of FNIP1.** *Science*. 2023;380:eabj5559. [DOI 10.1126/science.abj5559](https://doi.org/10.1126/science.abj5559)；[PMID 37079666](https://pubmed.ncbi.nlm.nih.gov/37079666/)。
15. **R15.** Zhao JF, et al. **Targeted dephosphorylation of TFEB promotes its nuclear translocation.** *iScience*. 2024;27:110432. [DOI 10.1016/j.isci.2024.110432](https://doi.org/10.1016/j.isci.2024.110432)；[PMID 39081292](https://pubmed.ncbi.nlm.nih.gov/39081292/)。
16. **R16.** Ruolo I, et al. **Investigation of dynamic regulation of TFEB nuclear shuttling by microfluidics and quantitative modelling.** *Commun Biol*. 2025;8:443. [DOI 10.1038/s42003-025-07870-x](https://doi.org/10.1038/s42003-025-07870-x)。
17. **R17.** Negoita F, et al. **AMPK promotes TFEB transcriptional activity through dephosphorylation at both MTORC1-dependent and -independent sites.** *Autophagy*. 2026;22:1386–1400. [DOI 10.1080/15548627.2026.2629720](https://doi.org/10.1080/15548627.2026.2629720)；[PMID 41661247](https://pubmed.ncbi.nlm.nih.gov/41661247/)。
18. **R18.** Xiong Y, et al. **TFEB degradation is regulated by an IKK/β-TrCP2 phosphorylation-ubiquitination cascade.** *Nat Commun*. 2026;17:4679. [DOI 10.1038/s41467-026-71001-1](https://doi.org/10.1038/s41467-026-71001-1)；[PMID 41922353](https://pubmed.ncbi.nlm.nih.gov/41922353/)。
19. **R19.** Esposito A, et al. **A TBK1/ULK1 signalling axis couples lysosomal stress to TFEB activation.** *Nature*. 在线发表 2026-09-30. [DOI 10.1038/s41586-026-11093-3](https://doi.org/10.1038/s41586-026-11093-3)；[PMID 42816627](https://pubmed.ncbi.nlm.nih.gov/42816627/)。本次检索未取得正式卷页，未补造。

**M01.** Klionsky DJ, et al. **Guidelines for the use and interpretation of assays for monitoring autophagy (4th edition).** *Autophagy*. 2021. [DOI 10.1080/15548627.2020.1797280](https://doi.org/10.1080/15548627.2020.1797280)；[PMID 33634751](https://pubmed.ncbi.nlm.nih.gov/33634751/)。本次阅读限摘要与多检测互补原则，不以此声称已核查具体条目。
