# Claude 独立审稿复查与措辞修订

日期：2026-10-25，America/New_York。起点：`main` 的 `5fcad00`（Codex 合并后的 32 页报告）。工作分支 `claude/reviewer-revision-20261025`。

## 范围与方法

- 用户要求：以审稿人角度判断实验设计是否合理、是否严谨；发现问题就补设计（必要时重新设计质粒）；删除 Codex 正文中的防御性写作；重画不符合学术风格的图；解决后 push。
- 全文 15 个 `sections/*.tex`、`main.tex` 与参考文献库逐段阅读。独立回算序列数字（HindIII 622、BamHI 2061、融合 ORF 2,169 nt、S142/S211/S466/S467/S469 的密码子与坐标、三个无标签构建 5,396 bp 且无 BamHI）：全部与报告一致，没有算术错误。
- 关键文献表述逐条回查 PubMed 摘要及可得的 PMC 全文，结果见下表 C14。
- 读取了用户与 Codex 的会话记录，以确认课程作业定位、工作量要求（“尽量不要让别人质疑”）和既有决定；未改变主线方案。
- 这是 Claude 的自审，不是外部同行评审，也不是作者签字。

## 发现的问题与处理

| ID | 审稿人会质疑的点 | 处理 | 仍待实验确认 |
|---|---|---|---|
| C01 | 主估计量 C 是“三因素交互”，而题目问题“3A 效应是否依赖融合架构”没有任何主要检验直接覆盖 | 新增 A_s = (3A−WT)_EGFP − (3A−WT)_untagged；主要估计量改为 A_M、A_T、C=A_M−A_T，Holm 校正；A_s 为 4 个细胞均值、C 为 8 个 | 方差与天数仍由 pilot 决定 |
| C02 | 文献中 Paquette（磷酸化激活）与 Negoita（去磷酸化激活）对 3A 的预测相反，报告没有写预期方向 | 研究延伸开头写出两种模型对 3A 的相反预测，并说明架构效应的表现形式；对比仍为双侧 | 无实测支持 |
| C03 | S467 同时是 Akt 位点（Palmieri 2017），两药差异不一定来自 AMPK | 引用 Palmieri 2017；两药下测 p-AKT(S473)；Akt 活性不同时在 WT/3A 上加药理学 Akt 抑制剂对比 | Akt 抑制剂剂量未定 |
| C04 | HeLa 缺 LKB1（Fogarty 2010），AMPK 结论的适用背景需要限定 | AMPK 模块和讨论明确 LKB1 缺失背景，用 p-ACC/p-AMPK(T172) 确认 MK-8722 作用；关键对比在表达 LKB1 的第二细胞系（如 HEK293 或 U2OS）重复 | 第二细胞系未设计 |
| C05 | 去掉连接肽和 EGFP 后无法归因于 EGFP 本身（报告多处承认但未设计） | **重新设计质粒**：TFEB-WT-linker 与 TFEB-3A-linker（TFEB 1–476 + GDPPVAT + TAA，5,417 bp、483 aa，HindIII/BamHI 保留，消化 1,439 + 3,978）。新增条件模块 R06（R02 出现架构效应后触发）：A_s 分解为 EGFP 分量和连接肽/终止位置分量，4 个对比一个 Holm 族；51 条件/天，最低 612 孔/天数 6 | 合成与测序未做 |
| C06 | 主终点 MCOLN1 在 HeLa 中可能没有 TFEB 依赖的药物诱导，主对比无信息 | pilot 只用 NT/WT 救援/KD 无救援孔，从 MCOLN1、LAMP1、CTSD 选一次最大 TFEB 依赖诱导者，不使用 3A/3T 数据，冻结；三者都不合格则 R02 不进行；MK-8722 不过门槛则只用 Torin1，主估计量退为 A_T，AMPK 模块取消 | 阈值待 pilot |
| C07 | C_d 由 8 个均值构成，六天的精度承诺没有方差依据 | 写出 Var(A)=4σ²、Var(C)=8σ²，天数 n ≥ (t·s_C/h)²，C 决定样本量 | s_C 由 pilot 估计 |
| C08 | 瞬转孔的 bulk RNA 混合了转染与未转染细胞，剂量难以匹配 | 每个构建测转染比例；匹配不达标时转为单拷贝同源整合细胞系（Flp-In 类）或低 MOI 慢病毒池 | 方案未设计 |
| C09 | TFE3 补偿只在“结果出现后”才触发 | pilot 前置比较 TFEB 单敲与 TFEB+TFE3 共敲；单敲仍保留至少一半诱导则采用共敲背景并冻结 | 一半为设计阈值 |
| C10 | 核心复现缺数值成功标准 | 三个主要对比：Holm 校正 p<0.05、日水平 log2(N/C) 变化至少 log2 1.5 且每个独立日同号，pilot 可修订阈值并记录；报告蛋白方向须与内源 IF 一致 | 阈值待 pilot |
| C11 | N/C 没有上限对照；报告蛋白与内源 TFEB 的相对丰度没有交代 | WT 加 60 min 莱普霉素 B（CRM1 抑制剂）定义核内上限，各基因型 N/C 同时报告为上限的比例；所选剂量与内源 TFEB 的比值 | 浓度待 pilot |
| C12 | 自噬流检测缺共识指南引用 | 引用 Klionsky 2021 | — |
| C13 | 转染用 DNA 为 miniprep，基因型与制备批次可能混杂 | 每个构建由两个独立克隆各制 midi 级去内毒素质粒，跨实验日交替使用；每批测 A260/280、A260/230、HindIII/BamHI 消化和连接处测序 | 试剂盒待定 |
| C14 | 引用与原文不符 | 见下节 | — |
| C15 | 图不符合学术风格或有渲染缺陷 | 见下节 | — |

所有阈值、剂量和样本量均为方案设计值，须在 pilot 后、确认数据收集前冻结；没有湿实验数据。

## 文献回查（PubMed 摘要及 PMC 全文）

- Palmieri 2017：摘要明确 Akt 在 Ser467 磷酸化 TFEB 并抑制其入核，且独立于 mTORC1。支持 C03。
- Fogarty 2010：摘要注明 HeLa 缺 LKB1。支持 C04。
- Klionsky 2021：自噬检测指南第 4 版。
- Negoita 2026（PMC13185439 全文）：确认 MK-8722 与 Torin1 使 C 端位点去磷酸化、AMPK 缺失取消 MK-8722 而非 Torin1 反应、GFP knock-in MEF 质谱；其引言转述 Paquette 的 GFP 标记 3A“不影响入核但取消转录活性”，支持报告原有叙述；Negoita 使用 S466T/S467T、S469T 作抗体特异性对照，并观察到 S469A 提高 p-S466/S467 信号。**全文中找不到“明确提出丙氨酸改变局部结构”**，已从正文删去该归因。
- Vega 2017：摘要只支持“S211A 提高核定位、TFEB 仍受 MTORC1 调控、S122 重要”，**不支持“S211A 保留 Torin1 响应”**，已改为摘要范围内的表述。
- Paquette 2021：只有摘要可得（PMC 全文为空）；“同时检验 S142A/S211A 组合”无法核实，已删去。
- Martina 2014（全文）：TFE3 过表达在 ARPE-19 中提高 MCOLN1 等基因，并在 TFEB 敲低 HeLa 中显著提高溶酶体基因表达。正文改为“过表达 TFE3 ……包括在 TFEB 敲低的 HeLa 细胞中”。
- 参考文献库新增 Palmieri2017、Klionsky2021、Fogarty2010；DOI 与卷页来自 PubMed 记录。

## 防御性写作清理

- 重写：摘要、背景（新增 Aims and Hypotheses）、克隆方法、对照构建、细胞方法、研究延伸、工作计划、预期结果、讨论、声明共 15 个正文文件中的限制性套话。
- 去掉的类型：“not a claim of / does not establish / cannot … ”式否定、反复出现的“设计选择而非实测”声明、“不是已验证克隆”式图注、“no claim of first discovery”。保留的是决定实验怎么做或怎么判读的规则，改写为肯定句（例如“独立培养日是生物学重复”，“结论以 …… 为准”）。
- 保留的诚实边界：报告是拟开展的设计，无湿实验数据；作者贡献与学号待核实；Working draft。它们在摘要、讨论和声明中各出现一次。
- 同时修正：AI 使用声明新增 Anthropic Claude（Claude Code）的实际用途。

## 图

- `control-schedule`、`research-design` 按期刊风格重画：去掉图内说明小字（移入图注）、增加图例/刻度/分组框线；研究设计图新增线性比例的构建示意、R06 条件臂和对比权重面板（权重由写明的定义推出）。
- 新增 `architecture-panel`：三种融合架构（untagged / linker-only / EGFP fusion）按序列文件读出的长度画出。
- 多面板质粒图原先把“Created by SnapGene”署名角标渲染成黑块（PDF→SVG→PDF 往返丢失透明蒙版），面板字母字体回退；现直接保存 PDF，字母嵌入 DejaVu Bold。
- 旧图备份于 `docs/备份/figures-before-academic-20261025/`（本地，不入 Git）。未使用 SnapGene 重新导出原生图；两张 linker-only 质粒的序列文件可直接导入 SnapGene，图集中只给架构示意图。

## 核验

- 新质粒：`scripts/B-figures/build_linker_controls.py` 生成，`scripts/B-audit/check_linker_controls.py` 独立重建期望序列并核对：9 项检查通过（整条 DNA、天然残基编号与翻译、单一 TAA 终止无内部终止、EGFP 完全移除、HindIII/BamHI 位点唯一、消化几何 1,439/3,978、FASTA/GenBank/蛋白一致、注释边界与骨架特征保留、输入输出哈希）。
- 条件计划：`check_research_condition_plan.py` 新增 C = A_M − A_T 的代数核对和 R06 条件/孔数（6 + 45 = 51；51×6×2 = 612），5 项通过。
- 报告一致性：`check_report_consistency.py` 新增 linker 构建与架构图绑定核对。
- 编译：Tectonic 0.17.0，34 页；最终日志无 Overfull、未定义引用或缺字。逐页目视记录见 `docs/协作/报告逐页核验.json`（34 页；其中 10 页为最终版直接查看，其余页的页图与此前直接查看的页图逐字节相同）。

## 没有验证的

- 没有湿实验；R06 构建、siRNA、qPCR 寡核苷酸、抗体性能、Akt 抑制剂剂量、莱普霉素 B 浓度和第二细胞系均未设计或验证。
- Paquette 2021 全文与 Vega 2017 全文未取得，对应表述只限摘要及 Negoita 引言转述的范围。
- 全球新颖性未证明；“未发现匹配比较”仅指本次定向检索。
- 本次复查不能保证没有其他质疑；作者署名、学号、实际贡献与课程提交要求仍待作者确认。
