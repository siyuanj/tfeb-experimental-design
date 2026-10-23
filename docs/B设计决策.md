# B 设计决策与证据边界

日期：2026-10-23；版本 B-v1。以下技术取舍为我的判断，按用户已授权的 TFEB/HeLa、基础克隆设计范围执行；不是用户逐项指定，也不是湿实验结论。

## 唯一参考、模板和构型

- TFEB：人 `NM_007162.3`，CDS 303–1733（1-based，+链），1431 nt含TGA；`NP_009093.1`，476 aa。与旧文献 `.2` 的CDS/蛋白相同；mRNA只有末端10A差异。未声称这是当前MANE Select。原始记录、下载URL及哈希在 `资料/序列/source_manifest.json`，翻译校验在 `设计/reference_verification.json`。
- PCR模板：拟订制与上述CDS一致的1431nt线性dsDNA、不做密码子优化；采购例为IDT gBlocks HiFi（产品使用序列确认的克隆中间体，仍不保证本项目每个分子无错误）。不预设HeLa表达哪个转录本，不声称已获得cDNA或已测浓度。下单可行性/QC需将来核验；最终克隆须全覆盖测序。
- 空载体：NCBI Clontech直接提交 `U55762.1`，4733bp环状pEGFP-N1。该记录是完整设计参考，厂家说明其由公开资料及部分测序汇编，不能冒充本次实物全长测序。使用前需核实空载体身份。
- 方向：TFEB–EGFP（C端EGFP），与Roczniak-Ferguson 2012的HeLa饥饿/Torin1定位先例最直接一致。Park2021有N端EGFP先例；Martina2013显示N端调控重要，但两者不证明N端标签一定失活，也不构成标签方向的直接优劣实验。
- 命名解释：pEGFP-N1 将目标 TFEB 接在 EGFP 的 N 侧，因此 EGFP 实际位于 TFEB 的 C 端；与参考选择线的命名说明一致。
- native TFEB ATG保留；TFEB原TGA移除；EGFP的ATG成为内部Met，最终由EGFP TAA终止。整个蛋白必须从最终质粒翻译核验。

## 酶对、接头及引物

- 使用HindIII-HF/BamHI-HF：两位点在空载体唯一、TFEB内不存在，均为不同的5′黏端（AGCT/GATC）。XhoI在TFEB内部出现，弃选。HindIII/KpnI虽可用且有历史先例，本设计用更靠近EGFP的BamHI缩短MCS残余片段；不声称重建Addgene #38119原序列。
- 两引物各6nt保护区，是参照NEB末端切割表的保守选择；表中3–5bp已达+++，其意义为50–100%，不是本PCR片段切割率实测。
- 正向尾部 `ACTGAC|AAGCTT|ACC`；反向尾部 `TATGAC|GGATCC|CC`；最后CC是coding-strand GG的反向互补。`GG`与载体残余19nt组成21nt linker `GGGGATCCACCGGTCGCCACC`，编码GDPPVAT。
- 曾检查GCCACC起始尾部，完整正向寡核苷酸模型发卡Tm为66.001°C，超过本设计55°C筛选阈值。改用ACC，保留−3 A及+4 G；并未改变TFEB蛋白。候选及完整结构结果在 `设计/primer_checks.json`。
- 主引物退火区各22nt。Tm模型为Primer3/SantaLucia（50mM单价盐假设、2mM Mg、总dNTP0.8mM、每引物500nM）；显式区分Q5厂商未公开的完整缓冲配方。65°C及62/65/68°C梯度为拟优化起点，不是实测条件。5′非退火尾巴不用于初期退火Tm。
- 特异性核查仅针对已定义合成模板的完整退火区和末端12mer两链精确匹配；没有做转录组/基因组级特异性声称。

## 片段、坐标与功能注释

- pydna执行带尾PCR、实际错位黏端酶切、兼容端连接、环化。独立核查路径不调用主线构建函数，检查字面重建、密码子翻译及故意错误。
- 胶图和连接质量换算采用标称cut-to-cut单链长度；pydna全跨度另列，不能混用：insert1439nt/1443跨度，backbone4695nt/4699跨度，最终6134bp。
- TFEB631–2058；linker2059–2079；EGFP2080–2799（含stop）；完整ORF631–2799。TFEB476 + linker7 + EGFP239 = 722aa，理论未修饰平均质量80.425kDa。Western迁移、剪切和磷酸化可能使表观条带不同。
- 元件来源为U55762.1 CDS与Clontech PT3027-5 p2，按插入造成的偏移映射并核对其序列未变。polyA标注是信号基序；SV40早期启动子列离散元件，非推测的完整连续边界。原MCS被切开，保留两部分分别注明。

## 自审与尚未验证

- 写后执行前审查了脚本的版本/哈希、正反引物、0/1-based边界、片段选择、环化原点及特征偏移。发现并修正：不可用会在`python -O`消失的assert作为来源校验；不为非方向性origin猜箭头；源模板FASTA必须与GenBank CDS再比较。
- 首次主路径实际成功，未出现构建失败；随后独立P5核查6组通过，4类错误注入按预期原因拒绝。输入/脚本哈希回查后已发布sequence-ready；不将此计算验收称为实验验证。
- 未做湿实验、供应商合成可行性审核、真实库存验证、实际Sanger或细胞表达。标签对定位的影响需细胞对照核验；定位改变不能证明自噬通量增加。

## 关键来源

- NCBI原始版本记录与来源索引：`资料/B来源索引.tsv`。
- Roczniak-Ferguson 2012：doi:10.1126/scisignal.2002790，Methods / Fig.3；Addgene38119标签/版本/克隆信息。
- Clontech PT3027-5 PR93634，1999-03-19，第2页坐标、第3页序列限制；原件仅本地保存 `docs/核查/B/sources/Clontech_PT3027-5_PR93634.pdf`。
- 厂家反应条件与服务指南：`docs/协作/方法来源.tsv`，由主线汇入总索引，未使用候选条件冒充定稿。
