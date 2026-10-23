# B 设计决策：TFEB、载体与融合方向

日期：2026-10-23。当前范围：基础方案的参考选择；不是完整克隆设计验收。

## 已选定的主方案

**人野生型 TFEB isoform 1（NM_007162.3 / NP_009093.1）＋ pEGFP-N1（U55762.1）＋ C 端 EGFP，即 TFEB–EGFP。**

用户决定保持基础方案，不增加必须完成的突变体或创新实验。上述具体构建选择由本轮文献与序列核对作出。后续计算读取 [construct_spec.json](../设计/construct_spec.json)；实测核查见 [reference_verification.json](../设计/reference_verification.json)。

| 项目 | 固定选择 | 核查结果或边界 |
|---|---|---|
| 物种、基因、蛋白 | Homo sapiens，TFEB，wild-type isoform 1 | 不使用小鼠、截短体、S211A 或双突变体 |
| 转录本 | NCBI RefSeq **NM_007162.3**，transcript variant 1 | 本轮下载完整 2,354 nt mRNA；不把整条 mRNA 当插入片段 |
| 蛋白 | **NP_009093.1，476 aa** | CDS 翻译与独立蛋白记录逐字一致 |
| CDS | 转录本正向 **303–1733**，1-based、包含两端 | 1,431 nt，含起始 ATG 和末端 TGA |
| TFEB 编码部分 | **1,428 nt，不含末端 TGA** | 保留完整 476 aa；不是最终带尾 PCR 产物长度 |
| 空载体参考 | **pEGFP-N1，GenBank U55762.1** | 厂商提交的完整 4,733 bp 环状参考序列 |
| 标签位置 | **TFEB–EGFP，EGFP 在 TFEB 的 C 端** | N→C 顺序为 TFEB、接头、EGFP |
| 表达及选择 | CMV 驱动；细菌 KanR；哺乳动物 Neo/G418 | G418 不用于细菌筛选；暂不为瞬时表达增加稳定筛选实验 |
| 文献构建 | **pEGFP-N1-TFEB，Addgene #38119** | 作为已发表架构参照，不声称已取得或逐碱基复刻成品质粒 |

序列原始记录：[NM_007162.3](https://www.ncbi.nlm.nih.gov/nuccore/NM_007162.3)、[NP_009093.1](https://www.ncbi.nlm.nih.gov/protein/NP_009093.1)、[U55762.1](https://www.ncbi.nlm.nih.gov/nuccore/U55762.1)。表达、标签及选择信息由 [Addgene #38119](https://www.addgene.org/38119/) 支持。

## 为什么选这个 TFEB 版本

Addgene #38119 指向历史 **NM_007162.2**。本轮实际取得 `.2`、`.3` 两版 GenBank 和独立 `NP_009093.1` 蛋白记录后比较：两版 CDS 坐标均为 303–1733，核酸编码区和 476 aa 蛋白完全相同；`.2` 的 2,364 nt 转录本等于 `.3` 的全部 2,354 nt 加末端 10 个 A。因此采用当前 `.3` 可保持与该历史编码序列一致，不引入蛋白改变。不是只根据 accession 名称或长度推定一致。

选择 variant 1 是为对齐既有构建的编码序列，不主张它是唯一合理的 TFEB 转录本，也不称它为 MANE Select。它与长 N 端 isoform 2 必须区分；新旧记录均明确标注 isoform 1。

## 为什么选 pEGFP-N1 和 C 端 EGFP

[Roczniak-Ferguson et al., 2012](https://doi.org/10.1126/scisignal.2002790) 的 *Cell Culture and Transfection* 方法记载 TFEB cDNA 亚克隆进 pEGFP-N1，使用 HindIII/KpnI；Figure 3 直接展示 HeLa 中 TFEB–GFP 对饥饿和 Torin1 的定位响应。[Addgene #38119](https://www.addgene.org/38119/) 明确写 EGFP 位于 insert 的 C 端，并注明删除 TFEB stop codon。

因此，本次沿用有直接 HeLa 定位依据的架构，便于解释和对照已有研究。原稿的 pEGFP-C1 / EGFP–TFEB 不再作为主方案；这不表示 N 端标签一定不可用，也没有证明 C 端标签完全不干扰 TFEB。过表达与标签本身仍是报告中的限制。

**名称容易混淆：pEGFP-N1 用于将目标蛋白接在 EGFP 的 N 侧，所以本项目中 EGFP 实际位于 TFEB 的 C 端。** 最终报告统一写 TFEB–EGFP，避免只写“GFP–TFEB”而让读者误解顺序。

## 参考空载体与成品质粒分开

- **U55762.1** 是 Clontech 于 1996 年提交的 pEGFP-N1 空载体参考，4,733 bp；原记录说明它综合数据库、文献及部分测序信息。它不是本实验室实物库存的测序结果，旧记录中的订购说明不代表今天仍有库存。
- **#38119** 已含 TFEB。官方序列页的下载/复制需要登录；本轮没有获取其完整序列。因此不把官方网页的成品质粒总长度当成本项目最终构建长度，也不承诺连接区与该成品完全相同。
- 当前作业保留“TFEB PCR → 空载体克隆 → 筛选 → Sanger”的设计内容；现成成品仅作文献参照。
- 当前空载体 GenBank 主要注释 EGFP 和 neomycin phosphotransferase 两个 CDS；完整 promoter、MCS、复制起点、poly(A) 等特征坐标仍须结合厂商资料补齐，不能把稀疏原注释称为完整质粒图。

## 模板与连接规则

为让后续 PCR 有一个明确输入，设计假定使用**按 NM_007162.3 CDS 合成、保留天然末端 TGA 的人 TFEB cDNA 模板克隆**。这是拟获取材料，不是已经采购或测序的样品；实际投入实验前须以供应商序列及测序核对。模板采用天然编码序列，不作密码子优化。

拟定表达顺序为 `CMV → TFEB 编码区（1,428 nt，含起始 ATG、不含末端 TGA）→ 待核查接头 → EGFP → 最终 stop`。TFEB 的 ATG 用作融合蛋白起始，不在 1,428 nt 之外额外添加一次；TFEB 末端 TGA 不进入融合构建。计划保留完整 EGFP 编码区及其末端 TAA，EGFP 的 ATG 在融合中编码内部 Met。翻译起始环境与接头碱基在引物设计阶段固定；不能仅凭此顺序就声称读框通过。

HindIII/KpnI 目前只是文献优先候选酶对；本轮未验证其内部位点、切割与连接后读框，也未设计引物。最终 TFEB–EGFP 的氨基酸数、理论分子量和质粒长度需在接头确定后计算，不能简单把参考长度相加。

## 本轮验证和下一步

本地 Python 3.12.15 / Biopython 1.88 实际检查：版本标识、CDS 翻译、终止、旧版一致性、独立蛋白一致性、载体环状和长度、载体 EGFP 翻译；均通过。脚本先自审并经独立静态审阅，无需 UF。两个参考的序列及文件 SHA-256 在 JSON 中分别保存，注释文件更新不等于序列改变。

下一步只推进基础克隆设计：补齐骨架特征，确定两端接头及酶对，设计引物，再重建和翻译完整构建。未核验实物模板、实物空载体、完整成品质粒、实际 PCR、表达或定位；本轮不能标为整个 B 已完成。
