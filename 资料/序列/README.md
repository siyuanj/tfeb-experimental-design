# 参考序列与复算

获取及核查日期：2026-10-23。只包含公开数据库参考与其衍生序列，不是实物测序结果。

原始 GenBank/GenPept 保留 NCBI 返回的全部字节，包括格式自带的行尾空格和末尾空行；`.gitattributes` 对这些记录关闭换行转换和空白告警。原创文档与 TSV 仍执行通常的空白检查。

| 文件 | 来源或生成方式 |
|---|---|
| `TFEB_reference.gb` | NCBI NM_007162.3，原始 GenBank |
| `TFEB_NM_007162.2.gb` | NCBI 历史版本，用于确认与 Addgene 页面引用的历史 RefSeq 编码区一致 |
| `TFEB_NP_009093.1.gp` | NCBI 独立蛋白 GenPept 记录 |
| `pEGFP-N1_U55762.1.gb` | NCBI 厂商提交空载体参考，原始 GenBank |
| `TFEB_cds.fasta` | 从 `.3` 注释提取 303–1733，1,431 nt，含 TGA |
| `TFEB_cds_no_stop.fasta` | 上述 CDS 去末端 TGA，1,428 nt，不含接头/酶切尾巴 |
| `TFEB_protein.fasta` | 翻译并核对 NP_009093.1，476 aa |
| `pEGFP-N1_U55762.1.fasta` | 空载体原序列，4,733 bp |

## 本轮下载命令

在项目根目录执行。重新获取时先保存为新文件并比较，不覆盖冻结参考：版本相同的注释记录仍可能更新。

```bash
curl --fail --location 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_007162.3&rettype=gb&retmode=text' -o 资料/序列/TFEB_reference.gb
curl --fail --location 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_007162.2&rettype=gb&retmode=text' -o 资料/序列/TFEB_NM_007162.2.gb
curl --fail --location 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_009093.1&rettype=gp&retmode=text' -o 资料/序列/TFEB_NP_009093.1.gp
curl --fail --location 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U55762.1&rettype=gb&retmode=text' -o 资料/序列/pEGFP-N1_U55762.1.gb
```

## 离线复算

现有课程环境：Python 3.12.15 / Biopython 1.88；无需网络、UF 或 GPU。

```bash
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/verify_reference_choices.py
```

其他设备可在含 Biopython 1.88 的 Python 3.12 环境运行同一脚本。脚本读取四个原始记录、生成四个 FASTA 和 `设计/reference_verification.json`；不会生成引物、接头或最终质粒。JSON 包含文件字节 SHA-256 与纯大写序列 SHA-256，两者不能混用。失败即报错，不把部分输出当成通过。
