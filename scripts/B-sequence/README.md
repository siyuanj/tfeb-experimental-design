# B 序列设计复算

日期：2026-10-23。所有命令在 TFEB 仓库根目录执行，使用既有本地环境。无需 UF、GPU、网页服务或 SnapGene 许可。

## 已核验环境

Python 3.12.15；pydna5.5.16；Biopython1.88；primer3-py2.4.0；DNA Features Viewer3.1.5；Matplotlib3.11.2；macOS27.0 ARM64。报告编译为Tectonic0.17.0。

解释器：`/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python`。原环境及安装脚本不修改。

## 运行顺序

1. `retrieve_sources.py`：从NCBI取明确版本记录，若已有原始文件字节发生变化立即拒绝覆盖。网络再取可因数据库注释更新改变文件字节，须人工比较；不应盲目接受新哈希。日常离线复算使用已保存原件，不必重新下载。
2. `prepare_spec.py`：核对四条原始记录并提取CDS/蛋白；写B-v1明确设计参数和模板方案。**会重写spec及派生参考FASTA**；修改输入后必须撤回ready并重验下游。
3. `build_construct.py`：Primer3模型筛查、pydna带尾PCR、双酶切、实际黏端连接、ORF翻译、元件映射；生成主路径核查。会重写本阶段派生序列/表/JSON，不修改原始GenBank。
4. `independent_verify.py`：标准库独立字面重建、密码子表翻译、质量求和及错误注入；不导入构建脚本。仅主路径成功不能发布ready。
5. `update_source_index.py`：把已复核主线/方法来源汇入统一台账；原说明书/论文全文位于Git忽略的核查目录。

示例：

```bash
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-sequence/prepare_spec.py
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-sequence/build_construct.py
'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-sequence/independent_verify.py
```

所有验收以真实输出的status、逐项checks、输入/脚本哈希为准。标准输出或退出码0本身不足以证明构建正确。最终图表及方法还需各自核查。

## 序列和长度定义

`construct_spec.json`为唯一设计选择入口；脚本中的参数默认是该版本的显式决策，变更方案需先修spec生成器与独立判据。原始CDS1431nt含TGA，实际TFEB编码1428nt不含stop；PCR1457bp含完整引物尾巴；切后1439nt标称插入不同于TFEB编码长。

表格坐标1-based inclusive，内部切片0-based half-open。pydna的4699/1443为包含错开黏端的双链跨度；连接及预测胶使用4695/1439标称单链长度。构建6134bp，融合蛋白722aa（含7aa linker），理论平均质量80.425kDa；不是Western实测迁移值。

`record_status.py`只更新交接文档：备份、读取后比较、前插并核验旧正文保留。它拒绝已观察到的并发冲突；发生写入后的冲突需重新读取合并，不能恢复旧备份覆盖同伴。
