# 工作进展

## 2026-10-23 12:43 EDT GitHub 协作与 ElegantPaper 初始化

### 这次做了什么

- 用户要求：建立公开 GitHub 仓库，以一份英文报告为最终交付；采用 ElegantPaper。作者为蒋思源、宋子研。宋子研尚未注册 GitHub，用户之后提供账号。
- 用户决定：蒋思源负责 B（分子克隆与质粒构建），宋子研负责 A（背景、细胞实验与课题叙事）。已按方法及预期结果拆分章节，并写入计划贡献；尚不代表两人的实际完成贡献。
- 建立 `报告/main.tex` 与 10 个章节文件、`references.bib`、编译配置；固定官方 ElegantPaper v0.12，原类文件未修改，保留 LPPL 许可证和来源。
- 原始 Word 提纲保留在 `资料/原始提纲/`。LaTeX 工作稿修正了终止密码子解释、Settembre 2012 引文和 Napolitano 2018 首作者；完善统计单位措辞。没有固定最终载体、融合方向或补造序列及实验结果。
- 新建 README、课程要求与待办、模板来源、同文 AGENTS/CLAUDE。公开仓库不收课程讲义、旧实验记录或学号。

### 现在真实状态

- 仓库地址：https://github.com/siyuanj/tfeb-experimental-design 。本次以 `--public` 创建，初始工作稿准备作为首个提交推送；推送与远端回查另行记录。
- 原提纲 SHA-256：`cd8df8b47ba31b6c3bf699bfd3f8e7b5253fba0cb564ad65a3d72ab582b78463`，与 MacBook 来源字节一致。
- 模板 upstream commit：`366d1e095a595b812b9a76bef1b6e245fbe417b3`；类文件 SHA-256：`29370455a0cdf0af94039a6483ba1433f5b067666771b0f96c4cca0159830173`。
- Codex 内置 LaTeX 编译器本次返回 success；Tectonic 0.17.0 本次导出 `报告/main.pdf`，4 页 A4、51942 字节，SHA-256 `895d05bc4fd24d9bb8e816d29ec10f67f01b104304a7fb6f1fb7f822d4f599c2`。
- 已查看 `docs/核查/page-1.png` 至 `page-4.png`：作者、正文、公式及参考文献完整，未见缺字、裁切或重叠。日志无 overfull/underfull、未定义引用或 missing character；有预期的 biblatex BibTeX fallback 提示。
- 自审范围：原提纲与迁移稿差别、作者分工、融合终止逻辑、引文元数据、图文状态、公开仓库范围与模板许可。AGENTS/CLAUDE 字节一致。

### 卡在哪

- 协作者权限等待宋子研注册 GitHub；下一步由用户提供准确账号后添加，不阻塞 B 的独立设计。
- 最终序列、载体与融合方向尚未固定；B 下一步核对来源并完成克隆设计，再把一致的构建信息交给 A。

### 还没验证的

- 当前 PDF 是模板与结构工作稿，未完成引物、Tm、完整质粒图、具体实验条件与全部文献核查，不能提交为最终报告。
- 尚未测试独立 XeLaTeX + BibTeX 或 Overleaf 构建；本次实际测试为内置编译器与 Tectonic。
- 课程截止日期及网络学堂当次提交要求、两人学号和正式英文署名未核实。没有提交作业或开展实验、科研计算。

### 要用户定的

- 宋子研注册后的 GitHub 用户名；其余身份信息在最终提交前补齐。当前仓库和模板工作无需等待。
