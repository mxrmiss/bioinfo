# bioinfo

一个面向科研日常的生信工具库与流程集合（workflow-first）。  
按子项目分舱管理：系统发育/比较基因组、RNA-seq、Swiss-Prot 注释、RBH 分歧景观、可视化模板与常用工具。

---

## What’s inside（你能用它做什么）

- **Phylogenomics**：OrthoFinder 同源聚类、严格单拷贝集合、MSA/树与下游接口（`phylo/`）
- **Codon-level downstream**：codon MSA 质控、PAML/codeml 等（`aphylo/`）
- **RNA-seq**：差异、富集、WGCNA 等（`rna-seq/`）
- **Swiss-Prot annotation**：DIAMOND best hit 注释与结构化汇总（`swissprot/`）
- **Divergence landscape (RBH)**：RBH 同源对驱动的 substitution/indel 统计与染色体/滑窗景观图（`snp_rbh/`）
- **Plotting kit**：顶刊风格作图脚本与模板（`magic/`）
- **Utilities**：常用小工具（`tools/`）

---

## Quick Start（从这里选一条路线）

按你的目标点进去读“子项目 README”，每个子项目都把入口、输入、输出讲得更细：

- 比较基因组 / 系统发育：[`phylo/`](phylo/)  
  入口说明：`phylo/README.phylo.txt`

- aphylo 下游（codon MSA / codeml）：[`aphylo/`](aphylo/)  
  入口说明：`aphylo/README.aphylo.txt`

- RNA-seq 分析：[`rna-seq/`](rna-seq/)  
  入口说明：`rna-seq/README.rna-seq.txt`

- Swiss-Prot 注释：[`swissprot/`](swissprot/)  
  入口说明：`swissprot/README/readme.md`

- RBH 分歧景观：[`snp_rbh/`](snp_rbh/)  
  入口说明：`snp_rbh/README/readme.md`

- 作图模板：[`magic/`](magic/)  
  入口说明：`magic/README_CN.txt`

---

## Repository layout（目录结构速览）

| Path | Purpose |
|------|---------|
| `phylo/` | 多物种 phylogenomics 主流程与发布包接口 |
| `aphylo/` | 基于发布包的 codon-level 深层分析 |
| `rna-seq/` | RNA-seq 分析工作流（差异/富集/WGCNA 等） |
| `swissprot/` | Swiss-Prot 注释（DIAMOND 比对 + best hit 解析） |
| `snp_rbh/` | RBH 同源对 → 对齐/回译 → 分歧统计 → 景观图 |
| `magic/` | R 作图“魔法包”（模板与脚本） |
| `tools/` | 辅助脚本集合 |
| `bak/` | 历史/备份（通常不作为主入口） |

---

## Conventions（通用约定）

大多数子项目遵循类似结构：

- `data/`：输入数据或索引
- `scripts/`：脚本与流程入口
- `results/`：主要产物（TSV/图/汇总）
- `logs/`：运行日志（便于复现与排错）
- `config.yaml`（如存在）：集中管理关键参数

---

## Dependencies（依赖概览）

不同子项目依赖不同，这里只列“最常见的一组”供你快速判断是否具备运行基础：

- **General**：Python 3.x，R
- **Alignment / Annotation**：DIAMOND，MAFFT，PAL2NAL（Perl）
- **Phylogenomics**：OrthoFinder（及其依赖）
- **Plotting (R)**：常用以 `ggplot2` 为核心（具体以子项目 README 为准）

建议按子项目分别维护环境，避免依赖互相干扰。

---

## Reproducibility notes（复现与写论文的小提示）

如果你在论文/报告中复用本仓库产物，建议在方法部分至少写清三件事：

1) 使用的子项目（如 `swissprot/`、`snp_rbh/`）  
2) 关键软件（如 DIAMOND/MAFFT/PAL2NAL/OrthoFinder 等）  
3) 核心阈值/参数来自哪里（对应 `config.yaml` 或子项目说明文档）

---

## License

本项目采用 **GNU GPL v3.0** 授权，详见根目录 [`LICENSE`](LICENSE)。
如需发布/分发修改后的版本，必须同样以 GPLv3 开源发布源码。


